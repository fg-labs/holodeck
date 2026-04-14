use std::collections::HashMap;
use std::ops::Index;

use noodles::fasta;

/// Metadata for a single reference sequence (contig).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceMetadata {
    /// 0-based positional index.
    index: usize,
    /// Contig name.
    name: String,
    /// Contig length in bases.
    length: usize,
}

impl SequenceMetadata {
    /// Create a new `SequenceMetadata` with the given index, name, and length.
    #[must_use]
    pub fn new(index: usize, name: String, length: usize) -> Self {
        Self { index, name, length }
    }

    /// Return the 0-based index of this sequence.
    #[must_use]
    pub fn index(&self) -> usize {
        self.index
    }

    /// Return the name of this sequence.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Return the length of this sequence in bases.
    #[must_use]
    pub fn length(&self) -> usize {
        self.length
    }
}

/// A lookup table mapping reference sequence names to their indices and lengths.
///
/// Constructed from a FASTA index, this type consolidates all name-index-length
/// mapping into a single shared structure.
#[derive(Debug, Clone)]
pub struct SequenceDictionary {
    /// Sequences in index order.
    sequences: Vec<SequenceMetadata>,
    /// Name to index mapping.
    name_to_index: HashMap<String, usize>,
}

impl SequenceDictionary {
    /// Return the number of sequences.
    #[must_use]
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    /// Return `true` if the dictionary contains no sequences.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    /// Look up a sequence by its 0-based index.
    #[must_use]
    pub fn get_by_index(&self, index: usize) -> Option<&SequenceMetadata> {
        self.sequences.get(index)
    }

    /// Look up a sequence by name.
    #[must_use]
    pub fn get_by_name(&self, name: &str) -> Option<&SequenceMetadata> {
        self.name_to_index.get(name).map(|&i| &self.sequences[i])
    }

    /// Iterate over all sequences in index order.
    pub fn iter(&self) -> impl Iterator<Item = &SequenceMetadata> {
        self.sequences.iter()
    }

    /// Return contig names in index order.
    #[must_use]
    pub fn names(&self) -> Vec<&str> {
        self.sequences.iter().map(|s| s.name.as_str()).collect()
    }

    /// Return the total length of all sequences combined.
    #[must_use]
    pub fn total_length(&self) -> u64 {
        self.sequences.iter().map(|s| s.length as u64).sum()
    }

    /// Build a dictionary from a pre-constructed list of entries.
    ///
    /// Useful for testing and for cases where the dictionary is constructed
    /// from sources other than a FASTA index. Entries are re-indexed by their
    /// position in the vector (the `index` field of each `SequenceMetadata`
    /// is ignored and replaced).
    #[must_use]
    pub fn from_entries(mut sequences: Vec<SequenceMetadata>) -> Self {
        for (i, seq) in sequences.iter_mut().enumerate() {
            seq.index = i;
        }
        let name_to_index =
            sequences.iter().enumerate().map(|(i, s)| (s.name.clone(), i)).collect();
        Self { sequences, name_to_index }
    }
}

impl Index<usize> for SequenceDictionary {
    type Output = SequenceMetadata;

    fn index(&self, index: usize) -> &SequenceMetadata {
        &self.sequences[index]
    }
}

impl Index<&str> for SequenceDictionary {
    type Output = SequenceMetadata;

    fn index(&self, name: &str) -> &SequenceMetadata {
        let i = self.name_to_index[name];
        &self.sequences[i]
    }
}

impl From<&fasta::fai::Index> for SequenceDictionary {
    fn from(index: &fasta::fai::Index) -> Self {
        let mut sequences = Vec::new();
        let mut name_to_index = HashMap::new();

        for (i, record) in index.as_ref().iter().enumerate() {
            let name = String::from_utf8_lossy(record.name().as_ref()).to_string();
            #[expect(clippy::cast_possible_truncation, reason = "FASTA index lengths fit in usize")]
            let length = record.length() as usize;
            name_to_index.insert(name.clone(), i);
            sequences.push(SequenceMetadata { index: i, name, length });
        }

        Self { sequences, name_to_index }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a dictionary from a list of (name, length) pairs for testing.
    fn make_dict(contigs: &[(&str, usize)]) -> SequenceDictionary {
        let sequences: Vec<SequenceMetadata> = contigs
            .iter()
            .enumerate()
            .map(|(i, &(name, length))| SequenceMetadata {
                index: i,
                name: name.to_string(),
                length,
            })
            .collect();
        let name_to_index = sequences.iter().map(|s| (s.name.clone(), s.index)).collect();
        SequenceDictionary { sequences, name_to_index }
    }

    #[test]
    fn test_len_and_is_empty() {
        let dict = make_dict(&[("chr1", 1000), ("chr2", 2000)]);
        assert_eq!(dict.len(), 2);
        assert!(!dict.is_empty());

        let empty = make_dict(&[]);
        assert!(empty.is_empty());
        assert_eq!(empty.len(), 0);
    }

    #[test]
    fn test_get_by_index() {
        let dict = make_dict(&[("chr1", 1000), ("chr2", 2000)]);
        let meta = dict.get_by_index(0).unwrap();
        assert_eq!(meta.name(), "chr1");
        assert_eq!(meta.length(), 1000);
        assert_eq!(meta.index(), 0);

        assert!(dict.get_by_index(2).is_none());
    }

    #[test]
    fn test_get_by_name() {
        let dict = make_dict(&[("chr1", 1000), ("chr2", 2000)]);
        let meta = dict.get_by_name("chr2").unwrap();
        assert_eq!(meta.index(), 1);
        assert_eq!(meta.length(), 2000);

        assert!(dict.get_by_name("chrZ").is_none());
    }

    #[test]
    fn test_index_by_usize() {
        let dict = make_dict(&[("chr1", 1000), ("chr2", 2000)]);
        assert_eq!(dict[0].name(), "chr1");
        assert_eq!(dict[1].name(), "chr2");
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn test_index_by_usize_out_of_bounds() {
        let dict = make_dict(&[("chr1", 1000)]);
        let _ = &dict[5];
    }

    #[test]
    fn test_index_by_str() {
        let dict = make_dict(&[("chr1", 1000), ("chr2", 2000)]);
        assert_eq!(dict["chr1"].length(), 1000);
    }

    #[test]
    #[should_panic(expected = "no entry found for key")]
    fn test_index_by_str_unknown() {
        let dict = make_dict(&[("chr1", 1000)]);
        let _ = &dict["nope"];
    }

    #[test]
    fn test_names() {
        let dict = make_dict(&[("chr1", 1000), ("chr2", 2000), ("chrX", 500)]);
        assert_eq!(dict.names(), vec!["chr1", "chr2", "chrX"]);
    }

    #[test]
    fn test_total_length() {
        let dict = make_dict(&[("chr1", 1000), ("chr2", 2000), ("chrX", 500)]);
        assert_eq!(dict.total_length(), 3500);
    }

    #[test]
    fn test_iter() {
        let dict = make_dict(&[("chr1", 1000), ("chr2", 2000)]);
        let names: Vec<&str> = dict.iter().map(SequenceMetadata::name).collect();
        assert_eq!(names, vec!["chr1", "chr2"]);
    }
}
