#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HOOKS_DIR="${REPO_ROOT}/scripts/hooks"
GIT_HOOKS_DIR="${REPO_ROOT}/.git/hooks"

echo "Installing git hooks from ${HOOKS_DIR} to ${GIT_HOOKS_DIR}..."

for hook in "${HOOKS_DIR}"/*; do
    hook_name="$(basename "${hook}")"
    target="${GIT_HOOKS_DIR}/${hook_name}"
    cp "${hook}" "${target}"
    chmod +x "${target}"
    echo "  Installed ${hook_name}"
done

echo "Done."
