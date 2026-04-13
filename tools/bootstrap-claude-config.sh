#!/usr/bin/env bash
# Bootstrap shared Claude Code config from ghobriallab/compbio-commons.
#
# Usage (run from the root of your analysis repo):
#   bash <(curl -fsSL https://raw.githubusercontent.com/ghobriallab/compbio-commons/main/tools/bootstrap-claude-config.sh)
#
# Or clone and run locally:
#   bash /path/to/common-base-scripts/tools/bootstrap-claude-config.sh

set -euo pipefail

COMMON_REPO="https://github.com/ghobriallab/compbio-commons.git"
BRANCH="main"
TMPDIR_BOOTSTRAP=$(mktemp -d)

cleanup() { rm -rf "$TMPDIR_BOOTSTRAP"; }
trap cleanup EXIT

# Detect repo root
REPO_ROOT=$(git rev-parse --show-toplevel 2>/dev/null || pwd)
echo "Target repo root: $REPO_ROOT"

# Confirm before proceeding
read -r -p "This will copy .claude/ and tools/ from common-base-scripts into $REPO_ROOT. Continue? [y/N] " confirm
case "$confirm" in
  [yY]*) ;;
  *) echo "Aborted."; exit 0 ;;
esac

echo "Fetching common-base-scripts..."
git clone --depth 1 --branch "$BRANCH" "$COMMON_REPO" "$TMPDIR_BOOTSTRAP/source"

# Copy .claude/
mkdir -p "$REPO_ROOT/.claude/commands"
cp -r "$TMPDIR_BOOTSTRAP/source/.claude/commands/" "$REPO_ROOT/.claude/commands/"
cp "$TMPDIR_BOOTSTRAP/source/.claude/settings.json" "$REPO_ROOT/.claude/settings.json"
echo "Copied .claude/"

# Copy tools/
mkdir -p "$REPO_ROOT/tools"
cp -r "$TMPDIR_BOOTSTRAP/source/tools/" "$REPO_ROOT/tools/"
echo "Copied tools/"

echo ""
echo "Done! Next steps:"
echo "  1. Review the copied files (git diff or git status)"
echo "  2. Commit: git add .claude/ tools/ && git commit -m 'chore: add shared Claude config from common-base-scripts'"
echo "  3. Add this repo to .github/sync-targets.txt in common-base-scripts to receive future updates automatically."
