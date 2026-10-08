#!/bin/zsh
set -euo pipefail

# launchd uses a minimal PATH; include Homebrew for Git, Git LFS, and GitHub CLI.
export PATH="/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"

readonly SOURCE_DIR="/Users/alhilaba/UG4_promesh_ogrid"
readonly BACKUP_DIR="/Users/alhilaba/UG4_promesh_ogrid_backup"
readonly BRANCH="ug4-build"
readonly GIT_BIN="/opt/homebrew/bin/git"

if [[ ! -d "$SOURCE_DIR" || ! -d "$BACKUP_DIR/.git" ]]; then
  print -u2 "UG4 source or backup repository is missing."
  exit 1
fi

if [[ ! -x "$GIT_BIN" ]]; then
  print -u2 "Homebrew Git is missing: $GIT_BIN"
  exit 1
fi

# Flatten the UG4 checkout without copying any of its nested Git databases.
# Preserve metadata and documentation maintained only by the backup repository.
rsync -a --delete \
  --exclude='.git' \
  --exclude='.gitattributes' \
  --exclude='backup_ug4.sh' \
  --exclude='BACKUP_DOCUMENTATION' \
  "$SOURCE_DIR/" "$BACKUP_DIR/"

"$GIT_BIN" -C "$BACKUP_DIR" add --all

if "$GIT_BIN" -C "$BACKUP_DIR" diff --cached --quiet; then
  print "No UG4 changes to back up."
  exit 0
fi

"$GIT_BIN" -C "$BACKUP_DIR" commit \
  -m "UG4 backup $(date '+%Y-%m-%d %H:%M:%S %z')"

if "$GIT_BIN" -C "$BACKUP_DIR" remote get-url origin >/dev/null 2>&1; then
  "$GIT_BIN" -C "$BACKUP_DIR" push origin "$BRANCH"
else
  print "Snapshot committed locally; no GitHub remote is configured yet."
fi
