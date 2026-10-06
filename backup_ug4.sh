#!/bin/zsh
set -euo pipefail

# launchd uses a minimal PATH; include Homebrew for git-lfs and GitHub CLI.
export PATH="/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"

readonly SOURCE_DIR="/Users/alhilaba/UG4_promesh_ogrid"
readonly BACKUP_DIR="/Users/alhilaba/UG4_promesh_ogrid_backup"
readonly BRANCH="ug4-build"

if [[ ! -d "$SOURCE_DIR" || ! -d "$BACKUP_DIR/.git" ]]; then
  print -u2 "UG4 source or backup repository is missing."
  exit 1
fi

# Flatten the UG4 checkout without copying any of its nested Git databases.
# Preserve the backup repository metadata and files maintained only by it.
rsync -a --delete \
  --exclude='.git' \
  --exclude='.gitattributes' \
  --exclude='backup_ug4.sh' \
  "$SOURCE_DIR/" "$BACKUP_DIR/"

git -C "$BACKUP_DIR" add --all

if git -C "$BACKUP_DIR" diff --cached --quiet; then
  print "No UG4 changes to back up."
  exit 0
fi

git -C "$BACKUP_DIR" commit -m "UG4 backup $(date '+%Y-%m-%d %H:%M:%S %z')"

if git -C "$BACKUP_DIR" remote get-url origin >/dev/null 2>&1; then
  git -C "$BACKUP_DIR" push origin "$BRANCH"
else
  print "Snapshot committed locally; no GitHub remote is configured yet."
fi
