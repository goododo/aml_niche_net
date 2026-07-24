#!/usr/bin/env bash
# daily_commit.sh ----
# Daily snapshot: commit the day's changes to code / small tables / figures / docs to the local git.
# INPUT  : current working-tree state
# OUTPUT : one commit (skipped if no change); log appended to 00_project/git_autocommit.log
# Usage  : bash scripts/99_admin/daily_commit.sh          # manual
#          crontab runs it daily at 23:30 (see note at end of file)
# [DECISION] Commit only, do not push, unless a remote named 'origin' is configured -- avoids
#            pushing unreviewed changes to a shared repo at midnight. Once a remote exists,
#            pushing starts automatically.

set -euo pipefail

PROJECT_ROOT="/FAST/gr10634/gaozy/aml_niche_net"
LOG="${PROJECT_ROOT}/00_project/git_autocommit.log"

cd "${PROJECT_ROOT}"

ts() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(ts)] $*" >> "${LOG}"; }

# Bail out loudly if this is not a git repo, do not fail silently
if [[ ! -d .git ]]; then
  log "ERROR: ${PROJECT_ROOT} is not a git repo, skipping"
  exit 1
fi

# Do not create an empty commit when there is nothing to commit
if git diff --quiet && git diff --cached --quiet && [[ -z "$(git status --porcelain)" ]]; then
  log "no changes, skipping"
  exit 0
fi

git add -A

# Summarize the change into the commit body for easier history browsing
N_FILES=$(git diff --cached --name-only | wc -l)
STAT=$(git diff --cached --shortstat)
CHANGED=$(git diff --cached --name-only | head -20)

git commit -q -F - <<EOF
auto: $(date '+%Y-%m-%d') daily snapshot (${N_FILES} files)

${STAT}

Changed files (up to 20):
${CHANGED}

Generated automatically by scripts/99_admin/daily_commit.sh.
EOF

log "committed ${N_FILES} files -- ${STAT}"

# Push only if a remote 'origin' is configured
if git remote get-url origin >/dev/null 2>&1; then
  BRANCH=$(git rev-parse --abbrev-ref HEAD)
  if git push -q origin "${BRANCH}" 2>>"${LOG}"; then
    log "pushed to origin/${BRANCH}"
  else
    log "WARNING: push failed (commit is kept locally)"
  fi
else
  log "no origin configured, committed locally only"
fi

# ── Install cron (run once by hand) ─────────────────────────────────
# crontab -e, then add this line (daily 23:30):
#   30 23 * * * bash /FAST/gr10634/gaozy/aml_niche_net/scripts/99_admin/daily_commit.sh
#
# View auto-commit history:
#   git log --oneline --grep='^auto:'
#   tail -20 /FAST/gr10634/gaozy/aml_niche_net/00_project/git_autocommit.log
