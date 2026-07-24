#!/usr/bin/env bash
# daily_commit.sh ----
# 每日自动快照：把当天对代码/小表/图/文档的改动提交到本地 git。
# INPUT  : 项目工作区当前状态
# OUTPUT : 一个 commit（无改动则跳过）；日志追加到 00_project/git_autocommit.log
# Usage  : bash scripts/99_admin/daily_commit.sh          # 手动跑
#          crontab 每日 23:30 自动跑（见文件末尾说明）
# [DECISION] 只提交不推送，除非配置了名为 origin 的远程 —— 避免半夜把
#            未审阅的改动推到共享仓库。远程一旦配置，自动开始推送。

set -euo pipefail

PROJECT_ROOT="/FAST/gr10634/gaozy/aml_niche_net"
LOG="${PROJECT_ROOT}/00_project/git_autocommit.log"

cd "${PROJECT_ROOT}"

ts() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(ts)] $*" >> "${LOG}"; }

# 不是 git 仓库就直接退出，不要静默失败
if [[ ! -d .git ]]; then
  log "ERROR: ${PROJECT_ROOT} 不是 git 仓库，跳过"
  exit 1
fi

# 无改动则不产生空 commit
if git diff --quiet && git diff --cached --quiet && [[ -z "$(git status --porcelain)" ]]; then
  log "无改动，跳过"
  exit 0
fi

git add -A

# 统计改动，写进 commit message 正文，便于日后翻历史
N_FILES=$(git diff --cached --name-only | wc -l)
STAT=$(git diff --cached --shortstat)
CHANGED=$(git diff --cached --name-only | head -20)

git commit -q -F - <<EOF
auto: $(date '+%Y-%m-%d') 每日快照 (${N_FILES} 个文件)

${STAT}

改动文件（最多列 20 个）：
${CHANGED}

由 scripts/99_admin/daily_commit.sh 自动生成。
EOF

log "已提交 ${N_FILES} 个文件 — ${STAT}"

# 只有配置了 origin 才推送
if git remote get-url origin >/dev/null 2>&1; then
  BRANCH=$(git rev-parse --abbrev-ref HEAD)
  if git push -q origin "${BRANCH}" 2>>"${LOG}"; then
    log "已推送到 origin/${BRANCH}"
  else
    log "WARNING: 推送失败（commit 已在本地保留）"
  fi
else
  log "未配置 origin，仅本地提交"
fi

# ── 安装 cron（手动执行一次）────────────────────────────────
# crontab -e 后加入这一行（每日 23:30）：
#   30 23 * * * bash /FAST/gr10634/gaozy/aml_niche_net/scripts/99_admin/daily_commit.sh
#
# 查看自动提交历史：
#   git log --oneline --grep='^auto:'
#   tail -20 /FAST/gr10634/gaozy/aml_niche_net/00_project/git_autocommit.log
