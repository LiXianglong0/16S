#!/usr/bin/env bash
# 将 origin 改为 Gitee，原 GitHub 保留为远程名 github
# 用法：bash scripts/setup-gitee-origin.sh <Gitee用户名> [仓库名]
set -euo pipefail

GITEE_USER="${1:-}"
REPO_NAME="${2:-16S}"
if [[ -z "$GITEE_USER" ]]; then
  echo "用法: bash scripts/setup-gitee-origin.sh <Gitee用户名> [仓库名]"
  exit 1
fi

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"

GITEE_URL="https://gitee.com/${GITEE_USER}/${REPO_NAME}.git"
GITHUB_URL="https://github.com/LiXianglong0/16S.git"

if git remote get-url origin >/dev/null 2>&1; then
  cur=$(git remote get-url origin)
  if [[ "$cur" == *"github.com"* ]]; then
    if git remote get-url github >/dev/null 2>&1; then
      git remote set-url github "$GITHUB_URL"
    else
      git remote rename origin github
      git remote set-url github "$GITHUB_URL"
    fi
  fi
fi

if git remote get-url origin >/dev/null 2>&1; then
  git remote set-url origin "$GITEE_URL"
else
  git remote add origin "$GITEE_URL"
fi

if ! git remote get-url github >/dev/null 2>&1; then
  git remote add github "$GITHUB_URL"
fi

echo "当前远程："
git remote -v
echo
echo "下一步："
echo "  1) 确认 Gitee 上已建空仓库: $GITEE_URL"
echo "  2) git push -u origin main"
echo "  3) 按 Gitee同步说明.md 配置 GitHub Secrets: GITEE_REPO / GITEE_TOKEN"
