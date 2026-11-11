#!/bin/bash
# 创建 16S 分析 conda 环境（最小依赖，避免完整 amplicon seed 的求解冲突）。
# 用法（在 WSL 中）：
#   bash /mnt/d/study/16S/pipeline/setup_env.sh

set -euo pipefail

ENV_NAME="16s"
PIPELINE_DIR=$(cd "$(dirname "$0")" && pwd)
ENV_YML="$PIPELINE_DIR/env/16s-minimal.yml"

if [[ ! -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
  echo "未找到 $HOME/miniconda3 ，请先在 WSL 中安装 Miniconda。"
  exit 1
fi

if [[ ! -f "$ENV_YML" ]]; then
  echo "缺少环境文件: $ENV_YML"
  exit 1
fi

# conda 激活脚本会读取未定义变量，这里暂时关闭 nounset
set +u
source "$HOME/miniconda3/etc/profile.d/conda.sh"
set -u

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "环境 ${ENV_NAME} 已存在，跳过创建。"
else
  echo "正在创建环境 ${ENV_NAME} ，体积较大，请等待求解与下载。"
  echo "使用本地文件: $ENV_YML"
  conda env create -n "$ENV_NAME" -f "$ENV_YML"
fi

echo "校验关键命令..."
set +u
conda activate "$ENV_NAME"
set -u
qiime --version
cutadapt --version
vsearch --version | head -n 1
trimmomatic -version
flash --version | head -n 1
biom --help >/dev/null

ROOT_DIR=$(cd "$PIPELINE_DIR/.." && pwd)
echo "环境 ${ENV_NAME} 已就绪。下一步："
echo "  conda activate ${ENV_NAME}"
echo "  # 原始数据放入 古菌项目/rawData 或 细菌项目/rawData 后："
echo "  bash ${ROOT_DIR}/古菌项目/run_16s.sh"
echo "  bash ${ROOT_DIR}/细菌项目/run_16s.sh"
echo "  # 若仍使用开发目录名 GJproject / XJproject，则改为对应路径。"
