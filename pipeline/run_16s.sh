#!/bin/bash
# 16S 扩增子分析主流程。
# 与模板交付的 valid 文件同一口径：质控不丢弃读段，切除引物，flash 合并。
# 有效序列统计在去除嵌合体之前。嵌合体由后面的 DADA2 denoise-single 去除。
# 用法：bash run_16s.sh /path/to/GJproject

set -euo pipefail

PROJECT_DIR="${1:-}"
if [[ -z "$PROJECT_DIR" ]]; then
  echo "用法: bash run_16s.sh 项目目录"
  exit 1
fi
PROJECT_DIR=$(cd "$PROJECT_DIR" && pwd)
PIPELINE_DIR=$(cd "$(dirname "$0")" && pwd)

if [[ ! -f "$PROJECT_DIR/config.sh" ]]; then
  echo "缺少配置文件: $PROJECT_DIR/config.sh"
  exit 1
fi
# shellcheck disable=SC1091
source "$PROJECT_DIR/config.sh"
PRIMER_F_EXTRA="${PRIMER_F_EXTRA:-}"
REF_SEQS="${REF_SEQS:-}"
REF_TAX="${REF_TAX:-}"
RUN_PICRUST2="${RUN_PICRUST2:-0}"
RUN_FUNGAL_FUNCTION="${RUN_FUNGAL_FUNCTION:-0}"
MARKER="${MARKER:-16S}"
THREADS="${THREADS:-8}"
MIN_OVERLAP="${MIN_OVERLAP:-8}"
MAX_OVERLAP="${MAX_OVERLAP:-300}"
FLASH_X="${FLASH_X:-0.6}"
CUTADAPT_E="${CUTADAPT_E:-0.45}"
CUTADAPT_O="${CUTADAPT_O:-8}"

RAW_DIR="$PROJECT_DIR/rawData"
VALID_DIR="$PROJECT_DIR/valid"
WORK_DIR="$PROJECT_DIR/work"
QIIME_DIR="$PROJECT_DIR/qiime2"

log() {
  echo "[$(date '+%F %T')] $*"
}

need_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "缺少命令: $1 。请先执行 bash $PIPELINE_DIR/setup_env.sh 并 conda activate 16s"
    exit 1
  fi
}

# 激活 conda 环境。已在目标环境中时不会重复创建。
activate_env() {
  if [[ "${CONDA_DEFAULT_ENV:-}" == "16s" ]]; then
    return
  fi
  if [[ ! -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    echo "未找到 miniconda3，无法激活环境 16s"
    exit 1
  fi
  set +u
  # shellcheck disable=SC1091
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
  conda activate 16s
  set -u
}

# 从 raw.{样本}.1.fq.gz 收集样本名，并确认双端文件都在。
collect_samples() {
  SAMPLES=()
  local name sample r2
  shopt -s nullglob
  for name in "$RAW_DIR"/raw.*.1.fq.gz; do
    name=$(basename "$name")
    sample=${name#raw.}
    sample=${sample%.1.fq.gz}
    r2="$RAW_DIR/raw.${sample}.2.fq.gz"
    if [[ ! -f "$r2" ]]; then
      echo "缺少配对文件: $r2"
      exit 1
    fi
    SAMPLES+=("$sample")
  done
  shopt -u nullglob
  if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "在 $RAW_DIR 中没有找到 raw.{样本}.1.fq.gz"
    exit 1
  fi
  local sorted=()
  mapfile -t sorted < <(printf '%s\n' "${SAMPLES[@]}" | sort -V)
  SAMPLES=("${sorted[@]}")
}

# 单个样本：质控、切引物、合并、去嵌合体，写出 valid/trim.{样本}.fq.gz
process_sample() {
  local sample="$1"
  local r1="$RAW_DIR/raw.${sample}.1.fq.gz"
  local r2="$RAW_DIR/raw.${sample}.2.fq.gz"
  local sd="$WORK_DIR/samples/$sample"
  mkdir -p "$sd" "$VALID_DIR"

  # 原始序列条数按 R1 计数，写入 rawdata.stat.xls
  local raw_n
  raw_n=$(( $(gzip -dc "$r1" | wc -l) / 4 ))
  printf '%s\t%s\n' "$sample" "$raw_n" >> "$WORK_DIR/raw_counts.tsv"
  log "样本 ${sample} 原始序列 ${raw_n}"

  # 第二步：只去掉两端低于 Q3 的碱基。模板 valid 条数等于原始读段对数，不能用滑窗把读段截短后丢弃
  log "样本 ${sample} 质控"
  trimmomatic PE -threads "$THREADS" -phred33 \
    "$r1" "$r2" \
    "$sd/trim.1.fq.gz" "$sd/trim.1.unpaired.fq.gz" \
    "$sd/trim.2.fq.gz" "$sd/trim.2.unpaired.fq.gz" \
    LEADING:3 TRAILING:3 MINLEN:50

  # 第三步：cutadapt 从 5' 端切除引物。错配放宽后，引物区有插入的读段仍能保留
  log "样本 ${sample} 匹配引物"
  if [[ -n "$PRIMER_F_EXTRA" ]]; then
    cutadapt -j "$THREADS" -e "$CUTADAPT_E" -O "$CUTADAPT_O" --discard-untrimmed --pair-filter=any \
      -g "^${PRIMER_F}" -g "^${PRIMER_F_EXTRA}" -G "^${PRIMER_R}" \
      -o "$sd/primer.1.fq" -p "$sd/primer.2.fq" \
      "$sd/trim.1.fq.gz" "$sd/trim.2.fq.gz"
  else
    cutadapt -j "$THREADS" -e "$CUTADAPT_E" -O "$CUTADAPT_O" --discard-untrimmed --pair-filter=any \
      -g "^${PRIMER_F}" -G "^${PRIMER_R}" \
      -o "$sd/primer.1.fq" -p "$sd/primer.2.fq" \
      "$sd/trim.1.fq.gz" "$sd/trim.2.fq.gz"
  fi

  # 第四步：flash 合并。扩增子约 380–410 bp，读长约 295 bp，重叠超过 flash 默认的 65 bp
  log "样本 ${sample} 合并双端"
  mkdir -p "$sd/flash"
  flash -t "$THREADS" -m "$MIN_OVERLAP" -M "$MAX_OVERLAP" -x "$FLASH_X" \
    -d "$sd/flash" -o merge \
    "$sd/primer.1.fq" "$sd/primer.2.fq"
  local merged="$sd/flash/merge.extendedFrags.fastq"
  if [[ ! -s "$merged" ]]; then
    echo "样本 ${sample} 没有合并成功的序列"
    exit 1
  fi

  # 模板 valid 文件是合并后的全部序列。嵌合体留给 DADA2，不在这里扣除
  gzip -c "$merged" > "$VALID_DIR/trim.${sample}.fq.gz"
}

# 写出 QIIME 2 单端 manifest 和分组表
write_manifest() {
  local manifest="$QIIME_DIR/manifest.tsv"
  local metadata="$QIIME_DIR/sample-metadata.tsv"
  printf 'sample-id\tabsolute-filepath\n' > "$manifest"
  printf 'sample-id\tgroup\n' > "$metadata"
  local sample group fq
  for sample in "${SAMPLES[@]}"; do
    fq=$(readlink -f "$VALID_DIR/trim.${sample}.fq.gz")
    printf '%s\t%s\n' "$sample" "$fq" >> "$manifest"
    group="${sample%-*}"
    printf '%s\t%s\n' "$sample" "$group" >> "$metadata"
  done
}

# QIIME 2 导入、DADA2 去噪、导出 ASV 表
run_qiime() {
  log "导入 valid 序列"
  qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path "$QIIME_DIR/manifest.tsv" \
    --output-path "$QIIME_DIR/data.qza" \
    --input-format SingleEndFastqManifestPhred33V2

  # trunc-len 0 表示不按固定长度截断，与报告命令一致
  log "DADA2 去噪"
  # 2025.10 比模板里的旧版多一个必填输出 base-transition-stats，去噪参数仍是 trunc-len 0
  qiime dada2 denoise-single \
    --i-demultiplexed-seqs "$QIIME_DIR/data.qza" \
    --p-trim-left 0 \
    --p-trunc-len 0 \
    --p-n-threads "$THREADS" \
    --o-table "$QIIME_DIR/table.qza" \
    --o-representative-sequences "$QIIME_DIR/rep-seqs.qza" \
    --o-denoising-stats "$QIIME_DIR/denoising-stats.qza" \
    --o-base-transition-stats "$QIIME_DIR/base-transition-stats.qza"

  log "导出特征表和代表序列"
  rm -rf "$QIIME_DIR/export"
  mkdir -p "$QIIME_DIR/export" "$QIIME_DIR/export-seqs" "$QIIME_DIR/export-stats"
  qiime tools export --input-path "$QIIME_DIR/table.qza" --output-path "$QIIME_DIR/export"
  qiime tools export --input-path "$QIIME_DIR/rep-seqs.qza" --output-path "$QIIME_DIR/export-seqs"
  qiime tools export --input-path "$QIIME_DIR/denoising-stats.qza" --output-path "$QIIME_DIR/export-stats"
  biom convert -i "$QIIME_DIR/export/feature-table.biom" -o "$QIIME_DIR/export/otu_table.tsv" --to-tsv
  if [[ -f "$QIIME_DIR/export-seqs/dna-sequences.fasta" ]]; then
    cp "$QIIME_DIR/export-seqs/dna-sequences.fasta" "$QIIME_DIR/export/dna-sequences.fasta"
  fi
}

# 物种注释。QIIME 2 的标准比对方法是 classify-consensus-vsearch：
# 相似度 0.8、覆盖度 0.8，共识比例 0.8，对应原报告里的置信度 0.8。
# 参考序列已是 qza 时直接使用；fasta 则先导入。
run_taxonomy() {
  local db_dir
  db_dir="$(cd "$PIPELINE_DIR/../database/silva138" && pwd)"
  if [[ -z "$REF_SEQS" && -f "$db_dir/silva-138-99-seqs.qza" ]]; then
    REF_SEQS="$db_dir/silva-138-99-seqs.qza"
    REF_TAX="$db_dir/silva-138-99-tax.qza"
  fi
  if [[ -z "$REF_SEQS" || -z "$REF_TAX" ]]; then
    log "未配置 REF_SEQS / REF_TAX，跳过分类学注释"
    return
  fi
  if [[ ! -f "$REF_SEQS" || ! -f "$REF_TAX" ]]; then
    echo "参考库文件不存在: $REF_SEQS 或 $REF_TAX"
    exit 1
  fi
  if [[ ! -f "$QIIME_DIR/rep-seqs.qza" || ! -f "$QIIME_DIR/export/otu_table.tsv" ]]; then
    echo "还没有 ASV 代表序列或特征表，无法做物种注释"
    exit 1
  fi
  local marker="${MARKER:-16S}"
  log "分类学注释（${marker}，VSEARCH 共识；参考库见 REF_SEQS）"
  local ref_seqs="$REF_SEQS" ref_tax="$REF_TAX"
  if [[ "$REF_SEQS" != *.qza ]]; then
    qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path "$REF_SEQS" \
      --output-path "$QIIME_DIR/ref-seqs.qza"
    ref_seqs="$QIIME_DIR/ref-seqs.qza"
  fi
  if [[ "$REF_TAX" != *.qza ]]; then
    local fmt="HeaderlessTSVTaxonomyFormat"
    if head -n 1 "$REF_TAX" | grep -qi 'taxon'; then
      fmt="TSVTaxonomyFormat"
    fi
    qiime tools import \
      --type 'FeatureData[Taxonomy]' \
      --input-format "$fmt" \
      --input-path "$REF_TAX" \
      --output-path "$QIIME_DIR/ref-tax.qza"
    ref_tax="$QIIME_DIR/ref-tax.qza"
  fi
  qiime feature-classifier classify-consensus-vsearch \
    --i-query "$QIIME_DIR/rep-seqs.qza" \
    --i-reference-reads "$ref_seqs" \
    --i-reference-taxonomy "$ref_tax" \
    --p-perc-identity 0.8 \
    --p-query-cov 0.8 \
    --p-min-consensus 0.8 \
    --p-strand both \
    --p-threads "$THREADS" \
    --o-classification "$QIIME_DIR/taxonomy.qza" \
    --o-search-results "$QIIME_DIR/taxonomy-search.qza"
  qiime tools export --input-path "$QIIME_DIR/taxonomy.qza" --output-path "$QIIME_DIR/export"
  log "写出原始丰度与抽平归一化丰度"
  python3 "$PIPELINE_DIR/taxa_abundance.py" "$PROJECT_DIR"
  log "构建 ASV 系统发育树"
  python3 "$PIPELINE_DIR/build_phylogeny.py" "$PROJECT_DIR"
  if [[ "${RUN_FUNGAL_FUNCTION:-0}" == "1" ]]; then
    log "真菌生态功能注释（FUNGuild + FungalTraits）"
    python3 "$PIPELINE_DIR/fungal_function_annotate.py" "$PROJECT_DIR" \
      || log "真菌功能注释跳过或失败（请检查 database/funguild 与 fungaltraits）"
  fi
}

# 可选：PICRUSt2 功能预测。参考库较大，默认不运行
run_picrust2() {
  if [[ "$RUN_PICRUST2" != "1" ]]; then
    log "RUN_PICRUST2 未开启，跳过功能预测"
    return
  fi
  if ! command -v picrust2_pipeline.py >/dev/null 2>&1; then
    echo "未找到 picrust2_pipeline.py，请先安装 PICRUSt2 后再把 RUN_PICRUST2 设为 1"
    exit 1
  fi
  log "PICRUSt2 功能预测"
  picrust2_pipeline.py \
    -s "$QIIME_DIR/export/dna-sequences.fasta" \
    -i "$QIIME_DIR/export/feature-table.biom" \
    -o "$QIIME_DIR/picrust2" \
    -p "$THREADS"
}

# 记录实际软件版本。单个命令失败时不影响主流程
write_versions() {
  one_version() {
    local text
    text=$("$@" 2>&1 | head -n 1 || true)
    echo "${text:-unknown}"
  }
  {
    echo -e "cutadapt\t$(one_version cutadapt --version)"
    echo -e "vsearch\t$(one_version vsearch --version)"
    echo -e "flash\t$(one_version flash --version)"
    echo -e "trimmomatic\t$(one_version trimmomatic -version)"
    echo -e "qiime2\t$(one_version qiime --version)"
  } > "$WORK_DIR/software_versions.txt" || true
}

main() {
  cd "$PROJECT_DIR"
  activate_env
  need_cmd trimmomatic
  need_cmd cutadapt
  need_cmd flash
  need_cmd vsearch
  need_cmd qiime
  need_cmd biom
  need_cmd python3

  mkdir -p "$WORK_DIR" "$VALID_DIR" "$QIIME_DIR"
  : > "$WORK_DIR/raw_counts.tsv"

  # 第一步：报告中的拆 barcode 使用 fastq-multx。当前原始数据已按样本拆分，直接进入质控
  log "原始数据已按样本拆分，跳过 fastq-multx"
  collect_samples
  log "共 ${#SAMPLES[@]} 个样本：${SAMPLES[*]}"

  local sample
  for sample in "${SAMPLES[@]}"; do
    process_sample "$sample"
  done

  log "写出序列统计表"
  python3 "$PIPELINE_DIR/fastq_stat.py" \
    --raw-counts "$WORK_DIR/raw_counts.tsv" \
    --valid-dir "$VALID_DIR" \
    --raw-out "$PROJECT_DIR/rawdata.stat.xls" \
    --valid-out "$PROJECT_DIR/valid.stat.xls"

  write_manifest
  run_qiime
  run_taxonomy
  run_picrust2
  write_versions

  log "生成报告"
  python3 "$PIPELINE_DIR/make_report.py" "$PROJECT_DIR"
  log "合并多样性分析到主报告"
  python3 "$PIPELINE_DIR/supplement_diversity.py" "$PROJECT_DIR"
  if [[ "${WRITE_PAPER:-1}" == "1" ]]; then
    log "调用智谱模型生成论文草稿"
    python3 "$PIPELINE_DIR/write_paper.py" "$PROJECT_DIR" || log "论文生成跳过或失败（可稍后单独运行 write_paper.py）"
  fi
  log "完成。报告: $PROJECT_DIR/amplicon-ngs.report.html"
  log "论文草稿: $PROJECT_DIR/paper/manuscript.md（若已生成）"
}

# 第二个参数为 taxonomy 时，只做物种注释和丰度表，不重跑质控与 DADA2。
if [[ "${2:-}" == "taxonomy" ]]; then
  cd "$PROJECT_DIR"
  activate_env
  run_taxonomy
else
  main
fi
