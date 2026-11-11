#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""根据本次分析结果，调用智谱免费模型生成一篇简短中文论文草稿。

用法:
  python3 write_paper.py 项目目录

API Key 读取顺序:
  1) 环境变量 ZHIPU_API_KEY
  2) pipeline/zhipu.env 中的 ZHIPU_API_KEY=...
该文件不要打进交付包。
"""

import json
import os
import ssl
import sys
import urllib.error
import urllib.request
from collections import defaultdict

PIPELINE_DIR = os.path.dirname(os.path.abspath(__file__))
API_URL = "https://open.bigmodel.cn/api/paas/v4/chat/completions"
# 智谱免费文本模型；若平台更名可在 zhipu.env 里设 ZHIPU_MODEL=
DEFAULT_MODEL = "glm-4-flash"


def load_env_file(path):
    cfg = {}
    if not os.path.isfile(path):
        return cfg
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            cfg[key.strip()] = value.strip().strip('"').strip("'")
    return cfg


def get_api_key():
    key = os.environ.get("ZHIPU_API_KEY", "").strip()
    if key:
        return key
    env = load_env_file(os.path.join(PIPELINE_DIR, "zhipu.env"))
    key = env.get("ZHIPU_API_KEY", "").strip()
    if not key:
        raise SystemExit(
            "未找到智谱 API Key。请设置环境变量 ZHIPU_API_KEY，"
            "或在 pipeline/zhipu.env 中写入 ZHIPU_API_KEY=你的密钥"
        )
    return key


def get_model():
    env = load_env_file(os.path.join(PIPELINE_DIR, "zhipu.env"))
    return os.environ.get("ZHIPU_MODEL") or env.get("ZHIPU_MODEL") or DEFAULT_MODEL


def load_config(project):
    cfg = {}
    path = os.path.join(project, "config.sh")
    if not os.path.isfile(path):
        return cfg
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            # 跳过 bash 赋值里的命令替换
            if "$(" in line or "${" in line or "BASH_SOURCE" in line:
                continue
            key, value = line.split("=", 1)
            cfg[key.strip()] = value.strip().strip('"').strip("'")
    return cfg


def read_tsv(path):
    if not os.path.isfile(path):
        return [], []
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if not lines:
        return [], []
    header = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:]]
    return header, rows


def phylum_means(project):
    """返回门水平平均相对丰度（%），取前 8。"""
    path = os.path.join(project, "01.ASV.tax", "01.original", "phylum.xls")
    if not os.path.isfile(path):
        return []
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        n = len(header) - 1
        data = []
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            data.append((parts[0], [float(x) for x in parts[1:n + 1]]))
    if not data or n <= 0:
        return []
    col_sum = [0.0] * n
    for _name, vals in data:
        for i, v in enumerate(vals):
            col_sum[i] += v
    means = []
    for name, vals in data:
        fr = [(vals[i] / col_sum[i] * 100.0 if col_sum[i] else 0.0) for i in range(n)]
        means.append((name, sum(fr) / n))
    means.sort(key=lambda x: -x[1])
    return means[:8]


def build_facts(project):
    cfg = load_config(project)
    facts = {
        "project_name": cfg.get("PROJECT_NAME", os.path.basename(project)),
        "organism": cfg.get("ORGANISM", ""),
        "region": cfg.get("REGION", ""),
        "primer_f": cfg.get("PRIMER_F_NAME", ""),
        "primer_r": cfg.get("PRIMER_R_NAME", ""),
    }
    _, valid = read_tsv(os.path.join(project, "valid.stat.xls"))
    facts["n_samples"] = len(valid)
    if valid:
        seqs = [int(float(r[1])) for r in valid if len(r) > 1]
        facts["valid_seq_min"] = min(seqs) if seqs else None
        facts["valid_seq_max"] = max(seqs) if seqs else None
        facts["valid_seq_mean"] = int(sum(seqs) / len(seqs)) if seqs else None

    asv_path = os.path.join(project, "01.ASV.tax", "01.original", "asv_taxa_table.xls")
    if os.path.isfile(asv_path):
        with open(asv_path, "r", encoding="utf-8") as handle:
            n = sum(1 for _ in handle) - 1
        facts["n_asv"] = n

    alpha_h, alpha_rows = read_tsv(
        os.path.join(project, "supplement", "alpha_diversity.tsv")
    )
    if alpha_rows:
        by_group = defaultdict(list)
        for row in alpha_rows:
            # Sample Group Depth Observed Shannon ...
            if len(row) < 5:
                continue
            by_group[row[1]].append(
                {
                    "observed": float(row[3]),
                    "shannon": float(row[4]),
                }
            )
        facts["alpha_by_group"] = {
            g: {
                "observed_mean": round(sum(x["observed"] for x in rows) / len(rows), 1),
                "shannon_mean": round(sum(x["shannon"] for x in rows) / len(rows), 3),
                "n": len(rows),
            }
            for g, rows in by_group.items()
        }

    _, perm = read_tsv(os.path.join(project, "supplement", "permanova.tsv"))
    if perm:
        facts["permanova"] = {
            "pseudoF": perm[0][1] if len(perm[0]) > 1 else "",
            "R2": perm[0][2] if len(perm[0]) > 2 else "",
            "p": perm[0][4] if len(perm[0]) > 4 else "",
        }

    top = phylum_means(project)
    if top:
        facts["top_phyla"] = [{"name": n, "mean_pct": round(v, 2)} for n, v in top]

    asv1 = ""
    if os.path.isfile(asv_path):
        with open(asv_path, "r", encoding="utf-8") as handle:
            handle.readline()
            line = handle.readline().rstrip("\n")
            if line:
                parts = line.split("\t")
                asv1 = parts[-1][:120]
    if asv1:
        facts["asv1_taxonomy"] = asv1
    return facts


def chat(api_key, model, system, user):
    payload = {
        "model": model,
        "messages": [
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ],
        "temperature": 0.4,
    }
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        API_URL,
        data=data,
        headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {api_key}",
        },
        method="POST",
    )
    ctx = ssl.create_default_context()
    try:
        with urllib.request.urlopen(req, context=ctx, timeout=180) as resp:
            body = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise SystemExit(f"智谱 API 调用失败 HTTP {exc.code}: {detail}") from exc
    try:
        return body["choices"][0]["message"]["content"]
    except (KeyError, IndexError, TypeError) as exc:
        raise SystemExit(f"智谱返回无法解析: {body}") from exc


def main():
    if len(sys.argv) != 2:
        raise SystemExit("用法: python3 write_paper.py 项目目录")
    project = os.path.abspath(sys.argv[1])
    api_key = get_api_key()
    model = get_model()
    facts = build_facts(project)
    facts_json = json.dumps(facts, ensure_ascii=False, indent=2)

    system = (
        "你是环境微生物学与生物信息学领域的科研写作助手。"
        "请根据给定的分析结果事实，写一篇结构完整、表述严谨的简短中文论文草稿。"
        "不要编造未提供的数据；不确定处用谨慎表述。"
        "输出使用 Markdown，包含：标题、摘要、关键词、引言、材料与方法、结果、讨论、结论。"
        "全文约 2500～4000 字为宜。"
    )
    user = (
        "以下是一次 16S 扩增子测序分析的关键事实（JSON）。"
        "请据此撰写简短论文草稿，可适当补充该研究领域的一般背景，但具体数字必须以 JSON 为准。\n\n"
        f"{facts_json}"
    )
    print(f"调用智谱模型 {model} …", flush=True)
    paper = chat(api_key, model, system, user)

    out_dir = os.path.join(project, "paper")
    os.makedirs(out_dir, exist_ok=True)
    md_path = os.path.join(out_dir, "manuscript.md")
    facts_path = os.path.join(out_dir, "facts.json")
    with open(facts_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(facts_json + "\n")
    with open(md_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(paper.strip() + "\n")
        handle.write("\n---\n")
        handle.write(f"\n> 本文由智谱模型 `{model}` 根据本次分析结果自动生成，仅供草稿参考，发表前请人工审校。\n")
    print(md_path)


if __name__ == "__main__":
    main()
