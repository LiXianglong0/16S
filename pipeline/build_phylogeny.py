#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""根据 ASV 代表序列生成比对与系统发育树，写入 01.ASV.tax/03.phylogeny_tree。

优先使用 taxonomy.aligned.tsv（模板对齐后），否则用 taxonomy.tsv。
ASV 编号与过滤规则与 taxa_abundance.py 一致。
"""

import os
import subprocess
import sys


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
            key, value = line.split("=", 1)
            cfg[key.strip()] = value.strip().strip('"').strip("'")
    return cfg


def load_fasta(path):
    seqs = {}
    name = None
    parts = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def read_taxonomy(path):
    taxa = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if parts[0] in ("Feature ID", "FeatureID"):
                continue
            taxon = parts[1].strip() if len(parts) > 1 else "Unassigned"
            taxa[parts[0]] = taxon.replace("; ", ";")
    return taxa


def keep_domain(taxon, domain):
    if not domain:
        return True
    compact = (taxon or "").replace("; ", ";")
    for prefix in (f"d__{domain}", f"k__{domain}"):
        if compact.startswith(prefix + ";") or compact == prefix:
            return True
    return False


def feature_order(otu_path, tax_path, domain):
    taxa = read_taxonomy(tax_path) if os.path.isfile(tax_path) else {}
    items = []
    with open(otu_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("# "):
                continue
            parts = line.split("\t")
            if parts[0] in ("#OTU ID", "OTU ID"):
                continue
            feature = parts[0]
            taxon = taxa.get(feature, "Unassigned")
            if not keep_domain(taxon, domain):
                continue
            total = sum(int(round(float(x))) for x in parts[1:])
            items.append((total, feature))
    items.sort(key=lambda item: (-item[0], item[1]))
    return [feature for _total, feature in items]


def write_fasta(path, pairs):
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for name, seq in pairs:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def build(project):
    project = os.path.abspath(project)
    cfg = load_config(project)
    domain = cfg.get("KEEP_DOMAIN", "").strip()
    fasta = os.path.join(project, "qiime2", "export", "dna-sequences.fasta")
    otu = os.path.join(project, "qiime2", "export", "otu_table.tsv")
    tax_aligned = os.path.join(project, "qiime2", "export", "taxonomy.aligned.tsv")
    tax = tax_aligned if os.path.isfile(tax_aligned) else os.path.join(
        project, "qiime2", "export", "taxonomy.tsv"
    )
    out_dir = os.path.join(project, "01.ASV.tax", "03.phylogeny_tree")
    if not os.path.isfile(fasta) or not os.path.isfile(otu):
        raise SystemExit("缺少代表序列或 ASV 表，无法建树")
    os.makedirs(out_dir, exist_ok=True)

    seqs = load_fasta(fasta)
    order = feature_order(otu, tax, domain)
    pairs = []
    for index, feature in enumerate(order, start=1):
        if feature in seqs:
            pairs.append((f"ASV{index}", seqs[feature]))
    if len(pairs) < 3:
        raise SystemExit("可用于建树的序列少于 3 条")

    unaligned = os.path.join(out_dir, "asv_rep.fasta")
    aligned = os.path.join(out_dir, "asv_rep_aligned.fasta")
    tree = os.path.join(out_dir, "phylogeny.tre")
    write_fasta(unaligned, pairs)

    print(f"比对 {len(pairs)} 条 ASV …", flush=True)
    with open(aligned, "w", encoding="utf-8", newline="\n") as handle:
        subprocess.run(
            ["mafft", "--quiet", "--thread", "-1", "--retree", "2", "--maxiterate", "2", unaligned],
            check=True,
            stdout=handle,
        )
    print("建树 …", flush=True)
    with open(tree, "w", encoding="utf-8", newline="\n") as handle:
        subprocess.run(
            ["FastTree", "-nt", "-gtr", aligned],
            check=True,
            stdout=handle,
            stderr=subprocess.DEVNULL,
        )
    with open(os.path.join(out_dir, "3.系统发育树说明.txt"), "w", encoding="utf-8", newline="\n") as handle:
        handle.write(
            "本文件夹：03.phylogeny_tree（系统发育树）\n\n"
            "文件说明\n"
            "- asv_rep.fasta：按总丰度排序后的 ASV 代表序列（ASV1、ASV2…）\n"
            "- asv_rep_aligned.fasta：MAFFT 多序列比对结果\n"
            "- phylogeny.tre：FastTree 构建的近似最大似然树（Newick）\n\n"
            "用途：Faith's PD、UniFrac 等依赖系统发育关系的多样性分析。\n"
            f"本项目共写入 {len(pairs)} 条 ASV"
            + (f"（KEEP_DOMAIN={domain}）" if domain else "")
            + "。\n"
        )
    print(out_dir)


def main():
    if len(sys.argv) != 2:
        raise SystemExit("用法: python3 build_phylogeny.py 项目目录")
    build(sys.argv[1])


if __name__ == "__main__":
    main()
