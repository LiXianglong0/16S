#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""真菌 ITS 生态功能注释：FUNGuild + FungalTraits（对照微信教程双库思路）。

输入：项目 01.ASV.tax/01.original/asv_taxa_table.xls
输出：01.ASV.tax/04.fungal_function/

数据库默认路径（可用环境变量或参数覆盖）：
  database/funguild/FUNGuild.txt
  database/fungaltraits/FungalTraits_genera.tsv

用法:
  python3 fungal_function_annotate.py 项目目录
"""

from __future__ import annotations

import csv
import os
import sys
from collections import OrderedDict, defaultdict

PIPELINE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(PIPELINE_DIR, ".."))


def db_paths():
    return {
        "funguild": os.environ.get(
            "FUNGUILD_DB",
            os.path.join(ROOT, "database", "funguild", "FUNGuild.txt"),
        ),
        "fungaltraits": os.environ.get(
            "FUNGALTRAITS_DB",
            next(
                (
                    p
                    for p in (
                        os.path.join(ROOT, "database", "fungaltraits", "FungalTraits_genera.tsv"),
                        os.path.join(ROOT, "database", "fungaltraits", "FungalTraits_genera.csv"),
                    )
                    if os.path.isfile(p)
                ),
                os.path.join(ROOT, "database", "fungaltraits", "FungalTraits_genera.tsv"),
            ),
        ),
    }


def read_asv_table(path):
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        samples = header[1:-1]
        rows = []
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            asv = parts[0]
            counts = [int(round(float(x))) for x in parts[1:-1]]
            tax = parts[-1].replace("; ", ";")
            rows.append({"id": asv, "counts": counts, "tax": tax, "total": sum(counts)})
    return samples, rows


def parse_taxonomy(tax):
    """解析 k__/p__/... 或 d__/... 为层级字典。"""
    levels = OrderedDict()
    if not tax or tax == "Unassigned":
        return levels
    for part in tax.split(";"):
        part = part.strip()
        if "__" in part:
            prefix, name = part.split("__", 1)
            key = {
                "k": "kingdom",
                "d": "domain",
                "p": "phylum",
                "c": "class",
                "o": "order",
                "f": "family",
                "g": "genus",
                "s": "species",
            }.get(prefix.lower())
            if key and name and name.lower() not in ("", "unidentified", "uncultured", "norank"):
                levels[key] = name
        elif part:
            # 无前缀时不猜测层级
            pass
    return levels


def load_funguild(path):
    """加载 FUNGuild 文本库（制表符或逗号分隔均可）。"""
    if not os.path.isfile(path):
        return None
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delim = "\t" if sample.count("\t") >= sample.count(",") else ","
        reader = csv.DictReader(handle, delimiter=delim)
        for row in reader:
            # 统一小写键
            norm = { (k or "").strip().lower(): (v or "").strip() for k, v in row.items() }
            taxon = norm.get("taxon") or norm.get("taxons") or ""
            if not taxon:
                continue
            rows.append(
                {
                    "taxon": taxon,
                    "trophicMode": norm.get("trophicmode") or norm.get("trophic_mode") or "",
                    "guild": norm.get("guild") or "",
                    "confidenceRanking": norm.get("confidenceranking")
                    or norm.get("confidence")
                    or "",
                    "taxonomicLevel": norm.get("taxonomiclevel") or norm.get("taxonomic_level") or "",
                    "growthForm": norm.get("growthform") or "",
                    "trait": norm.get("trait") or "",
                }
            )
    # 按分类名小写建索引
    by_name = defaultdict(list)
    for row in rows:
        by_name[row["taxon"].lower()].append(row)
    return {"rows": rows, "by_name": by_name}


def load_fungaltraits(path):
    if not os.path.isfile(path):
        return None
    by_genus = defaultdict(list)
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delim = "\t" if sample.count("\t") >= sample.count(",") else ","
        reader = csv.DictReader(handle, delimiter=delim)
        for row in reader:
            norm = {}
            for k, v in row.items():
                if k is None or not isinstance(v, str):
                    continue
                norm[k.strip().lower()] = v.strip()
            genus = norm.get("genus") or norm.get("genera") or ""
            if not genus:
                continue
            by_genus[genus.lower()].append(
                {
                    "GENUS": genus,
                    "primary_lifestyle": norm.get("primary_lifestyle")
                    or norm.get("primarylifestyle")
                    or "",
                    "Secondary_lifestyle": norm.get("secondary_lifestyle")
                    or norm.get("secondarylifestyle")
                    or "",
                }
            )
    return by_genus


def match_funguild(levels, fg_db):
    """自种→属→科…逐级尝试匹配 FUNGuild。"""
    if not fg_db:
        return None
    order = ["species", "genus", "family", "order", "class", "phylum", "kingdom", "domain"]
    for key in order:
        name = levels.get(key)
        if not name:
            continue
        hits = fg_db["by_name"].get(name.lower())
        if hits:
            hit = hits[0]
            hit = dict(hit)
            hit["matched_on"] = f"{key}:{name}"
            return hit
    return None


def match_fungaltraits(levels, ft_db):
    if not ft_db:
        return None
    genus = levels.get("genus")
    if not genus:
        return None
    hits = ft_db.get(genus.lower())
    if not hits:
        return None
    lifestyles = sorted({h["primary_lifestyle"] for h in hits if h["primary_lifestyle"]})
    if len(lifestyles) > 1:
        return {
            "GENUS": genus,
            "primary_lifestyle": "|".join(lifestyles),
            "Secondary_lifestyle": hits[0].get("Secondary_lifestyle", ""),
            "note": "multiple_primary_lifestyle",
        }
    return hits[0]


def write_tsv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(x) for x in row) + "\n")


def composition(samples, annotated, field, n):
    """按功能字段汇总绝对计数与相对丰度（分母=全部 reads）。"""
    abs_map = OrderedDict()
    depth = [0] * n
    for row in annotated:
        for i, v in enumerate(row["counts"]):
            depth[i] += v
        label = row.get(field) or "Unassigned"
        if label not in abs_map:
            abs_map[label] = [0] * n
        for i, v in enumerate(row["counts"]):
            abs_map[label][i] += v
    abs_rows = []
    rel_rows = []
    for label, vals in sorted(abs_map.items(), key=lambda x: -sum(x[1])):
        abs_rows.append([label] + vals)
        rel_rows.append(
            [label]
            + [
                f"{(vals[i] / depth[i] if depth[i] else 0.0):.6f}"
                for i in range(n)
            ]
        )
    return abs_rows, rel_rows, depth


def missing_db_message(paths):
    lines = [
        "缺少真菌功能注释数据库，请先下载到下列路径后再运行本脚本：",
        "",
        f"1) FUNGuild → {paths['funguild']}",
        "   推荐：https://github.com/UMNFuN/FUNGuild",
        "   或将 FUNGuild.txt / funguild_db 导出为制表符文本，需含列：",
        "   taxon, trophicMode, guild, confidenceRanking",
        "",
        f"2) FungalTraits → {paths['fungaltraits']}",
        "   来源：Põlme et al. FungalTraits 论文补充材料",
        "   整理为 TSV/CSV，至少含列：GENUS, primary_lifestyle",
        "",
        "详细说明见 ZJproject/下载清单.txt 与 README.md。",
    ]
    return "\n".join(lines)


def main():
    if len(sys.argv) != 2:
        raise SystemExit("用法: python3 fungal_function_annotate.py 项目目录")
    project = os.path.abspath(sys.argv[1])
    asv_path = os.path.join(project, "01.ASV.tax", "01.original", "asv_taxa_table.xls")
    if not os.path.isfile(asv_path):
        raise SystemExit(f"缺少 ASV 表: {asv_path}（请先完成 ITS 注释与丰度）")

    paths = db_paths()
    fg_db = load_funguild(paths["funguild"])
    ft_db = load_fungaltraits(paths["fungaltraits"])
    if fg_db is None and ft_db is None:
        raise SystemExit(missing_db_message(paths))

    samples, rows = read_asv_table(asv_path)
    n = len(samples)
    annotated = []
    fg_hit = ft_hit = 0
    fg_reads = ft_reads = 0
    total_reads = 0

    for row in rows:
        levels = parse_taxonomy(row["tax"])
        fg = match_funguild(levels, fg_db) if fg_db else None
        ft = match_fungaltraits(levels, ft_db) if ft_db else None
        if fg:
            fg_hit += 1
            fg_reads += row["total"]
        if ft:
            ft_hit += 1
            ft_reads += row["total"]
        total_reads += row["total"]
        annotated.append(
            {
                **row,
                "genus": levels.get("genus", ""),
                "FG_trophicMode": (fg or {}).get("trophicMode", ""),
                "FG_guild": (fg or {}).get("guild", ""),
                "FG_confidence": (fg or {}).get("confidenceRanking", ""),
                "FG_matched_on": (fg or {}).get("matched_on", ""),
                "FT_primary_lifestyle": (ft or {}).get("primary_lifestyle", ""),
                "FT_secondary_lifestyle": (ft or {}).get("Secondary_lifestyle", ""),
                "FT_note": (ft or {}).get("note", ""),
            }
        )

    out = os.path.join(project, "01.ASV.tax", "04.fungal_function")
    os.makedirs(out, exist_ok=True)

    ann_header = [
        "ASV ID",
        *samples,
        "taxonomy",
        "genus",
        "FG_trophicMode",
        "FG_guild",
        "FG_confidence",
        "FG_matched_on",
        "FT_primary_lifestyle",
        "FT_secondary_lifestyle",
        "FT_note",
    ]
    ann_rows = []
    for r in annotated:
        ann_rows.append(
            [
                r["id"],
                *r["counts"],
                r["tax"],
                r["genus"],
                r["FG_trophicMode"],
                r["FG_guild"],
                r["FG_confidence"],
                r["FG_matched_on"],
                r["FT_primary_lifestyle"],
                r["FT_secondary_lifestyle"],
                r["FT_note"],
            ]
        )
    write_tsv(os.path.join(out, "ASV_annotated.xls"), ann_header, ann_rows)

    cov_rows = [
        ["metric", "value"],
        ["n_asv", len(annotated)],
        ["total_reads", total_reads],
        ["FUNGuild_asv_hit", fg_hit],
        ["FUNGuild_asv_rate", f"{(fg_hit / len(annotated) if annotated else 0):.4f}"],
        ["FUNGuild_read_hit", fg_reads],
        ["FUNGuild_read_rate", f"{(fg_reads / total_reads if total_reads else 0):.4f}"],
        ["FungalTraits_asv_hit", ft_hit],
        ["FungalTraits_asv_rate", f"{(ft_hit / len(annotated) if annotated else 0):.4f}"],
        ["FungalTraits_read_hit", ft_reads],
        ["FungalTraits_read_rate", f"{(ft_reads / total_reads if total_reads else 0):.4f}"],
        ["FUNGuild_db", paths["funguild"] if fg_db else "MISSING"],
        ["FungalTraits_db", paths["fungaltraits"] if ft_db else "MISSING"],
    ]
    write_tsv(os.path.join(out, "coverage.xls"), ["metric", "value"], cov_rows[1:])

    for field, abs_name, rel_name in (
        ("FG_trophicMode", "FG_trophicMode.absolute.xls", "FG_trophicMode.percent.xls"),
        ("FG_guild", "FG_guild.absolute.xls", "FG_guild.percent.xls"),
        ("FT_primary_lifestyle", "FT_primary_lifestyle.absolute.xls", "FT_primary_lifestyle.percent.xls"),
    ):
        abs_rows, rel_rows, _depth = composition(samples, annotated, field, n)
        write_tsv(os.path.join(out, abs_name), ["Function"] + samples, abs_rows)
        write_tsv(os.path.join(out, rel_name), ["Function"] + samples, rel_rows)

    with open(os.path.join(out, "4.真菌功能注释说明.txt"), "w", encoding="utf-8", newline="\n") as handle:
        handle.write(
            "本文件夹：04.fungal_function（真菌生态功能注释）\n\n"
            "方法对照：FUNGuild + FungalTraits 双库独立注释（不强制合并）。\n"
            "输入：01.original/asv_taxa_table.xls\n"
            "相对丰度分母：各样本全部 ITS reads（all_reads）。\n\n"
            "主要文件\n"
            "- ASV_annotated.xls：ASV 追加双库生态属性\n"
            "- coverage.xls：注释覆盖率\n"
            "- FG_trophicMode.* / FG_guild.*：FUNGuild 组成\n"
            "- FT_primary_lifestyle.*：FungalTraits 组成\n\n"
            "说明：此为基于分类学的生态属性推断，不是功能基因预测。\n"
            f"FUNGuild 命中 ASV {fg_hit}/{len(annotated)}；"
            f"FungalTraits 命中 ASV {ft_hit}/{len(annotated)}。\n"
        )
    if fg_db is None:
        print("警告: 未加载 FUNGuild，仅输出 FungalTraits（若可用）", flush=True)
    if ft_db is None:
        print("警告: 未加载 FungalTraits，仅输出 FUNGuild（若可用）", flush=True)
    print(out)


if __name__ == "__main__":
    main()
