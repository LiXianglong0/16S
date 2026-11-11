#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""在已有 ASV 表上补充多样性、样本间距离和 ASV 分布，不改动上游质控与 DADA2。"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

FONT = "/mnt/c/Windows/Fonts/msyh.ttc"
if os.path.isfile(FONT):
    font_manager.fontManager.addfont(FONT)
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei"]
plt.rcParams["axes.unicode_minus"] = False

PREFERRED_GROUP_ORDER = ["CK", "GRA", "S200-100", "S400-100", "S700-100"]
PREFERRED_COLORS = {
    "CK": "#4C78A8",
    "GRA": "#F58518",
    "S200-100": "#54A24B",
    "S400-100": "#E45756",
    "S700-100": "#72B7B2",
}
PALETTE = [
    "#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2",
    "#B279A2", "#FF9DA6", "#9D755D", "#BAB0AC", "#4C9A8B",
    "#F2B5D4", "#86BCB6", "#F6C85F", "#8B6B4C", "#6B5B95",
    "#88B04B", "#F7CAC9", "#92A8D1", "#955251", "#B565A7",
    "#009B77", "#DD4124", "#D65076", "#45B8AC",
]
N_PERM = 999
RARE_REPS = 8
TOP_N = 20


def group_of(sample):
    """样本名去掉末尾重复号（如 CK-1 → CK）；无连字符则整名作为分组。"""
    if "-" in sample:
        head, tail = sample.rsplit("-", 1)
        if tail.isdigit():
            return head
    return sample


def ordered_groups(samples):
    """优先沿用已知分组顺序，其余按名称排序追加。"""
    seen = []
    for s in samples:
        g = group_of(s)
        if g not in seen:
            seen.append(g)
    preferred = [g for g in PREFERRED_GROUP_ORDER if g in seen]
    rest = sorted(g for g in seen if g not in preferred)
    return preferred + rest


def color_of(group, group_names=None):
    if group in PREFERRED_COLORS:
        return PREFERRED_COLORS[group]
    names = group_names or []
    if group in names:
        return PALETTE[names.index(group) % len(PALETTE)]
    return PALETTE[hash(group) % len(PALETTE)]


def load_otu(path):
    """读取 biom 转出的 otu_table.tsv。第一列是 ASV，其余列是样本计数。"""
    with open(path, "r", encoding="utf-8") as handle:
        lines = []
        for line in handle:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("# "):
                continue
            lines.append(line)
    header = lines[0].split("\t")
    if header[0] not in ("#OTU ID", "OTU ID"):
        raise SystemExit(f"无法识别特征表表头: {path}")
    samples = header[1:]
    ids = []
    rows = []
    for line in lines[1:]:
        parts = line.split("\t")
        ids.append(parts[0])
        rows.append([int(round(float(x))) for x in parts[1:]])
    table = np.asarray(rows, dtype=np.int64)
    return ids, samples, table


def shannon(counts):
    total = counts.sum()
    if total <= 0:
        return 0.0
    p = counts[counts > 0] / total
    return float(-(p * np.log(p)).sum())


def simpson(counts):
    total = counts.sum()
    if total <= 0:
        return 0.0
    p = counts[counts > 0] / total
    return float(1.0 - (p * p).sum())


def chao1(counts):
    obs = int((counts > 0).sum())
    f1 = int((counts == 1).sum())
    f2 = int((counts == 2).sum())
    if f2 > 0:
        return obs + (f1 * f1) / (2.0 * f2)
    if f1 > 1:
        return obs + f1 * (f1 - 1) / 2.0
    return float(obs)


def coverage(counts):
    total = int(counts.sum())
    if total <= 0:
        return 0.0
    f1 = int((counts == 1).sum())
    return 1.0 - f1 / total


def alpha_table(samples, table):
    rows = []
    for i, sample in enumerate(samples):
        counts = table[:, i]
        rows.append({
            "sample": sample,
            "group": group_of(sample),
            "depth": int(counts.sum()),
            "observed": int((counts > 0).sum()),
            "shannon": shannon(counts),
            "simpson": simpson(counts),
            "chao1": chao1(counts),
            "coverage": coverage(counts),
        })
    return rows


def rarefy_curves(samples, table, rng):
    """每个样本按自身测序深度做稀释曲线，重复取平均。"""
    curves = {}
    for i, sample in enumerate(samples):
        counts = table[:, i]
        depth = int(counts.sum())
        if depth < 1000:
            continue
        points = np.unique(np.linspace(1000, depth, 10, dtype=int))
        means = []
        for point in points:
            observed = []
            for _ in range(RARE_REPS):
                draw = rng.multivariate_hypergeometric(counts, int(point))
                observed.append(int(np.count_nonzero(draw)))
            means.append(float(np.mean(observed)))
        curves[sample] = (points, means)
    return curves


def bray_curtis(table):
    """相对丰度上的 Bray-Curtis 距离。"""
    totals = table.sum(axis=0).astype(float)
    totals[totals == 0] = 1.0
    rel = table / totals
    n = rel.shape[1]
    dist = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            denom = rel[:, i].sum() + rel[:, j].sum()
            value = np.abs(rel[:, i] - rel[:, j]).sum() / denom if denom else 0.0
            dist[i, j] = dist[j, i] = value
    return dist


def pcoa(dist):
    """经典主坐标分析。负特征值不进入坐标。"""
    n = dist.shape[0]
    d2 = dist * dist
    cent = np.eye(n) - np.ones((n, n)) / n
    matrix = -0.5 * cent @ d2 @ cent
    eigvals, eigvecs = np.linalg.eigh(matrix)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    positive = eigvals > 1e-10
    pos_sum = eigvals[positive].sum()
    coords = np.zeros((n, 2))
    variance = [0.0, 0.0]
    used = 0
    for k in range(len(eigvals)):
        if eigvals[k] <= 1e-10:
            continue
        if used < 2:
            coords[:, used] = eigvecs[:, k] * np.sqrt(eigvals[k])
            variance[used] = 100.0 * eigvals[k] / pos_sum if pos_sum else 0.0
        used += 1
    return coords, variance


def permanova(dist, groups, rng):
    """Anderson PERMANOVA。置换分组标签，不置换距离。"""
    labels = np.asarray(groups)
    uniq = []
    for g in labels:
        if g not in uniq:
            uniq.append(g)

    def pseudo_f(lab):
        n = len(lab)
        ss_t = np.sum(dist[np.triu_indices(n, 1)] ** 2) / n
        ss_w = 0.0
        for name in uniq:
            idx = np.where(lab == name)[0]
            ng = len(idx)
            if ng < 2:
                continue
            sub = dist[np.ix_(idx, idx)]
            ss_w += np.sum(sub[np.triu_indices(ng, 1)] ** 2) / ng
        ss_a = ss_t - ss_w
        a = len(uniq)
        df_a = a - 1
        df_w = n - a
        if ss_w <= 0 or df_a <= 0 or df_w <= 0:
            return 0.0, ss_a, ss_t
        return (ss_a / df_a) / (ss_w / df_w), ss_a, ss_t

    f_obs, ss_a, ss_t = pseudo_f(labels)
    greater = 0
    for _ in range(N_PERM):
        perm = rng.permutation(labels)
        f_perm, _, _ = pseudo_f(perm)
        if f_perm >= f_obs:
            greater += 1
    p_value = (greater + 1) / (N_PERM + 1)
    r2 = ss_a / ss_t if ss_t else 0.0
    return f_obs, r2, p_value


def presence_by_group(samples, table):
    """组内任一重复计数大于 0 即视为该组存在。"""
    groups = [group_of(s) for s in samples]
    names = ordered_groups(samples)
    present = {}
    for name in names:
        idx = [i for i, g in enumerate(groups) if g == name]
        present[name] = table[:, idx].sum(axis=1) > 0
    return names, present


def write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(x) for x in row) + "\n")


def setup_font_ok():
    return os.path.isfile(FONT)


def plot_rarefaction(curves, samples, out_path):
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    group_names = ordered_groups(samples)
    for sample in samples:
        if sample not in curves:
            continue
        xs, ys = curves[sample]
        group = group_of(sample)
        ax.plot(
            xs, ys,
            color=color_of(group, group_names),
            linewidth=1.4,
            label=group,
        )
    handles, labels = ax.get_legend_handles_labels()
    seen = {}
    for handle, label in zip(handles, labels):
        seen.setdefault(label, handle)
    # 样本即分组时图例过多，只保留前若干项
    items = list(seen.items())  # (label, handle)
    if len(items) > 12:
        ax.legend(
            [h for _, h in items[:12]],
            [n for n, _ in items[:12]],
            frameon=False, title="分组(部分)", fontsize=7,
        )
    else:
        ax.legend([h for _, h in items], [n for n, _ in items], frameon=False, title="分组")
    ax.set_xlabel("抽取序列数")
    ax.set_ylabel("观察到的 ASV 数")
    ax.set_title("稀释曲线")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def plot_pcoa(samples, coords, variance, out_path):
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    group_names = ordered_groups(samples)
    for group in group_names:
        idx = [i for i, s in enumerate(samples) if group_of(s) == group]
        if not idx:
            continue
        ax.scatter(
            coords[idx, 0], coords[idx, 1],
            s=46, color=color_of(group, group_names), label=group, zorder=3,
        )
        for i in idx:
            ax.text(coords[i, 0], coords[i, 1], " " + samples[i], fontsize=7, color="#444444")
    ax.axhline(0, color="#dddddd", linewidth=0.8)
    ax.axvline(0, color="#dddddd", linewidth=0.8)
    ax.set_xlabel(f"PCoA1 ({variance[0]:.1f}%)")
    ax.set_ylabel(f"PCoA2 ({variance[1]:.1f}%)")
    ax.set_title("Bray-Curtis 主坐标分析")
    if len(group_names) <= 12:
        ax.legend(frameon=False, title="分组")
    else:
        ax.legend(frameon=False, title="分组", fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def plot_heatmap(ids, samples, table, out_path):
    totals = table.sum(axis=1)
    order = np.argsort(totals)[::-1][:TOP_N]
    sub_ids = [ids[i][:8] for i in order]
    sub = table[order, :].astype(float)
    depth = sub.sum(axis=0)
    depth[depth == 0] = 1.0
    rel = sub / depth * 100.0
    if rel.shape[0] > 2:
        row_order = leaves_list(linkage(rel, method="average"))
        rel = rel[row_order, :]
        sub_ids = [sub_ids[i] for i in row_order]
    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    image = ax.imshow(rel, aspect="auto", cmap="YlGnBu")
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(sub_ids)))
    ax.set_yticklabels(sub_ids, fontsize=8)
    ax.set_title("丰度前 20 的 ASV（相对丰度 %）")
    fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02, label="%")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return sub_ids, rel


def html_table(headers, rows):
    head = "".join(f"<th>{h}</th>" for h in headers)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{c}</td>" for c in row) + "</tr>")
    return (
        '<table class="table table-bordered table-striped" style="width: auto;">'
        f"<thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"
    )


def build_section(alpha, f_stat, r2, p_value, groups, present, min_depth):
    alpha_rows = []
    for row in alpha:
        alpha_rows.append([
            row["sample"], row["group"], row["depth"], row["observed"],
            f"{row['shannon']:.2f}", f"{row['simpson']:.3f}",
            f"{row['chao1']:.1f}", f"{row['coverage']:.3f}",
        ])
    means = []
    for name in groups:
        part = [r for r in alpha if r["group"] == name]
        means.append([
            name, len(part),
            f"{np.mean([r['observed'] for r in part]):.1f}",
            f"{np.mean([r['shannon'] for r in part]):.2f}",
            f"{np.mean([r['simpson'] for r in part]):.3f}",
        ])
    uniq_rows = []
    for name in groups:
        mask = present[name]
        others = np.zeros(mask.shape, dtype=bool)
        for other in groups:
            if other != name:
                others |= present[other]
        uniq_rows.append([name, int(mask.sum()), int((mask & ~others).sum())])
    pair_header = ["分组"] + groups
    pair_rows = []
    for left in groups:
        line = [left]
        for right in groups:
            line.append(int((present[left] & present[right]).sum()))
        pair_rows.append(line)
    all_mask = np.ones(next(iter(present.values())).shape, dtype=bool)
    for name in groups:
        all_mask &= present[name]
    p_text = "<0.001" if p_value < 0.001 else f"{p_value:.3f}"
    return f"""<!-- supplement-start -->
        <h3 id="chapter1_3" class="sub-title">1.3 多样性分析</h3>
        <p>本节直接使用 1.2.3 导出的 ASV 表计算多样性，不重新做质控、合并或 DADA2。</p>
        <p>分组规则：样本名若为「组名-重复号」（如 CK-1），则去掉末尾数字作为分组；否则整名作为分组。当前分组：{"、".join(groups)}（共 {len(groups)} 组）。多样性按 ASV 计数计算。Shannon 为自然对数，Simpson 为 1 − Σp²。稀释曲线在每个样本自身深度内重复抽取 {RARE_REPS} 次。Bray-Curtis 使用相对丰度。PERMANOVA 置换 {N_PERM} 次。</p>

        <h4 id="chapter1_3_1" class="sub-title">1.3.1 Alpha 多样性</h4>
        <p>Depth 是该样本 ASV 表中的序列数，Observed 是观察到的 ASV 数，Chao1 是丰富度估计，Coverage 是 Good's coverage。各组测序深度接近，最低深度为 {min_depth}。</p>
        {html_table(["Sample", "Group", "Depth", "Observed", "Shannon", "Simpson", "Chao1", "Coverage"], alpha_rows)}
        <p>各组平均值：</p>
        {html_table(["Group", "N", "Observed", "Shannon", "Simpson"], means)}
        <div align="center"><img src="supplement/rarefaction.png" style="max-width: 820px; width: 100%;"></div>
        <p>稀释曲线趋于平缓时，说明继续增加测序深度也不会明显增加 ASV 数。</p>

        <h4 id="chapter1_3_2" class="sub-title">1.3.2 Beta 多样性</h4>
        <p>样本间距离为 Bray-Curtis。主坐标分析展示样本远近。PERMANOVA 检验组间距离是否大于组内距离。pseudo-F = {f_stat:.3f}，R² = {r2:.3f}，p = {p_text}。若每组重复较少或样本即分组，p 值仅作描述参考。</p>
        <div align="center"><img src="supplement/pcoa.png" style="max-width: 720px; width: 100%;"></div>

        <h4 id="chapter1_3_3" class="sub-title">1.3.3 ASV 分布</h4>
        <p>一个 ASV 在组内任一重复中计数大于 0，即记为该组存在。各组共同存在的 ASV 有 {int(all_mask.sum())} 个。下表对角线是该组 ASV 数，非对角线是两组共有 ASV 数。</p>
        {html_table(["Group", "ASV数", "本组特有"], uniq_rows)}
        {html_table(pair_header, pair_rows)}
        <p>热图给出总丰度最高的 {TOP_N} 个 ASV 在各样本中的相对丰度，行按丰度模式聚类，列保持样本顺序。ASV 名称取特征编号前 8 位。</p>
        <div align="center"><img src="supplement/heatmap.png" style="max-width: 920px; width: 100%;"></div>
        <br>
<!-- supplement-end -->
"""


def strip_markers(text):
    """去掉旧的补充章节与导航，便于重新插入。"""
    for start_mark, end_mark in (
        ("<!-- supplement-start -->", "<!-- supplement-end -->"),
        ("<!-- supplement-nav-start -->", "<!-- supplement-nav-end -->"),
    ):
        while True:
            start = text.find(start_mark)
            end = text.find(end_mark)
            if start == -1 or end == -1 or end <= start:
                break
            text = text[:start] + text[end + len(end_mark) :]
    return text


def inject_into_report(report_path, section):
    """把多样性章节并入 amplicon-ngs.report.html，不再单独生成 diversity.report.html。"""
    if not os.path.isfile(report_path):
        raise SystemExit(f"还没有主报告，请先运行 make_report.py: {report_path}")
    with open(report_path, "r", encoding="utf-8") as handle:
        text = handle.read()
    text = strip_markers(text)

    nav = """<!-- supplement-nav-start -->
                    <li><a href="#chapter1_3">&nbsp;&nbsp;1.3 多样性分析</a></li>
                    <li><a href="#chapter1_3_1">&nbsp;&nbsp;&nbsp;&nbsp;1.3.1 Alpha多样性</a></li>
                    <li><a href="#chapter1_3_2">&nbsp;&nbsp;&nbsp;&nbsp;1.3.2 Beta多样性</a></li>
                    <li><a href="#chapter1_3_3">&nbsp;&nbsp;&nbsp;&nbsp;1.3.3 ASV分布</a></li>
<!-- supplement-nav-end -->
"""
    nav_anchor = '<li><a href="#chapter1_2_3">&nbsp;&nbsp;1.2.3 qiime2导入valid数据</a></li>'
    if nav_anchor in text and "<!-- supplement-nav-start -->" not in text:
        text = text.replace(nav_anchor, nav_anchor + "\n" + nav, 1)

    # 插在「附录」之前
    appendix_mark = '<h2 id="chapter2" class="main-title">附录</h2>'
    if appendix_mark not in text:
        raise SystemExit(f"主报告找不到附录锚点，无法合并多样性章节: {report_path}")
    text = text.replace(appendix_mark, section + "\n        " + appendix_mark, 1)

    with open(report_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)

    # 旧的独立多样性报告不再维护
    old = os.path.join(os.path.dirname(report_path), "diversity.report.html")
    if os.path.isfile(old):
        os.remove(old)
    return report_path


def main():
    if len(sys.argv) != 2:
        raise SystemExit("用法: python3 supplement_diversity.py 项目目录")
    project = os.path.abspath(sys.argv[1])
    otu = os.path.join(project, "qiime2", "export", "otu_table.tsv")
    if not os.path.isfile(otu):
        raise SystemExit(f"还没有 ASV 表: {otu}")
    report = os.path.join(project, "amplicon-ngs.report.html")
    ids, samples, table = load_otu(otu)
    rng = np.random.default_rng(1)
    alpha = alpha_table(samples, table)
    curves = rarefy_curves(samples, table, rng)
    dist = bray_curtis(table)
    coords, variance = pcoa(dist)
    groups = [group_of(s) for s in samples]
    f_stat, r2, p_value = permanova(dist, groups, rng)
    group_names, present = presence_by_group(samples, table)

    out_dir = os.path.join(project, "supplement")
    os.makedirs(out_dir, exist_ok=True)
    write_tsv(
        os.path.join(out_dir, "alpha_diversity.tsv"),
        ["Sample", "Group", "Depth", "Observed", "Shannon", "Simpson", "Chao1", "Coverage"],
        [[r["sample"], r["group"], r["depth"], r["observed"], f"{r['shannon']:.4f}",
          f"{r['simpson']:.4f}", f"{r['chao1']:.2f}", f"{r['coverage']:.4f}"] for r in alpha],
    )
    write_tsv(
        os.path.join(out_dir, "bray_curtis.tsv"),
        ["Sample"] + samples,
        [[samples[i]] + [f"{dist[i, j]:.6f}" for j in range(len(samples))] for i in range(len(samples))],
    )
    write_tsv(
        os.path.join(out_dir, "permanova.tsv"),
        ["Test", "pseudoF", "R2", "permutations", "p"],
        [["group", f"{f_stat:.4f}", f"{r2:.4f}", N_PERM, f"{p_value:.4f}"]],
    )
    plot_rarefaction(curves, samples, os.path.join(out_dir, "rarefaction.png"))
    plot_pcoa(samples, coords, variance, os.path.join(out_dir, "pcoa.png"))
    plot_heatmap(ids, samples, table, os.path.join(out_dir, "heatmap.png"))
    section = build_section(
        alpha, f_stat, r2, p_value, group_names, present, min(r["depth"] for r in alpha)
    )
    out = inject_into_report(report, section)
    print(f"pseudo-F={f_stat:.3f} R2={r2:.3f} p={p_value:.3f}")
    print(out)


if __name__ == "__main__":
    main()
