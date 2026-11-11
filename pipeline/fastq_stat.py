#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""统计 FASTQ 条数、碱基数和平均长度，并写出与报告一致的 xls 表。"""

import argparse
import gzip
import os


def open_text(path):
    """按是否 gzip 打开文本。"""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8", errors="replace")


def stat_fastq(path):
    """统计一个 FASTQ 的序列条数和碱基数。"""
    n = 0
    bases = 0
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip()
            handle.readline()
            handle.readline()
            n += 1
            bases += len(seq)
    # 模板平均长度是向零截断到两位小数，末尾的 0 不写出，例如 408.40 写成 408.4
    avg = format_avg(bases, n)
    return n, bases, avg


def format_avg(bases, n):
    """碱基数除以条数，向零截断到两位小数。"""
    if not n:
        return "0"
    scaled = (bases * 100) // n
    whole, frac = divmod(scaled, 100)
    return f"{whole}.{frac:02d}".rstrip("0").rstrip(".")


def write_tables(raw_counts, valid_dir, raw_out, valid_out):
    """根据原始条数和 valid 目录中的合并序列写出两张统计表。"""
    order = []
    with open(raw_counts, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            sample = line.split("\t", 1)[0]
            order.append(sample)

    # 模板的 rawdata.stat.xls 与 valid.stat.xls 条数相同，统计的都是合并后的有效序列
    valid_rows = []
    for sample in order:
        fq = os.path.join(valid_dir, f"trim.{sample}.fq.gz")
        n, bases, avg = stat_fastq(fq)
        valid_rows.append((sample, n, bases, avg))

    with open(raw_out, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("Sample\tSequences\n")
        for sample, n, bases, avg in valid_rows:
            handle.write(f"{sample}\t{n}\n")

    with open(valid_out, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("Sample\tSequences\tBases(bp)\tAverageLength(bp)\n")
        for sample, n, bases, avg in valid_rows:
            handle.write(f"{sample}\t{n}\t{bases}\t{avg}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-counts", required=True)
    parser.add_argument("--valid-dir", required=True)
    parser.add_argument("--raw-out", required=True)
    parser.add_argument("--valid-out", required=True)
    args = parser.parse_args()
    write_tables(args.raw_counts, args.valid_dir, args.raw_out, args.valid_out)


if __name__ == "__main__":
    main()
