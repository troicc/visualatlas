"""Audit RMechDB with the formal checker: replay every step's arrows and compare with the written product.

【标注】本文件由 AI（Claude）编写，作者待精读，见 docs/STUDY_checker.md。

    python scripts/audit_rmechdb.py <path to rmechdb_data/all.csv> [out.csv]

Prints counts per issue code and per step type; optionally writes every failing row (line number, step type,
issue codes, messages) to out.csv. The data themselves are CC BY-NC-ND and stay outside the repository.
"""

import csv
import sys
from collections import Counter

from rdkit import RDLogger

from visualatlas import rmechdb
from visualatlas.check import check_step

RDLogger.DisableLog("rdApp.*")


def main(path, out=None):
    total, parse_fail, ok = 0, 0, 0
    by_code, by_kind, fail_by_kind = Counter(), Counter(), Counter()
    rows = []
    for line, rec in rmechdb.read(path):
        total += 1
        if isinstance(rec, Exception):
            parse_fail += 1
            rows.append((line, "", "parse_error", repr(rec)[:300]))
            continue
        by_kind[rec.kind] += 1
        try:
            issues = check_step(rec.step, rec.product)
        except Exception as e:  # noqa: BLE001
            issues = None
            rows.append((line, rec.kind, "checker_error", repr(e)[:300]))
            by_code["checker_error"] += 1
            fail_by_kind[rec.kind] += 1
            continue
        if not issues:
            ok += 1
            continue
        fail_by_kind[rec.kind] += 1
        for c in {i.code for i in issues}:
            by_code[c] += 1
        rows.append(
            (line, rec.kind, ";".join(sorted({i.code for i in issues})), " | ".join(i.message for i in issues)[:500])
        )
    print(f"steps: {total}   parsed: {total - parse_fail}   pass: {ok}   fail: {total - parse_fail - ok}")
    print("failing steps per issue code (a step can have several):")
    for c, n in by_code.most_common():
        print(f"  {c:16s} {n}")
    print("fail / total per step type:")
    for k, n in by_kind.most_common():
        print(f"  {k:16s} {fail_by_kind[k]} / {n}")
    if out:
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["line", "kind", "codes", "messages"])
            w.writerows(rows)


if __name__ == "__main__":
    main(*sys.argv[1:3])
