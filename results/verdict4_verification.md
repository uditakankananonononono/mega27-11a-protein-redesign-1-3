# Verdict #4 factual re-derivation - 11a-owned items (addendum 2026-09-27)

Re-derived from committed results, 2026-09-27. Underlying analyses live in the
shared 11a/11b workspace (11b repo paths noted); verdict text archived verbatim
in docs/JUDGE_VERDICT_USER_2026-09-27.md.

| Verdict claim | Committed evidence | Status |
|---|---|---|
| MYOGLOBIN r=0.455 contaminated by 41 same-protein rows; drops to 0.24-0.29 | 11b results/myoglobin_eval.md: full r=0.455, 41 training rows same protein (1BVC), 53/134 exact duplicates; decontaminated r=0.288 (excl. exact), 0.241 (excl. same-position) | CONFIRMED |
| 2VUK carries superstable quadruple-mutant background | 11b results/rescue_scan_2vuk_v2.md lines 25/65: M133L/V203A/N239Y/N268D superstable background already present in 2VUK | CONFIRMED |
| Greedy multi-mutant sums implausible (+46 kcal/mol) | 11b results/rescue_scan_2vuk_v2.md line 13: greedy sum +46.0 kcal/mol, flagged implausible, ranks only; v2 recomputation +42.6 under independence, same caveat | CONFIRMED |
| v1 aromatic identity bias up to +1.04 kcal/mol | 11b results/model_v1_bias_analysis.md: ->TYR +1.04, ->TRP +0.89, ->ILE +0.85; identity shortcut, drove v2 redesign | CONFIRMED |

Disposition: all four 11a-owned claims stand on committed evidence, and all
four were ALREADY disclosed in our own audits before the verdict - the verdict
independently reached the same numbers. Locked negatives, no claim softened.
