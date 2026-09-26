# Pre-registration addendum - 2026-09-26 (rules 1-8 revival)
Locked BEFORE any new outcome is scored. Original registrations unchanged.

## A1 - Rescue-finding pivot (rule 6: stuck negative -> ChatGPT redirection first)
Verified state (live clone 0f0b4104, 19/19 tests green with PYTHONPATH=src 4:24 PM):
no experimentally validated rescue mutation; S669 PDB-group holdout GNN MAE 1.422 vs
GeoDDG-Seq 1.310 and DDMut 1.290 kcal/mol (honest non-beat, comparator training
overlap unknown). Locked next steps, in order:
1. Rule-6 ChatGPT redirection session on the stuck rescue finding (verbatim logged);
   pivot options independently verified before adoption.
2. Model-improvement iteration on the SAME PDB-group holdout (locked split
   unchanged): feature/architecture changes chosen BEFORE outcomes; each iteration a
   new addendum. BEAT gate = holdout MAE below 1.290 (DDMut) with paired bootstrap CI
   excluding the comparator value, comparator-overlap caveat documented.
3. Rescue-discovery criterion (locked): a rescue candidate counts only with (a)
   top-decile predicted stabilization on the disease structure, (b) independent
   literature or database support (documented suppressor or homolog evidence), (c)
   no train/test leakage through the PDB grouping. Zero-pass = documented negative +
   further pivot, never terminal.

## Judge rounds
Minimum 10, each producing a concrete novelty improvement (rule 8), verbatim in
JUDGE_ROUNDS.md.
