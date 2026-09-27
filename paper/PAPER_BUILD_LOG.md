# Paper-build evidence log, 2026-09-27

Branch: `paper-build`, paper directory only.

- Added a locked evidence audit to the existing manuscript, tracing the corrected S669 669-row mapping and 66-variant held-out comparison to `results/s669_group_holdout.json`. The GNN MAE 1.4219 remains worse than stored GeoDDG-Seq 1.3101 on the same rows. The row-versus-protein weighting reversal and bootstrap interval come from `results/group_sensitivity.json`. Disease-target material is framed as context, not a measured rescue.
- The original manuscript specified `fontspec` and Times New Roman, but this sandbox lacks `fontspec` and licensed TNR. To rebuild, substituted `mathptmx` (Nimbus Roman, **not genuine TNR**) for that font declaration. pdfLaTeX output `paper/manuscript.pdf` has 21 pages by `pdfinfo`, from 20 previously. Visual check of pages 20-21 found readable audit text and no clipping. The 50-page and genuine-TNR gates remain open.

## Follow-up expansion, 2026-09-27

- Added `supplement_audit.tex` from committed `results/group_sensitivity.json`, `sparse_group_sensitivity.json`, `error_strata.json`, `leave_one_pdb_out.json`, `fireprot_discordance.json` and `fireprot_mapping_audit.json`: the 24 held-out PDB errors, threshold sensitivity, four observed-label bins, 94 per-PDB ridge errors, and partial-export discordance. It explicitly distinguishes a leave-one-PDB-out ridge rerun from the fixed 66-row GNN evaluation.
- Rebuilt twice with pdfLaTeX: 25 pages, up from 21, Times-style Nimbus Roman rather than licensed Times New Roman. Visual inspection of pages 22-23 found readable tables and corrected unresolved references/inequality glyphs. Still well short of 50 pages; no extra experiment or rescue claim was manufactured.

## Morning target-evidence expansion, 2026-09-27

- Added `target_audit.tex` from committed disease-site numbering, six wild/mutant structure overlays and TP53 function-report counts. Canonical-to-PDB offsets remain explicit (SOD1 5 to PDB 4; SOD1 94 to PDB 93). Structure and public functional reports are not a new GNN rescue or a matched double-mutant assay.
- pdfLaTeX twice now produces 27 pages, up from 25. Visual check of pages 26-27 confirms readable tables and no clipped content. Times-style Nimbus Roman substitute; 50-page and genuine-TNR gates remain open.
