# Paper-build evidence log, 2026-09-27

Branch: `paper-build`, paper directory only.

- Added a locked evidence audit to the existing manuscript, tracing the corrected S669 669-row mapping and 66-variant held-out comparison to `results/s669_group_holdout.json`. The GNN MAE 1.4219 remains worse than stored GeoDDG-Seq 1.3101 on the same rows. The row-versus-protein weighting reversal and bootstrap interval come from `results/group_sensitivity.json`. Disease-target material is framed as context, not a measured rescue.
- The original manuscript specified `fontspec` and Times New Roman, but this sandbox lacks `fontspec` and licensed TNR. To rebuild, substituted `mathptmx` (Nimbus Roman, **not genuine TNR**) for that font declaration. pdfLaTeX output `paper/manuscript.pdf` has 21 pages by `pdfinfo`, from 20 previously. Visual check of pages 20-21 found readable audit text and no clipping. The 50-page and genuine-TNR gates remain open.
