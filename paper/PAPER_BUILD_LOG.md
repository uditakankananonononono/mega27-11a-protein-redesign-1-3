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

## Structure-input reproducibility appendix, 2026-09-27

- Added `accession_manifest.tex`, a compressed human-readable index of the 163 PDB accessions actually used, their analysis use category and a 12-character SHA-256 prefix drawn from committed `results/accession_count.json`. Full hashes remain in that JSON. Added the narrow exact-sequence-overlap audit (69 train, 24 test, zero exact overlap); this does not rule out homologs or the third-party comparator's historical training overlap. The table is a source ledger, not independent experiments or rescue results.
- Rebuilt twice: 31 pages from 27. Visually checked new pages 28-31; repeated headings, table columns and hashes read cleanly without clipping. Nimbus Roman substitute. Page 31 is sparse after the 163-entry ledger ends; 50 substantive pages and genuine Times New Roman remain open.

## Judge requirement amended, 2026-09-27 10:00 IST

The owner changed the counted ChatGPT check from ten rounds to **one round per project, supplied by her** (authenticated WhatsApp message `wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMDJCMTZGRTVEMkQwMTFBQzc4MQA=`, 10:00:07 IST). Earlier ten-round language in scientific preregistration and root judge ledgers is retained as historical text pending a canonical science-branch update. The paper branch does not claim the one-round gate met merely because a ChatGPT conversation exists: the user-provided verdict and its integration must be linked to this project in the canonical ledger. Supplementary Gemini or other LLM consults remain logged but do not count. This amendment changes a process gate, not any scientific result, page or font gate.
