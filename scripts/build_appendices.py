"""Create auditable full-data appendices from committed results; no fabricated examples."""
import json
from pathlib import Path
R=Path('results'); p=Path('paper/manuscript.tex'); s=p.read_text(); j=lambda n:json.loads((R/n).read_text())
loo=j('leave_one_pdb_out.json');group=j('group_sensitivity.json');count=j('accession_count.json');strata=j('error_strata.json')
head='\\clearpage\n\\section{Supplementary audit tables}\nThe tables below expose measured groups and their uneven row counts. They are not new cohorts, independent replicates, or prospective rescue screens. Their machine-readable records and predictions are in the result JSON files. Numeric precision in tables is rounded for display; all summary calculations use stored full-precision values.\n'
t=head
# 94 PDB leave-one-out groups
t+='\\subsection{Leave-one-PDB-out ridge baseline}\nEach row fits a train-only ridge model excluding the named structure. The GNN was not retrained in this table.\\\n\\begin{longtable}{lrrr}\\hline PDB & Variants & Ridge MAE & Train-mean MAE\\\\\\hline\\endfirsthead\\hline PDB & Variants & Ridge MAE & Train-mean MAE\\\\\\hline\\endhead\n'
for v in sorted(loo['groups'],key=lambda z:z['pdb']): t+=f"{v['pdb']} & {v['n']} & {v['ridge_mae']:.3f} & {v['train_mean_mae']:.3f}\\\\\n"
t+='\\hline\\end{longtable}\n'
t+='\\clearpage\\subsection{Fixed test-set protein-group error}\nThe signed difference is GNN MAE minus stored GeoDDG-Seq MAE on that PDB; positive means higher graph-model error. A one-variant row is not a reliable protein-wide comparison.\n\\begin{longtable}{lrrrr}\\hline PDB & Variants & GNN MAE & GeoDDG-Seq MAE & Difference\\\\\\hline\\endfirsthead\\hline PDB & Variants & GNN MAE & GeoDDG-Seq MAE & Difference\\\\\\hline\\endhead\n'
for v in sorted(group['grouped'],key=lambda z:z['pdb']):t+=f"{v['pdb']} & {v['n']} & {v['gnn_mae']:.3f} & {v['geo_mae']:.3f} & {v['delta']:+.3f}\\\\\n"
t+='\\hline\\end{longtable}\n'
t+='\\clearpage\\subsection{Accession-backed use ledger}\nA structure enters this list only if an S669 graph used it, the partial FireProtDB audit directly matched mutation numbering, or a separate observed disease-mutant comparison used it. Multiple rows do not represent independent source studies. Full SHA-256 digests appear in the machine-readable ledger; short prefixes below are audit aids, not a substitute for full checksums.\n\\begin{longtable}{llll}\\hline PDB & S669 graph & FireProtDB direct & SHA-256 prefix\\\\\\hline\\endfirsthead\\hline PDB & S669 graph & FireProtDB direct & SHA-256 prefix\\\\\\hline\\endhead\n'
for v in count['accessions']:
 uses=v['used_in'];t+=f"{v['pdb']} & {'yes' if 'S669 graph' in uses else 'no'} & {'yes' if 'FireProtDB direct mutation-number match' in uses else 'no'} & \\texttt{{{v['sha256'][:12]}}}\\\\\n"
t+='\\hline\\end{longtable}\n'
# Annotate panel and short methodological discussion to make audit readable rather than pure pages.
t+='\\subsection{What these counts cannot establish}\nThe accession total counts structures used by distinct computational analyses, not distinct biological experiments. A PDB can contain variants with shared sequence background, assay labels can derive from related publications, and S669 predictor columns may have unknown historical training overlap with these proteins. FireProtDB is a partial download with a direct-numbering screen, not an exhaustively harmonized stability dataset. The fixed protein holdout cannot by itself justify a disease-rescue nomination. Published p53 suppressor evidence shows that functional rescue may happen with little stability gain, so a future follow-up needs disease-mutant-specific double-mutant measurements and functional endpoints. The guarded rank-check command deliberately abstains in their absence.\n'
marker='\\section{Source ledger}'
assert s.count(marker)==1;s=s.replace(marker,t+'\n'+marker);p.write_text(s);print('appended',len(loo['groups']),len(group['grouped']),len(count['accessions']))
