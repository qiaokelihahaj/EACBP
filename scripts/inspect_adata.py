import sys
from pathlib import Path
import anndata as ad
import numpy as np

p = sys.argv[1] if len(sys.argv) > 1 else "/public/home/qiaoke/eacbp_project/test_kb_out/counts_unfiltered/adata.h5ad"
adata = ad.read_h5ad(p)
print(f"=== Real AnnData Quantified from FASTQ ===")
print(f"AnnData shape: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"Total non-zero entries: {adata.X.nnz if hasattr(adata.X, 'nnz') else np.count_nonzero(adata.X)}")
print(f"First 5 barcodes: {list(adata.obs_names[:5])}")
print(f"First 10 genes: {list(adata.var_names[:10])}")
if "Kat8" in adata.var_names:
    print(f"Target gene 'Kat8' found in dataset!")
elif "gene_name" in adata.var and "Kat8" in adata.var["gene_name"].values:
    print(f"Target gene 'Kat8' found in var['gene_name']!")
