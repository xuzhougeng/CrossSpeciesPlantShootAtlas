Each species should include the following input data:

* **species.h5ad**: h5ad file saved by Scanpy
* **species.pep**: Protein sequences

If your data is in Seurat format, you can convert it using our provided script:

```bash
Rscript ./rds2h5ad.r input/species.rds data/mtx/species
```

The script will generate the following files:

- data/mtx/species.mtx
- data/mtx/species.rows
- data/mtx/speices.cols
- data/mtx/species.metadata.csv

You can load these files and save them as `.h5ad` with the following code:

```python
prefix = 'species'
mtx_file = "data/mtx/" + prefix + '.mtx'
row_file = "data/mtx/" + prefix + '.rows'
col_file = "data/mtx/" + prefix + '.cols'
metadata_file = "data/mtx/" + prefix + ".metadata.csv"

adata = sc.read_mtx(mtx_file)

genes = pd.read_csv(row_file, header=None)
cells = pd.read_csv(col_file, header=None)
metadata = pd.read_csv(metadata_file, index_col=0)

adata = adata.T
adata.var_names = genes[0]
adata.obs_names = cells[0]
adata.obs = metadata

# **Additional step:** Ensure that a `cell_type` column is included. For example
# adata.obs['cell_type'] = adata.obs['seurat_high']
adata.write_h5ad("data/h5ad/" + f'{prefix}.h5ad')
```

The SAMap/SATURN pipeline is an exploratory analytical process that supports multiple species. However, for convenience, we consider only two species in this context:

**species1**: 

- data/h5ad/sp1.h5ad
- data/pep/sp1.pep

**species2**: 

- data/h5ad/sp2.h5ad
- data/pep/sp2.pep


Meanwhile, we assume that there is one sample from sp2, with the batch information recorded using orig.ident.