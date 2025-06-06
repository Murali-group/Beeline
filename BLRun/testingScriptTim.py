import numpy as np
import pandas as pd

store = pd.HDFStore('rank_total_gene_rpkm.h5')
# store = pd.HDFStore('bone_marrow_cell.h5')
df = pd.read_hdf(store).reset_index()
print(df.shape)