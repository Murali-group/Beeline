import numpy as np
import pandas as pd

data = np.load("inputs/human-scRNA/hESC/hESC-500-ChIP/CNNC/NEPDF_data/Nxdata_tf0.npy")
print(data.shape)
ydata = np.load("inputs/human-scRNA/hESC/hESC-500-ChIP/CNNC/NEPDF_data/ydata_tf0.npy")
print(ydata.size)
zdata = np.load("inputs/human-scRNA/hESC/hESC-500-ChIP/CNNC/NEPDF_data/zdata_tf0.npy")
print(zdata.size)
