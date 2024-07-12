import pandas as pd
from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names
from distributed import Client, LocalCluster
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, default='')
item = parser.parse_args().filename
dt = pd.read_csv(item,index_col=0)
dt.columns = [str(x) for x in dt.columns]
print(dt)
TFname= ['KLF6', 'TCEB3', 'LYL1', 'SMARCC1', 'TCOF1', 'ZNF267', 'ZEB2', 'MNDA', 'ETS2', 'BAZ2B', 'POU2F2', 'MEF2C', 'KDM5A', 'PDLIM7', 'HDGF', 'ZBTB16', 'ZNF350', 'STAT3', 'TAF1B', 'HIST2H2BE', 'DHX38', 'TP53', 'SMAD3', 'MXD4', 'ARID5B', 'USF2', 'KDM2A', 'HIVEP3', 'MYBL1', 'HIST1H1E', 'ZNF593', 'BATF', 'TAX1BP3', 'TRIM28', 'CBFB', 'CHD4', 'ZBTB38', 'PBX2', 'CTNNBIP1', 'SERTAD2', 'ZMYND11', 'NCOA4', 'PER1', 'ID3', 'POLR2A', 'CDKN1A', 'TGFB1', 'ZNF277', 'MAPK1', 'NEAT1', 'SP3', 'MAX', 'SMARCA2', 'REL', 'SIN3A', 'NR4A1', 'ASCL2', 'JUND', 'TFDP2', 'BHLHE40', 'NFKBIA', 'HTT', 'SOX4', 'SPI1', 'FOS', 'CITED2', 'CREM', 'PURA', 'HEXIM1', 'PKNOX1', 'CEBPB', 'HHEX', 'BRD8', 'RUNX3', 'MAFB', 'EOMES', 'SERTAD3', 'ZNF143', 'ZNF467', 'AKT1', 'ATF6', 'PTTG1', 'TBX21', 'UIMC1', 'IRF5', 'EED', 'ID1', 'IRF8', 'HOPX', 'SUGP2', 'JUN', 'TAF6L', 'PDLIM1', 'SPIB', 'HIST1H1C', 'RNF19A', 'CREBBP', 'IRF1', 'SUZ12', 'CHD8', 'HDAC5', 'BLZF1', 'SHPRH', 'CUX1', 'RELB', 'GTF3C1', 'FOSB', 'MLXIP', 'NFIC', 'IRF7', 'BBC3', 'GTF2I', 'MKL1', 'POLR1C', 'CEBPD', 'SMARCD2', 'IKZF3', 'SLA2']
client = Client(processes = False)
gene_name = list(dt.columns)
print(gene_name)
TFname = list(set(TFname) & set(gene_name))
print(dt)
network = grnboost2(dt, client_or_address = client,gene_names=list(dt.columns),tf_names=TFname)
network.to_csv(item+'TF_grnboost2.csv')
