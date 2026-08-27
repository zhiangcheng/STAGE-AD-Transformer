import os
import pandas as pd
import glob


def get_tiss_name(path):
    fname = os.path.basename(path) 
    parts = fname.split(".")
    return parts[0]
 
GOT_BASE = False
og_phenotype_groups = glob.glob("/media/desk15/iy2120/TWAS/Breast-cancer-Example/data/subtype_gwas/input/by_tiss_phenotype_groups/*.gz") # Download the data from https://predictdb.org/post/2023/02/22/how-to-map-introns-to-genes/ or https://uchicago.box.com/s/xy71r0su6refrfggivrwc7kpvsxysmmn

for p_group in og_phenotype_groups:
    tissue = get_tiss_name(p_group)
    df = pd.read_table(p_group)
    df['intron_id'] = "{}.".format(tissue) + df['intron_id']
 
    if not GOT_BASE:
        base_df = df
        GOT_BASE = True
    else:
        base_df = base_df.append(df)

base_df = base_df.drop_duplicates()
base_df.to_csv("/media/desk15/iy2120/TWAS/Breast-cancer-Example/data/subtype_gwas/input/combine_phenotype_groups.txt.gz", sep="\t", index=False,compression="gzip")
