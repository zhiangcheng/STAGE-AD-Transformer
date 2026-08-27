# Make sure all .db files have been closed in sqlite3!
import sqlite3
import pandas as pd
import gzip
import os
import shutil

# Configure paths
db_dir = "/media/desk15/iy2120/TWAS/Breast-cancer-Example/pipeline/05_focus/input/patched/sqtl_dbs/before_patched/"
map_dir = "/media/desk15/iy2120/TWAS/Breast-cancer-Example/data/subtype_gwas/input/by_tiss_phenotype_groups/"  # Directory containing .leafcutter.phenotype_groups.txt.gz files
annotation_file = "/media/desk15/iy2120/TWAS/Breast-cancer-Example/pipeline/05_focus/input/gencode_v26_annotation.txt"  # Replace with the actual path

tissues = [
    'Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial',
    'Brain_Amygdala', 'Brain_Anterior_cingulate_cortex_BA24', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere',
    'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus', 'Brain_Hypothalamus',
    'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia', 'Brain_Spinal_cord_cervical_c-1', 'Brain_Substantia_nigra',
    'Breast_Mammary_Tissue', 'Cells_Cultured_fibroblasts', 'Cells_EBV-transformed_lymphocytes', 'Colon_Sigmoid', 'Colon_Transverse',
    'Esophagus_Gastroesophageal_Junction', 'Esophagus_Mucosa', 'Esophagus_Muscularis', 'Heart_Atrial_Appendage', 'Heart_Left_Ventricle',
    'Kidney_Cortex', 'Liver', 'Lung', 'Minor_Salivary_Gland', 'Muscle_Skeletal', 'Nerve_Tibial', 'Ovary', 'Pancreas', 'Pituitary',
    'Prostate', 'Skin_Not_Sun_Exposed_Suprapubic', 'Skin_Sun_Exposed_Lower_leg', 'Small_Intestine_Terminal_Ileum', 'Spleen', 'Stomach',
    'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood'
]

# ------------------------- Load the gene annotation file -------------------------
# Annotation file format: gene_id  gene_symbol  chrom  start  end  gene_type
# No header; tab-delimited; gene IDs may include version numbers
gene_annot = {}
if os.path.exists(annotation_file):
    with open(annotation_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                gene_id, symbol, chrom, start, end, gtype = parts[:6]
                gene_annot[gene_id] = {'symbol': symbol, 'type': gtype}
    print(f"Loaded the gene annotation file: {len(gene_annot)} records")
else:
    print(f"Warning: annotation file {annotation_file} does not exist; the extra table will lack symbol and type information")
    gene_annot = {}

# ------------------------- Process each tissue -------------------------
for tissue in tissues:
    db_path = os.path.join(db_dir, f"mashr_{tissue}.db")
    if not os.path.exists(db_path):
        print(f"Skipping {tissue}: database does not exist")
        continue

    map_file = os.path.join(map_dir, f"{tissue}.leafcutter.phenotype_groups.txt.gz")
    if not os.path.exists(map_file):
        print(f"Skipping {tissue}: mapping file does not exist")
        continue

    # 1. Read the mapping file: third column (intron_id) -> second column (gene_id)
    mapping = {}
    with gzip.open(map_file, 'rt') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            # Column order: gtx_intron_id, gene_id, intron_id, cluster_id
            gtx_id, gene_id, intron_id, cluster_id = parts[0], parts[1], parts[2], parts[3]
            mapping[intron_id] = gene_id   # Map simplified intron IDs to gene IDs

    print(f"{tissue}: read {len(mapping)} mappings")

    # 2. Back up the database
    backup_path = db_path + ".bak"
    if not os.path.exists(backup_path):
        shutil.copy2(db_path, backup_path)
        print(f"Backed up to {backup_path}")

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # ==================== Process the weights table ====================
    try:
        df = pd.read_sql_query("SELECT * FROM weights", conn)
    except Exception as e:
        print(f"Failed to read the weights table: {e}")
        conn.close()
        continue

    if df.empty:
        print(f"{tissue}: weights table is empty; skipping")
        conn.close()
        continue

    # The first column contains simplified intron IDs
    intron_col = df.columns[0]
    print(f"Intron column name: {intron_col}")

    # Add the gene column
    df['gene'] = df[intron_col].map(mapping)

    # Remove rows that cannot be mapped
    before = len(df)
    df = df.dropna(subset=['gene'])
    after = len(df)
    print(f"Successfully mapped {after}/{before} rows")

    if after == 0:
        print(f"{tissue}: no records were mapped successfully; skipping")
        conn.close()
        continue

    # Remove duplicates after adding the gene column and dropping unmapped rows
    # Assume that the variant ID is in the second column (adjust to the actual column name)
    variant_col = df.columns[1]  # Adjust to the actual column name, such as 'rsid' or 'variant_id'
    
    # Group by gene and variant ID and average the weight column
    df = df.groupby(['gene', variant_col], as_index=False).agg({
    # For other columns, decide which value to retain if values may differ across rows
    'varID': 'first',
    # If ref_allele and eff_allele are present, assume their values match and retain the first
    'ref_allele': 'first',
    'eff_allele': 'first',
    # Average only the weight column here; retain the first value for other columns, or handle them as appropriate
    'weight': 'mean',
    # Additional column-handling rules can be added here
    })
    
    print(f"Rows after deduplication: {len(df)}")
    
    # Build the new table: replace the original intron column with gene and keep the other columns unchanged
    other_cols = [col for col in df.columns if col not in [intron_col, 'gene']]
    # Rename the 'gene_id' column to 'gene'
    # df_unique = df_unique.rename(columns={'gene_id': 'gene'})
    new_cols = ['gene'] + other_cols
    df_new = df[new_cols]
    
    # Replace the original table
    cursor.execute("ALTER TABLE weights RENAME TO weights_old")
    df_new.to_sql('weights', conn, index=False, if_exists='replace')
    cursor.execute("DROP TABLE weights_old")
    conn.commit()
    print(f"{tissue}: finished processing the weights table")
    
    # ==================== Process the extra table ====================
    try:
        df_extra = pd.read_sql_query("SELECT * FROM extra", conn)
    except Exception as e:
        print(f"Failed to read the extra table: {e}; skipping")
        conn.commit()
        conn.close()
        continue

    if df_extra.empty:
        print(f"{tissue}: extra table is empty; skipping")
    else:
        # Map intron IDs to gene IDs
        extra_intron_col = 'gene'
        df_extra['mapped_gene'] = df_extra[extra_intron_col].map(mapping)

        before_extra = len(df_extra)
        df_extra = df_extra.dropna(subset=['mapped_gene'])
        after_extra = len(df_extra)
        print(f"Successfully mapped {after_extra}/{before_extra} rows in the extra table")

        if after_extra == 0:
            print(f"{tissue}: no records in the extra table were mapped successfully; skipping")
        else:
            # Replace the original gene column
            df_extra['gene'] = df_extra['mapped_gene']
            df_extra.drop(columns=['mapped_gene'], inplace=True)

            # Get gene symbols and types from the annotation file using the gene column; version numbers do not need to be removed
            df_extra['genename'] = df_extra['gene'].map(lambda x: gene_annot.get(x, {}).get('symbol', x))
            df_extra['gene_type'] = df_extra['gene'].map(lambda x: gene_annot.get(x, {}).get('type', 'intron'))

            # Process performance-metric columns by converting them to numeric values
            perf_cols = ['n.snps.in.model', 'pred.perf.R2', 'pred.perf.pval', 'pred.perf.qval']
            for col in perf_cols:
                if col in df_extra.columns:
                    df_extra[col] = pd.to_numeric(df_extra[col], errors='coerce')

            # Merge records by retaining the one with the largest R-squared for each gene; if all R-squared values are NaN, retain the first record
            if 'pred.perf.R2' in df_extra.columns and df_extra['pred.perf.R2'].notna().any():
                df_extra_sorted = df_extra.sort_values('pred.perf.R2', ascending=False)
            else:
                df_extra_sorted = df_extra
            df_extra_unique = df_extra_sorted.groupby('gene').first().reset_index()

            # Reorder the columns to match the standard extra table
            final_cols = ['gene', 'genename', 'gene_type', 'n.snps.in.model', 'pred.perf.R2', 'pred.perf.pval', 'pred.perf.qval']
            final_cols = [c for c in final_cols if c in df_extra_unique.columns]
            df_extra_final = df_extra_unique[final_cols]

            # Replace the original table
            cursor.execute("ALTER TABLE extra RENAME TO extra_old")
            df_extra_final.to_sql('extra', conn, index=False, if_exists='replace')
            cursor.execute("DROP TABLE extra_old")
            print(f"{tissue}: finished processing the extra table (gene symbols and types added)")

    conn.commit()
    conn.close()
    print(f"Finished processing {tissue}")
