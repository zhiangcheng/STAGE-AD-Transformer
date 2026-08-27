CLASS_NAMES = ["null", "shared", "female_biased", "male_biased", "sex_interaction"]
CLASS_TO_ID = {v:i for i,v in enumerate(CLASS_NAMES)}
ID_TO_CLASS = {i:v for v,i in CLASS_TO_ID.items()}

TISSUES = ["Brain cortex","Frontal cortex","Hippocampus","Whole blood","Monocyte","T cell","Artery","Liver","Adipose"]
CELL_TYPES = ["Microglia","Astrocyte","Oligodendrocyte","Endothelial","Neuron","OPC","Peripheral immune"]
MECHANISMS = ["Microglial lipid sensing","Astrocyte inflammatory reactivity","Vascular endothelial dysfunction","Oligodendrocyte myelin regulation","Peripheral immune regulation","APOE lipid tau axis","Neuronal synaptic vulnerability","X chromosome dosage or XCI"]

KNOWN_AD_GENES = ["APOE","TOMM40","APOC1","TREM2","BIN1","CLU","ABCA7","CR1","PICALM","CD2AP","INPP5D","MS4A6A","MS4A4A","PLCG2","ABI3","SORL1","EPHA1","HLA-DRB1","ACE","ADAM10","FERMT2","CASS4","PTK2B","SLC24A4","RIN3"]
FEMALE_GENES = ["APOE","CLU","INPP5D","MS4A6A","ESR1","ESR2","PGR","CYP19A1","IL1RAP"]
MALE_GENES = ["ACE","PTK2B","SLC24A4","VCAM1","COL4A1","MOG","MBP","NDUFA9"]
X_GENES = ["XIST","KDM6A","AR","TSPAN6","OPHN1","IL1RAPL1","PLP1","TREM2L_X"]

VARIANT_FEATURES = ["MAF","INFO","LD_score","distance_to_tss_kb","Z_total","Z_female","Z_male","Z_interaction","neglog10P_total","neglog10P_female","neglog10P_male","neglog10P_interaction","Beta_total","Beta_female","Beta_male","Beta_interaction","sex_delta_beta","abs_sex_delta_beta","PIP_total","PIP_female","PIP_male","PIP_interaction","coloc_brain_pph4","coloc_blood_pph4","coloc_immune_pph4","deepsea_delta","enformer_brain_delta","sei_regulatory_score","cadd_score","conservation_score","x_is_chrX","x_is_nonPAR","xci_escape","xci_subject","xci_variable"]
GENE_FEATURES = ["distance_kb","eqtl_z","sqtl_z","pqtl_z","twas_z","coloc_pph4","brain_expression","ad_de_z","constraint_score"]
TISSUE_FEATURES = ["eqtl_z","sqtl_z","coloc_pph4","enhancer","ldsc_enrichment","expression"]
CELL_FEATURES = ["marker_score","ad_de_score","enhancer_score","expression_specificity"]
