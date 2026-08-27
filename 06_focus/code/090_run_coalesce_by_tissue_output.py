import os
import pandas
from glob import glob
import re
from pandas.errors import EmptyDataError

STUDY="final_metal_ad_kunkle_pgcalz_ukb"
OUT_DIR="../output/final_output/"

EG_FNAME="eqtl_focus_14tiss_final_metal_ad_kunkle_pgcalz_ukb_Chr1__PM__Brain_Amygdala.focus.tsv"


TISSUES = [
'Brain_Amygdala',
'Brain_Anterior_cingulate_cortex_BA24',
'Brain_Caudate_basal_ganglia',
'Brain_Cerebellar_Hemisphere',
'Brain_Cerebellum',
'Brain_Cortex',
'Brain_Frontal_Cortex_BA9',
'Brain_Hippocampus',
'Brain_Hypothalamus',
'Brain_Nucleus_accumbens_basal_ganglia',
'Brain_Putamen_basal_ganglia',
'Brain_Spinal_cord_cervical_c-1',
'Brain_Substantia_nigra',
'Nerve_Tibial'
]

def main(eqtl_or_sqtl, tiss_tag="_14tiss_"):
    if eqtl_or_sqtl == "eqtl":
        pass
    elif eqtl_or_sqtl == "sqtl":
        pass
    else:
        raise Exception(f"Variable <eqtl_or_sqtl> should be one of 'eqtl' or 'sqtl', not: {eqtl_or_sqtl}")

    F=f"../output/final_output/focus_{eqtl_or_sqtl}{tiss_tag}{STUDY}_by_tissue.focus.tsv"
    d=[]
    count = 0
    for i in range(1,23): # (1-22)
        for cur_tiss in TISSUES:
            search_name = f"../output/by_tissue/{eqtl_or_sqtl}_focus{tiss_tag}{STUDY}_Chr{i}__PM__{cur_tiss}.focus.tsv" 
            #if i == 1: print(search_name)
            files_for_chrom = glob(search_name)
            if len(files_for_chrom) == 0:
                print("{}: {} missing for chr {}".format(eqtl_or_sqtl, cur_tiss, i))
            for result_file in files_for_chrom:
                fname_only = os.path.basename(result_file)
                # extract the tissue name after "__PM__" in the file
                tissue_target = re.search("__PM__(.*).focus.tsv", fname_only)[1]
                count += 1
                #print(result_file)
                try:
                	d_ = pandas.read_table(result_file)
                except EmptyDataError:
                	print(f"empty file skipped: {result_file}")
                	continue
                d_["tissue_target"] = tissue_target
                d.append(d_)
                
    print("")
    print(F)
    print(count)
    print("Expect count of {} if there are 14 tissues".format(242))
    if len(d) == 0:
    	print(f"No valid files found for {eqtl_or_sqtl}")
    	return
    final = pandas.concat(d).drop_duplicates()
    os.makedirs(os.path.dirname(F), exist_ok=True)
    final.to_csv(F, sep="\t", index=False)


if __name__ == "__main__":
   main("eqtl")
   main("sqtl")
   print("File one executed when ran directly")
else:
   print("File one imported")
