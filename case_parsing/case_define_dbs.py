import os
import csv
import json
import pandas as pd
import numpy as np
project_dir = os.path.dirname(os.getcwd())

#####################################################################################
# THIS Module is responsible for the generating of the Mutation Profiles of the Cases
#####################################################################################

def database_parsing_with_related_info_CCLE(maf):
    maf_dbs = pd.read_csv(maf,delimiter=",")
    del maf_dbs["ExAC_AF"],maf_dbs["Variant_annotation"],maf_dbs["CGA_WES_AC"],\
        maf_dbs["HC_AC"],maf_dbs["RD_AC"],maf_dbs["RNAseq_AC"],\
        maf_dbs["SangerWES_AC"],maf_dbs["WGS_AC"]
    case_ls = list(set(maf_dbs["DepMap_ID"].to_list())) #number of my cases is 1771get from CCLE

    for case in case_ls:
        contain_values = maf_dbs[maf_dbs["DepMap_ID"].str.contains(case)]

        defined_set_to_case = [[mutation_line["Hugo_Symbol"], mutation_line["Protein_Change"], mutation_line["Variant_Classification"], mutation_line["Variant_Type"],
                                mutation_line["isDeleterious"]]
                               for index,mutation_line in contain_values.iterrows()]

        file = open(project_dir+"/dbs_case_related_info/prepared_cases/"+case+".json","w")
        json.dump(defined_set_to_case,file)
        file.close()

def database_parsing_with_related_info_TCGA_GENIE(maf):
    maf_dbs = pd.read_csv(maf,delimiter="\t")
    case_ls = list(set(maf_dbs["Tumor_Sample_Barcode"].to_list()))
    print(len(case_ls))
    for case in case_ls:
        contain_values = maf_dbs[maf_dbs["Tumor_Sample_Barcode"].str.contains(case)]

        defined_set_to_case = [[mutation_line["Hugo_Symbol"], mutation_line["HGVSp_Short"], mutation_line["Variant_Classification"]]
                               for index,mutation_line in contain_values.iterrows()]

        file = open(project_dir+"/dbs_case_related_info/prepared_cases/"+case+".json","w")
        json.dump(defined_set_to_case,file)
        file.close()

if __name__ != '__main__':
    database_parsing_with_related_info_CCLE(project_dir+"/dbs_mutations/CCLE_mutations.csv")
    database_parsing_with_related_info_TCGA_GENIE(project_dir+"/dbs_mutations/TCGA_GENIE_MAF_File.tsv")
    #System has run one time to generate dbs.