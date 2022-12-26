import os
import math
import warnings
import pandas as pd

warnings.filterwarnings("ignore")
project_dir = os.path.dirname(os.getcwd())


def rating_characteristic_of_case_lise_from_given_mutation(case_ls, case_status_dict):
    from collections import Counter

    characteristic_ls = list()
    for status_of_case in case_ls:
        if "TCGA" not in status_of_case:
            characteristic_ls.append(case_status_dict[status_of_case])
        else:
            characteristic_ls.append("Primary")
    characteristic_dict = dict(Counter(characteristic_ls))

    return characteristic_dict


def calculate_characteristic_rate_from_chr_dict(chr_dict):
    try:
        if "Primary" not in chr_dict.keys() and "Metastasis" in chr_dict.keys():
            chr_rate_primary, primary_count = 0, 0
            chr_rate_metastatic, metastasis_count  = 1, chr_dict["Metastasis"]
        elif "Metastasis" not in chr_dict.keys() and "Primary" in chr_dict.keys():
            chr_rate_metastatic, metastasis_count = 0, 0
            chr_rate_primary, primary_count = 1, chr_dict["Primary"]

        elif "Metastasis" not in chr_dict.keys() and "Primary" not in chr_dict.keys():
            chr_rate_primary, primary_count = 0, 0
            chr_rate_metastatic, metastasis_count = 0, 0
        else:
            chr_rate_primary, primary_count = characteristic_dict["Primary"] / (
                    characteristic_dict["Primary"] + characteristic_dict["Metastasis"]), characteristic_dict["Primary"]
            chr_rate_metastatic, metastasis_count = characteristic_dict["Metastasis"] / (
                    characteristic_dict["Primary"] + characteristic_dict["Metastasis"]), characteristic_dict["Metastasis"]

        return chr_rate_primary, chr_rate_metastatic, primary_count ,metastasis_count
    except:
        return 0


def calculate_overall_span_of_cohort(chr_dict, counted_status_overall_info):
    try:
        if "Primary" not in chr_dict.keys():
            primary_overall_span = 0
            try:
                metastasis_overall_span = math.pow(10, 5) * chr_dict["Metastasis"] / counted_status_overall_info["Metastasis"]
            except:
                metastasis_overall_span = 0

        elif "Metastasis" not in chr_dict.keys():
            metastasis_overall_span = 0
            try:
                primary_overall_span = math.pow(10, 5) * chr_dict["Primary"] / counted_status_overall_info["Primary"]
            except:
                primary_overall_span = 0

        else:
            primary_overall_span = math.pow(10, 5) * chr_dict["Primary"] / counted_status_overall_info["Primary"]
            metastasis_overall_span = math.pow(10, 5) * chr_dict["Metastasis"] / counted_status_overall_info["Metastasis"]
        return primary_overall_span, metastasis_overall_span
    except:
        return 0

def record_pandas_dbs_and_write(PM_Chr_Span_DBS):
    characteristics_and_overall_span_dbs = pd.DataFrame(PM_Chr_Span_DBS, columns=["Hugo_Sym_Position_Implemented","Primary_Count_of_MT" ,"Metastatic_Count_of_MT","Primary_Characteristics_Rate", "Metastatic_Characteristics_Rate","Primary_Tissue_Span_Cohort", "Metastatic_Tissue_Span_Cohort"])
    characteristics_and_overall_span_dbs.to_csv(project_dir+"/prepared_dbs/PM_PAN_Analysis_Results.csv",index=False)

############################################################################################################
from case_parsing.frequency_calculation import *

# - Functions were defined in the frequency_calculation scripts are going to use in this script to find
# relationships of genes being metastatic/primary
# - All required functions from previous scripts being present in
# the __main__ stage as a pipeline, this script is continue of it.
# - At the end of the script all primary/metastatic
# analysis are going to be done.
# ###############################################################################################################

if __name__ == '__main__':  ##In this main block evaluation just be done for the GENIE and TCGA databases !!!
    # 1 - In this part mutation_list is going to be created and patient dictionary is requested
    mutation_ls = obtain_mutated_gene(project_dir + "/dbs_mutations/"
                                                    "TCGA_GENIE_MAF_File.tsv", positions=True,
                                      altered=False)
    mutations_Counted_sorted_by_Freq = find_frequency_of_mutations_in_cohort(mutation_ls)

    case_status_dict, counted_status_overall_info = patient_status_define(project_dir + "/dbs_case_related_info/"
                                                                                        "GENIE_Primary_Metastasis_Info_new.txt")

    maf_dbs = pd.read_csv(project_dir + "/dbs_mutations/TCGA_GENIE_MAF_File.tsv", delimiter="\t")

    # 2 - This is the part of scoring mutations by its P/M characterization.

    PM_Chr_Span_DBS = list()
    # ----------MUTATION ITERATION---------- #
    for mutation_given in mutations_Counted_sorted_by_Freq.keys().__reversed__():  # The main reason why list is reversed is most freqs are being end of the list.
        ### DENOTES THE MUTATION BELONG BLOCK FOR AN ITERATION ###
        print("Objected Mutation : " + str(mutation_given))

        Hugo_symbol_given, Position = mutation_given.split("_")[0], mutation_given.split("_")[1]
        mutation_relevant_dbs = maf_dbs[(maf_dbs[
                                             "Hugo_Symbol"] == Hugo_symbol_given)  # FILTER the overall databases to find the relavent patients with mutation
                                        & (maf_dbs["WT_Residue"] == Position)]

        related_patients_with_given_mutation = mutation_relevant_dbs[
            "Tumor_Sample_Barcode"].to_list()  # Relevant patients are presence at this list.

        characteristic_dict = rating_characteristic_of_case_lise_from_given_mutation(
            related_patients_with_given_mutation, case_status_dict)

        # And Giving the info of given between occuring mutations how much percent metastatic/primary in line below
        primary_rate, metastasis_rate, primary_count, metastasis_count = calculate_characteristic_rate_from_chr_dict(characteristic_dict)
        #print("The chr. rate of Primary | Metastasis  ---> " + str(primary_rate) + " | " + str(metastasis_rate))

        # And how much they were spanned to plot
        span_rate_primary, span_rate_metastatic = calculate_overall_span_of_cohort(characteristic_dict,
                                                                                   counted_status_overall_info)
        #print("The chr. rate of Primary | Metastasis  ---> " + str(span_rate_primary) + " 10^-3 | " + str(span_rate_metastatic) + " 10^-3")
        #print("############################################################################################################")

        ################### THIS IS THE PART OF EVERYTHING RECORDING A LIST ############################
        PM_Chr_Span_DBS.append([mutation_given, round(primary_rate,3), round(metastasis_rate,3), span_rate_primary, span_rate_metastatic, primary_count, metastasis_count])

    ## END OF THE PIPELINE COLUMNS SEND TO PANDAS TO GEN DATAFRAME !!!
    ## HEADLINES ---->
    # "Hugo_Sym_Position_Implemented",
    # "Primary_Count_of_MT" , "Metastatic_Count_of_MT",
    # "Primary_Characteristics_Rate", "Metastatic_Characteristics_Rate",
    # "Primary_Tissue_Span_Cohort", "Metastatic_Tissue_Span_Cohort"])

    record_pandas_dbs_and_write(PM_Chr_Span_DBS)
    #------------------------------------------------------ END ------------------------------------------------------#