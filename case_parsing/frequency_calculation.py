import os
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
project_dir = os.path.dirname(os.getcwd())

def patient_status_define(file):
    from collections import Counter

    dbs = pd.read_csv(file,delimiter="\t")
    case, stage = dbs["SAMPLE_ID"].to_list(), dbs["SAMPLE_TYPE"].to_list()
    list_to_dict_for_patient_status = tuple(zip(case,stage))
    case_state_dbs_dict, status_overall_calculation =dict(), list()
    for x in list_to_dict_for_patient_status:
        case_state_dbs_dict[x[0]] = x[1]
        status_overall_calculation.append(x[1])

    status_overall_calculation_counted = dict(Counter(status_overall_calculation))

    return case_state_dbs_dict, status_overall_calculation_counted

def clean_out_MT(mutation_position,altered_residue=False):
    try:
        if "p." in mutation_position:
            mutation_position = mutation_position.split("p.")[1]

        if "fs" in mutation_position:
            return mutation_position.split("fs")[0]

        else:
            if altered_residue:
                return mutation_position
            else:
                return mutation_position[:-1]
    except:
        return "NaN"


def obtain_mutated_gene(main_maf, positions=True, altered=True):
    if "GENIE" in main_maf or "TCGA" in main_maf:
        maf_dbs = pd.read_csv(main_maf, delimiter="\t")
        mutation_column = list(maf_dbs["Hugo_Symbol"].to_list())
        WT_column = list(maf_dbs["WildTypeResidue"].to_list())
        positions_column = list(maf_dbs["ResidueNumber"].to_list())
        MT_column = list(maf_dbs["MutantResidue"].to_list())

        if positions and (altered is not True):
            mutation_w_residue = [element[0]+"_"+element[1]+str(element[2]) for element in tuple(zip(mutation_column,WT_column,positions_column))]
            return mutation_w_residue

        elif positions and altered:
            mutation_w_residue = [element[0] + "_" + element[1] + str(element[2]) + element[3] for element in
                                  tuple(zip(mutation_column, WT_column, positions_column,MT_column))]
            return mutation_w_residue

        else:
            return mutation_column

    if "CCLE" in main_maf:
        maf_dbs = pd.read_csv(main_maf, delimiter=",")
        mutation_column = list(maf_dbs["Hugo_Symbol"].to_list())

        if positions and (altered is not True):
            WT_position_column = [clean_out_MT(element, altered_residue=False) for element in
                                  maf_dbs["Protein_Change"].to_list()]
            mutation_w_residue = [element[0] + "_" + str(element[1]) for element in
                                  tuple(zip(mutation_column, WT_position_column))]
            return mutation_w_residue

        if positions and altered:
            WT_position_column = [clean_out_MT(element, altered_residue=True) for element in
                                  maf_dbs["Protein_Change"].to_list()]

            mutation_w_residue = [element[0] + "_" + str(element[1]) for element in
                                  tuple(zip(mutation_column, WT_position_column))]
            return mutation_w_residue

        else:
            return mutation_column


def find_frequency_of_mutations_in_cohort(mutation_ls):
    from collections import Counter

    x = dict(Counter(mutation_ls))
    sorted_dict = {k: v for k, v in sorted(x.items(), key=lambda item: item[1])}
    return sorted_dict


if __name__ != '__main__':
    mutation_ls = obtain_mutated_gene(project_dir+"/dbs_mutations/TCGA_GENIE_MAF_File.tsv", positions=False, altered=False)
    mutation_ls.extend(obtain_mutated_gene(project_dir+"/dbs_mutations/CCLE_mutations.csv", positions=False, altered=False))
