import os
import csv
import plotly.express as px
import pandas as pd
import numpy as np

""" This Block Reads dataset as Analysis Results!"""
PM_Pan_Analysis_Result_Indexed=pd.read_csv("prepared_dbs/PM_PAN_Analysis_Results.csv")
PM_Pan_Analysis_Result_Indexed_Filtered = PM_Pan_Analysis_Result_Indexed[(PM_Pan_Analysis_Result_Indexed["Metastatic_Tissue_Span_Cohort"] > 30) | (PM_Pan_Analysis_Result_Indexed["Primary_Tissue_Span_Cohort"] > 30)]
PM_Pan_Analysis_Result=pd.read_csv("prepared_dbs/PM_PAN_Analysis_Results.csv").set_index("Hugo_Sym_Position_Implemented")
PM_Pan_Analysis_Result_Filtered = PM_Pan_Analysis_Result[(PM_Pan_Analysis_Result["Metastatic_Tissue_Span_Cohort"] > 30) | (PM_Pan_Analysis_Result["Primary_Tissue_Span_Cohort"] > 30)]

if __name__ != '__main__': ##Plotting by the characteristics on Scatter [NAIVE]
    fig = px.scatter(PM_Pan_Analysis_Result_Filtered)
    fig.show()

if __name__ != '__main__':
    fig = px.scatter_3d(PM_Pan_Analysis_Result_Filtered, x="Metastatic_Characteristics_Rate",y="Primary_Tissue_Span_Cohort", z="Metastatic_Tissue_Span_Cohort",text=PM_Pan_Analysis_Result_Filtered.index)
    fig.show()

if __name__ != '__main__':
    fig = px.scatter(PM_Pan_Analysis_Result_Filtered, x="Metastatic_Characteristics_Rate",y="Primary_Tissue_Span_Cohort",text=PM_Pan_Analysis_Result_Filtered.index)
    fig.show()

if __name__ == '__main__':
    overall_mutations_list=PM_Pan_Analysis_Result_Indexed_Filtered["Hugo_Sym_Position_Implemented"].to_list()
    overall_genes_list=[mutation.split("_")[0] for mutation in overall_mutations_list]
    for mutation in overall_genes_list:
        bar_df_list = list()
        objected_df=PM_Pan_Analysis_Result_Indexed_Filtered[PM_Pan_Analysis_Result_Indexed_Filtered["Hugo_Sym_Position_Implemented"].str.contains(mutation)]
        for index,row in objected_df.iterrows():
            bar_df_list.append([row["Hugo_Sym_Position_Implemented"],"Primary",row["Primary_Characteristics_Rate"]])
            bar_df_list.append([row["Hugo_Sym_Position_Implemented"], "Metastasic", row["Metastatic_Characteristics_Rate"]])
        bar_plot_df_pandas = pd.DataFrame(bar_df_list, columns=["Hugo_Sym_Position_Implemented","Status","Rate"])
        fig=px.bar(bar_plot_df_pandas, x="Status",y="Rate",color="Hugo_Sym_Position_Implemented")
        fig.show()

