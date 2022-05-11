import pandas as pd
import numpy as np
import csv
import ast
import networkx as nx

variant_df = pd.read_excel("TOF_WES_filtered_variants.xlsx")
protein_df = pd.read_csv("Subnetworks of biological function.csv")
all_prots = []
for i in protein_df.iterrows():
    all_prots.extend(ast.literal_eval(i[1]["Proteins"]))
all_prots = list(sorted(set(all_prots)))

variant_df = variant_df.loc[variant_df.impact_severity == "HIGH"]
variant_df = variant_df.loc[variant_df.type != "indel"]
variant_df = variant_df[["start", "gene",  "type",
                         "impact_severity", "variant_samples", "het_samples", "hom_alt_samples", "case_count"]]


candidates = {}
patients = {}
for index, row in variant_df.iterrows():
    if row["gene"] in all_prots and row["gene"] not in candidates:
        occurence = 1
        position = row["start"]
        candidates[row["gene"]] = occurence
        if type(row["hom_alt_samples"]) == str:
            for case in row["hom_alt_samples"].split(","):
                if case not in patients:
                    patients[case] = [row["gene"]]
                else:
                    patients[case].append(row["gene"])
        if type(row["het_samples"]) == str:
            for case2 in row["het_samples"].split(","):
                if case2 not in patients:
                    patients[case2] = [row["gene"]]
                else:
                    patients[case2].append(row["gene"])
    elif row["gene"] in all_prots and row["gene"] in candidates:
        occurence += 1
        candidates[row["gene"]] = occurence
        if type(row["hom_alt_samples"]) == str:
            for case in row["hom_alt_samples"].split(","):
                if case not in patients:
                    patients[case] = [row["gene"]]
                else:
                    patients[case].append(
                        row["gene"])
        if type(row["het_samples"]) == str:
            for case2 in row["het_samples"].split(","):
                if case2 not in patients:
                    patients[case2] = [row["gene"]]
                else:
                    patients[case2].append(
                        row["gene"])
patients = {patient: list(set(case)) for (patient, case) in patients.items()}


# Dict {Protein:[patients]}
patients_rev = {}
for index, row in variant_df.iterrows():
    if row["gene"] in all_prots and row["gene"] not in patients_rev:
        patients_rev[row["gene"]] = row["variant_samples"].split(",")
    elif row["gene"] in all_prots and row["gene"] in patients_rev:
        patients_rev[row["gene"]].extend(row["variant_samples"].split(","))
    elif row["gene"] not in all_prots and row["gene"] not in patients_rev:
        patients_rev[row["gene"]] = []
patients_rev = {patient: list(set(case))
                for (patient, case) in patients_rev.items()}


output = []
writer2 = csv.writer(
    open('Variants per biological function_COUNT.csv', 'w', newline=""))
for index, row in protein_df.iterrows():
    intersection = set(candidates) & set(ast.literal_eval(row["Proteins"]))
    if len(intersection) > 0:
        count = []
        for protein in list(intersection):
            count.extend(patients_rev[protein])
        count = set(count)
        output.append([row["Biological function"],
                      ast.literal_eval(row["Proteins"]), len(count)])
writer2.writerow(["Biological function", "Proteins", "Cases"])
for i in output:
    writer2.writerow(i)


# Write it into a csv file
output = []
writer = csv.writer(
    open('Variants per biological function.csv', 'w', newline=""))
for index, row in protein_df.iterrows():
    intersection = set(candidates) & set(ast.literal_eval(row["Proteins"]))
    if len(intersection) == 1:
        output.append([row["Biological function"], str(intersection).strip("{'}"),
                       candidates.get(str(intersection).strip("{'}"))])
    elif len(intersection) > 1:
        result = [row["Biological function"], list(intersection)[
            0].strip("{'}"), str(candidates.get(list(intersection)[0].strip("{'}")))]
        for protein in list(intersection)[1:]:
            result.append(protein)
            result.append(str(candidates.get(protein)))
        output.append(result)

output = sorted(output, key=len, reverse=True)
for i in output:
    writer.writerow(i)
