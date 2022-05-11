import csv
import ast
import numpy as np
import pandas as pd
import statsmodels.stats.api as sms

# Reading all biological functions and proteins and create a Dict out of it
biological_functions_file = pd.read_csv(
    "Variants per biological function_COUNT.csv")
protein_collection = {}
for index, line in biological_functions_file.iterrows():
    frame = []
    for protein in ast.literal_eval(line[1]):
        frame.append(protein)
    protein_collection[line[0]] = frame
foo = sorted(protein_collection.values(), key=len, reverse=True)
for index, i in enumerate(foo):
    foo[index] = len(i)
foo = list(sorted(set(foo)))

numbering = []
dict_of_models = {}
for size in foo:
    for j in protein_collection.items():
        if size == len(j[1]):
            dict_of_models[j[0]] = j[1]
            break


prots_df = pd.read_csv("Subnetworks of biological function.csv")
original_proteins = []
for index, line in prots_df.iterrows():
    original_proteins.extend(ast.literal_eval(line["Proteins"]))
original_proteins = list(sorted(set(original_proteins)))


#  Reading variants file, extracting High impact variants and trimming the dataset
variant_df = pd.read_excel("TOF_WES_filtered_variants.xlsx")
variant_df = variant_df.loc[variant_df.impact_severity == "HIGH"]
variant_df = variant_df.loc[variant_df.type != "indel"]
variant_df = variant_df[["start", "gene",  "type",
                         "impact_severity", "variant_samples", "case_count"]]


# Function for permutation and mapping
def permutate_map(reference):
    permutated_proteins = (np.random.permutation(reference))
    mapped_proteins = {}
    for index, protein in enumerate(reference):
        mapped_proteins[protein] = permutated_proteins[index]
    return mapped_proteins


# Dict {Protein:[cases]}
cases = {}
for index, row in variant_df.iterrows():
    if row["gene"] in original_proteins and row["gene"] not in cases:
        cases[row["gene"]] = row["variant_samples"].split(",")
    elif row["gene"] in original_proteins and row["gene"] in cases:
        cases[row["gene"]].extend(row["variant_samples"].split(","))
cases = {patient: list(set(case))
         for (patient, case) in cases.items()}
for remainder in original_proteins:
    if remainder not in cases:
        cases[remainder] = []


# Function for changing protein order
def changing_order(BOs):
    new_order = {"Cases #0": []}
    new_map = permutate_map(original_proteins)
    for biological_function, proteins in BOs.items():
        counting = []
        for protein in proteins:
            counting.extend(cases[new_map[protein]])
        counting = set(counting)
        counting = len(counting)
        new_order["Cases #0"].append(counting)
    return new_order


# Run n permutation and write each result in a new column, then save it into a .csv
BOs = {"Biological function": list(
    (biological_functions_file["Biological function"]))}
results_by_chance = pd.DataFrame(BOs)
results_by_chance["Proteins"] = biological_functions_file["Proteins"]
results_by_chance["Observed cases"] = biological_functions_file["Cases"]

results_by_chance = pd.DataFrame(changing_order(protein_collection))


for size in dict_of_models.items():
    BOs = {"Biological function": [size[0]]}
    results_by_chance = pd.DataFrame(BOs)
    results_by_chance["Proteins"] = [size[1]]
    n_cases = biological_functions_file.loc[biological_functions_file["Biological function"] == size[0]]
    results_by_chance["Observed cases"] = list(n_cases["Cases"])

    for i in range(10000):
        series = "Cases #" + str(i+1)
        results_by_chance[series] = pd.DataFrame(
            changing_order(dict_of_models))["Cases #0"]
    model_string = "Variants per biological function_BY_CHANCE_" + \
        str(len(size[1])) + ".csv"
    results_by_chance.to_csv(
        model_string, index=False)
    print(len(size[1]))


result = []
file_name = "Variants per biological function_BY_CHANCE_.csv"
for model_size in foo:
    file_name = "Variants per biological function_BY_CHANCE_" + \
        str(model_size) + ".csv"
    model_file = pd.read_csv(file_name)
    permutations = model_file.iloc[0][3:].quantile(
        [((0.01/len(foo))/2), 1-((0.01/len(foo))/2)])
    conf_int = list(permutations)
    frame2 = []
    for index, line in biological_functions_file.iterrows():
        if model_size == len(ast.literal_eval(line[1])):
            if line[2] < conf_int[0] or line[2] > conf_int[1]:
                frame2.append(line[0])
                frame2.append(line[1])
                frame2.append(line[2])
                frame2.append(conf_int)
                frame2.append("SIGNIFICANT")
            else:
                frame2.append(line[0])
                frame2.append(line[1])
                frame2.append(line[2])
                frame2.append(conf_int)
                frame2.append("NOT SIGNIFICANT")
            result.append(frame2)
            frame2 = []
    print(model_size)
writer = csv.writer(
    open('Significance.csv', 'w', newline=""))
for i in result:
    writer.writerow(i)
