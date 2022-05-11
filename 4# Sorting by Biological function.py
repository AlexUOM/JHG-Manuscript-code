import csv
import ast
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms.tree.mst import prim_mst_edges
import numpy as np

source_file = csv.reader(
    open("P-P and GO annotations and functions.csv"))
content = []
for line in source_file:
    frame = [line[0], line[1], ast.literal_eval(line[5])]
    content.append(frame)

# Parsing all up-to-date ontologies that are biological funtions
biological_Collection = {}
with open('go.obo', "r") as biological_functions_file:
    for line in biological_functions_file:
        if '[Term]' in line:
            frame = []
            frame.append(biological_functions_file.readline().strip())
            frame.append(biological_functions_file.readline().strip())
            frame.append(biological_functions_file.readline().strip())
            if "obsolete" in frame[1] or "biological_process" not in frame[2]:
                frame = []
            else:
                biological_Collection[frame[0][4:]] = frame[1][6:]

# Make a dictionary that has {biological function: "protein_A", "protein_B", ...}
collection = {}
networks = {}
for line in content:
    biological_process_network = nx.Graph()
    for term in line[2]:
        if term not in collection:
            collection[term] = [line[0], line[1]]
            networks[term] = [(line[0], line[1])]
        else:
            collection[term].append(line[0])
            collection[term].append(line[1])
            networks[term].append((line[0], line[1]))


# Extracting all subnetworks from each biological process
subnetworks = {}
for biological_process, proteins in networks.items():
    G = nx.Graph()
    G.add_edges_from(proteins)
    subgraphs = list(nx.connected_components(G))
    if len(subgraphs) == 1:
        subnetworks[biological_process] = subgraphs[0]
    else:
        count = 1
        for subgraph in subgraphs:
            dict_key = biological_process + "_" + str(count)
            subnetworks[dict_key] = subgraph
            count += 1


# If subnetworks share the same proteinsm consider them as one
temp = {}
for key, value in subnetworks.items():
    if tuple(value) in temp:
        temp[tuple(value)].append(key)
    else:
        temp[tuple(value)] = [key]
e = {}
for value, keys in temp.items():
    e[';'.join(keys)] = value
subnetworks = e

# remove dumplicated proteins from each key
collection = {a: list(set(b)) for a, b in collection.items()}

# remove keys with identical values (SAME ORDER)
d2 = {tuple(v): k for k, v in collection.items()}
collection = {v: list(k) for k, v in d2.items()}

# Writing whole biological processes into a .csv file
writer = csv.writer(
    open('Sorting by biological function.csv', 'w', newline=""))
for row in collection:
    writer.writerow([row, collection.get(row)])

# Writing biological process subnetworks into a .csv file
writer2 = csv.writer(
    open('Subnetworks of biological function.csv', 'w', newline=""))
for row in subnetworks:
    writer2.writerow([row, subnetworks.get(row)])

biological_functions_file.close()
