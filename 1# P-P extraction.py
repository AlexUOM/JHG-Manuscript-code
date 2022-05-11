import csv
from collections import OrderedDict

csv_file = open("BIOGRID-ORGANISM-Homo_sapiens-4.3.195.tab3.txt")
csv_reader = csv.reader(csv_file, delimiter="\t")
rows = [row for row in csv_reader]
shortlisted_rows = []
# extracting all row that have physical interraction (no RNA) only if both proteins are in humans
for row in rows:
    if "Affinity Capture-Luminescence" in row[11] or "Affinity Capture-MS" in row[11] or "Affinity Capture-Western" in row[11] or "Biochemical Activity" in row[11] or "Co-crystal Structure" in row[11] or "Co-fractionation" in row[11] or "Co-localization" in row[11] or "Co-purification" in row[11] or "Far Western" in row[11] or "FRET" in row[11] or "PCA" in row[11] or "Protein-peptide" in row[11] or "Proximity Label-MS" in row[11] or "Reconstituted Complex" in row[11] or "Two-hybrid" in row[11]:
        if "9606" in row[15] and "9606" in row[16]:
            shortlisted_rows.append(row)

# extact protein A and protein B, then remove duplicates and identical reversed pairs
proteins = []
for a in shortlisted_rows:
    proteins.append(tuple(a[7:9]))
proteins = list(set(tuple(sorted(unique)) for unique in proteins))
proteins = [tuple((Pa, Pb)) for Pa, Pb in proteins if Pa != Pb]
proteins = sorted(proteins, key=lambda element: (element[0], element[1]))

# write every protein pair as a row in a CSV file
f = open("Shortlisted proteins.csv", "w", newline='')
thewriter = csv.writer(f)
for i in proteins:
    thewriter.writerow([i[0], i[1]])

csv_file.close()
f.close()
