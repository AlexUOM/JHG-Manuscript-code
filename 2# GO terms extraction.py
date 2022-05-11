import csv

csv_reader = csv.reader(open("goa_human.gaf"), delimiter="\t")
shortlisted_rows = []
for i in enumerate(csv_reader):
    if i[0] >= 40:
        if i[1][8] == "P" and i[1][6] != "IEA" and i[1][6] != "NAS" and i[1][6] != "TAS" and i[1][6] != "IC" and i[1][6] != "ND":
            shortlisted_rows.append(i[1])

reader = csv.reader(open('Shortlisted proteins.csv', 'r'))
Prots = []
P_pairs = []
for a, b in reader:
    P_pairs.append((a, b))
    Prots.append(a)
    Prots.append(b)
Prots = list(dict.fromkeys(sorted(Prots)))
GOcollection = {}
for candidate in shortlisted_rows:
    if candidate[2] in Prots:
        GOcollection.setdefault(candidate[2], [])
        GOcollection[candidate[2]].append(candidate[4])
GOcollection = {a: list(set(b)) for a, b in GOcollection.items()}


# intersections
intersections = []
output = []
for Pa, Pb in P_pairs:
    if Pa in GOcollection and Pb in GOcollection:
        for mutual1 in GOcollection.get(Pa):
            for mutual2 in GOcollection.get(Pb):
                if mutual1 == mutual2:
                    intersections.append(mutual1)
        row = (Pa, Pb, GOcollection.get(Pa),
               GOcollection.get(Pb), intersections)
        intersections = []
        output.append(row)
print(output[:10])

writer = csv.writer(open('P-P and GO annotations.csv', 'w', newline=""))
for row in output:
    writer.writerow([row[0], row[1], row[2], row[3], row[4]])
