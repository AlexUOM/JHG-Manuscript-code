import csv
import ast

# Parsing all intersections from previous files
annotations_input_file = csv.reader(open("P-P and GO annotations.csv"))
output = []
intersections = []
for row in annotations_input_file:
    output.append(row)
    res = ast.literal_eval(row[4])
    intersections.append(res)

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

# For each GO term assign the respective biological function
c = 0
for GOterm in intersections:
    if GOterm:
        processes = []
        for k in GOterm:
            processes.append((biological_Collection.get(k)))
        output[c].append(processes)
    else:
        output[c].append("[]")
    c += 1


output = [i for i in output if len(i[4]) > 2]


# Write the output in a csv file in rows
writer = csv.writer(
    open('P-P and GO annotations and functions.csv', 'w', newline=""))
for row in output:
    writer.writerow([row[0], row[1], row[2], row[3],
                    row[4], row[5]])

biological_functions_file.close()
