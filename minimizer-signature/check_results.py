import sys

res_files = [
    ("./acc_results/acc_test_" + sys.argv[1] + "_B1427", "B.1.427"),
    ("./acc_results/acc_test_" + sys.argv[1] + "_B1526", "B.1.526"),
    ("./acc_results/acc_test_" + sys.argv[1] + "_P1", "P.1"),
]


acc_vals = {}
for filepath, key in res_files:
    with open(filepath) as fin:
        for line in fin.readlines():
            line = line.strip().split()
            name = line[0]
            acc  = float(line[1])

            if name not in acc_vals:
                acc_vals[name] = {}
            acc_vals[name][key] = acc


correct = 0
count = 0
for line, vals in acc_vals.items():
    name = line.split("_")[0][1:]
    if name == max(vals, key=vals.get):
        correct += 1
    count += 1

print("Accuracy:", correct*1.0/count)