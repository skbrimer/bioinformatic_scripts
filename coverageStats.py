import glob
import csv

statsFiles = [i for i in glob.glob("*.tsv")]

for item in statsFiles:
    with open(item, "r") as depthFile:
        csv_reader = csv.reader(depthFile, delimiter="\t")
        zero_count = 0
        line_count = 0
        for row in csv_reader:
            if int(row[-1]) == 0:
                zero_count += 1
            else:
                line_count += 1
    coverage =round(1-(zero_count/(zero_count + line_count)),2)
    stats = item.split("_")[0]
    print(stats + "\t" + str(coverage))
