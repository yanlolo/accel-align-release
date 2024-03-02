from __future__ import division
import pandas as pd, swifter, sys, csv

ref_path = '/home/yiqing/yan/input/simulate/10m-se-1/sv-10m-100-se-align.sam'
align_path = '/home/yiqing/yan/output/yiqing_hash.sam'
output_path = '/home/yiqing/yan/output//accalign.csv'
# ref_path = str(sys.argv[1])
# align_path = str(sys.argv[2])
# output_path = str(sys.argv[3])

with open(ref_path, 'r') as file:
    for ref_skip, line in enumerate(file):
        if 'NM' in line:
            rlen = len(line.split('\t')[9])
            match = rlen / 10
            break

print('Skipping the first ' + str(ref_skip) + ' lines in the reference file')
print('Read length: ' + str(rlen) + ', at least match' + str(match))

with open(align_path, 'r') as file:
    for align_skip, line in enumerate(file):
        if 'AS' in line and 'NM' in line:
            break

print('Skipping the first ' + str(align_skip) + ' lines in the aligned SAM')


colName = ["QNAME", "FLAG", "RNAME", "POS"]
colIndex = [0, 1, 2, 3]

def get_distance(x):
    return abs(x["POS_x"] - x["POS_y"])

ref = pd.read_csv(ref_path, sep = "\t", header = None, skiprows = ref_skip, usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
align = pd.read_csv(align_path, sep = "\t", header = None, skiprows = align_skip, usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)

align["QNAME"] = align["QNAME"].str.split('/').str[0]

## specially for minimap2, because it will produce multiple aligned places in sam file for one read
if len(ref) != len(align):
    align = align.groupby('QNAME').head(2).reset_index(drop=True)

aligned = align[align["RNAME"]!="*"]

merged = pd.merge(ref, aligned, how='inner', left_on=['QNAME', 'FLAG', 'RNAME'], right_on=['QNAME', 'FLAG', 'RNAME'])
merged['diff'] = merged.swifter.apply(get_distance, axis=1)

total = len(ref)
counts = merged.groupby("diff").size()
df = counts.to_frame("count").reset_index()
df["percent"] = df["count"] / total * 100

res = df[df["diff"] < match]
other = df[df["diff"] >= match]
match_cnt = res["count"].sum()
unmatch_cnt = other["count"].sum()
not_align_cnt = len(align[align["RNAME"]=="*"])

cols = ['diff', 'count', 'percent']
df.to_csv(output_path, index = False)


exact_match = df[df["diff"] == 0]
exact_match_cnt = exact_match["count"].sum()
print("Total strands: ", total)
print("Not aligned(%):", not_align_cnt/total* 100)
print("Correctly aligned(%):", match_cnt/total* 100)
print("Exact aligned(%):", exact_match_cnt/total* 100)
