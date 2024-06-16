## why minimap2 mapped, AA not?

from __future__ import division
import pandas as pd, swifter, sys, csv

ref_path = "/home/yiqing/yan/input/simulate/1m/sv-1m-100-pe-align.sam"
aa_path = "/home/yiqing/yan/output/accalign-pair.sam"
other_path = "/home/yiqing/yan/output/accalign.sam"
# ref_path = str(sys.argv[1])
# align_path = str(sys.argv[2])
# output_path = str(sys.argv[3])

with open(ref_path, 'r') as file:
    for ref_skip, line in enumerate(file):
        if 'NM' in line:
            rlen = len(line.split('\t')[9])
            match = rlen / 10
            break

print('Skipping the header: first ' + str(ref_skip) + ' lines in the reference file')
print('Read length: ' + str(rlen) + 'bp, at least match: ' + str(match))

with open(aa_path, 'r') as file:
    for aa_skip, line in enumerate(file):
        if 'AS' in line and 'NM' in line:
            break

with open(other_path, 'r') as file:
    for other_skip, line in enumerate(file):
        if 'AS' in line and 'NM' in line:
            break

def get_aa_unmap_other_map(x):
    return (x["FLAG-aa"] != x["FLAG"] or (abs(x["POS-aa"] - x["POS"]) > match))

def get_distance(x):
    return x["RNAME"] == x["RNAME-o"] and x["POS"] == x["POS-o"] and (x["RNAME"] != x["RNAME-aa"] or x["POS"] != x["POS-aa"])

colName = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ"]
colIndex = [0, 1, 2, 3, 4]

aa = pd.read_csv(aa_path, sep = "\t", header = None, skiprows = aa_skip, usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
aa['qname'] = aa['QNAME'].str.split('.').str[1].astype('int')
aa = aa.sort_values(by = ['qname', 'FLAG']).reset_index()
aa = aa.drop(['qname', 'index'], axis=1)

other = pd.read_csv(other_path, sep = "\t", header = None, skiprows = other_skip, usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
other['qname'] = other['QNAME'].str.split('.').str[1].astype('int')
other = other.sort_values(by = ['qname', 'FLAG']).reset_index()
other = other.drop(['qname', 'index'], axis=1)

ref = pd.read_csv(ref_path, sep = "\t", header = None, skiprows = ref_skip, usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
ref['qname'] = ref['QNAME'].str.split('.').str[1].astype('int')
ref = ref.sort_values(by = ['qname', 'FLAG']).reset_index()
ref = ref.drop(['qname', 'index'], axis=1)

aa = aa.rename(columns={"QNAME": "QNAME-aa", "FLAG": "FLAG-aa", "RNAME": "RNAME-aa", "POS": "POS-aa", "MAPQ": "MAPQ-aa"})
other = other.rename(columns={"QNAME": "QNAME-o", "FLAG": "FLAG-o", "RNAME": "RNAME-o", "POS": "POS-o", "MAPQ": "MAPQ-o"})

# other = other.groupby('QNAME').head(2).reset_index(drop=True)

merged = pd.concat([other, ref, aa], axis=1)
merged['diff'] = merged.swifter.apply(get_distance, axis=1)
merged[merged['diff']]


pd.set_option('display.float_format', '{:.0f}'.format)
diff = merged[merged['diff']]
haha=diff[['QNAME-o', 'FLAG-o', "RNAME-o", "POS-o", 'QNAME-aa', 'FLAG-aa', "RNAME-aa", "POS-aa"]]

haha[haha['RNAME-aa']=='*']
haha = haha[haha['RNAME-aa']!='*']
chr_pos = {'chr1': 0, 'chr2': 249250621, 'chr3': 492449994, 'chr4': 690472424, 'chr5': 881626700,
           'chr6': 1062541960, 'chr7': 1233657027, 'chr8': 1392795690, 'chr9': 1539159712, 'chr10': 1680373143,
           'chr11': 1815907890, 'chr12': 1950914406, 'chr13': 2084766301, 'chr14': 2199936179, 'chr15': 2307285719,
           'chr16': 2409817111, 'chr17': 2500171864, 'chr18': 2581367074, 'chr19': 2659444322, 'chr20': 2718573305,
           'chr21': 2781598825, 'chr22': 2829728720, 'chrX': 2881033286, 'chrY': 3036303846, 'chrM': 3095677412}
haha["p_o"] =  haha["POS-o"] + haha["RNAME-o"].map(chr_pos)
haha["p_aa"] =  haha["POS-aa"] + haha["RNAME-aa"].map(chr_pos)
haha['cmp_pos']= haha["p_o"]>haha["p_aa"]

len(haha[(haha["FLAG-o"]==haha["FLAG-aa"])])
len(haha[(haha["FLAG-o"]!=haha["FLAG-aa"])])

len(haha[(haha["FLAG-o"]==haha["FLAG-aa"]) & haha['cmp_pos']])
len(haha[(haha["FLAG-o"]==haha["FLAG-aa"]) & (~haha['cmp_pos'])])

len(haha[(haha["FLAG-o"]!=haha["FLAG-aa"]) & (haha["FLAG-o"]==0)])
len(haha[(haha["FLAG-o"]!=haha["FLAG-aa"]) & (haha["FLAG-o"]==16)])

len(haha[(haha["FLAG-o"]!=haha["FLAG-aa"]) & (haha["FLAG-o"]==0) & haha['cmp_pos']])
len(haha[(haha["FLAG-o"]!=haha["FLAG-aa"]) & (haha["FLAG-o"]==0) & ~haha['cmp_pos']])
len(haha[(haha["FLAG-o"]!=haha["FLAG-aa"]) & (haha["FLAG-o"]==16) & haha['cmp_pos']])
len(haha[(haha["FLAG-o"]!=haha["FLAG-aa"]) & (haha["FLAG-o"]==16) & ~haha['cmp_pos']])

