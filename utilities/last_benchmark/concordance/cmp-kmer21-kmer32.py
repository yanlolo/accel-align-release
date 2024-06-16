from __future__ import division

import pandas as pd
import swifter
import sys
import csv

workspace = '/home/yiqing/'
output_path = workspace + 'yan/output-bio150/'
aa21_path = output_path + 'accalign-kmer21-final.sam'
aa32_path = output_path + 'accalign-kmer32-final.sam'
bwa_path = output_path + 'bwa.sam'

aligner = 'accalign'
match = 15 #10%, 15
base_skip = 26


colName = ["QNAME", "FLAG", "RNAME", "POS"]
colIndex = [0, 1, 2, 3]

def get_distance(x):
    return x["RNAME"] == x["RNAME-bwa"] and (abs(x["POS-bwa"] - x["POS"]) < match or abs(x["POS-bowtie2"] - x["POS"]) < match)

def get_other_distance(x):
    return x["RNAME-other"] == x["RNAME-bwa"] and (abs(x["POS-bwa"] - x["POS-other"]) < match or abs(x["POS-bowtie2"] - x["POS-other"]) < match)


def get_skiprow(aligner):
    if aligner == "bwa" or aligner == "minimap2" or aligner == "mem2":
        return base_skip
    if aligner.startswith("acc") or aligner == "bowtie2" or aligner == "mrFast" or aligner == "fsva" or aligner == "subread" or aligner == "strobealign":
        return base_skip+1
    if aligner == "snap":
        return base_skip+2

def get_concordance(x):
    return x["RNAME-bwa"] == x["RNAME-bowtie2"] and abs(x["POS-bwa"] - x["POS-bowtie2"]) < match


diff32 = pd.read_csv(output_path + 'diff-32-final.csv')
diff21 = pd.read_csv(output_path + 'diff-21-final.csv')
diff32['FLAG'] = (diff32['FLAG'] & 0x40) | (diff32['FLAG'] & 0x80)
diff21['FLAG'] = (diff21['FLAG'] & 0x40) | (diff21['FLAG'] & 0x80)

merged = pd.merge(diff32, diff21, indicator=True, how='outer', left_on=['QNAME', 'FLAG'], right_on=['QNAME', 'FLAG'])
rows_only_in_diff21= merged[merged['_merge'] == 'left_only']

rows_only_in_diff21.to_csv(output_path + 'diff-21ok-21ko.csv', index = False)
distinct_values = rows_only_in_diff21['QNAME'].unique()
diff_force_align_qname = pd.DataFrame({'QNAME': distinct_values})

colName = ["QNAME", "FLAG", "RNAME", "POS", "CIGAR", "NM"]
colIndex = [0, 1, 2, 3, 5, 11]
bwa = pd.read_csv(bwa_path, sep = "\t", header = None, skiprows = get_skiprow("bwa"), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
diff_force_align_bwa = pd.merge(diff_force_align_qname, bwa, how='left', left_on=['QNAME'], right_on=['QNAME'])
diff_force_align_bwa = diff_force_align_bwa.groupby('QNAME').head(2).reset_index(drop=True)

res = diff_force_align_bwa[diff_force_align_bwa['CIGAR']!='*']
# 45210
# - 392 cigar is * --> bwa force align
# --> there is one mate has more than 4 errors --> kmer32 not find the correct one

aa32 = pd.read_csv(aa32_path, sep = "\t", header = None, skiprows = get_skiprow("accalign"), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
diff_aa32 = pd.merge(diff_force_align_qname, aa32, how='left', left_on=['QNAME'], right_on=['QNAME'])
diff_aa32 = diff_aa32.groupby('QNAME').head(2).reset_index(drop=True)



###check the cigar, NM
diff_force_align = pd.read_csv(output_path + 'diff-forcealign.csv')
distinct_values = diff_force_align['QNAME'].unique()
diff_force_align_qname = pd.DataFrame({'QNAME': distinct_values})
colName = ["QNAME", "FLAG", "RNAME", "POS", "CIGAR", "NM"]
colIndex = [0, 1, 2, 3, 5, 11]
bwa = pd.read_csv(bwa_path, sep = "\t", header = None, skiprows = get_skiprow("bwa"), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
diff_force_align_bwa = pd.merge(diff_force_align_qname, bwa, how='left', left_on=['QNAME'], right_on=['QNAME'])
diff_force_align_bwa = diff_force_align_bwa.groupby('QNAME').head(2).reset_index(drop=True)
res = diff_force_align_bwa[diff_force_align_bwa['CIGAR']!='*']
res['CIGAR'].value_counts()
res['NM'].value_counts()

bwa = pd.read_csv(bwa_path, sep = "\t", header = None, skiprows = get_skiprow("bwa"), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
bowtie2 = pd.read_csv(bowtie2_path, sep = "\t", header = None, skiprows = get_skiprow("bowtie2"), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
bwa = bwa.rename(columns={"QNAME": "QNAME-bwa", "FLAG": "FLAG-bwa", "RNAME": "RNAME-bwa", "POS": "POS-bwa", "MAPQ": "MAPQ-bwa"})
bowtie2 = bowtie2.rename(columns={"QNAME": "QNAME-bowtie2", "FLAG": "FLAG-bowtie2", "RNAME": "RNAME-bowtie2", "POS": "POS-bowtie2", "MAPQ": "MAPQ-bowtie2"})
bwa = bwa.groupby('QNAME-bwa').head(2).reset_index(drop=True)
bwa = bwa.sort_values(by = ['QNAME-bwa', 'FLAG-bwa']).reset_index()
bwa = bwa.drop(['index'], axis=1)
bowtie2 = bowtie2.sort_values(by = ['QNAME-bowtie2', 'FLAG-bowtie2']).reset_index()
bowtie2 = bowtie2.drop(['index'], axis=1)

ref = pd.concat([bwa, bowtie2], axis=1)
ref["same"] = ref.swifter.apply(get_concordance, axis=1)
concordance = ref[ref["same"]]

other = pd.read_csv(other_path, sep ="\t", header = None, skiprows = get_skiprow(aligner), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
other = other.rename(columns={"QNAME": "QNAME-other", "FLAG": "FLAG-other", "RNAME": "RNAME-other", "POS": "POS-other", "MAPQ": "MAPQ-other"})
other = other.sort_values(by = ['QNAME-other', 'FLAG-other']).reset_index()
other = other.drop(['index'], axis=1)

merged = pd.concat([other, concordance], axis=1)
merged['other-same'] = merged.swifter.apply(get_other_distance, axis=1)
merged[merged['other-same']]

aa = pd.read_csv(aa_path, sep ="\t", header = None, skiprows = get_skiprow(aligner), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
aa = aa.sort_values(by=['QNAME', 'FLAG']).reset_index()
aa = aa.drop(['index'], axis=1)

aa_merged = pd.concat([aa, merged], axis=1)
aa_merged['aa_same'] = aa_merged.swifter.apply(get_distance, axis=1)
aa_merged[aa_merged['aa_same']]

not_same = aa_merged[merged['other-same'] & (~aa_merged['aa_same'])]
not_same = not_same[['QNAME', 'FLAG', 'RNAME', 'POS', 'FLAG-bwa', 'RNAME-bwa', 'POS-bwa']]
not_same.to_csv(output_path + 'diff-21-final.csv', index = False)

not_same = pd.read_csv(output_path + 'diff-unalign21.csv')
not_same_bwa = pd.merge(not_same, bwa, how='left', left_on=['QNAME', 'FLAG-bwa'], right_on=['QNAME', 'FLAG'])
not_same_bwa['AS'] = not_same_bwa['AS']

########### analysis the reads not aligned by bwa, but AA aligned [Force align]--- START ###########
diff_force_align = pd.read_csv(output_path + 'diff-unalign21.csv')

colName = ["QNAME", "FLAG", "RNAME", "POS", "CIGAR"]
colIndex = [0, 1, 2, 3, 5]
bwa = pd.read_csv(bwa_path, sep = "\t", header = None, skiprows = get_skiprow("bwa"), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)

diff_force_align_bwa = pd.merge(diff_force_align, bwa, how='left', left_on=['QNAME', 'FLAG-bwa'], right_on=['QNAME', 'FLAG'])
res = diff_force_align_bwa[(diff_force_align_bwa['POS-bwa']!=0) & (diff_force_align_bwa['CIGAR']=='*')]
res.to_csv(output_path + 'diff-forcealign.csv', index = False)

###check the cigar, NM
diff_force_align = pd.read_csv(output_path + 'diff-forcealign.csv')
distinct_values = diff_force_align['QNAME'].unique()
diff_force_align_qname = pd.DataFrame({'QNAME': distinct_values})
colName = ["QNAME", "FLAG", "RNAME", "POS", "CIGAR", "NM"]
colIndex = [0, 1, 2, 3, 5, 11]
bwa = pd.read_csv(bwa_path, sep = "\t", header = None, skiprows = get_skiprow("bwa"), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
diff_force_align_bwa = pd.merge(diff_force_align_qname, bwa, how='left', left_on=['QNAME'], right_on=['QNAME'])
diff_force_align_bwa = diff_force_align_bwa.groupby('QNAME').head(2).reset_index(drop=True)
res = diff_force_align_bwa[diff_force_align_bwa['CIGAR']!='*']
res['CIGAR'].value_counts()
res['NM'].value_counts()

########### analysis the reads not aligned by bwa, but AA aligned [Force align]--- START ###########




########### analysis the reads not aligned by bwa, but AA aligned [UNAIGN some reads]--- START ###########
diff_not_align = not_same[not_same['POS-bwa']==0]
diff_not_align.to_csv(output_path + 'diff-notalign.csv', index = False)

colName = ["QNAME", "FLAG", "RNAME", "POS", "CIGAR", "NM"]
colIndex = [0, 1, 2, 3, 5, 11]
aa = pd.read_csv(aa_path, sep ="\t", header = None, skiprows = get_skiprow(aligner), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
diff_not_align_aa = pd.merge(diff_not_align, aa, how='left', left_on=['QNAME', 'FLAG'], right_on=['QNAME', 'FLAG'])

import re
# Function to extract 'M' values
def extract_m_value(column_value):
    match = re.search(r'(\d+)M', column_value)
    return int(match.group(1)) if match else None

diff_not_align_aa['M'] = diff_not_align_aa["CIGAR"].apply(lambda x: extract_m_value(x))
diff_not_align_aa['NM'] = diff_not_align_aa['NM'].str.split(':').str[2].astype('int')
diff_not_align_aa['matched'] = diff_not_align_aa['M'] - diff_not_align_aa['NM']
sorted_counts = diff_not_align_aa['matched'].value_counts().sort_values(ascending=False)
##[?M - NM, nb reads]

########### analysis the reads not aligned by bwa, but AA aligned [UNAIGN some reads]--- END ###########



print("Total strands: ", len(bwa))
print("BWA and Bowtie aligned same(concordance): ", len(concordance))
print("Other same as concordance: ", len(merged[merged['other-same']]))
print("AA same as concordance: ", len(aa_merged[aa_merged['aa_same']]))

print("Not aligned BWA(%):", len(bwa[bwa["RNAME-bwa"] == "*"])/len(bwa)* 100)
print("Not aligned Bowtie(%):", len(bowtie2[bowtie2["RNAME-bowtie2"] == "*"])/len(bowtie2)* 100)
print("Not aligned Other(%):", len(other[other["RNAME-other"] == "*"])/len(other)* 100)
print("Not aligned AA(%):", len(aa[aa["RNAME"] == "*"])/len(aa)* 100)


