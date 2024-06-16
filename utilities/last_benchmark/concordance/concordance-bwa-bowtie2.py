# ref_path = str(sys.argv[1])
# align_path = str(sys.argv[2])
# output_path = str(sys.argv[3])
# aligner = str(sys.argv[4])
# match = int(sys.argv[5])
# base_skip = int(sys.argv[6])

from __future__ import division
import pandas as pd, swifter, sys, csv

workspace = '/home/yiqing/'
output_path = workspace + 'yan/output-bio250/'
bwa_path = output_path + 'bwa.sam'
bowtie2_path = output_path + 'bowtie2.sam'

colName = ["QNAME", "FLAG", "RNAME", "POS"]
colIndex = [0, 1, 2, 3]

with open(bwa_path, 'r') as file:
    for  line in file:
        if 'NM' in line:
            match = len(line.split('\t')[9]) * 0.9
            break


def get_skiprow(align_path):
    align_skip = 0
    with open(align_path, 'r') as file:
        for align_skip, line in enumerate(file):
            if 'AS' in line and 'NM' in line:
                break
        return align_skip


def get_distance(x):
    return x["RNAME"] == x["RNAME-bwa"] and (abs(x["POS-bwa"] - x["POS"]) < match or abs(x["POS-bowtie2"] - x["POS"]) < match)

def get_other_distance(x):
    return x["RNAME-other"] == x["RNAME-bwa"] and (abs(x["POS-bwa"] - x["POS-other"]) < match or abs(x["POS-bowtie2"] - x["POS-other"]) < match)

def get_concordance(x):
    return x["RNAME-bwa"] == x["RNAME-bowtie2"] and abs(x["POS-bwa"] - x["POS-bowtie2"]) < match

bwa = pd.read_csv(bwa_path, sep = "\t", header = None, skiprows = get_skiprow(bwa_path), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
bwa = bwa.rename(columns={"QNAME": "QNAME-bwa", "FLAG": "FLAG-bwa", "RNAME": "RNAME-bwa", "POS": "POS-bwa", "MAPQ": "MAPQ-bwa"})
bwa = bwa.groupby('QNAME-bwa').head(2).reset_index(drop=True)
bwa = bwa.sort_values(by = ['QNAME-bwa', 'FLAG-bwa']).reset_index()
bwa = bwa.drop(['index'], axis=1)

bowtie2 = pd.read_csv(bowtie2_path, sep = "\t", header = None, skiprows = get_skiprow(bowtie2_path), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
bowtie2 = bowtie2.rename(columns={"QNAME": "QNAME-bowtie2", "FLAG": "FLAG-bowtie2", "RNAME": "RNAME-bowtie2", "POS": "POS-bowtie2", "MAPQ": "MAPQ-bowtie2"})
bowtie2 = bowtie2.groupby('QNAME-bowtie2').head(2).reset_index(drop=True)
bowtie2 = bowtie2.sort_values(by = ['QNAME-bowtie2', 'FLAG-bowtie2']).reset_index()
bowtie2 = bowtie2.drop(['index'], axis=1)

ref = pd.concat([bwa, bowtie2], axis=1)
ref["same"] = ref.swifter.apply(get_concordance, axis=1)

ref.to_csv(output_path + 'ref.sam')

###########
output_path = workspace + 'yan/output-bio250/'
other_path = output_path + 'accalign-cano-32-250.sam'
ref = pd.read_csv(output_path + 'ref.sam')
other = pd.read_csv(other_path, sep ="\t", header = None, skiprows = get_skiprow(other_path), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
other = other.rename(columns={"QNAME": "QNAME-other", "FLAG": "FLAG-other", "RNAME": "RNAME-other", "POS": "POS-other", "MAPQ": "MAPQ-other"})
other = other.groupby('QNAME-other').head(2).reset_index(drop=True)
other = other.sort_values(by = ['QNAME-other', 'FLAG-other']).reset_index()
other = other.drop(['index'], axis=1)

other_merged = pd.concat([other, ref], axis=1)
other_merged['other-same'] = other_merged.swifter.apply(get_other_distance, axis=1)
other_merged[other_merged['same'] & other_merged['other-same']]


###########

other_path = output_path + 'accalign-cano.sam'
cano = pd.read_csv(other_path, sep ="\t", header = None, skiprows = get_skiprow(other_path), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
cano = cano.rename(columns={"QNAME": "QNAME-c", "FLAG": "FLAG-c", "RNAME": "RNAME-c", "POS": "POS-c", "MAPQ": "MAPQ-c"})
cano = cano.groupby('QNAME-c').head(2).reset_index(drop=True)
cano = cano.sort_values(by = ['QNAME-c', 'FLAG-c']).reset_index()
cano = cano.drop(['index'], axis=1)

other_path = output_path + 'accalign-non-cano.sam'
non_cano = pd.read_csv(other_path, sep ="\t", header = None, skiprows = get_skiprow(other_path), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
non_cano = non_cano.rename(columns={"QNAME": "QNAME-n", "FLAG": "FLAG-n", "RNAME": "RNAME-n", "POS": "POS-n", "MAPQ": "MAPQ-n"})
non_cano = non_cano.groupby('QNAME-n').head(2).reset_index(drop=True)
non_cano = non_cano.sort_values(by = ['QNAME-n', 'FLAG-n']).reset_index()
non_cano = non_cano.drop(['index'], axis=1)

def get_distance(x):
    return x["RNAME-c"] == x["RNAME-n"] and x["POS-c"] == x["POS-n"]

other_merged = pd.concat([cano, non_cano], axis=1)
other_merged['other-same'] = other_merged.swifter.apply(get_distance, axis=1)
diff=other_merged[~other_merged['other-same']]
diff.to_csv(output_path + 'diff-cano-noncano.sam')

28   HISEQ1:93:H2YHMBCXX:1:1101:10000:40611      83    chrX  58578846  HISEQ1:93:H2YHMBCXX:1:1101:10000:40611      83    chrX  61692090       False
29   HISEQ1:93:H2YHMBCXX:1:1101:10000:40611     163    chrX  58578760  HISEQ1:93:H2YHMBCXX:1:1101:10000:40611     163    chrX  61692004       False
31   HISEQ1:93:H2YHMBCXX:1:1101:10000:43311     147    chr1  13139461  HISEQ1:93:H2YHMBCXX:1:1101:10000:43311     147    chr1  13139495       False
33   HISEQ1:93:H2YHMBCXX:1:1101:10000:43362     147   chr19   3145182  HISEQ1:93:H2YHMBCXX:1:1101:10000:43362     147   chr19   3145079       False
70   HISEQ1:93:H2YHMBCXX:1:1101:10000:90499      99    chr2  92276678  HISEQ1:93:H2YHMBCXX:1:1101:10000:90499      99   chr11  48859465       False
71   HISEQ1:93:H2YHMBCXX:1:1101:10000:90499     147    chr2  92275801  HISEQ1:93:H2YHMBCXX:1:1101:10000:90499     147   chr11  48858969       False
78   HISEQ1:93:H2YHMBCXX:1:1101:10000:98515      83    chr7  23794138  HISEQ1:93:H2YHMBCXX:1:1101:10000:98515      83    chr7  23794137       False
82   HISEQ1:93:H2YHMBCXX:1:1101:10001:19719      99    chr2  92320754  HISEQ1:93:H2YHMBCXX:1:1101:10001:19719      99    chr2  92325168       False
83   HISEQ1:93:H2YHMBCXX:1:1101:10001:19719     147    chr2  92321326  HISEQ1:93:H2YHMBCXX:1:1101:10001:19719     147    chr2  92325168       False
168  HISEQ1:93:H2YHMBCXX:1:1101:10001:98077      83    chr8  43794715  HISEQ1:93:H2YHMBCXX:1:1101:10001:98077      83    chr2  92312698       False


aa_path = output_path + 'accalign-kmer32-kmer21.sam'

aa = pd.read_csv(aa_path, sep ="\t", header = None, skiprows = get_skiprow(aligner), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
aa = aa.sort_values(by=['QNAME', 'FLAG']).reset_index()
aa = aa.drop(['index'], axis=1)

aa_merged = pd.concat([aa, ref], axis=1)
aa_merged['aa_same'] = aa_merged.swifter.apply(get_distance, axis=1)
aa_same = aa_merged[aa_merged['same'] & aa_merged['aa_same']]
print("Not aligned AA(%):", len(aa[aa["RNAME"] == "*"])/len(aa)* 100)
print("AA same as concordance: ", len(aa_same))

aa_diff = aa_merged[aa_merged['same'] & (~aa_merged['aa_same'])]
aa_diff = aa_diff[['QNAME', 'FLAG', 'RNAME', 'POS', 'FLAG-bwa', 'RNAME-bwa', 'POS-bwa']]
aa_diff.to_csv(output_path + 'diff-two-index.csv', index = False)

print("Total strands: ", len(bwa))
print("BWA and Bowtie aligned same(concordance): ", len(concordance))
print("Other same as concordance: ", len(merged[merged['other-same']]))

print("Not aligned BWA(%):", len(bwa[bwa["RNAME-bwa"] == "*"])/len(bwa)* 100)
print("Not aligned Bowtie(%):", len(bowtie2[bowtie2["RNAME-bowtie2"] == "*"])/len(bowtie2)* 100)
print("Not aligned Other(%):", len(other[other["RNAME-other"] == "*"])/len(other)* 100)
print("Not aligned AA(%):", len(aa[aa["RNAME"] == "*"])/len(aa)* 100)


########### two index--- START ###########
colName = ["QNAME", "FLAG", "RNAME", "POS", "AS"]
colIndex = [0, 1, 2, 3, 12]

aa32_path = output_path + 'accalign-kmer32.sam'
aa32 = pd.read_csv(aa32_path, sep ="\t| ", header = None, skiprows = get_skiprow(aligner), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
aa32 = aa32.sort_values(by=['QNAME', 'FLAG']).reset_index()
aa32 = aa32.drop(['index'], axis=1)
aa32['AS'] = aa32['AS'].str.split(':').str[2].astype('int')
aa32_merged = pd.concat([aa32, ref], axis=1)
aa32_merged['aa_same'] = aa32_merged.swifter.apply(get_distance, axis=1)

aa21_path = output_path + 'accalign-kmer21.sam'
aa21 = pd.read_csv(aa21_path, sep ="\t| ", header = None, skiprows = get_skiprow(aligner), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
aa21 = aa21.sort_values(by=['QNAME', 'FLAG']).reset_index()
aa21 = aa21.drop(['index'], axis=1)
aa21['AS'] = aa21['AS'].str.split(':').str[2].astype('int')
aa21_merged = pd.concat([aa21, ref], axis=1)
aa21_merged['aa_same'] = aa21_merged.swifter.apply(get_distance, axis=1)

bad_thres = 248
bad32 = aa32[aa32['AS'] < bad_thres] # 2 mismatch
distinct_values = bad32['QNAME'].unique()
bad32_qname = pd.DataFrame({'QNAME': distinct_values})

good32 = aa32_merged[~aa32_merged['QNAME'].isin(bad32_qname['QNAME'])]
good32_same = good32[good32['same'] & good32['aa_same']]
print("High AS:", len(good32))
print("High AS aligned:", len(good32_same))

bad32_21 = pd.merge(bad32_qname, aa21_merged, how='left', left_on=['QNAME'], right_on=['QNAME'])
bad32_21_same = bad32_21[bad32_21['same'] & bad32_21['aa_same']]
print("Low AS:", len(bad32_21))
print("Low AS aligned:", len(bad32_21_same))

res = pd.concat([good32_same, bad32_21_same], ignore_index=True)
res.to_csv(output_path + 'two-index-manually-good.csv', index = False)

print("Not aligned AA(%):", len(aa[aa["RNAME"] == "*"])/len(aa)* 100)


########### two index--- END ###########


########### analysis the reads diff because AA not aligned [STEP]--- START ###########
# not_same = pd.read_csv(output_path + 'diff-step16-cheap.csv')
aa_diff_aa = pd.merge(aa_diff, aa, how='left', left_on=['QNAME', 'FLAG'], right_on=['QNAME', 'FLAG'])
aa_diff_aa_unalign = aa_diff_aa[aa_diff_aa['POS_x']==0]

distinct_values = aa_diff_aa_unalign['QNAME'].unique()
diff_unalign_qname = pd.DataFrame({'QNAME': distinct_values})
colName = ["QNAME", "FLAG", "RNAME", "POS", "CIGAR", "NM"]
colIndex = [0, 1, 2, 3, 5, 11]
bwa = pd.read_csv(bwa_path, sep = "\t", header = None, skiprows = get_skiprow("bwa"), usecols= colIndex, names = colName, quoting=csv.QUOTE_NONE)
diff_unalign_bwa = pd.merge(diff_unalign_qname, bwa, how='left', left_on=['QNAME'], right_on=['QNAME'])
diff_unalign_bwa = diff_unalign_bwa.groupby('QNAME').head(2).reset_index(drop=True)

diff_unalign_bwa.to_csv(output_path + 'diff-step16-cheap-unalign.csv', index = False)


########### analysis the reads diff because AA not aligned [STEP]--- END ###########


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





