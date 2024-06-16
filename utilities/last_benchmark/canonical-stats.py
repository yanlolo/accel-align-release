
from __future__ import division
import pandas as pd, swifter, sys, csv

cano_path = "/home/yiqing/code/accel-align-release/log-non-cano-32-250"
cano = pd.read_csv(cano_path, sep = ";", header = None)
cano = cano.rename(columns={0: "name", 1: "cnt-cano", 2: "cnt-pair-cano"})
cano['cnt-cano'].sum()
cano['cnt-pair-cano'].sum()


non_cano_path = "/home/yiqing/code/accel-align-release/log-non-cano-150"

cano = pd.read_csv(cano_path, sep = ";", header = None)
non_cano = pd.read_csv(non_cano_path, sep = ";", header = None)
cano = cano.rename(columns={0: "name", 1: "cnt-cano"})
non_cano = non_cano.rename(columns={0: "name", 1: "non-cnt-cano"})

cano['cnt-cano'].sum()
non_cano['non-cnt-cano'].sum()



cano_embed = pd.read_csv(cano_embed_path, sep = ";", header = None)
non_cano_embed = pd.read_csv(non_cano_embed_path, sep = ";", header = None)
cano_embed = cano_embed.rename(columns={0: "name", 1: "cnt-cano-embed"})
non_cano_embed = non_cano_embed.rename(columns={0: "name", 1: "non-cnt-cano-embed"})

cano_grp = cano.groupby('name').sum()
non_cano_grp = non_cano.groupby('name').sum()
cano_embed_grp = cano_embed.groupby('name').sum()
non_cano_embed_grp = non_cano_embed.groupby('name').sum()

cano_grp = cano_grp.sort_values(by = ['name']).reset_index()
non_cano_grp = non_cano_grp.sort_values(by = ['name']).reset_index()
cano_embed_grp = cano_embed_grp.sort_values(by = ['name']).reset_index()
non_cano_embed_grp = non_cano_embed_grp.sort_values(by = ['name']).reset_index()



cano_grp['cnt-cano'].sum()
non_cano_grp['non-cnt-cano'].sum()
cano_embed_grp['cnt-cano-embed'].sum()
non_cano_embed_grp['non-cnt-cano-embed'].sum()

merged = pd.concat([cano_grp, non_cano_grp, cano_embed_grp, non_cano_embed_grp], axis=1)
merged[merged['cnt-cano']<merged['non-cnt-cano']]

merged['same'] = merged[(merged['cnt-cano']<=merged['non-cnt-cano']) & (merged['cnt-cano-embed']>merged['non-cnt-cano-embed'])]
name='HISEQ1:93:H2YHMBCXX:1:1101:10001:77506'
cano_grp[cano_grp['name']==name]
non_cano_grp[non_cano_grp['name']==name]
cano_embed_grp[cano_embed_grp['name']==name]
non_cano_embed_grp[non_cano_embed_grp['name']==name]


>>> merged[(merged['cnt-cano']<=merged['non-cnt-cano']) & (merged['cnt-cano-embed']>merged['non-cnt-cano-embed'])].head(20)
name  cnt-cano                                    name  ...  cnt-cano-embed                                    name  non-cnt-cano-embed
72    HISEQ1:93:H2YHMBCXX:1:1101:10001:77506      1424  HISEQ1:93:H2YHMBCXX:1:1101:10001:77506  ...              72  HISEQ1:93:H2YHMBCXX:1:1101:10001:77506                   0
180   HISEQ1:93:H2YHMBCXX:1:1101:10005:15220       664  HISEQ1:93:H2YHMBCXX:1:1101:10005:15220  ...              88  HISEQ1:93:H2YHMBCXX:1:1101:10005:15220                   0
301   HISEQ1:93:H2YHMBCXX:1:1101:10008:59684        62  HISEQ1:93:H2YHMBCXX:1:1101:10008:59684  ...              62  HISEQ1:93:H2YHMBCXX:1:1101:10008:59684                   0
323   HISEQ1:93:H2YHMBCXX:1:1101:10009:32117      2065  HISEQ1:93:H2YHMBCXX:1:1101:10009:32117  ...              62  HISEQ1:93:H2YHMBCXX:1:1101:10009:32117                   0
330   HISEQ1:93:H2YHMBCXX:1:1101:10009:50505      1148  HISEQ1:93:H2YHMBCXX:1:1101:10009:50505  ...              85  HISEQ1:93:H2YHMBCXX:1:1101:10009:50505                   0
476   HISEQ1:93:H2YHMBCXX:1:1101:10013:66153       644  HISEQ1:93:H2YHMBCXX:1:1101:10013:66153  ...              75  HISEQ1:93:H2YHMBCXX:1:1101:10013:66153                   0
492   HISEQ1:93:H2YHMBCXX:1:1101:10013:97713       436  HISEQ1:93:H2YHMBCXX:1:1101:10013:97713  ...             195  HISEQ1:93:H2YHMBCXX:1:1101:10013:97713                   0
627   HISEQ1:93:H2YHMBCXX:1:1101:10018:18995      1870  HISEQ1:93:H2YHMBCXX:1:1101:10018:18995  ...             119  HISEQ1:93:H2YHMBCXX:1:1101:10018:18995                 111
670   HISEQ1:93:H2YHMBCXX:1:1101:10019:12969       166  HISEQ1:93:H2YHMBCXX:1:1101:10019:12969  ...              85  HISEQ1:93:H2YHMBCXX:1:1101:10019:12969                   0
783   HISEQ1:93:H2YHMBCXX:1:1101:10022:43009       403  HISEQ1:93:H2YHMBCXX:1:1101:10022:43009  ...              76  HISEQ1:93:H2YHMBCXX:1:1101:10022:43009                   0
949   HISEQ1:93:H2YHMBCXX:1:1101:10027:30502       271  HISEQ1:93:H2YHMBCXX:1:1101:10027:30502  ...             271  HISEQ1:93:H2YHMBCXX:1:1101:10027:30502                 201
984   HISEQ1:93:H2YHMBCXX:1:1101:10027:97573      7395  HISEQ1:93:H2YHMBCXX:1:1101:10027:97573  ...             406  HISEQ1:93:H2YHMBCXX:1:1101:10027:97573                   0
1213  HISEQ1:93:H2YHMBCXX:1:1101:10035:60265       224  HISEQ1:93:H2YHMBCXX:1:1101:10035:60265  ...              62  HISEQ1:93:H2YHMBCXX:1:1101:10035:60265                   0
1397  HISEQ1:93:H2YHMBCXX:1:1101:10040:63939       626  HISEQ1:93:H2YHMBCXX:1:1101:10040:63939  ...              72  HISEQ1:93:H2YHMBCXX:1:1101:10040:63939                   0
1506  HISEQ1:93:H2YHMBCXX:1:1101:10043:98508       865  HISEQ1:93:H2YHMBCXX:1:1101:10043:98508  ...             128  HISEQ1:93:H2YHMBCXX:1:1101:10043:98508                   0
1521  HISEQ1:93:H2YHMBCXX:1:1101:10044:51745       113  HISEQ1:93:H2YHMBCXX:1:1101:10044:51745  ...              56  HISEQ1:93:H2YHMBCXX:1:1101:10044:51745                   0
1551  HISEQ1:93:H2YHMBCXX:1:1101:10045:23351       293  HISEQ1:93:H2YHMBCXX:1:1101:10045:23351  ...              80  HISEQ1:93:H2YHMBCXX:1:1101:10045:23351                   0
1575  HISEQ1:93:H2YHMBCXX:1:1101:10045:93429       108  HISEQ1:93:H2YHMBCXX:1:1101:10045:93429  ...             108  HISEQ1:93:H2YHMBCXX:1:1101:10045:93429                   0
1733  HISEQ1:93:H2YHMBCXX:1:1101:10050:30788      1231  HISEQ1:93:H2YHMBCXX:1:1101:10050:30788  ...             265  HISEQ1:93:H2YHMBCXX:1:1101:10050:30788                 194
1882  HISEQ1:93:H2YHMBCXX:1:1101:10054:92124      2040  HISEQ1:93:H2YHMBCXX:1:1101:10054:92124  ...              84  HISEQ1:93:H2YHMBCXX:1:1101:10054:92124                   0



merged[~merged['same']]

cano_grp['cnt-cano'].sum()
non_cano_grp['non-cnt-cano'].sum()

## disable the slide
>>> cano_grp['cnt-cano'].sum()
252385251
>>> non_cano_grp['non-cnt-cano'].sum()
307862888

## candidates pass to embed (slide)
>>> cano_grp['cnt-cano'].sum()
146515010
>>> non_cano_grp['non-cnt-cano'].sum()
137128390


>>> non_cano[non_cano['name']=='HISEQ1:93:H2YHMBCXX:1:1101:10001:19719']
name  non-cnt-cano
487320  HISEQ1:93:H2YHMBCXX:1:1101:10001:19719             3
487321  HISEQ1:93:H2YHMBCXX:1:1101:10001:19719             0
487322  HISEQ1:93:H2YHMBCXX:1:1101:10001:19719             0
487323  HISEQ1:93:H2YHMBCXX:1:1101:10001:19719             4
>>> cano[cano['name']=='HISEQ1:93:H2YHMBCXX:1:1101:10001:19719']
name  cnt-cano
487320  HISEQ1:93:H2YHMBCXX:1:1101:10001:19719         3
487321  HISEQ1:93:H2YHMBCXX:1:1101:10001:19719         0
487322  HISEQ1:93:H2YHMBCXX:1:1101:10001:19719         9
487323  HISEQ1:93:H2YHMBCXX:1:1101:10001:19719         0

>>> non_cano[non_cano['name']=='HISEQ1:93:H2YHMBCXX:1:1102:9998:31193']
name  non-cnt-cano
3580040  HISEQ1:93:H2YHMBCXX:1:1102:9998:31193             0
3580041  HISEQ1:93:H2YHMBCXX:1:1102:9998:31193             1
3580042  HISEQ1:93:H2YHMBCXX:1:1102:9998:31193             1
3580043  HISEQ1:93:H2YHMBCXX:1:1102:9998:31193             0
>>> cano[cano['name']=='HISEQ1:93:H2YHMBCXX:1:1102:9998:31193']
name  cnt-cano
3580040  HISEQ1:93:H2YHMBCXX:1:1102:9998:31193         1
3580041  HISEQ1:93:H2YHMBCXX:1:1102:9998:31193         0
3580042  HISEQ1:93:H2YHMBCXX:1:1102:9998:31193       610
3580043  HISEQ1:93:H2YHMBCXX:1:1102:9998:31193         0


merged[merged['cnt-cano']>merged['non-cnt-cano']]

merged[merged['cnt-cano']<merged['non-cnt-cano']]


>>> merged[merged['cnt-cano']>merged['non-cnt-cano']]
name  cnt-cano                                    name  non-cnt-cano   same
14      HISEQ1:93:H2YHMBCXX:1:1101:10000:40611        60  HISEQ1:93:H2YHMBCXX:1:1101:10000:40611            59  False
41      HISEQ1:93:H2YHMBCXX:1:1101:10001:19719        12  HISEQ1:93:H2YHMBCXX:1:1101:10001:19719             7  False
54      HISEQ1:93:H2YHMBCXX:1:1101:10001:37547         3  HISEQ1:93:H2YHMBCXX:1:1101:10001:37547             2  False
55      HISEQ1:93:H2YHMBCXX:1:1101:10001:38124         3  HISEQ1:93:H2YHMBCXX:1:1101:10001:38124             2  False
72      HISEQ1:93:H2YHMBCXX:1:1101:10001:77506        72  HISEQ1:93:H2YHMBCXX:1:1101:10001:77506             0  False
...                                        ...       ...                                     ...           ...    ...
999952   HISEQ1:93:H2YHMBCXX:1:1102:9996:18275         9   HISEQ1:93:H2YHMBCXX:1:1102:9996:18275             7  False
999974   HISEQ1:93:H2YHMBCXX:1:1102:9997:28079         3   HISEQ1:93:H2YHMBCXX:1:1102:9997:28079             2  False
999975   HISEQ1:93:H2YHMBCXX:1:1102:9997:29913      7380   HISEQ1:93:H2YHMBCXX:1:1102:9997:29913             0  False
999976   HISEQ1:93:H2YHMBCXX:1:1102:9997:33355       130   HISEQ1:93:H2YHMBCXX:1:1102:9997:33355           127  False
999987   HISEQ1:93:H2YHMBCXX:1:1102:9998:31193       611   HISEQ1:93:H2YHMBCXX:1:1102:9998:31193             2  False




