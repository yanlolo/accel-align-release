from __future__ import division
import pandas as pd, sys, csv, swifter
import matplotlib.pyplot as plt

path = '/home/yiqing/data/fsva-hg37/hg37.fna.hash21'
out_dir = '/home/yiqing/yan/output/'

data = pd.read_csv(path, sep = ",", header = None, skiprows = 2, quoting=csv.QUOTE_NONE)


plt.yscale('log')
# plt.bar(data[0], data[1], color='blue', alpha=0.7)
plt.hist(data[1], bins=500, color='blue', edgecolor='black')
plt.xlabel('key value')
plt.ylabel('nb of positions per key')
plt.savefig(out_dir + "kmer21-non-cano.png")

ehit = data[1].apply(lambda x: x**2).sum()
print('[yyan-log] nb of keys: ' + str(len(data)))
print('[yyan-log] nb of positions: ' + str(data[1].sum()))
print('[yyan-log] Average (nb pos per key): ', data[1].mean())
print('[yyan-log] Variance (nb pos per key): ', data[1].var())
print('[yyan-log] E-hits (nb pos per key*nb pos per key/total nb of positions): ', ehit/data[1].sum())








name = 'mini'

# path = '/Users/yan/Desktop/eurecom_code/accel-align-release/data/hg37-mini.fna.hash'
# name = "mod2_30_prime_HXX0"
output = "/home/yiqing/yan/output/"
# name = str(sys.argv[1])
# output = str(sys.argv[2])

high_freq = 1000

data = pd.read_csv(path, sep = "\t", header = None, skiprows = 2, quoting=csv.QUOTE_NONE)
ehit = data[0].apply(lambda x: x**2).sum()
high_freq_data = data[data[0] > high_freq]


sorted_data = data.sort_values(by=0)

value_counts = data[0].value_counts().reindex(range(20), fill_value=0)

low_freq_data = data[data[0] <= high_freq]

plt.yscale('log')
plt.hist(low_freq_data[0], bins=50, color='blue', edgecolor='black')
plt.xlabel('nb of positions per key')
plt.ylabel('nb of keys')
plt.title('Histogram')
plt.savefig(output + name + ".png")

print('[yyan-log] mode:', name)
print('[yyan-log] nb of keys: ' + str(len(data)))
print('[yyan-log] nb of positions: ' + str(data[0].sum()))
print('[yyan-log] Average (nb pos per key): ', data[0].mean())
print('[yyan-log] Variance (nb pos per key): ', data[0].var())
print('[yyan-log] E-hits (nb pos per key*nb pos per key/total nb of positions): ', ehit/data[0].sum())
print('[yyan-log] nb of high freq seed (key level): ', len(high_freq_data))
print('[yyan-log] fraction of high freq seed (key level) (%): ', len(high_freq_data)/len(data)*100)
print('[yyan-log] nb of high freq seed (position level): ', high_freq_data[0].sum())
print('[yyan-log] fraction of high freq seed (position level) (%): ', high_freq_data[0].sum()/data[0].sum())



