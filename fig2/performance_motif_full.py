#!/usr/bin/python3
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
from matplotlib.patches import Rectangle
import numpy as np

cutoffs = [1,2,3,4,5]

with open(sys.argv[1], 'r') as fh:
    header1 = fh.readline().split()
    results1 = fh.readlines()

with open(sys.argv[2], 'r') as fh:
    header2 = fh.readline().split()
    results2 = fh.readlines()

rms_measure = sys.argv[3]


def collect_data(header, data, rms):
    complexes = [l.split()[0] for l in data]
    columns = []
    for i, h in enumerate(header):
        if h.startswith(rms):
            columns.append(i)
    methods = {}
    for c in columns:
        method = header[c].split(':')[1].strip()
        method_values = []
        for line in data:
            method_values.append(float(line.split()[c]))
        methods[method] = method_values
    return complexes, methods


complexes1, methods1 = collect_data(header1, results1, rms_measure)
complexes2, methods2 = collect_data(header2, results2, rms_measure)


def under_cutoff(values):
    cumulative_count = []
    n = 0
    for c in cutoffs:
        for v in values:
            if v < c:
                n += 1
        cumulative_count.append(n)
        n = 0
    return cumulative_count

counts1 = {}
for a in methods1.keys():
    counts1[a] = under_cutoff(methods1[a])
counts2 = {}
for b in methods2.keys():
    counts2[b] = under_cutoff(methods2[b])

fig = plt.figure(figsize=(5,4), tight_layout=1)

colormap = cm.get_cmap('tab10')
for count in counts1.keys():
    if count == 'PatchMAN':
        pm1, = plt.plot(cutoffs, [v/len(complexes1) for v in counts1[count]], '^-', linewidth=3, label=count, color='firebrick', markersize=12)
    elif count == 'PFPD':
        pfpd1, = plt.plot(cutoffs, [v/len(complexes1) for v in counts1[count]], '.-', linewidth=3, label=count, color='cornflowerblue', markersize=12)
    else:
        continue
        # plt.plot(cutoffs, [v/len(complexes1) for v in counts1[count]], '.-', linewidth=3, label=count, color='darkblue', markersize=12)

for count in counts2.keys():
    if count == 'PatchMAN':
        pm2, = plt.plot(cutoffs, [v/len(complexes2) for v in counts2[count]], linestyle='dashed', label='PatchMAN motif', color='firebrick', linewidth=3, markersize=12)
    elif count == 'PFPD':
        pfpd2, = plt.plot(cutoffs, [v/len(complexes2) for v in counts2[count]], linestyle='dashed', linewidth=3, label='PFPD motif', color='cornflowerblue', markersize=12)


# create blank rectangle
extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)

#Create organized list containing all handles for table. Extra represent empty space
legend_handle = [extra, extra, extra, extra, pm1, pm2, extra, pfpd1,pfpd2]
#Define the labels
label_row_1 = ["   ", "Full peptide", "Motif only"]
label_j_1 = ["PatchMAN"]
label_j_2 = ["PFPD"]
label_empty = [""]

#organize labels for table construction
# legend_labels = np.concatenate([label_row_1, label_j_1, label_j_2, label_empty * 3])
legend_labels = np.concatenate([label_row_1, label_j_1, label_empty*2, label_j_2,  label_empty*2])
plt.legend(legend_handle, legend_labels,
           ncol = 3, shadow = True, handletextpad = -2)


# plt.legend(fontsize=12)
plt.grid(b=1, which='both', axis='y', linewidth=0.3)
plt.ylabel('Complexes within cutoff', fontsize=16)
plt.xlabel('RMSD Cutoff, $\AA$', fontsize=16)
plt.xticks(cutoffs, fontsize=14)
plt.yticks([0,0.2,0.4,0.6,0.8], fontsize=14)

# plt.title('Overall performance', fontweight='bold', fontsize=22)
#plt.minorticks_on()
plt.rc('axes', linewidth=1)
# plt.show()
plt.savefig('motif_full.png', dpi=300)
