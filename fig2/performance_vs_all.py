#!/usr/bin/python3
import matplotlib.pyplot as plt
from matplotlib import cm
import sys

cutoffs = [1,2,3,4,5]

with open(sys.argv[1], 'r') as fh:
    header1 = fh.readline().split()
    results1 = fh.readlines()

# with open(sys.argv[2], 'r') as fh:
#     header2 = fh.readline().split()
#     results2 = fh.readlines()

rms_measure = sys.argv[2]


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
            try:
                method_values.append(float(line.split()[c]))
            except ValueError:
                continue
        methods[method] = method_values
    return complexes, methods


complexes1, methods1 = collect_data(header1, results1, rms_measure)
# complexes2, methods2 = collect_data(header2, results2, rms_measure)


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
# counts2 = {}
# for b in methods2.keys():
#     counts2[b] = under_cutoff(methods2[b])

fig = plt.figure(figsize=(5,4), tight_layout=1)

colormap = cm.get_cmap('Set1')
i=2
for count in counts1.keys():
    if count == 'PatchMAN':
        plt.plot(cutoffs, [v/len(complexes1) for v in counts1[count]], '^-', linewidth=3, label=count, color='firebrick', markersize=12)
    elif count == 'AF2':
        continue
        # plt.plot(cutoffs, [v/len(complexes1) for v in counts1[count]], '*-', linewidth=3, label=count, color='gold', markersize=12)
    elif count == 'PFPD':
        continue
        # plt.plot(cutoffs, [v/len(complexes1) for v in counts1[count]], '.-', linewidth=3, label=count, color='dodgerblue', markersize=12)
    else:
        plt.plot(cutoffs, [v/len(complexes1) for v in counts1[count]], '.-', linewidth=3, label=count,
                 color=colormap.colors[i], markersize=12)
        i +=1

plt.legend(fontsize=12)
plt.grid(b=1, which='both', axis='y', linewidth=0.3)
plt.ylabel('Complexes within cutoff', fontsize=16)
plt.xlabel('RMSD Cutoff, $\AA$', fontsize=16)
plt.xticks(cutoffs, fontsize=14)
plt.yticks([0,0.2,0.4,0.6,0.8], fontsize=14)

# plt.title('Overall performance', fontweight='bold', fontsize=22)
#plt.minorticks_on()
plt.rc('axes', linewidth=1)
plt.show()
# plt.savefig('patchman_vs_all.png', dpi=600)
