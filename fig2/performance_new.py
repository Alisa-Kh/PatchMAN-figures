#!/usr/bin/python3
import matplotlib.pyplot as plt

import sys

cutoffs = [1,2,3,4,5]

with open(sys.argv[1], 'r') as fh:
    header = fh.readline().split()
    all_results = fh.readlines()

with open(sys.argv[2], 'r') as fh:
    n_cmplxes = [c.strip() for c in fh.readlines()]

complexes = []
pm = []
af = []

new_only = [l for l in all_results if l.split()[0].upper() in n_cmplxes]
for line in new_only:
    complexes.append(line.split()[0])
    pm.append(float(line.split()[header.index('rmsBB_if:PatchMAN')]))
pm_all = []
for line in all_results:
    complexes.append(line.split()[0])
    pm_all.append(float(line.split()[header.index('rmsBB_if:PatchMAN')]))

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

pm_all_count = under_cutoff(pm_all)
pm_count = under_cutoff(pm)

fig = plt.figure(figsize=(5,4), tight_layout=1)

plt.plot(cutoffs, [v/len(pm) for v in pm_count], '*-', label='PatchMAN new targets', color='coral', linewidth=3, markersize=12)
plt.plot(cutoffs, [v/len(pm_all) for v in pm_all_count], '^-', linewidth=3, label='PatchMAN all', color='firebrick', markersize=12)

plt.legend(fontsize=12)
plt.grid(b=1, which='both', axis='y', linewidth=0.3)
plt.ylabel('Complexes within cutoff', fontsize=16)
plt.xlabel('RMSD Cutoff, $\AA$', fontsize=16)
plt.xticks(cutoffs, fontsize=14)
plt.yticks([0,0.2,0.4,0.6,0.8], fontsize=14)

# plt.title('Overall performance', fontweight='bold', fontsize=18)
#plt.minorticks_on()
plt.rc('axes', linewidth=1)
plt.show()
# plt.savefig('patchman_performance_onNEW.png', dpi=600)
