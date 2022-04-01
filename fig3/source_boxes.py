import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


Full = 14
Motif = 12
overall = Full + Motif

motif_dataset = ['1CZY', '1EG4', '1ELW', '1JD5', '1JWG', '1MFG', '1NTV', '1RXZ', '1SSH', '1X2R', '2A3I', '2CCH']
nonmotif = ['1AWR', '1ER8', '1NVR', '1NX1', '1OU8', '1U00', '2B9H', '2C3I', '2DS8', '2FMF', '2H9M', '2HPL', '2O02', '3D1E']

all_cmplxs = motif_dataset + nonmotif

with open('all_top1per_under5A_src', 'r') as f:
    all_top_models = f.readlines()

with open('best_models_sources', 'r') as f:
    tmp = f.readlines()
best_models = {}
for l in tmp:
    best_models[l.split()[0]] = l.split()[1].strip()

names = []
monomers = 0
interfaces = 0
sources = []

for i,l in enumerate(all_top_models):
    if len(l.strip()) == 4 and i == 0:
        names.append(l.strip())
    elif len(l.strip()) == 4:
        names.append(l.strip())
        sources.append((monomers, interfaces))
        monomers = 0
        interfaces = 0
    else:
        if l.strip().split()[1] == 'monomer':
            monomers += 1
        else:
            interfaces += 1

sources.append((monomers, interfaces))
for n in names:
    print(n)

with open('allSorted_by_lowest_rmsd_filters', 'r') as s:
    sorted_names = [l.strip() for l in s.readlines() if l.strip() in names]
    sorted_names_annotated = ['* ' + l if l in motif_dataset else l for l in sorted_names]

# The position of the bars on the x-axis
xaxis = list(range(len(names)))

barWidth = 0.8

fig = plt.figure(figsize=(5,4))
ax = fig.add_subplot(111)

for i,n in enumerate(sorted_names):
    idx = names.index(n)
    if best_models[n] == 'interface':
        b1 = ax.bar(i, sources[idx][0], color='orange', width=barWidth)
        b2 = ax.bar(i, sources[idx][1], bottom=sources[idx][0], color='royalblue', edgecolor='mediumblue', linewidth=1.5, width=barWidth)
    else:
        b1 = ax.bar(i, sources[idx][0], color='orange', edgecolor='darkorange', linewidth=1.5, width=barWidth)
        b2 = ax.bar(i, sources[idx][1], bottom=sources[idx][0], color='royalblue', width=barWidth)


orange_patch = mpatches.Patch(color='orange', label='monomers')
blue_patch = mpatches.Patch(color='royalblue', label='interfaces')
ax.legend(handles=[blue_patch, orange_patch])

# # Custom X axis
plt.xticks(xaxis, sorted_names_annotated, fontweight='bold', rotation='vertical')
plt.xlabel("Complex", fontsize=12)
plt.ylabel("Top scoring models under 5$\AA$", fontsize=12)

idx_1nvr = names.index('1NVR')
plt.ylim(0,sources[idx_1nvr][0] + sources[idx_1nvr][1] + 1)

plt.gcf().subplots_adjust(bottom=0.2)

# Show graphic
# plt.show()
plt.savefig('template_sources.png', dpi=600, bbox_inches='tight')
