import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


motif_dataset = ['1CZY', '1EG4', '1ELW', '1JD5', '1JWG', '1MFG', '1NTV', '1RXZ', '1SSH', '1X2R', '2A3I', '2CCH']
nonmotif = ['1AWR', '1ER8', '1NVR', '1NX1', '1OU8', '1U00', '2B9H', '2C3I', '2DS8', '2FMF', '2H9M', '2HPL', '2O02', '3D1E']

fig = plt.figure(figsize=(5,4))

with open('all_top1per_under5A_pepalign', 'r') as fh:
    all_top_models = fh.readlines()

with open('best_models_pepids', 'r') as s:
    tmp = s.readlines()
top_names = [l.split()[0].strip() for l in tmp]
top_ids = [float(l.split()[1].strip()) for l in tmp]
names_motif_indicated = ['* ' + l if l in motif_dataset else l for l in top_names]

names = []
values = []
for i,l in enumerate(all_top_models):
    if len(l.strip()) == 4 and i ==0:
        names.append(l.strip())
        values_for_complex = []
    elif len(l.strip()) == 4:
        names.append(l.strip())
        values.append(values_for_complex)
        values_for_complex = []
    else:
        values_for_complex.append(float(l.split()[3]))
values.append(values_for_complex)


with open('allSorted_by_lowest_rmsd_filters', 'r') as s:
    sorted_names = [l.strip() for l in s.readlines() if l.strip() in names]
    names_motif_indicated = ['* ' + l if l in motif_dataset else l for l in sorted_names]

# The position of the bars on the x-axis
xaxis = list(range(len(top_names)))

barWidth = 0.8

ax = fig.add_subplot(111)

for i,n in enumerate(sorted_names):
    idx = names.index(n)
    id0 = len([c for c in values[idx] if c == 0])
    id30 = len([c for c in values[idx] if 0<c <=0.3])
    id60 = len([c for c in values[idx] if 0.3<c <=0.6])
    id90 = len([c for c in values[idx] if 0.6<c <=0.9])
    id100 = len([c for c in values[idx] if c > 0.9])
    if top_ids[idx] == 0:
        b1 = ax.bar(i, id0, color='royalblue', edgecolor='mediumblue', linewidth=1.5, width=barWidth)
    else:
        b1 = ax.bar(i, id0, color='royalblue', width=barWidth)
    if 0 < top_ids[idx] <= 0.3:
        b2 = ax.bar(i, id30, bottom=id0, color='limegreen', edgecolor='darkgreen', linewidth=1.5, width=barWidth)
    else:
        b2 = ax.bar(i, id30, bottom=id0, color='limegreen', width=barWidth)
    if 0.3 < top_ids[idx] <= 0.6:
        b3 = ax.bar(i, id60, bottom=id0+id30,color='gold', edgecolor='orange', linewidth=1.5, width=barWidth)
    else:
        b3 = ax.bar(i, id60, bottom=id0+id30,color='gold', width=barWidth)
    if 0.6 < top_ids[idx] <= 0.9:
        b4 = ax.bar(i, id90, bottom=id0+id30+id60, color='orange', edgecolor='darkorange', linewidth=1.5, width=barWidth)
    else:
        b4 = ax.bar(i, id90, bottom=id0+id30+id60, color='orange', width=barWidth)
    if top_ids[idx] > 0.9:
        b5 = ax.bar(i, id100, bottom=id0+id30+id60+id90, color='indianred', edgecolor='darkred', linewidth=1.5, width=barWidth)
    else:
        b5 = ax.bar(i, id100, bottom=id0+id30+id60+id90, color='indianred', width=barWidth)

    blue_patch = mpatches.Patch(color='royalblue', label='0% ID')
    green_patch = mpatches.Patch(color='limegreen', label='<30% ID')
    yellow_patch = mpatches.Patch(color='gold', label='30%-60% ID')
    orange_patch = mpatches.Patch(color='orange', label='60%-90% ID')
    red_patch = mpatches.Patch(color='indianred', label='>90% ID')
    ax.legend(handles=[blue_patch, green_patch, yellow_patch, orange_patch, red_patch])

# # Custom X axis
plt.xticks(xaxis, names_motif_indicated, fontweight='bold', rotation='vertical')
plt.xlabel("Complex", fontsize=12)
plt.ylabel("Top scoring models under 5$\AA$", fontsize=12)

idx_1nvr = names.index('1NVR')
plt.ylim(0,len(values[idx_1nvr]) + 1)

plt.gcf().subplots_adjust(bottom=0.2)

# Show graphic
# plt.show()
plt.savefig('pepseq_ids.png', dpi=600, bbox_inches='tight')
