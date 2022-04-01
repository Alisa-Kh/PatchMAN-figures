#!/usr/bin/python3
import math
import sys

import pandas as pd
import matplotlib.pyplot as plt
import os


def header_index(headers: list):
    output = ""
    cnt = 0
    while cnt <= len(headers) - 1:
        for i in range(3):
            output += '{:2d} {:25s}\t'.format(cnt, headers[cnt])
            cnt += 1
            if cnt >= len(headers):
                break
        output += '\n'
    return output


if __name__ == '__main__':
    filename = sys.argv[1]
    top10 = None
    if len(sys.argv) > 2:
        color = sys.argv[2]
        plot_name = sys.argv[3]
        score_term = sys.argv[4]
    if len(sys.argv) == 6:
        top10 = sys.argv[5]
    elif len(sys.argv) == 2:
        color = 'b'
        plot_name = 'score_vs_rmsbb'
        score_term = 'reweighted_sc'

    rm_percent = 0.51 # default remove the worst 50%

    with open(sys.argv[5], 'r') as fh:
        tmplist = fh.readlines()
    match_rmsds = {}
    for l in tmplist:
        match_rmsds[l.split()[0]] = l.split()[1]

    dfs = []
    headers = None
    raw_data = []
    with open(filename, 'r') as f:
        next(f)  # skip sequence line
        headers = next(f).split()[1:]
        for row in f:
            raw_data.append(row.split()[1:])

    raw_data = pd.DataFrame(raw_data, columns=headers)
    for head in headers[:-1]:
        raw_data[head] = raw_data[head].astype(float)

    y = score_term
    if score_term == 'reweighted_sc':
        y_label = 'Rosetta Score'
    elif score_term == 'I_sc':
        y_label = 'Rosetta Interface Score'
    else:
        y_label = 'Rosetta Score'
    x = 'rmsBB_if'
    x_label = 'RMSD, Å'

    # for df, name in dfs:
    df = raw_data.sort_values(by=y, ascending=False)

    # match_rmsds_sorted = []
    # for name in df.description:
    #     match_rmsds_sorted.append(float(match_rmsds['_'.join(name.split('_')[:4])]))
    #
    # df['Match_rms'] = match_rmsds_sorted

    import matplotlib as mpl
    import matplotlib.cm as cm

    norm = mpl.colors.Normalize(vmin=0, vmax=2.5)
    cmap = cm.summer

    # top1_percent = int(len(df) * 0.01)
    df = df[int(len(df) * rm_percent):].reset_index()

    match_rmsds_sorted = []
    for name in df.description:
        match_rmsds_sorted.append(float(match_rmsds['_'.join(name.split('_')[:4])]))

    df['Match_rms'] = match_rmsds_sorted

    fig = plt.figure(figsize=(5,4), tight_layout=1)
    ax = fig.add_subplot(111)

    plt.plot(df[x].values, df[y].values, 'o', color=cmap(df['Match_rms'].values), alpha=0.6)
    # For FERM min run
    # plt.plot(df[x].values[-1], df[y].values[-1], '*', color='red', markersize=12, alpha=0.8) # 118_1j19_319_326
    # plt.plot(df[x].values[-2], df[y].values[-2], 's', color='blue', markersize=8, alpha=0.8) # 125_2zpy_296_303
    # plt.plot(1.778, -1035.935, '^', color='purple', markersize=12, alpha=0.8) # 121_2i2a_178_185

    # For FERM nomin run
    # plt.plot(4.072, -1014.494, '*', color='red', markersize=12, alpha=0.8) # 118_1j19_319_326
    # plt.plot(1.379, -996.877, 's', color='blue', markersize=8, alpha=0.8) # 125_2zpy_296_303
    # plt.plot(2.045, -934.767, '^', color='purple', markersize=8, alpha=0.8) # 121_2i2a_178_185

    # For 5jiu 2.5
    # plt.plot(4.571, -702.584, '*', color='red', markersize=12, alpha=0.8)

    if top10:
        with open(top10, 'r') as tc:
            top_clusters = [l.split()[0] for l in tc.readlines()][1:11]
        df_top10 = df.loc[df['description'].isin(top_clusters)]
        plt.plot(df_top10[x].values, df_top10[y].values, 'o', color='r', alpha=1)
    
    plt.vlines(5, min(df[y].values)-100, max(df[y].values), color='red')
    plt.xlim([0, 35])
    if min(df[y].values)+100 < max(df[y].values):
        plt.ylim([min(df[y].values)-5, min(df[y].values)+80])
    else:
        plt.ylim([min(df[y].values)-5, max(df[y].values)])

    left, right = plt.xlim()
    plt.ylabel(y_label, fontsize=14, fontweight='bold')
    plt.xlabel(x_label, fontsize=14, fontweight='bold')

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    ax.set_title('2.5Å cutoff', y=1.0, pad=-14, loc='right', color='w', fontweight='bold', fontsize=14)

    # ax.axhspan(min(df[y].values) - 10, df[y].values[-top1_percent],  facecolor='lightgreen', alpha=0.5)

    plt.show()
    # plt.savefig(plot_name + '.png', dpi=600)
