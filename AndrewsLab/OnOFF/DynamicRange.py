# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 15:35:55 2023

@author: zengm
"""
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.ticker as ticker
from pathlib import Path
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
plt.rcParams["errorbar.capsize"] = 2.0
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['axes.xmargin'] = 1
plt.rcParams['axes.ymargin'] = 1

plt.tick_params(
    axis='both',          # changes apply to the y-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelright=False,
    direction='in',
    colors='k',
    labelcolor='k') # labels along the bottom edge are off

##For controlling the appearance of ticklabels
params = {'axes.titlesize': 8, 'xtick.labelsize': 7, 'ytick.labelsize':7, 'axes.labelsize': 8}
plt.rcParams.update(params)
plt.rcParams['axes.titlepad'] = 5

fig,ax = plt.subplots(nrows=1,ncols=3, figsize=(6, 1.44), sharey=True, gridspec_kw={'wspace':0, 'width_ratios': [5,7,0]})
# The width of graph = (# of x-axis item) * 0.5

ax1 = ax[0]
ax2 = ax[1]
ax3 = ax[2]

ax2.spines['left'].set_visible(False)
ax3.spines['left'].set_visible(False)

x_minor = mpl.ticker.LogLocator(base=10.0, subs="auto", numticks=999)

for axi in ax.flat:
    axi.spines['right'].set_visible(False)
    axi.spines['top'].set_visible(False)
    axi.set_yscale('log')
    axi.set_ylim(0.1,100)
    axi.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10,numticks=15)) 
    axi.yaxis.set_minor_locator(x_minor) 
    axi.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    axi.tick_params(which='both', width=0.5)
    axi.tick_params(which='major', length=4) 
    axi.tick_params(which='minor', length=2)
    axi.tick_params(which='both', direction='in')    
    axi.tick_params(axis='both', which='major', pad=5)

def rbs(input_file, rbs_num, ax):
    df1 =  pd.read_csv(f'{input_file}_{rbs_num}_dr.csv', converters={'Version': lambda x: str(x)[1:]})
    df_rbs = df1[df1['RBS'] == rbs_num]
    errors_dr = df_rbs.groupby(['Version'], sort=True)['dr'].std().reset_index()
    avg_dr = df_rbs.groupby(['Version'], sort=True)['dr'].mean().reset_index()
    
    In=avg_dr['Version']
    x=np.arange(len(In))

    width=0.4
    ax.bar(x, avg_dr['dr'], yerr=errors_dr['dr'], width=width, bottom=0.001, align='center', log=False, color='black', 
        error_kw=dict(lw=0.5, capthick=0.5), edgecolor='None', alpha=0.65, linewidth = 0.5)
    # ax.legend(prop={"size":8}, loc = 1, fancybox = True, framealpha=0, ncol=3)
    ax.set_xticks(x)
    ax.set_xticklabels(In, rotation=90, ha = 'center')
    
rbs('Rhl', 1, ax1)
# rbs(2, ax2)
rbs('Rhl', 3, ax2)

ax1.set_xbound(-0.5, 4.5)  # lower = -0.5, upper = number of items - 0.5
ax2.set_xbound(-0.5, 6.5)
# ax3.set_xbound(-0.5, 6.5)

ax2.tick_params(
    axis='y',          # changes apply to the y-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the left edge are off
    length = 0)

ax3.tick_params(
    axis='y',          # changes apply to the y-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the left edge are off
    length = 0)

ax1.set_ylabel('Dynamic range')
ax1.set_xlabel('$\mathregular{P_{Rhl}}$ variants')
ax2.set_ylabel('')
# ax3.set_ylabel('')

fig.set_size_inches(3,0.72)
fig.autofmt_xdate(rotation=90, ha='center')

plt.savefig('Dynamic_range_example.svg',transparent=True)


