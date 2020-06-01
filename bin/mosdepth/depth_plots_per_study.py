from numpy import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from glob import glob
import ptitprince as pt
import os 

di = os.getcwd()
len(glob(di + '/*.region.dist.txt'))
savefigs = True
figs_dir = di + '/'

def export_fig(axis,text, fname):
    if savefigs:
        axis.text()
        axis.savefig(fname, bbox_inches='tight')

lst_df = []
for file in glob(di + '/*.region.dist.txt'):
    #print(file)
    df_tmp = pd.read_table(file, header=None, names = ['chr', 'depth', 'percent'])
    df_tmp = df_tmp.query("chr == 'total'")
    df_tmp['sample'] = file
    lst_df.append(df_tmp)
df = pd.concat(lst_df)
df.to_csv(di + '/depth_overlap.csv', index=None)
#df = df.query('depth < 60')
depths = [10, 15, 20, 30, 50] #, 100]
depths.reverse()
df = df[df['depth'].isin(depths)]
df.head()
len(unique(df['sample']))

# make plot
mpl.rcParams['figure.dpi'] = 300
dx = "depth"; dy = "percent"; ort = "h"; pal = "colorblind"; sigma = .2
f, ax = plt.subplots(figsize=(7, 5))
pt.RainCloud(x = dx, y = dy, data = df, width_viol = .9, ax = ax, palette=pal, scale="area", jitter=1, pointplot = False)
ax.set_ylabel('Fraction of exome covered at min. depth')
ax.set_xlabel('Coverage at least')
# labels = [item.get_text() for item in ax.get_xticklabels()]
# labels = [">=" + x + "X" for x in labels]
# ax.set_xticklabels(labels)
sns.despine()
ax.set_xlim(-0.7)

if savefigs:
    plt.savefig(di + '/figure1.png', bbox_inches='tight')

df.query('depth == 20').percent.hist(bins=50)

if savefigs:
    plt.savefig(di + '/figure2.png', bbox_inches='tight')

f, ax = plt.subplots(figsize=(7, 5))

sns.violinplot(x = dx, y = dy, data = df, ax = ax, palette=pal, scale="area", jitter=1, pointplot = False)

ax.set_ylabel('Proportion of exome covered (across all samples)')
ax.set_xlabel('Minimal coverage')

# labels = [item.get_text() for item in ax.get_xticklabels()]
# labels = [">=" + x + "X" for x in labels]

# ax.set_xticklabels(labels)
sns.despine()

if savefigs:
    plt.savefig(di + '/figure3.png', bbox_inches='tight')
