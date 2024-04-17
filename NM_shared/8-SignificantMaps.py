import os
import pandas as pd
import numpy as np
import mne
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import mannwhitneyu
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import pairwise_distances
from sklearn.cluster import AgglomerativeClustering
from matplotlib import rcParams
from scipy.stats import mannwhitneyu
import matplotlib as mpl
import matplotlib.colors as mcolors

mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20

import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from pyS3.utils import *

zs_path_root = "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Results/Zs"
zs_folders = [file for file in os.listdir(zs_path_root) if file.startswith("Allch")]
zs_file = "zsAllCh.csv"
result_path = "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Results/DeviationMaps"
nch = 19

BUCKET = 'nm-neurodeg'
channels_file = 'NKEEG19_channels.sfp'
TMP_DOWNLOAD_DIR = os.path.expanduser('~/tmp/nm-tmp/')
os.makedirs(TMP_DOWNLOAD_DIR, exist_ok=True)
tmp_bool = download_file_from_s3(BUCKET, channels_file, TMP_DOWNLOAD_DIR)
if tmp_bool:
    montage19 = mne.channels.read_custom_montage(os.path.join(TMP_DOWNLOAD_DIR, channels_file))
    os.remove(os.path.join(TMP_DOWNLOAD_DIR, channels_file))


mont1020 = mne.channels.make_standard_montage('standard_1020')
kept_channels = montage19.ch_names
ind = [mont1020.ch_names.index(value) for value in kept_channels]
mont1020_19 = mont1020.copy()
mont1020_19.ch_names = [mont1020.ch_names[x] for x in ind]
kept_channel_info = [mont1020.dig[x+3] for x in ind]
mont1020_19.dig = mont1020.dig[0:3]+kept_channel_info
print(mont1020_19.ch_names)
info = mne.create_info(mont1020_19.ch_names, sfreq=200, ch_types='eeg')
info.set_montage(mont1020_19)
print(info)

column_names = [f"zs_X{i}" for i in range(1, nch+1)]

# Red map
darkest_color = '#730128' 
lightest_color = 'white'
adjusted_start_color = mcolors.rgb_to_hsv(mcolors.hex2color(darkest_color))
adjusted_start_color[1] = min(1.0, adjusted_start_color[1] * 1.2)  # Increase saturation
start_color = mcolors.rgb2hex(mcolors.hsv_to_rgb(adjusted_start_color))
adjusted_end_color = mcolors.rgb_to_hsv(mcolors.hex2color(lightest_color))
adjusted_end_color[1] = min(1.0, adjusted_end_color[1] * 1.2)  # Increase saturation
end_color = mcolors.rgb2hex(mcolors.hsv_to_rgb(adjusted_end_color))
cmap_r = mcolors.LinearSegmentedColormap.from_list("red_map", [start_color, end_color])
cmap_r = cmap_r.reversed()

# Blue map
darkest_color = '#0A2B6B'
lightest_color = 'white'
adjusted_start_color = mcolors.rgb_to_hsv(mcolors.hex2color(darkest_color))
adjusted_start_color[1] = min(1.0, adjusted_start_color[1] * 1.2)  # Increase saturation
start_color = mcolors.rgb2hex(mcolors.hsv_to_rgb(adjusted_start_color))
adjusted_end_color = mcolors.rgb_to_hsv(mcolors.hex2color(lightest_color))
adjusted_end_color[1] = min(1.0, adjusted_end_color[1] * 1.2)  # Increase saturation
end_color = mcolors.rgb2hex(mcolors.hsv_to_rgb(adjusted_end_color))
cmap_b = mcolors.LinearSegmentedColormap.from_list("blue_map", [start_color, end_color])
cmap_b = cmap_b.reversed()


# Red-Blue map
start_color = '#730128'
middle_color = 'white' 
end_color = '#0A2B6B'
adjusted_start_color = mcolors.rgb_to_hsv(mcolors.hex2color(start_color))
adjusted_start_color[1] = min(1.0, adjusted_start_color[1] * 1.2)  # Increase saturation
start_color = mcolors.rgb2hex(mcolors.hsv_to_rgb(adjusted_start_color))
adjusted_middle_color = mcolors.rgb_to_hsv(mcolors.hex2color(middle_color))
adjusted_middle_color[1] = min(1.0, adjusted_middle_color[1] * 1.2)  # Increase saturation
middle_color = mcolors.rgb2hex(mcolors.hsv_to_rgb(adjusted_middle_color))
adjusted_end_color = mcolors.rgb_to_hsv(mcolors.hex2color(end_color))
adjusted_end_color[1] = min(1.0, adjusted_end_color[1] * 1.2)  # Increase saturation
end_color = mcolors.rgb2hex(mcolors.hsv_to_rgb(adjusted_end_color))
cmap_br = mcolors.LinearSegmentedColormap.from_list("custom_map", [start_color, middle_color, end_color])
        
        
def zs_IO(zs_df, column_names, group):
    zs_subset = zs_df[zs_df['group']==group]
    zs_subset = zs_subset[column_names]
    IO_pos= zs_subset.apply(lambda x: x.map(lambda y: 1 if y > 2 else 0))
    IO_neg = zs_subset.apply(lambda x: x.map(lambda y: 1 if y < -2 else 0))
    IO_posneg = IO_pos + IO_neg
    zs_subset['global_id'] = zs_df[zs_df['group']==group]['global_id']
    zs_subset['dataset'] = zs_df[zs_df['group']==group]['dataset']
    return IO_pos, IO_neg, IO_posneg, zs_subset

def zs_counts (zs_df,column_names, group, zs_path):

        count_greater_2, count_lower_2, _, zs_subset = zs_IO(zs_df, column_names, group)
        if not os.path.exists(zs_path):
            os.makedirs(zs_path)

        temp_p = pd.concat([zs_df[zs_df['group']==group]['global_id'], count_greater_2], axis=1)
        temp_p.to_csv(zs_path+f'count_greater_2.csv', index = False)
        temp_n = pd.concat([zs_df[zs_df['group']==group]['global_id'], count_lower_2], axis=1)
        temp_n.to_csv(zs_path+f'count_lower_2.csv', index = False)
        
        temp_p = count_greater_2
        temp_n = count_lower_2
        count_greater_2 = count_greater_2[count_greater_2.sum(axis=1) != 0]
        count_lower_2 = count_lower_2[count_lower_2.sum(axis=1) != 0]
        x = 100/nch
        l1 = 100/(count_greater_2.shape[0]) #subjects
        l2 = 100/(count_lower_2.shape[0]) #subjects
        return count_greater_2.sum(axis = 0)*l1, temp_p.sum(axis = 1), count_lower_2.sum(axis =0)*l2, temp_n.sum(axis =1) , zs_subset['global_id'], zs_subset['dataset']
    
    
def table2fig (df, name):    
    f_bands = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma']
    # Group the DataFrame by 'f_band'
    df['group'] = df['group'].replace(group_map)
    grouped = df.groupby(pd.Categorical(df['f_band'], categories=f_bands, ordered=True))

    # Create subplots
    fig, axes = plt.subplots(nrows=5, figsize=(10, 5 * len(grouped)))

    # Custom column labels
    custom_col_labels = ['Group', "% at least one + ED","median [range] +", "% at least one - ED", "median [range] -"]
    rcParams['font.family'] = 'monospace'
    rcParams['font.monospace'] = 'DejaVu Sans Mono'
    color_dict = {'HC_tr': '#D4EBF2',  # Soft Blue
                'HC(test)': '#FADBD8',  # Soft Red 
                'PD': '#D4EDD2',   # Soft Green
                'AD': '#E6E6E6',  # Light Gray
                }

    # Iterate over groups and create tables in subplots
    for i, (fb_name, group_data) in enumerate(grouped):
        table_data = group_data.drop('f_band', axis=1).values
        colors = [color_dict[val] for val in group_data['group']]
        repeated_colors = np.array([colors]).T.repeat(len(group_data.columns)-1, axis=1)
    
        # Add a table to the subplot with custom column labels
        table = axes[i].table(cellText=table_data,
                            colLabels=custom_col_labels,
                            cellLoc='center',
                            loc='center', 
                            cellColours=repeated_colors)
        table.set_fontsize(40)
        title = axes[i].set_title(f'Table for {fb_name}', pad=0, loc='center')

        axes[i].axis('off')
        
    plt.tight_layout()
    plt.savefig(f'{result_path}/{name}.png')


def add_significance_stars(group_data, reference_group, ax):
    for variable in ['posChED', 'negChED']:
        if variable in group_data.columns and variable in reference_group.columns:
            for group_name in group_data['group'].unique():
                if group_name != 'HC(test)':
                    stat, p_value = mannwhitneyu(reference_group[variable], group_data[group_data['group'] == group_name][variable])
                    if p_value < 0.05:
                        print(group_name, variable)
                        y_position = -0.7
                        x_position = group_data['group'].unique().tolist().index(group_name)  # Get the x-coordinate for the group
                        if variable=='posChED':
                            if p_value <0.01:
                                ax.text(x_position-0.22, y_position, '*', ha='center', va='center', fontsize=18, fontweight='bold')
                                ax.text(x_position-0.18, y_position, '*', ha='center', va='center', fontsize=18, fontweight='bold')
                            else:
                                ax.text(x_position-0.2, y_position, '*', ha='center', va='center', fontsize=18, fontweight='bold')
                        else:
                            if p_value <0.01:
                                ax.text(x_position+0.22, y_position, '*', ha='center', va='center', fontsize=18, fontweight='bold')
                                ax.text(x_position+0.18, y_position, '*', ha='center', va='center', fontsize=18, fontweight='bold')
                            else:
                                ax.text(x_position+0.2, y_position, '*', ha='center', va='center', fontsize=18, fontweight='bold')
                        
                        stat_g, p_value_g = mannwhitneyu(reference_group[variable], group_data[group_data['group'] == group_name][variable], alternative = 'greater' )
                        if p_value_g < 0.05:
                            print(f'HC is greater than {group_name}  at {variable}')
                        stat_l, p_value_l = mannwhitneyu(reference_group[variable], group_data[group_data['group'] == group_name][variable], alternative = 'less' )
                        if p_value_l < 0.05:
                            print(f'HC is less than {group_name} at {variable}')
                        


vmin = 0
vmax = 100
max_p = float('-inf')
max_n = float('-inf')

group_list = ["HC_te", "PD", "AD"]
group_map = {
    "HC_te": "HC(test)",
    "PD": "PD",
    "AD": "AD",
}
atleastOneED = []
EDsubj_ch_list = []
chPsubj_list = []
for folder in (zs_folders):
    print(folder)

    zs_path = f'{zs_path_root}/{folder}/{zs_file}'
    zs_df = pd.read_csv(zs_path)
    zs_df = zs_df[zs_df['group'] != "HC_tr"]
    fig, axes = plt.subplots(nrows=len(group_list), ncols=2, figsize=(8, 24))

    for i, group in enumerate(group_list):
        print(group)

        IO_pos,IO_neg,_,_ = zs_IO(zs_df, column_names, group)
        temp_pos = round((1 - (IO_pos.sum(axis=1) == 0).sum() / IO_pos.shape[0]) * 100, 2)
        temp_neg = round((1 - (IO_neg.sum(axis=1) == 0).sum() / IO_neg.shape[0]) * 100, 2)
        
        atleastOneED.append(pd.DataFrame({"group": [group], "f_band": [folder[6:]], "pos_atleastOne": [temp_pos], "pos_range": f'{IO_pos.sum(axis=1).median()} [{IO_pos.sum(axis=1).min()}-{IO_pos.sum(axis=1).max()}]', "neg_atleastOne": [temp_neg],  "neg_range": f'{IO_neg.sum(axis=1).median()} [{IO_neg.sum(axis=1).min()}-{IO_neg.sum(axis=1).max()}]'}))
        
        pos_subj_ch, pos_ch_subj, neg_subj_ch, neg_ch_subj, id, dataset = zs_counts(zs_df, column_names, group, f'{result_path}/{folder}/{group}/')

        temp_df = pd.DataFrame(pos_subj_ch).T
        temp_df = temp_df.assign(group=group, f_band=folder[6:], posneg = 'pos')
        EDsubj_ch_list.append(temp_df)
        
        temp_df = pd.DataFrame(neg_subj_ch).T
        temp_df = temp_df.assign(group=group, f_band=folder[6:], posneg = 'neg')
        EDsubj_ch_list.append(temp_df)
        
        temp_chPsubj = pd.DataFrame({'global_id':id,'posChED':pos_ch_subj, 'negChED':neg_ch_subj, 'group': group, 'dataset':dataset, "f_band": folder[6:]})

        chPsubj_list.append(temp_chPsubj)

        print(f'max ch + : {round(pos_subj_ch.max(),2)}, max ch - {round(neg_subj_ch.max(),2)}, group: {group}, fb:{folder[6:]}')

        max_p = max(round(pos_subj_ch.max(),2), max_p)
        max_n = max(round(neg_subj_ch.max(),2), max_n)

        im1, cm1 = mne.viz.plot_topomap(pos_subj_ch, pos=info, axes=axes[i, 0], show=False, vlim=[vmin, vmax], cmap=cmap_b)
        divider1 = make_axes_locatable(axes[i, 0])
        cax1 = divider1.append_axes('right', size='5%')
        cbar1 = plt.colorbar(im1, cax=cax1)
        axes[i, 0].set_title(f"% > 2 {group_map[group]}-{folder[5:]}")

        im2, cm2 = mne.viz.plot_topomap(neg_subj_ch, pos=info, axes=axes[i, 1], show=False, vlim=[vmin, vmax], cmap=cmap_r)
        divider2 = make_axes_locatable(axes[i, 1])
        cax2 = divider2.append_axes('right', size='5%')
        cbar2 = plt.colorbar(im2, cax=cax2)
        axes[i, 1].set_title(f"% < -2 {group_map[group]}-{folder[5:]}")

    fig.savefig(f'{result_path}/{folder}/topo.png', transparent=True)
    plt.close(fig)

        
print(f'max_p= {max_p}, max_n={max_n}')

df_atleastchOneED = pd.concat(atleastOneED, ignore_index=True)
df_atleastchOneED.to_csv(os.path.join(result_path, 'atleastOne.csv'), index = False)

df_EDsubj_ch = pd.concat(EDsubj_ch_list, ignore_index=True)
df_EDsubj_ch.to_csv(os.path.join(result_path, 'EDsubjectsPch.csv'), index = False)

chPsubj = pd.concat(chPsubj_list)
chPsubj.to_csv(f'{result_path}/EDchPsubj.csv')

#------------------------------- Perc of at least one ED Tables
table2fig(df_atleastchOneED, 'atleastone-tables')

#------------------------------ number of ED channels P subject - violin plot

fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(15, 36))
fig.suptitle('Violin Plot for number of +/- ED for each subject')

df = chPsubj
# Group the DataFrame by 'f_band'
df = df[df['group'] != "HC_tr"]
df['group'] = df['group'].replace(group_map)

f_bands = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma']
grouped_df = df.groupby(pd.Categorical(df['f_band'], categories=f_bands, ordered=True))

# Define colors for positive and negative data
colors = {'posChED': '#0A2B6B', 'negChED': '#730128'}


# Iterate through the frequency bands and create violin plots with scatter for each group
for (f_band, group), ax in zip(grouped_df, axes.flatten()):
    melted_data = pd.melt(group[['posChED', 'negChED', 'group']].reset_index(),
                         id_vars=['index', 'group'], var_name='variable', value_name='value')

    sns.violinplot(x='group', y='value', hue='variable', split=False, inner='quart', cut = 0, palette=colors, alpha=0.3,
                   data=melted_data, ax=ax)

    strip = sns.stripplot(x='group', y='value', hue='variable', data=melted_data, palette=colors, jitter=True, dodge=True, ax=ax)
    
    print(f'***{f_band}')
    add_significance_stars(group, df[df['group'] == 'HC(test)'], ax)

    # Customize legend labels
    handles, labels = strip.get_legend_handles_labels()
    labels = ['+', '-']
    ax.legend(handles, labels)
    
    ax.set_title(f'Frequency Band: {f_band}')
    ax.set_xlabel('Group')
    ax.set_ylabel('Count')

# Adjust layout
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(f'{result_path}/violinbChPsubj.png')

fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(15, 36))
fig.suptitle('Violin Plot for number of +/- ED for each subject')


#------------------------------ Density plot for overlap %

colors2 = {'pos': '#0A2B6B', 'neg': '#730128'}

fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(25, 36))
fig.suptitle('Density Plot for Overlap %') 

df = df_EDsubj_ch
df = df[df['group'] != "HC_tr"]
# Group the DataFrame by 'f_band'
grouped_df = df.groupby('f_band')
df_long = pd.melt(df, id_vars=['f_band', 'group', 'posneg'], var_name='column', value_name='value')
print(df_long.head())
# Create separate plots for pos and neg groups
g = sns.FacetGrid(df_long, col="f_band", row="group", hue="posneg", margin_titles=True, palette=colors2)
g.map(sns.kdeplot, "value", shade=True)

# Add legend
g.add_legend()

# Adjust layout
g.fig.tight_layout()
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(f'{result_path}/OverlapDensity.png')
