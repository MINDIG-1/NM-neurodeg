import joblib
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as plt
import shutil
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import time
from statsmodels.stats.multitest import multipletests  # Importing the multipletests function

from pyS3.utils import *
import random
random.seed(42)
# np.random.seed(42)

zs_path_root = "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Results/Zs"
zs_folders = [file for file in os.listdir(zs_path_root) if file.startswith("Allch")]
zs_file = "zsAllCh.csv"
devmaps_path = "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Results/DeviationMaps"
result_path = "/home/ubuntu/Documents/codes/MEEG-normative-modeling/NM_shared/Results/Sig2"

os.makedirs(result_path, exist_ok=True)


def zs_IO(zs_df, column_names, group):
    zs_subset = zs_df[zs_df['group']==group][column_names].values
    IO_pos = np.where(zs_subset > 2, 1, 0)
    IO_neg = np.where(zs_subset < -2, 1, 0)
    return IO_pos, IO_neg


def zs_counts (IO):
    IO = IO[IO.sum(axis=1) != 0]
    if IO.shape[0] != 0:
        l = 100 / IO.shape[0] #subjects
    else:
        l = 0
    subj_ch = IO.sum(axis = 0)*l
    return subj_ch


def permutation_test(ip, zs_df, column_names, group_list, folder, col_dict, perm_test):
    fb = folder.split('-')[1]
    print(f'Band {fb} --- Permutation {ip}')
    
    nullval_pos_dict = {group: [] for group in group_list}
    nullval_neg_dict = {group: [] for group in group_list}

    if perm_test == 'group':
        labels_shuffled = zs_df['group'].sample(frac=1).reset_index(drop=True)
        zs_df_shuffled = zs_df.copy()
        zs_df_shuffled['group'] = np.array(labels_shuffled)
    elif perm_test == 'spatial':
        zs_df_shuffled = zs_df.copy()
        zs_df_shuffled[column_names] = np.random.permutation(zs_df[column_names].values.T).T

    IO_pos_HC, IO_neg_HC = zs_IO(zs_df_shuffled, column_names, "HC_te")
    HC_temp_pos = zs_counts(IO_pos_HC)
    HC_temp_neg = zs_counts(IO_neg_HC)

    for i, group in enumerate(group_list):
        IO_pos, IO_neg = zs_IO(zs_df_shuffled, column_names, group)

        if level == 'region' or level == 'spectral':
            grp_temp = zs_counts(IO_pos)
        grp_temp = grp_temp - HC_temp_pos
        nullval_pos_dict[group].append(grp_temp)

        if level == 'region' or level == 'spectral':
            grp_temp = zs_counts(IO_neg)
        grp_temp = grp_temp - HC_temp_neg
        nullval_neg_dict[group].append(grp_temp.reshape(1, -1))

    return nullval_pos_dict, nullval_neg_dict


def apply_fdr_correction(p_values):
    _, p_values_corrected, _, _ = multipletests(p_values, method='fdr_bh')
    return p_values_corrected

def find_rank(sorted_arr_pos, true_val_pos):
    indices = []
    for col_idx in range(sorted_arr_pos.shape[1]):
        # Find the index of the first occurrence of true_val_pos in the current column
        row_idx = np.argmax(sorted_arr_pos[:, col_idx] == true_val_pos[col_idx])
        # Append the row index to the list of indices
        indices.append(row_idx)
    return indices



nperm = 5000
perm_test = 'group' #group' # 'group' or 'spatial'
level = 'region' # 'network' or 'region' or 'spectral'

group_list = ["HC_te", "PD", "AD"]


#######################
zs_path_base = f'{zs_path_root}/{zs_folders[1]}/{zs_file}'
zs_df_base = pd.read_csv(zs_path_base)
column_names = [col for col in zs_df_base.columns if col.startswith('zs_X')]
if  level == 'region':
    col_dict = {f"zs_X{i}": f"zs_X{i}" for i in range(0, len(column_names))}
    col_list = [f"zs_X{i}" for i in range(1, len(column_names)+1)]
    true_net = pd.read_csv(f'{devmaps_path}/EDsubjectsPch.csv')
    true_net_pos = true_net[true_net['posneg']=='pos']
    true_net_neg = true_net[true_net['posneg']=='neg']
    true_net_pos = true_net_pos[column_names + ['group', 'f_band']]
    true_net_neg = true_net_neg[column_names + ['group', 'f_band']]
elif level == 'spectral':
    col_dict = {f"zs_X{i}": f"zs_X{i}" for i in range(1, len(column_names)+1)}
    col_list = [f"zs_X{i}" for i in range(1, len(column_names)+1)]
    true_net = pd.read_csv(f'{devmaps_path}/EDsubjectsPch.csv')
    true_net_pos = true_net[true_net['posneg']=='pos']
    true_net_neg = true_net[true_net['posneg']=='neg']
    true_net_pos = true_net_pos[column_names + ['group', 'f_band']]
    true_net_neg = true_net_neg[column_names + ['group', 'f_band']]

        
zs_folders = ['Allch-alpha']
for folder in zs_folders:

    fb = folder.split('-')[1]
    print(fb)

    nullval_pos_dict = {gp: [] for gp in group_list}
    nullval_neg_dict = {gp: [] for gp in group_list}
    all_nullval_pos = []
    all_nullval_neg = []
    
    df_pval_pos_all = pd.DataFrame(columns=list(col_dict.keys())+['group', 'f_band'])
    df_pval_neg_all = pd.DataFrame(columns=list(col_dict.keys())+['group', 'f_band'])
    df_pval_pos_fdr_all = pd.DataFrame(columns=list(col_dict.keys())+['group', 'f_band'])
    df_pval_neg_fdr_all = pd.DataFrame(columns=list(col_dict.keys())+['group', 'f_band'])

    true_net_pos_fb = np.array(true_net_pos[true_net_pos['f_band']==fb])
    true_net_neg_fb = np.array(true_net_neg[true_net_neg['f_band']==fb])

    true_HC_pos = true_net_pos_fb[np.where(true_net_pos_fb == "HC_te")[0],:]
    true_HC_neg = true_net_neg_fb[np.where(true_net_neg_fb == "HC_te")[0],:]

    true_net_pos_fb[:, :-2] = true_net_pos_fb[:, :-2] - true_HC_pos[:, :-2]
    true_net_neg_fb[:, :-2] = true_net_neg_fb[:, :-2] - true_HC_neg[:, :-2]
    
    zs_path = f'{zs_path_root}/{folder}/{zs_file}'
    zs_df = pd.read_csv(zs_path, index_col=False)
    zs_df = zs_df[zs_df['group']!= "HC_tr"]


    # # Use joblib to parallelize permutation tests
    parallel_result = joblib.Parallel(n_jobs=-1)(
        joblib.delayed(permutation_test)(ip, zs_df, column_names, group_list, folder, col_dict, perm_test)
        for ip in range(2, nperm + 1)
    )

    all_nullval_pos.extend(result[0] for result in parallel_result)
    all_nullval_neg.extend(result[1] for result in parallel_result)

    for group in group_list:
        nullval_pos_dict[group] = [d[group][0].reshape(1, -1) for d in all_nullval_pos]        
        nullval_pos_dict[group].append(true_net_pos_fb[np.where(true_net_pos_fb == group)[0],:-2])
        nullval_neg_dict[group] = [d[group][0] for d in all_nullval_neg]
        nullval_neg_dict[group].append(true_net_neg_fb[np.where(true_net_neg_fb == group)[0],:-2])

        df_nullval_pos = np.concatenate(nullval_pos_dict[group])
        df_nullval_neg = np.concatenate(nullval_neg_dict[group])
        

        # Create an empty DataFrame with specific column names as networks for pvalues
        df_pval_pos = pd.DataFrame(columns=list(col_dict.keys())+['group', 'f_band'])
        df_pval_neg = pd.DataFrame(columns=list(col_dict.keys())+['group', 'f_band'])
        df_pval_pos.loc[0, 'group'], df_pval_neg.loc[0, 'group'] = group, group
        df_pval_pos.loc[0, 'f_band'], df_pval_neg.loc[0, 'f_band'] = fb, fb
        
        true_val_pos, true_val_neg = true_net_pos_fb[np.where(true_net_pos_fb == group)[0],:-2][0], true_net_neg_fb[np.where(true_net_neg_fb == group)[0],:-2][0]
        sorted_arr_pos_des, sorted_arr_neg_des = np.sort(df_nullval_pos, axis=0)[::-1],  np.sort(df_nullval_neg, axis=0)[::-1]
        sorted_arr_pos_asc, sorted_arr_neg_asc = np.sort(df_nullval_pos, axis=0)[::1],  np.sort(df_nullval_neg, axis=0)[::1]
        
        rank_pos_des = find_rank(sorted_arr_pos_des,true_val_pos)
        rank_neg_des = find_rank(sorted_arr_neg_des,true_val_neg)
        rank_pos_asc = find_rank(sorted_arr_pos_asc,true_val_pos)
        rank_neg_asc = find_rank(sorted_arr_neg_asc,true_val_neg)
        min_rank_pos = []
        min_rank_neg = []

        # Compare ranks for each column and keep the lower rank
        for col_index in range(sorted_arr_pos_des.shape[1]):
            min_rank_pos.append(min(rank_pos_des[col_index], rank_pos_asc[col_index]))
            min_rank_neg.append(min(rank_neg_des[col_index], rank_neg_asc[col_index]))

        # Convert lists to NumPy arrays
        min_rank_pos = np.array(min_rank_pos)
        min_rank_neg = np.array(min_rank_neg)


        pval_pos_list = [round(value / nperm, 4) for value in min_rank_pos]
        pval_neg_list = [round(value / nperm, 4) for value in min_rank_neg]
        df_pval_pos.iloc[0,:-2] = pval_pos_list
        df_pval_neg.iloc[0,:-2] = pval_neg_list

        
        print(f"{group} - pos: {np.sum((np.array(pval_pos_list) <0.05))} - neg: {np.sum(np.array(pval_neg_list) <0.05)}")

        df_pval_pos_all = pd.concat([df_pval_pos_all, df_pval_pos], ignore_index=True)
        df_pval_neg_all = pd.concat([df_pval_neg_all, df_pval_neg], ignore_index=True)
        
        # Apply FDR correction to p-values
        df_pval_pos_fdr = df_pval_pos.copy()
        df_pval_neg_fdr = df_pval_neg.copy()
        df_pval_pos_fdr[column_names] = apply_fdr_correction(pval_pos_list)
        df_pval_neg_fdr[column_names] = apply_fdr_correction(pval_neg_list)
        df_pval_pos_fdr_all = pd.concat([df_pval_pos_fdr_all, df_pval_pos_fdr], ignore_index=True)
        df_pval_neg_fdr_all = pd.concat([df_pval_neg_fdr_all, df_pval_neg_fdr], ignore_index=True)


    pval_save = os.path.join(result_path, folder)
    if not os.path.exists(pval_save):
        os.makedirs(pval_save)
        
    df_pval_pos_all.to_csv(os.path.join(pval_save, f'pval_{level}_{perm_test}perm_pos.csv'), index = False)
    df_pval_neg_all.to_csv(os.path.join(pval_save, f'pval_{level}_{perm_test}perm_neg.csv'), index = False)

    df_pval_pos_fdr_all.to_csv(os.path.join(pval_save, f'pval_{level}_{perm_test}perm_pos_fdr.csv'), index = False)
    df_pval_neg_fdr_all.to_csv(os.path.join(pval_save, f'pval_{level}_{perm_test}perm_neg_fdr.csv'), index = False)
