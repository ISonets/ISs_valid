#!/usr/bin/python3
import numpy as np
import pandas as pd
import os
import sys
# important hardcoded variables
isescan_cols = [0,2,3,4,5]
blast_cols = [0,2,3,4,7]
correct_order = ['seqID', 'Title','isBegin', 'isEnd', 'Length']
filename_blast = '.txt'
coverage_threshold = 1.25
length_threshold = 0.5
#  reading .depth and ISEscan files
depth_file = pd.read_csv(sys.argv[1], sep = '\t', header = None, names = ['Contig', 'Coord', 'Depth'])
folder_path = os.path.dirname(sys.argv[1])
is_scan_res = pd.read_csv(sys.argv[2], usecols = isescan_cols, skiprows = 1, names = ['seqID', 'Title','isBegin', 'isEnd', 'Length'])
#  coverage median across assembly
median_cov = np.median(depth_file['Depth'])
#  coverage median and length for each contig
depth_by_contig_df = depth_file.groupby('Contig')['Depth'].median().reset_index()
length_by_contig_df = depth_file.groupby('Contig')['Coord'].size().reset_index()
length_by_contig_df.columns = ['Contig', 'Length']
# reading BLAST files from folder and storing as one df
blast_res = pd.DataFrame()
for file_name in os.listdir(folder_path):
    if file_name.endswith(filename_blast):
        file_path = os.path.join(folder_path, file_name)
        #  checking if file is empty
        if os.path.getsize(file_path) > 0:
            try:
                temp_df = pd.read_csv(file_path, sep = '\t', usecols=blast_cols, names = ['seqID', 'Length', 'isBegin', 'isEnd', 'Title'])
                blast_res = pd.concat([blast_res, temp_df])
            except:
                pass
#  merging files in 1 dataframe. NEW 3/05/24
is_scan_res = is_scan_res[correct_order]
blast_res = blast_res[correct_order]
comb_file = pd.concat([is_scan_res, blast_res],ignore_index=True)
#  calculating median coverage for each IS/determinant according to contig using .depth file
median_depth_calc_df =[]
for _, row in comb_file.iterrows():
    start = row['isBegin']
    stop = row['isEnd']
    name = row['seqID']
    depth_of_contig = depth_file[depth_file['Contig'] == name]
    window_data = depth_of_contig[(depth_of_contig['Coord'] >= start) & (depth_of_contig['Coord'] <= stop)]
    median_window = window_data['Depth'].median()
    median_depth_calc_df.append(median_window)
comb_file['Median_Depth'] = median_depth_calc_df
#  checking if determinant found by ISEscan/BLAST is over-covered or excceds 50% of contig's length
comb_file_check_result_df = pd.DataFrame(columns=comb_file.columns)
for _, row in comb_file.iterrows():
    name = row['seqID']
    median_of_determt = row['Median_Depth']
    length_of_determt = row['Length']
    median_by_contig = depth_by_contig_df.loc[depth_by_contig_df['Contig'] == name, 'Depth'].values[0]
    length_by_contig = length_by_contig_df.loc[length_by_contig_df['Contig'] == name, 'Length'].values[0]
    try:
        if length_of_determt >= length_threshold * length_by_contig :
            if median_by_contig >= length_threshold * median_cov:
                comb_file_check_result_df = pd.concat([comb_file_check_result_df, row.to_frame().T], ignore_index=True)
        else:
            if median_of_determt >= length_threshold * median_by_contig:
                comb_file_check_result_df = pd.concat([comb_file_check_result_df, row.to_frame().T], ignore_index=True)
    except:
        pass
filename_comb_file_check_res = os.path.join(folder_path, 'determinants_verified' + '.csv')
comb_file_check_result_df.to_csv(filename_comb_file_check_res, sep ='\t', header = True, index = False)
#  NEW 24/03/24: saving individual contig as separate IS
overcovered_whole_contig = pd.DataFrame(columns=depth_by_contig_df.columns)
for _, row in depth_by_contig_df.iterrows():
    if row['Depth'] >= median_cov:
        overcovered_whole_contig = pd.concat([overcovered_whole_contig,row.to_frame().T],ignore_index=True)
filename_overcovered_contigs = os.path.join(folder_path, 'overcovered_contigs' + '.csv')
overcovered_whole_contig.to_csv(filename_overcovered_contigs, sep='\t', header=True, index=False)
