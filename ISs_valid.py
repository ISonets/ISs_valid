#!/usr/bin/python
import numpy as np
import pandas as pd
import os
import sys
import argparse
from tqdm import tqdm

isecan_cols = [0,2,3,4]
blast_cols = [0,2,3,4,7]
#  thresholds also saved here, in case I will get rid of argparse
#threshold_cov = 1.5
#threshold_len = 0.5

def extract_info_from_depth_file(res_folder):
    '''This functon reads .depth file of assembly made from SPAdes/Flye assembly(depends on data type) 
    and calculating median coverage for each contig, also calculating length of each contig and stores it in
    separate dataframes for ease of access.'''
    #  reading .depth file
    for file_name in os.listdir(res_folder):
        if file_name.endswith('.depth'):
            depth_file = pd.read_csv(file_name, sep = '\t', header = None, names = ['Contig', 'Coord', 'Depth'])
            depth_by_contig_df = depth_file.groupby('Contig')['Depth'].median().reset_index()
            length_by_contig_df = depth_file.groupby('Contig')['Coord'].size().reset_index()
            length_by_contig_df.columns = ['Contig', 'Length']
    #  AFAIK I should return, but maybe not?
    return depth_by_contig_df, length_by_contig_df
#  it seems I should add call to this function into the main(), and provied calls to other functions here. So it would be like main()==>extract_info()==>other functions. But I'm not sure.

def median_calculator(res_folder):
    '''This function calculates median coverage for each determinant stored in vairous files from ISEscan or BLAST. 
    Data formats are different,so we nedd first read each file, and later using start and stop coordinates, 
    we calculate median coverage to compare it with contig coverage. 
    Column names for BLAST file are renamed using ISEscan as an example to simplify coding.'''
    #  reading files, either BLAST or ISEscan. They have different extension and structure, so we read them differently.
    for file_name in os.listdir(res_folder):
        if file_name.endswith('.csv'):
            data = pd.read_csv(file_name, usecols = isescan_cols)
            #  ISEscan file doesn't have length column, like a BLAST, so we create it.
            data['Length'] = int(data['isEnd']) - int(data['isBegin'])
        elif file_name.endswith('.txt'):
            data = pd.read_csv(filename, sep = '\t', header = None, usecols = blast_cols, names = ['seqID', 'Length', 'isBegin', 'isEnd', 'Title'])
    #  calculating median
    #  empty list to store loop results, later add them as an end column for data.
    median_cov_calc_output_res = []
    for _, row in data.iterrows():
        start = row['isBegin']
        stop = row['isEnd']
        determinant_window_data = depth_file[(depth_file['Coord'] >= start) & (depth_file['Coord'] <= stop)]
        median_cov_in_window = determinant_window_data['Depth'].median()
        median_cov_calc_output_res = median_cov_calc_output_res.concat([median_cov_calc_output_res,median_cov_in_window])
        data['Median_Depth'] = median_cov_calc_output_res
        data_mod = data.dropna(axis = 0)
    return data_mod

def compare_coverage(depth_by_contig, data_mod):
    '''This function compares coverage of each determinant with median coverage of each contig. 
    If it exceeds 1.5 times or more than contig coverage, store it in separate file.'''
    sel_determinants_cov = pd.DataFrame(columns=data_mod.columns)
    for _, row in data_mod.iterrows():
        name = row['seqID']
        median_cov_of_determinant = row['Median_Depth']
        print(name, median_cov_of_determinant)
        # try/except is here, because BLAST sometimes shits its pants. If smth wrong, skip the fucking row entirely
        try:
            contig_median_cov = depth_by_contig_df.loc[depth_by_contig_df['Contig'] == name, 'Depth'].values[0]
            if median_cov_of_determinant >= threshold_cov * contig_median_cov:
                sel_determinants_cov = sel_determinants_cov.concat([sel_determinants_cov,row])
        except:
            pass
    #  dunno how to pass filename here(((
    sel_determinants_cov.to_csv()

def compare_length(length_by_contig, data_mod):
    ''' This function compares length of determinant and length of contig. 
    If len(determninant) 50% or more of contig length, store it in separate file.'''
    sel_determinants_length = pd.DataFrame(columns=data.columns)
    for _, row in data.iterrows():
        name = row['seqID']
        length_of_determinant = row['Length']
        #  same logic as for coverage, thanks BLAST sukablyat'
        try:
            contig_length = length_by_contig_df.loc[length_by_contig_df['Contig'] == name, 'Length'].values[0]
            if length_of_determinant >= threshold_len * contig_length:
                sel_determinants_length = sel_determinants_length.concat([sel_determinants_length,row])
        except:
            pass
    #  same shit here, dunno how to save files this way(
    sel_determinants_length.to_csv()

#  main function do not working, didn't add inside it any claas to next functions(I'm just dumbass). Fix ASAP
def main():
    '''Main function to filter vairous files by length and coverage'''
    parser = argparse.ArgumentParser(description='ISs_valid is a script to filter results of ISEScan and BLAST for rRNA, ISs and phage genes by length and covarage.')
    group1 = parser.add_argument_group("Path to BLAST/ISEscan folder")
    group2 = parser.add_argument_group("Thresholds")
    group3 = parser.add_argument_group("Path to store filtered BLAST/ISEscan files")
    group1.add_argument('-input_dir',
                        type=str,
                        help='Path to folder where ISEscan and BLAST res are.',
                        default=None) 
    group2.add_argument('-threshold_cov',
                        type=int,
                        help='Magic number for filtration by coverage',
                        default=1.5)
    group2.add_argument('-threshold_len',
                        type=int,
                        help='Magic number for filtration by length',
                        default=0.5)
    group3.add_argument('-output_dir',
                        type=str,
                        help='Path to output directory',
                        default='./') 
    args = parser.parse_args()
    res_folder = args.input_dir
    threshold_cov = args.threshold_cov
    threshold_len = args.threshold_len
    out_dir = args.output_dir
    if args.input_dir == None:
        print('Error: Provide a path fo folder with depth file, ISEscan and BLAST results. Exiting!')
        sys.exit()

if __name__ == "__main__":
    main()
