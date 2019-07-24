#!/usr/bin/env python3

import h5py
import os
import sys
import pandas as pd
import numpy as np
from dtaidistance import dtw
import itertools

## input data (first commandline argument)
input_file_name = sys.argv[1]
#data_all = pd.read_csv(input_directory,header=0)

## penalty
penalty = int(sys.argv[2])

## k-mer list
kmer_file = pd.read_csv(str(sys.argv[3]), header=None)
kmer_subset =list(kmer_file[0])

## output folder
out_folder = str(sys.argv[4])

############## Process fast5 file
hf = h5py.File(input_file_name, 'r')


# hf.keys() will have the respective folders: [u'Analyses', u'PreviousReadInfo', u'Raw', u'UniqueGlobalKey']
# read id extraction
read_id = list(hf['/Raw/Reads/'].keys())[0]   ## Use list() for python 3

## extracting signal values
# signal location in fast5: '/Raw/Reads/Read_746'
# adding read id to the path
signal_path = '/Raw/Reads/' + read_id + '/Signal'
# index and signal
signal_array = list(hf.get(signal_path))
# events location in the fast5: '/Analyses/Basecall_1D_001/BaseCalled_template/Events/'
# Columns in events index: object 2: start, 4: 5mer, 5: move value
events_table = hf.get('/Analyses/Basecall_1D_001/BaseCalled_template/Events/')
kmer_signal_index_array = []
move_val = 0
combo_str = ''

for line in events_table[()]:
    ele = str(line).split(", ")
    move_val = int(ele[5])
    if move_val != 0:
        kmer = (ele[4]).replace('\'', '')[1:6]
        start = int(ele[2])
        # extract start and stop into an array
        combo_str = combo_str + str(start-1)      # this will introduce a -1 in the first line of the array
        kmer_signal_index_array.append(combo_str)
        # print combo_str
        combo_str = (kmer + "\t" + str(start) + "\t")
        
readID_list = []
kmer_list = []
signal_list = []

# extracting signal based on the start and stop
for coordinates in kmer_signal_index_array[1:]:
    ele2 = coordinates.split("\t")
    kmer = ele2[0]
    start = int(ele2[1])
    stop = int(ele2[2])
    signal_values = "_".join(map(str, signal_array[start:stop+1]))
    readID_list.append(read_id)
    kmer_list.append(kmer)
    signal_list.append(signal_values)
    
data_dict = {'read_ID':readID_list, 'kmer':kmer_list, 'values': signal_list}
data_all = pd.DataFrame(data=data_dict)
data_all = data_all.sort_values(by=['kmer'])

###############

############### subsetting 5-mers
for kmer in kmer_subset:
    if '*' in kmer:
        pos = kmer.index('*')
        if pos == 0:
            kmer_subset.append('A'+kmer[pos+1:])
            kmer_subset.append('C'+kmer[pos+1:])
            kmer_subset.append('G'+kmer[pos+1:])
            kmer_subset.append('T'+kmer[pos+1:])
        elif pos == len(kmer)-1:
            kmer_subset.append(kmer[0:pos]+'A')
            kmer_subset.append(kmer[0:pos]+'C')
            kmer_subset.append(kmer[0:pos]+'G')
            kmer_subset.append(kmer[0:pos]+'T')
        else:
            kmer_subset.append(kmer[0:pos]+'A'+kmer[pos+1:])
            kmer_subset.append(kmer[0:pos]+'C'+kmer[pos+1:])
            kmer_subset.append(kmer[0:pos]+'G'+kmer[pos+1:])
            kmer_subset.append(kmer[0:pos]+'T'+kmer[pos+1:])
             
new_kmer_subset = [kmer for kmer in kmer_subset if ('*' in kmer) == False]
data_all = data_all[data_all['kmer'].isin(new_kmer_subset)]
print(new_kmer_subset)
###############


## For storing boxplots data
median_dist = []
min_dist = []
max_dist = []
q1_dist = []
q3_dist = []
## For storing violin plot data
mean_dist_all = []


## calculating distance matrices for t-SNE
series_data = []
for i in data_all.index:
    val = map(int,data_all.loc[i,'values'].split("_"))
    series_data.append(np.array(list(val), dtype=np.double))

ds = dtw.distance_matrix_fast(series_data,show_progress=True, penalty=penalty)

ds[np.tril_indices(ds.shape[0],k=-1)] = ds.T[np.tril_indices(ds.shape[0],k=-1)]
np.fill_diagonal(ds,0)
ds = pd.DataFrame(ds)
ds.index = data_all['kmer']
ds.columns = data_all['kmer']

os.makedirs(out_folder+'/distance_matrices/')
os.makedirs(out_folder+'/raw_signal/')
for kmer_row in data_all['kmer'].unique():
    for kmer_column in data_all['kmer'].unique():
				
        current_ds = ds.loc[[kmer_row], [kmer_column]]
        current_ds.columns = range(0,len(current_ds.columns))

        if kmer_row == kmer_column:
            y = current_ds.mean(axis = 1)
            median_dist.append(np.median(y))
            min_dist.append(np.min(y))
            max_dist.append(np.max(y))
            q1_dist.append(np.percentile(y,25))
            q3_dist.append(np.percentile(y,75))
            mean_dist_all.append(y)		
						
		## output location for distance matrices
        outdir = out_folder+'/distance_matrices/'+ kmer_row + '/'					
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        filename = kmer_row + '_' + kmer_column + '.csv'
        fullname = os.path.join(outdir, filename)
        
        current_ds.to_csv(fullname, index=False)

## creating a dataset for violinplot
mean_dist_all = list(itertools.chain.from_iterable(mean_dist_all))  ## combine list of lists into one list
violin_dict = dict([('kmer',list(data_all.kmer)), ('ave_distance',mean_dist_all)])
violin_df = pd.DataFrame.from_dict(violin_dict)
violin_df.to_csv(out_folder+'/violin_data.csv', index=False)   ## output directory for violinplot data


## creating a dataset for boxplot
kmer_list = data_all['kmer'].unique().tolist()
boxplot_dict = dict([('kmer', kmer_list),
                          ('min', min_dist),
                          ('max', max_dist),
                          ('median', median_dist),
                          ('q1', q1_dist),
                          ('q3', q3_dist)    
                    ])

boxplot_df = pd.DataFrame.from_dict(boxplot_dict)
boxplot_df.to_csv(out_folder+'/boxplot_data.csv', index=False)   ## output directory for boxplot data


## creating datasets for drawing raw signals in d3
data_all.kmer.unique()
for kmer in data_all.kmer.unique():
    data_kmer = data_all.loc[data_all['kmer']==kmer,:]
    filename = out_folder+'/raw_signal/'+kmer+'_signal.csv'    ## output directory for raw signal data
    data_kmer.to_csv(filename, index=False)






