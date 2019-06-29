#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
from dtaidistance import dtw
import itertools

## input data (first commandline argument)
input_directory = str(sys.argv[1])
data_all = pd.read_csv(input_directory,header=0)

## penalty
penalty = int(sys.argv[2])

## output folder
out_folder = str(sys.argv[3])


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
				
        current_ds = ds.loc[kmer_row, kmer_column]
        current_ds.columns = range(0,len(current_ds))

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






