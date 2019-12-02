# Sequoia

A web interface tool for visualizing similarities of nanopore sequencing data

In order to use the tool, two separated steps must be executed because the tool comprises of two major parts. The first part is a pipeline of python scripts that extract signal data form the input data, compute similarity matrix across all signals using the dynamic time warping algorithm, and perform summary statistics for visualization. The second part is the visualization based on the data generated from the first part.	

## Install

``` git clone https://<i></i>github.com/dnonatar/Sequoia.git ```

## Backend Computation 

Python 3 and the following libraries are required.
* pandas
* numpy
* h5py
* dtaidistance 

Please refer to https://dtaidistance.readthedocs.io/en/latest/ for further details on the dtaidistance library.

#### Usage

``` python backend_computation.py arg1 arg2 arg3 arg4 ```

Parameters 
* arg1 : directory of a fast5 input file 
* arg2 : dynamic time warping penalty (for the simplest case, set equal to 0)
* arg3 : csv file with a list of 5-mer of interest
* arg4 : output directory

After running backend_computation.py, a new folder containing subfolders and files necessary to generate the visualization will be created.

Example:

``` python ./backend_computation.py data.csv 0 kmer_list.csv ./myoutput ```

kmer_list.csv would contains a list of 5-mers in a single column (no header). A symbol (*) is used for either A, C, T, G.

For example, if the list below is provided, the code would  extract all 5-mers with AAC as the first three letters, and also CCCCC. 

|------|
|AAC** |  
|CCCCC |  
  

## Visualization 

To start the web interface for visualization, investigate to the output directory where the output folder is placed and use the following command to start a local web server.

``` python -m http.server ```

Then, open the url and select index.html.

Below is how the landing page looks like. In the 'Type Input Folder' textbox, type in the name of the output folder you would like to explore. For instance, type 'myoutput' if you would like to use the output generated in the example above. Once click 'Choose', the boxplots for each of the 5-mers will show up. 

![image 1](/images/landing_page.png)

## Demo
This is the visualization demo with pre-loaded data.

http://research.vis.ninja/sequoia/

