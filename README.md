# Sequoia

A web interface tool for visualizing similarities of nanopore sequencing data

To use the tool, two separated steps must be executed. The first step is to run a python script via command line which will extract signal data from the input data and perform necessary computations. The second step is to open the visualization interface using the computational output from the first step.	

The visualization demo with pre-loaded m5C data can be found at https://khreda.com/sequoia_demo/index.html.

## Install

``` git clone https://github.com/dnonatar/Sequoia.git ```

## Backend Computation 

Python 3 and the following libraries are required.
* pandas (0.25 or newer)
* numpy (1.17 or newer)
* h5py (2.9 or newer)
* dtaidistance (1.2.3 or newer)

Please refer to https://dtaidistance.readthedocs.io/en/latest/ for further details on the dtaidistance library.

#### Usage

``` python backend_computation.py [--f1 F1] [--f2 F2] [--p P] [--k K] [--s S] [--o O] ```

Parameters 
* F1 : Fast5 file for modified signals (required)
* F2 : Fast5 file for unmodified signals (optional)
* P : Dynamic time warping penalty (optional, default = 0)
* K : CSV file with a list of 5-mer of interest (required)
* S : Sample size for each 5-mer (optional, default = 100)
* O : Output directory (required)


Example:

``` python backend_computation.py --f1 test.fast5 --k kmer_list.csv --o myputput ```

kmer_list.csv would contains a list of 5-mers in a single column (no header). A symbol (*) is used for either A, C, T, G.

For example, if the list below is provided, the code would  extract all 5-mers with AAC as the first three letters, and also CCCCC. 


AAC**   
CCCCC 


After running backend_computation.py, a new folder containing subfolders and files to be used for the visualization will be created.

## Visualization 

To start the web interface for visualization, investigate to the output directory where the output folder is placed and use the following command to start a local web server.

``` python -m http.server ```

Then, open the url and select index.html.

Below is how the landing page looks like. In the 'Type Input Folder' textbox, type in the name of the output folder you would like to explore. For instance, type 'myoutput' and click 'Choose' if you would like to use the output generated in the example above. 

![image 1](/images/landing_page.png)



