# Sequoia

A web interface tool that can visualize the similarity of nanopore sequencing data

In order to use the tool, two separated steps must be executed.  This is because the tool comprises of two major parts. The first part is the dynamic time warping algorithm that measures the similarities between signals. The second part is the visualization based on the data generated from the first part.	

## Overview

1.) Backend Computation (Python)

required python library
* pandas
* numpy
* os 
* sys
* dtaidistance (python3 only)

please refer to https://dtaidistance.readthedocs.io/en/latest/ for further details on the dtaidistance library.

Input format (csv file)

read_ID | kmer | values
------------ | ------------- | ------------
...... | AACAA | 628_671_629_658_673_718_717_698_691_700....
...... | AACCC | 628_671_629_658_673_718_717_698_691_700....
...... | ...... | ......
 
##### Usage

``` python  backend_computation.py arg1 arg2 ```

Parameters 
* arg1 : input file (need to show / explain the input format)
* arg2 :  DTW penalty

Run backend_computation.py in the directory where you want data stored. 

example:

``` python backend_computation.py data.csv 0 ```

After the script is run, a new folder named data_0 will be created in your directory. The ‘0’ in data_0 comes from the penalty choice. 

2.) Visualization (D3.js)

First, download the web server for Chrome from here
https://chrome.google.com/webstore/detail/web-server-for-chrome/ofhbbkphhbklhfoeikjpcbhemlocgigb?hl=en

To start the visualization:
1. Place current.html into the same directory as the data folder from the previous step
2. Open the Web Server from Chrome
3. In the Web Server, click at CHOOSE FOLDER and select the directory that you place the data and current.html
4. open the Web Server URL and click curent.html


![image 1]('./image/first_screen.png')

