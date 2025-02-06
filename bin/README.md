# Overview
This script to run convergent variant finder (ConVarFinder) analysis to idenify convergent single sequence variants of target species which share a convergent trait you are interested in. Based on this script and demo data, you can find ConSAVs of avian vocal learners mutually exclusive to amino acids of vocal non-learners.
## Input Files
It needs multiple peptide sequence aglinments as input files (as [fasta](https://en.wikipedia.org/wiki/FASTA_format) format).
## Variables in Script
-l=    strings of list of target species as comma seperated (csv)

-ttf=  input file of templete tree

-tfmt= file format of binary tree with ancestral nodes

-sfmt= file format of input with ancestral reconstructions

-ip=   input path

-op=   output path  

## Output Files
It generates a text file (.txt) as a output with a summary of CSAV analysis and the full amino acids of whole species list at identified CSAV sites.
- - -

# System Requirments
## Hardware requirements
This script requires only a standard computer with enough RAM to support the in-memory operations. It was developed and tested in a standard computer (AMD Ryzen 5 2400G 3.6GHz and 16GB RAM)

## Software requirements
### OS Requirements
This script is supported for *Windows OS* and *Linux*. The script has been tested on the following systems:
* Windows 10 Education
* Linux: Ubuntu 18.04.3 LTS

## Python Dependencies
This script mainly depends on the Python scientific stack.

    numpy
    scipy

## Running time with the demo data
* under 1 sec
- - -

# Running the script
* (option1) If you want to clone the github link:
<pre>
<code>
git clone https://github.com/chulbioinfo/ConVarFinder/
cd ConVarFinder/bin/
Python ConVarFinder.py [options]
</code>
</pre>


- - -



# License
This project is covered under the Apache 2.0 License.
