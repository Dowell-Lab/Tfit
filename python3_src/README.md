## Python3 Tfit

This is the specific README for the python3 version of Tfit.
For the general/c++ README, that can be found in the base directory.

# Usage:

To use Tfit with a slurm system, a sample script has been included as python3_tfit_sample
This script essentially loads all dependant modules needed, and runs using example files. 
Edits would need to be made to the source directories and the files that are being used.

To run this code direcectly in the command line:

` $python3 main.py runModel -i /path/to/BedGraph -wo /path/to/out/dir -sc hg19 `

This runs "model across" on the BedGraph file, with refrence genome hg19


