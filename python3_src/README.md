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

If you wish to, there is a "harcoded" version written in to the main.py file. If you edit this to include
all relevant files, and paths. All that is needed to run is:

` $python3 main.py`

This will as if the flags formatData, and FstitchSingleIsoform where given.

So it is equal to this:

``` 
$python3 main.py formatData FStitchSingleIsoform -ref annotation_files/hg19_TSS.bed \\
-ffs examples/SRR1105737.pos.sorted.BedGraph -rfs examples/SRR1105737.minus.sorted.BedGraph \\
-fbg examples/SRR1105737.pos.sorted.BedGraph -rbg examples/SRR1105737.minus.sorted.BedGraph - \\
wo examples/single_isoform_FStitch.tsv
```
