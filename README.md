# Tfit
Transcription fit (or Tfit) implements a finite mixture model to identify signatures of RNA Polymerase II (RNAPII) activity through the identification of bidirectional or divergent transcription in nascent transcription assays such as Global Run-On (GRO) and Precision Run-on (PRO) followed by sequencing data. Tfit is separated into two modules:

1. `bidir` which generates calls for preliminary regions of interest 
2. `model` which attempts to model RNAPII behavior for the preliminary outputs and make refined predictions for activity

Output from these modules can be imported into any genome browser. An example bed file is shown below.  

![Screenshot of data output](images/Example_Snapshot.png)

The `bidir` module will compute local likelihood statistics given a fixed template mixture model (estimated either from promoter associated transcription) or specified explicitly by the user (the former is encouraged). This method is fast and will finish in about 30 minutes on a single node single CPU machine. The output will be a bed file corresponding to areas of possible bidirectional transcription. This output is discussed more in later sections.

The `model` module will compute full MLE estimates of the mixture model at user specified regions of the genome provided as a bed file. This bed file may be the output from the bidir module. It is recommended that users fit MLE estimates to the output of the bidir module as this will decrease false positives. Three files will output from this module:

* Full MLE output description (\<jobName\>.tsv) containing information about pausing, spreading, and strand probability of each predicted region
* A BED7 file contianing predicted regions of activity (\<outFile\>.bed); the first three columns are the standard chr, start, end, columns 4-7 contain model information for each region
* Stanard output with detailed job run information (\<jobName\>.log)

This module is much more computationally expensive than `bidir` and can take up to 10-12 hours depending on sample size.

### System Requirements and Makefile
Transcription Fit (TFit) is written in the C/C++ programming language that requires GNU compilers >4.7.3 and =< 7.1.0. Tfit uses the popular openMPI framework to perform massive parallelization via multithreading on multiple core, single node systems or multiple core, multiple node compute clusters and therefore also has an MPI dependency requiring a verion >3.0. After cloning this repo, please change directory into /where/you/clone/this/repo/Tfit/src/ and run make clean (removing any existing binaries) followed by make.


```
$ cd  /where/you/clone/this/repo/Tfit/src/

$ make clean

Successfully removed binaries!

$ make
```

If the program compiles successfully you should see the following output.

```
=========================================
GCC version: 7.1.0
main                  : done
load                  : done
split                 : done
model                 : done
across_segments       : done
template_matching     : done
read_in_parameters    : done
model_selection       : done
error_stdo_logging    : done
MPI_comm              : done
density_profiler      : done
bootstrap             : done
biPWD_main            : done
model_main            : done
select_main           : done
FDR                   : done
BIC                   : done
ParamWrapper          : done
old_template_matching : done
linking               : done
=========================================
Successfully compiled!
```

If your program, did not compile properly it is likely that you do not have the correct dependencies. The three significant dependencies are listed below. 

1.**c++11** (this ships with most recent versions of GCC, please visit https://gcc.gnu.org/install/)

2. **openmp** (this ships with most recent versions of GCC, please visit https://gcc.gnu.org/install/)

3. **MPI** (this needs to installed and configured and serves as a wrapper for GCC, please visit https://www.open-mpi.org/faq/)

In short, the make file requires the path to mpic++ (install and config openMPI) to be in your PATH.

### Utilizing openMP and MPI
Tfit is written using openMP and MPI to perform massive parallelization. If your institution has a large compute cluster, than Tfit will operate well across multiple cores and nodes. The `bidir` module runs relatively quickly, so this is only recommended for the `model` module. To invoke 4 MPI processes (i.e. run across 4 nodes) run:

```
$ mpirun -np 4 Tfit model [arguments]
```

This highly system dependent, but openMPI is a well maintained library with lots of resources to help aid in getting Tfit up and running on your compute cluster, please consult https://www.open-mpi.org/faq/ for further reference. 


### Running and Installing Through Docker Container
Running Tfit via a Docker container requires that Docker is installed and running. In Windows and Mac OSX you must run the Docker/Tfit script from the Docker Quickstart Terminal (https://www.docker.com/products/docker-toolbox).

To run Tfit via docker, simply run the Docker/Tfit script with the same arguments and options as you would use when running the standard Tfit program. This script will check if you have the biofrontiers/fstitch_tfit:latest image and download it if you do not. It will then pass the arguments you provide to the Tfit program inside the container and output the results and logs to the location you specify.

Advanced users may want to build their own Docker images using the provided Dockerfile. To do so, run the following command from within the Docker directory:

docker build -t [image name] .

And alter the Docker/Tfit script, changing the REPOSITORY and TAG variables so that they match your built image.

## Running Tfit
The basic argument to run TFit is as follows:

```
$ export OMP_NUM_THREADS=16

$ Tfit [module] [arguments]
```
***In order to properly allocate the number of cores used in MPI, you must include the first command above in bash***! 

To view the usage statement, use following standard arguments (*must be compiled*):
```
$ /src/Tfit -h 

$ /src/Tfit --help
```

### Coverage File (bedGraph) Format Requirements

The provided coverage file must be in bedGraph format (See the UCSC description<sup>2</sup>) which is a BED4 file where the fourth column represents coverage over the annotated start/end positions. For example:

```
chr    start    end    coverage
1      0        100     3
1      107      117     1
```

***IMPORTANT:*** Your bedGraph file **should not contain 0 values** and should be **non-normalized**. FStitch performs an internal normalization.

There are two main tools for generating bedGraph coverage files from BAM files, deepTools and BEDTools. By default, deepTools bamCoverage will have discrete bins (50bp) and will therefore calculate average coverage over regions, rather that contigs of regions with equal read coverage, and "smooth" the data. While this is not a problem for visualizaiton at smaller bins, it will conflict with normalization. Therefore, we recommend using default BEDtools<sup>3</sup> genomecov settings:
    
```
$ bedtools genomecov -ibam <file.bam> -g <file.bedGraph> -bg -s <+/->
```

Specifying five prime (-5 argument) in the “genomecov” may allow for cleaner annotations however unspecified five prime bed works just fine as well. 

Often times, it is useful to merge the positive and negative strands for conversion to TDF if you are using IGV as your genome browser so you can view both strands on one track. While Tfit will accept separated positive/negative strand files, you can provide this concatenated pos/neg bedGraph file as long as it is formatted with appropriate integers in the fourth column (i.e. -\<integer> for a negative strand coverage region). If two separate pos/neg strand files are provided, Tfit will simply concatenate them internally. An example of how to generate an appropriately formatted concatenated pos/neg bedGraph file is as follows:

```
awk 'BEGIN{FS=OFS="\t"} {$4=-$4}1' ROOTNAME.neg.bedGraph \
 > ROOTNAME.neg.formatted.bedGraph

cat \
 ROOTNAME.pos.bedGraph \
 <(grep -v '^@' ROOTNAME.neg.formatted.bedGraph) \
 | sortBed \
 > ROOTNAME.sort.cat.bedGraph
```

Note that the last command, sortBed, will require that BEDTools be specified in your PATH. Sorting is not required for Tfit, but is good practice as it will expedite daa processing and is required for conversion to TDF (for visualization in IGV) if so desired.

## Tfit bidir
The bidir module scans across the genome for areas resembling bidirectional transcription by comparing a fixed template mixture model (user provided parameters or parameters estimated from promoter regions) to a noise model (uniform distribution) by a Likelihood ratio score (LLR).

In short, the `bidir` will output putitive regions of RNAPII transcriptional activity. ***This is only intended as a pre-filter to the `model` module!***

**Required Arguments:**

| Flag                         | Type | Description |
|------------------------------|------|-------------|
| -f     --forward             | \<FILE.pos.bedGraph>      | Forward (positive) strand bedGraph file (***Required only if -bg not specified***)
| -r     --reverse             | \<FILE.neg.bedGraph>      | Reverse (negative) strand bedGraph file (***Required only if -bg not specified***)
| -bg    --bedgraph            | \<FILE.cat.bedGraph>      | Concatenated pos/neg bedGraph file (***Required only if -f and -r not specified***) 
| -N     --jobName             | \<MYJOB>                  | Job name for log output (no extensions).
| -o     --output              | \<PRELIMHITS.bed>         | User specified output name and directory for BED file (.bed extension)
| -l     --logOut              | \</project/logs/>         | Log file output directory. Will contain stdout (.log) with run specifications with --jobName as the rootname.


**Optional Arguments:**

| Flag                | Type | Description |
|---------------------|------|-------------|
| -s     --segment    | \<SEGFILE.bed>            | BED file that specifies sample regions of interest (e.g. FStitch output)
| -chr   --chromosome | \<chrX>                   | Run bidir only on the specified chromosome. Default = all
| -n     --threads    | \<integer>                | Number of threads to run the job on; 16 recommended. Default=1

The recommended segment file above is the output obtained through FStitch (see FStitch @ https://github.com/Dowell-Lab/FStitch). This will filter regions of interest to only the areas of active transcription. Be aware that this file should not be too segmented (i.e. "active" regions will likely contain some 0 values). For more details, see the FStich README. Providing this segment file will signficantly reduce runtime.

Putting these arguments together, an example command is as follows:

```
$ export OMP_NUM_THREADS=16

$ Tfit bidir \
    -bg <FILE.cat.bedGraph> \
    -N <MYJOB> \
    -o <PRELIMHITS.bed> \
    -l </project/logs> \
    -s <SEGFILE.bed> \
    -n 16
```

The output BED file will contain 5 columns: chr, start, end, identifier (ME_X), and the MLE estimates for the predicted region.

```
chr1    11179   18104   PRELIM_0
chr1    18254   22704   PRELIM_1
chr1    26704   32079   PRELIM_2
chr1    33329   37529   PRELIM_3
chr1    97904   98879   PRELIM_4
chr1    135479  141229  PRELIM_5
chr1    179504  193204  PRELIM_6
chr1    193354  204254  PRELIM_7
chr1    257854  264954  PRELIM_8
chr1    347779  352229  PRELIM_9
chr1    358204  360304  PRELIM_10
```

## Tfit model
Unlike the `bidir` module which utilizes an average or template version of the mixture model to scan the entire genome quickly, the `model` module will attempt to find (by maximum likelihood estimation, MLE) the best set of parameters (sigma,lambda, pi, w) on a per region of interest basis. Such a routine is especially valuable if it is believed that pausing or strand bias is changing following experimental perturbation. In addition, running the model module on the PRELIMHITS.bed file will greatly decrease the number of false positives as the MLE estimates will more accurately reflect the underlying structure of the region of interest rather than a static template model.  In short, MLE estimates are computed by the EM algorithm which is a convergent optimization method found commonly in gaussian mixture modeling and cluster analysis. 

The `model` module is therefore meant as an extension of the `bidir` module and as such should always be run in succession.

**Required Arguments:**

| Flag                         | Type | Description |
|------------------------------|------|-------------|
| -f     --forward             | \<FILE.pos.bedGraph>      | Forward (positive) strand bedGraph file (***Required only if -bg not specified***)
| -r     --reverse             | \<FILE.neg.bedGraph>      | Reverse (negative) strand bedGraph file (***Required only if -bg not specified***)
| -bg    --bedgraph            | \<FILE.cat.bedGraph>      | Concatenated pos/neg bedGraph file (***Required only if -f and -r not specified***) 
| -N     --jobName             | \<MYJOB>                  | Job name for log output (no extensions).
| -o     --output              | \<FINALHITS.bed>          | User specified output name and directory for BED file (.bed extension)
| -l     --logOut              | \</project/logs/>         | Log file output directory. Will contain stdout (.log) with run specifications with --jobName as the rootname.

**Optional Arguments:**

| Flag                | Type | Description |
|---------------------|------|-------------|
| -bd    --bidirs     | \<BIDIRS.bed>             | BED file with either TSS's or annotated bidirectionals from your sample
| -s     --segment    | \<PRELIMHITS.bed>         | BED file that specifies sample regions of interest (e.g. FStitch output)
| -chr   --chromosome | \<chrX>                   | Run bidir only on the specified chromosome. Default = all
| -n     --threads    | \<integer>                | Number of threads to run the job on; 16 recommended. Default=1

***IMPORTANT***: The -s --segment <PRELIMHITS.bed> is listed as an optional argument, but is highly recommended to reduce false positives and reduce overall runtime for bidirectional modeling. That said, Tfit `model` will run without this argument specified across the entire genome and may be a useful diagnostic tool.

Alternatively, you could specify FStitch segment output (see https://github.com/Dowell-Lab/FStitch for more details) to model over entire gene bodies.

Another important optional options above is the -bd || --bidirs flag. This flag will attempt to fine-tune the default modeling parameters to the regions of interest that it is provided. This can be either a RefSeq type file of transcriptional start sites (which will bias your results towards robust bidirectionals) or a user-provided annotation BED file similar to that used in FStitch (see https://github.com/Dowell-Lab/FStitch for more details). Only a BED3 is required in Tfit (chr, start, end). This is useful if your data is of lower complexity than what is ideal for Tfit modeling.

Putting these arguments together, the following is an example run using MPI on 4 nodes:

```
$ export OMP_NUM_THREADS=16

$ mpirun -np 4 Tfit model \
    -bg <FILE.cat.bedGraph> \
    -N <MYJOB> \
    -o <FINALHITS.bed> \
    -l </project/logs> \
    -s <PRELIMHITS.bed> \
    -bd <BIDIRS.bed> \
    -n 16
```

which will attempt to model the preliminary calls output by invoking the `bidir` module and adjust the modeling parameters to the user-specified bidirectional file. The -np flag specifies the number of nodes while the number provided OMP_NUM_THREADS and -n should match. Together, this run will use 64 threads/cores across 4 nodes.

After the model module has finished, Tfit will output two files in the user specified output directories: 
1. Full model estimates: \<jobName>\_models_MLE.tsv
2. Predicted regions of bidirectionals: FINALHITS.bed

The first file, \<jobName>\_models_MLE.tsv, gives a detailed account for each interval of interest from -k input parameters, a list for each finite mixture model fits -mink to -maxk, the log-likelihood score, and MLE estimates for the center of the bidirectional transcript (mu), variance in mu (sigma), entry length (lambda), strand bias (pi), distance between bidirectional peaks (foot print) and all the associated mixture weights (basically w). In addition, stats on the number of reads over that interval etc. This can be used for manual Bayesian Information Criterion (BIC) calculations. An example output of this file is below with outputs using 4 model fits (0, 1, 2, or 3 bidirectionals for each input region given in the PRELIMHITS.bed):

```
>ME_0|chr1:9178-20105|1700.000000,38403.000000
~0,-195273.185480
~1,-194868.534602       9178.002509     70.710678       100.000098      0.500000        0.000000        0.000697,0.996517,0.000697      20105.000000    9178.000000
~2,-195006.021454       9178.002507,9178.002507 70.710678,70.710678     100.000098,100.000098   0.500000,0.500000       0.000000,0.000000       0.000694,0.033837,0.000694|0.000694,0.959224,0.000694   20105.000000,20105.000000       9178.000000,9178.000000
~3,-195142.583693       9178.002488,9178.002488,9178.002488     70.710678,70.710678,70.710678   100.000097,100.000097,100.000097        0.500000,0.500000,0.500000      0.000000,0.000000,0.000000      0.000691,0.216681,0.000691|0.000691,0.556277,0.000691|0.000691,0.216681,0.000691        20105.000000,20105.000000,20105.000000  9178.000000,9178.000000,9178.000000

```

The second file, FINALHITS.bed, provides a new bed file, where the center of each bed interval corresponds to the center of the bidirectional peak and the width corresponds the estimated standard deviation around that estimate following BIC model comparison. This again can be used for downstream analysis and comprises the most accurate set of bidirectional prediction centers that Tfit can currently offer. An example output of this file is given below.

```
chr1    28764   29839   BIDIR_2|1.191337   29302.207103    0.362018        425.818173
chr1    96496   96833   BIDIR_4|1.697595   96665.182490    0.678465        91.942768
chr1    197298  201275  BIDIR_7|1.152818   199287.083530   0.094496        1935.048257
chr1    199335  200287  BIDIR_7|1.152818   199811.534585   0.248295        379.689945
chr1    259282  260470  BIDIR_8|1.171639   259876.693769   0.907025        571.777800
chr1    368477  368553  BIDIR_12|4.583669  368515.378599   0.487298        20.000000
chr1    376808  377436  BIDIR_13|1.433726  377122.645102   0.380428        275.047612
chr1    377247  377484  BIDIR_13|1.433726  377366.182346   0.436560        59.999043
```

The first three columns are the same as most standard BED files (chr, start, stop) for regions of interest. The fourth column is the annotation for the predicted bidirectional. The integer following BIDIR_ tells you which prediction region from PRELIMHITS.bed that that bidirectional was called from. The floating point following BIDIR_X|\<FLOAT> is the BIC model estimate value. Columns 5, 6, and 7 are mu, w, and lambda respectively (polymerase loading, pausing ratio, and exponential decay).


## Questions and Comments
For questions, issues, and updates, see Dowell Lab @ https://github.com/Dowell-Lab/Tfit



## References
If you find Tfit useful in your analysis please cite:

J. Azofeifa and R. Dowell A generative model for the behavior of RNA polymerase Bioinformatics 2016 

http://bioinformatics.oxfordjournals.org/content/early/2016/09/23/bioinformatics.btw599.abstract


