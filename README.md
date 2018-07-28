# GMIM: Genome-wide multi-track integrator using minimum Bayes' Factor
GMIM is a Java-based commandline tool that integrates multiple genome-wide NGS data via the minimum Bayes' Factor from density of sequence reads mapped on the genome. To integrate multiple NGS data, homo-technique replicates or hetero-technique NGS data can be integrated into the single data track based on the complements of the minimum Bayes' factor (*cMBF*) (*cMBF = 1 - exp(-Z^2/2)*) (range: 0-1) calculated from signal-to-background ratios as a function of genomic position. 

This is **NOT** a peak-calling algorithm as it does not specifically locate peaks (though you can use a peak-calling tools on the output!). 

### Installation
You can download the .jar file from the Epithelial Systems Biology Laboratory website [here](https://esbl.nhlbi.nih.gov/Bioinformatic%20Tools.htm). Move the .jar file to your working folder, and follow the run steps below!

**Requires: Java 1.7**

## Part 1: Calculation of the complement of minimum Bayes' Factor
The first GMIM program **calculates** the cMBF for each position for later integration. It first splits a given genome/chromosome file into a directory of separate chromosomes to facilitate parallelizing computation. For each chromosome file, the program uses a "moving window" to calculate the cMBF. 

### Input file
GMIM takes as input a *sorted* BED file of sequence read counts (generated standard tools, such as *samtools* (depth), *bedtools* (coverage, genomecov), etc.)
The input file (bed file format) should contain the following (whitespace/tab separated) columns:
`chr  start end read_count`

The commandline options for changing the parameters are as below:
```
-h,--help
-i,--input <arg>         [req] input file path, must be a bed file
-m,--medMult <arg>       [opt] median multiple, cannot be 0                     default: 1
-o,--output <arg>        [opt] output file base name                            default: "out"
-w,--windowSize <arg>    [req] window size for calculating cMBF, must be        default: 10000 (bp)
                          a multiple of interval size
-z,--defaultZero <arg>   [opt] default number to replace zero                   default: 0.1
```

### Run
`java -jar GMIM.jar -i [Input.bed] -w [10000] [other options]`

### Output
This part produces a directory of chromosome bedGraph files with the calculated cMBFs. 
The files have the following tab-separated format, `chr start end cMBF`, and include headers identifying the file as a bedGraph file (see the [UCSC bedGraph Format](https://genome.ucsc.edu/goldenpath/help/bedgraph.html)).
The directory also includes a final bedGraph file that concatenates all chromosomes in the order of the given file. 



## Part Two: Integration
The second part of the GMIM program **integrates** the calculated cMBFs for the different tracks by multiplying the cMBF at each position. 

### Input file
Inputs should be Part One generated files, or have the same format (ordered as `chr start end cMBF`). 
An output filename also needs to included as the **first** input. 

### Run
`java -cp GMIM.jar Integration [Integration_Output_File_Name.bedGraph] [file1.bedGraph] [file2.bedGraph]`

### Output
The output will be a bedGraph file of the same format (`chr start end cMBF`) and of the specified name. 



## Example Run(s)


## Things to Consider
- Window size: Larger is more accurate, but will take longer - 10K bp window takes about an hour with 16 threads and 16 GB of space (*should I put this here*)
- Median multiple: Generally, a larger median multiple will reduce the background level, but may also reduce signal. (*Find a browser to show*)
- Default zero: (*same as above?*)
- If you have multiple separated regions of the same chromosome, please place them into separate files to run the Driver (Part One). 
- Running GMIM when low on space may cause errors such as improper file reading (`Improper file formatting`).

