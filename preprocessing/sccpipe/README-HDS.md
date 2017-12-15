# sccpipe
A pipeline for scRNA sequence analysis -- for c1 only 

```shell
# install sccpipe 
git clone http://git.mclab.me:10080/Yincong/sccpipe

Usage
sccpipe -I <Input File Directory> -O <Output File Directory> 
```
# Input File Directory contain Fastq files; Output File Directory should be a new created directory

> For the first time, the sccpipe.conf should be filled

```shell
### sccpipe.conf
################################ scc.conf #######################################
#																				#
#	 			Configure the basic parameters for sccpipe						#
#																				#
#################################################################################

#The full path of your genome index
Indexpath=<>

#The full Path of your GTF file 
GTFFile=<>

#the CPU Thread  default = (maximum CPU cores -4)
Core_use=20

#The full path of Barcode file in the format like "ATATAATA	COL", split by TAB
barcodefile=/home/ggj/sccpipe/barcode.txt

#the CPU Thread  default = (maximum CPU cores -4)
Core_use=20

#number of mismatch (default 1) -- to 20 cols
Missm=1

#Mode for polyA trimming(default 1) -- (only 3 mode)
Mode=1

# define the polyA length(default 6)
polyAcounts=6

#################################################################################

## Dependencies (make sure the dependencies are in your enviroment)
* [STAR](https://github.com/alexdobin/STAR)    -- Aligment Tools 
* [featureCounts](http://subread.sourceforge.net/)  --for quantify the gene expression(pip install
* [scater](http://bioconductor.org/packages/release/bioc/html/scater.html)   -- (Bioconductor R package -for constructing SCESET)

## Installation (Test on Ubuntu)

# install STAR 

# install featureCounts

# install scatter
enter R
source("https://bioconductor.org/biocLite.R")
biocLite("scater")

```

### Add to environment
```shell
#open bashrc file by ` vim ~/.bashrc` 
add the following to the bottom of bashrc
export PATH=<the full path where you put the pipeline>/sccpipe:$PATH
export PATH=<the full path where you put the pipeline>/sccpipe/scripts:$PATH
export PATH=<the full path where you put the pipeline>/STAR:$PATH
export PATH=<the full path where you put the pipeline>/featurecouts:$PATH

#then write and quit: wq
source ~/.bashrc
```
