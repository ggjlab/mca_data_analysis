# sccpipe
A pipeline for scRNA sequence analysis -- for c1 only 

```shell
Usage
sccpipe -I <Input File Directory> -O <Output File Directory> -B <Barcode file in text>
```
> For the first time, the scc.conf should be filled 

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
Core_use=<>
#number of mismatch (default 1)
Missm=1
#Mode for polyA trimming(default 1)
Mode=1
```

## Dependencies (make sure the dependencies are in your enviroment)
* ~~[paralell](https://www.gnu.org/software/parallel/) -- (for parallel computing)~~
* [STAR](https://github.com/alexdobin/STAR)    -- Aligment Tools 
* ~~[cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html) -- A python lib for adapter and ploy A/T trimming(pip install)~~
* ~~[BBmap](https://sourceforge.net/projects/bbmap/)    -- we use the subfunction bbduk2 to filter the barcode~~
* ~~[HTSeq](http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)    -- for quantify the gene expression(pip install)~~
* [featureCounts](http://subread.sourceforge.net/)  --for quantify the gene expression(pip install
* [scater](http://bioconductor.org/packages/release/bioc/html/scater.html)   -- (Bioconductor R package -for constructing SCESET)

## Installation (Test on Ubuntu)
```shell
# install parallel
# sudo apt-get install parallel  

# install STAR 

#install cutadpat (if pip is not in place, you need to install pip first)
#sudo pip install cutadpat

# install bbmap (download )
#tar zxvf BBMap_36.84.tar.gz
#cd bbmap/jni
#make -f makefile.linux
#cd ../../

# install HTSeq
#sudo pip install HTSeq

# install featureCounts



# install scatter
# enter R
source("https://bioconductor.org/biocLite.R")
biocLite("scater")
# last step 
git clone https://github.com/bioeauty/sccpipe.git
or downlaod from git.mclab.me
```

### Add to environment

```shell
#open bashrc file by ` vim ~/.bashrc` and add the following to the bottom of bashrc
export PATH=<the full path where you put the pipeline>/sccpipe:$PATH
export PATH=<the full path where you put the pipeline>/sccpipe/scripts:$PATH
#then write and quit
source ~/.bashrc
```
