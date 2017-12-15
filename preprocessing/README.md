## Data preprocess protocal
---

### an example

```shell


#use the index sequence to get raw.data
 ./sccpipe -I INPUT/  -O OUTPUT/

#cut off the sccpipe when the first step is finished 
#filiter the sequence including the CGACTCACTACAGGG,TCGGTGACACGATCG,TTTTTTTTTTTT 
./bbmap/bbduk2.sh in=INPUT/Undetermined_S0_L007_R1_001.fastq in2=INPUT/Undetermined_S0_L007_R2_001.fastq outm=INPUT/Undetermined_S0_L007_c1_R1_001.fastq outm2=INPUT/Undetermined_S0_L007_c1_R2_001.fastq fliteral=CGACTCACTACAGGG k=15 skipr2=t hdist=3 -Xmx58g

ggj@ggj-HP-Z840-Workstation:~$ /home/ggj/Documents/Genome-tools/bbmap/bbduk2.sh in=INPUT/Undetermined_S0_L007_c1_R1_001.fastq in2=INPUT/Undetermined_S0_L007_c1_R2_001.fastq outm=INPUT/Undetermined_S0_L007_c2_R1_001.fastq outm2=INPUT/Undetermined_S0_L007_c2_R2_001.fastq fliteral=TCGGTGACACGATCG k=15 skipr2=t hdist=3 -Xmx58g

./bbmap/bbduk2.sh in=INPUT/Undetermined_S0_L007_c2_R1_001.fastq in2=INPUT/Undetermined_S0_L007_c2_R2_001.fastq outm=INPUT/Undetermined_S0_L007_c3_R1_001.fastq outm2=INPUT/Undetermined_S0_L007_c3_R2_001.fastq fliteral=TTTTTTTTTTTT k=12 skipr2=t hdist=3 -Xmx58g

#get the bam file
java -Xmx58g -jar Drop-seq_tools-1.12/Drop-seq_tools-1.12/3rdParty/picard/picard.jar  FastqToSam F1=INPUT/Undetermined_S0_L007_c3_R1_001.fastq F2=INPUT/Undetermined_S0_L007_c3_R2_001.fastq  O=INPUT/20171020Lung.bam QUALITY_FORMAT=Standard SAMPLE_NAME=sample_name
#get dge
./Drop-seq_tools-1.12/Drop-seq_alignment_DS.sh  -g STAR_Reference_Mouse/genomeDir -r /home/ggj/Documents/STAR_Reference_Mouse/Mus_musculus.GRCm38.88.fa -d  Drop-seq_tools-1.12/Drop-seq_tools-1.12/ -o INPUT/ -t INPUT/ -s STAR-2.5.2a/source/STAR  INPUT/20171020Lung.bam



```
