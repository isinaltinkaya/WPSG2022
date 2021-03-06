The dataset for this exercise consists of two [bamfiles](ancient_data/) these are mapped fastq files.

**Aims**:

  - Get experience with navigating files associated with high throughput sequencing data

  - Figure out which sample looks ancient and which one looks modern




### Set some environment variables and paths
```
cd workshop_materials/a07_angsd/
```
First set some paths   

    # NB this must be done every time you open a new terminal
    SAMTOOLS=/home/wpsg/software/samtools-1.15.1/samtools
    
    # Set path to folder containing data
    DATAFOLDER=github_repos/ancient_data/
    
    # Set folder to assembly
    HG38=hg38/
    
    #MAPDAMAGE?
    MAPDAMAGE=/home/wpsg/miniconda3/bin/mapDamage


Let us first validate that we setup our variables correctly


```
ls ${SAMTOOLS} ${HG38} ${DATAFOLDER} ${MAPDAMAGE}
```

### Summary statistic and conversion

We have two files called holmes.sam.gz and sherlock.sam.gz. Let us copy these files to our working directory.

```
cp ${DATAFOLDER}/* .
```

Let us try to view the uncompressed version of one of the files.

```
gunzip -c sherlock.sam.gz|less -NS
```

We are piping the output from the uncompression into the less commands and are supplying less with the arguments `-NS`, this will give us line number and we can use the arrows to go up and down left and right. You exit the commands by pressing `q`.

 1. Identify the header
 2. How many lines does the header span?

Let us see how many reads we have in these two files

```
gunzip -c sherlock.sam.gz |grep -P "^@" -v |wc -l
gunzip -c holmes.sam.gz |grep -P "^@" -v |wc -l
```
### Converting SAM to BAM

Now let us convert the SAM files to BAM files

```
${SAMTOOLS} view -b sherlock.sam.gz -o sherlock.bam
${SAMTOOLS} view -b holmes.sam.gz -o holmes.bam
```

Try to look at one of the BAM files using `samtools view`, is there are difference from uncompressing the `.sam.gz` file?

```
${SAMTOOLS} view sherlock.bam|less -SN
```

Are the reads ordered?

### Sorting BAM files

We can sort the file using `samtools sort`

```
${SAMTOOLS} sort sherlock.bam -OBAM -o sherlock.sorted.bam
${SAMTOOLS} sort holmes.bam -OBAM -o holmes.sorted.bam
```

Visually inspect that they are sorted by piping `samtools view` to `less` command

```
${SAMTOOLS} view sherlock.sorted.bam|less -SN
```

### Endogenous content

Let us find the proportion of mapped reads (endogenous content) for each file

```
${SAMTOOLS} flagstat sherlock.sorted.bam
${SAMTOOLS} flagstat holmes.sorted.bam
```
What is reported? what is the percentage of mapped reads? Which file has the largest proportion of mapped reads?

### Readlength distribution

Let us look at the distribution of read lengths (for the mapped reads)

```
${SAMTOOLS} view -F4 sherlock.sorted.bam |awk '{print length($10)}'|sort -n |uniq -c >sherlock.readlength
${SAMTOOLS} view -F4 holmes.sorted.bam |awk '{print length($10)}'|sort -n |uniq -c >holmes.readlength
```

Let us load the data into R and plot it. Alternatively you can login to Rstudio through your browser and navigate to the directory with 'setwd("/home/wpsg/workshop_materials/a07_angsd/")

```
sherlock <-read.table("sherlock.readlength")
holmes <-read.table("holmes.readlength")
par(mfrow=c(1,2))
plot(sherlock[,2],sherlock[,1],main="Sherlock",type='l',col=1,lwd=2)
plot(holmes[,2],holmes[,1],main="Holmes",type='l',col=2,lwd=2)

```

If you had problems with the R commands the result can be found [here](results/sherlock.holmes.rlen.pdf)

### Clonality
Let us also look at the proportion of duplicates (clonality)

```
${SAMTOOLS} rmdup -s sherlock.sorted.bam sherlock.sorted.rmdup.bam
${SAMTOOLS} rmdup -s holmes.sorted.bam holmes.sorted.rmdup.bam 
```

What is reported? Which ones contains the most duplicates?

### Damage signal

First we need to index the bam files
```
${SAMTOOLS} index sherlock.sorted.rmdup.bam 
${SAMTOOLS} index holmes.sorted.rmdup.bam 
```

Then we run the MapDamage software for visualizing the damage patterns

```
${MAPDAMAGE} -i sherlock.sorted.rmdup.bam -r ${HG38}/genome.fa --no-stats
${MAPDAMAGE} -i holmes.sorted.rmdup.bam -r ${HG38}/genome.fa --no-stats
```
Where are the output files located?
If you had problems generating the files we have also precomputed the nucleotide misincorporation for each cycle of the read. A plot of these can be found at [sherlock](results/sherlock.nmis.pdf) [holmes](results/holmes.nmis.pdf)

Which one looks ancient? 

Hint: Consider checking

1. Read lengths
2. Damage signals
3. Clonality
4. Endogenous content
