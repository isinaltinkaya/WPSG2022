The dataset for this exercise consists of two [bamfiles](ancient_data/) these are mapped fastq files.

**Aims**:

  - Get experience with navigating files associated with high throughput sequencing data

  - Figure out which sample looks ancient and which one looks modern




# Set some environment variables and paths

First set some paths   

    # NB this must be done every time you open a new terminal

    
    # Set path to samtools program
    SAMTOOLS=samtools/samtools
    
    # Set path to folder containing data
    DATAFOLDER=dt/ancient_data/


Let us first validate that we setup our variables correctly


```
ls ${SAMTOOLS} ${DATAFOLDER}
```

# Summary statistic and convertion

We have two files called holmes.sam.gz and sherlock.sam.gz. Let us copy these files to our working directory.

```
cp ${DATAFOLDER}/* .
```

Let us try to view the uncompressed version of one of the files.

```
gunzip -c sherlock.sam.gz|less -NS
```
We are piping the output from the uncompression into the less commands and are supplying less with the arguments -NS, this will give us line number and we can use the arrows to got up and down left and right. You exit the commands by pressing q.

 1. Identify the header
 2. How many lines does the header span

Let us see how many reads we have in these two files

```
gunzip -c sherlock.sam.gz |grep -P "^@" -v |wc -l
gunzip -c holmes.sam.gz |grep -P "^@" -v |wc -l
```

Now let us convert the SAM files to BAM files

```
${SAMTOOLS} view -b sherlock.sam.gz -o sherlock.bam
${SAMTOOLS} view -b holmes.sam.gz -o holmes.bam
```

Try to look at one of the BAM files using samtools view, is there are difference from uncompressing the .sam.gz file?

```
${SAMTOOLS} view sherlock.bam|less -SN
```
Are the reads ordered?

We can sort the file using samtools sort

```
${SAMTOOLS} sort sherlock.bam -OBAM -o sherlock.sorted.bam
${SAMTOOLS} sort holmes.bam -OBAM -o holmes.sorted.bam
```

Visually inspect that they are sorted using previous view pipe less command
```
${SAMTOOLS} view sherlock.sorted.bam|less -SN
```

Let us find the proportion of reads the maps (endogenous content) for each file

```
${SAMTOOLS} flagstat sherlock.sorted.bam
${SAMTOOLS} flagstat holmes.sorted.bam
```

Let us look at the distribution of read lengths (for the read that maps)

```
${SAMTOOLS} view -F4 sherlock.sorted.bam |awk '{print length($10)}'|sort -n |uniq -c >sherlock.readlength
${SAMTOOLS} view -F4 holmes.sorted.bam |awk '{print length($10)}'|sort -n |uniq -c >holmes.readlength
```

Let us load the data into R and plot it

```
sherlock <-read.table("sherlock.readlength")
holmes <-read.table("holmes.readlength")
par(mfrow=c(1,2))
plot(sherlock[,2],sherlock[,1],main="Sherlock",type='l',col=1,lwd=2)
plot(holmes[,2],holmes[,1],main="Holmes",type='l',col=2,lwd=2)

```
If you had problems with the R commands the result can be found [here](results/sherlock.holmes.rlen.pdf)


Let us also look at the proportion of duplicates(clonality)

```
${SAMTOOLS} rmdup -s sherlock.sorted.bam sherlock.sorted.rmdup
${SAMTOOLS} rmdup -s holmes.sorted.bam holmes.sorted.rmdup
```

What is reported, which ones contains the most duplicates?

From these files we have also precomputed the nucleotide misincorporaton for each cycle of the read. A plot of these can be found [sherlock](results/sherlock.nmis.pdf) [holmes](results/holmes.nmis.pdf)

First lets set some filter to remove the worst reads (minMapQ), remove
the worst of the bases (minQ).

    FILTERS="-minMapQ 30 -minQ 20"

Lets set some options that means we will calculate genotype likelihoods
using the GATK model (gl) and calculate the site allele frequency
likelihoods (saf)

    OPT=" -dosaf 1 -gl 2"

Generate site frequency likelihoods using ANGSD

    $ANGSD -b  YRI.filelist  -anc $ANC -out yri $FILTERS $OPT -ref $REF &
    $ANGSD -b  JPT.filelist  -anc $ANC -out jpt $FILTERS $OPT -ref $REF &
    $ANGSD -b  CEU.filelist  -anc $ANC -out ceu $FILTERS $OPT -ref $REF

The run time is a couple of minutes

If it talks to long then you can copy the results using this command:

    cp dt/yri.saf* .
    cp dt/ceu.saf* .
    cp dt/jpt.saf* .

Estimate the site frequency spectrum for each of the 3 populations
without having to call genotypes or variable sites directly from the
site frequency likelihoods

    #calculate the 1 dimensional SFS
    $REALSFS yri.saf.idx > yri.sfs
    $REALSFS jpt.saf.idx > jpt.sfs
    $REALSFS ceu.saf.idx > ceu.sfs

In order to plot the results open R and make a barplot

``` 
 ##run in R                      
#plot the results
nnorm <- function(x) x/sum(x)
#expected number of sites with 1:20 derived alleles
res <- rbind(
  YRI=scan("yri.sfs")[-1],
  JPI=scan("jpt.sfs")[-1],
  CEU=scan("ceu.sfs")[-1]
)
colnames(res) <- 1:20

# density instead of expected counts
res <- t(apply(res,1,nnorm))

#plot the none ancestral sites
barplot(res,beside=T,legend=c("YRI","JPT","CEU"),names=1:20,main="realSFS non ancestral sites")

#plot the polymorphic sites. 
resPoly <- t(apply(res[,-20],1,nnorm))
barplot(resPoly,beside=T,legend=c("YRI","JPT","CEU"),names=1:19,main="realSFS polymorphic sites")

#due the very limited amount of sites
#downsample to 5 individuals (10 chromosome) and exclude fixed derived
downsampleSFS <- function(x,chr){ #x 1:2n , chr < 2n
    n<-length(x)
    mat <- sapply(1:chr,function(i) choose(1:n,i)*choose(n- (1:n),chr-i)/choose(n,chr))
    nnorm( as.vector(t(mat) %*% x)[-chr] )
}
resDown <- t(apply(res,1,downsampleSFS,chr=10))
barplot(resDown,beside=T,legend=c("YRI","JPT","CEU"),names=1:9,main="realSFS downsampled polymorphic sites")
```
If you had problems with the above commands the plots can also be found [her](results/yri.jpt.ceu.1dsfs.pdf)
  - Which population has the largest population size?

  - The data is a small subset of the genome (2Mb). If you had analysed
    6Mb it sould have looked like
    [this](results/realSFS4.pdf)

  - The analysed whole chromosome for the 1000G individual look [like
    this](results/full1dsfs.pdf)

lets use the sfs to calculate some statistics for the population

``` 

 ##run in R                      
## read sfs
yri<-scan("yri.sfs");
jpt<-scan("jpt.sfs");
ceu<-scan("ceu.sfs");

x<-ceu #change this one to try one of the other populations 

nSites<-sum(x)   #Number of sites where we have data
nSeg<-sum(x[c(-1,-21)])    #Number of segregating sites
an <- function(n) sum(1/1:(n-1)) 
thetaW <- nSeg/an(20) # Wattersons Theta
thetaW / 2.5e-8 / nSites / 4 # effective population size
```

The above example is for the African population. Try to run it for all
three populations.

  - which has the largest populations size

  - which has the largest variability (fraction of
    polymorphic/segregating sites)

## Fst and PBS In order to estimate Fst between two population we will need to estimate the 2-dimensional frequency spectrum from the site allele frequency likelihoods

    #calculate the 2D SFS 
    $REALSFS yri.saf.idx ceu.saf.idx >yri.ceu.ml &
    $REALSFS yri.saf.idx jpt.saf.idx >yri.jpt.ml &
    $REALSFS jpt.saf.idx ceu.saf.idx >jpt.ceu.ml

Plot the results in R

``` 
 ##run in R                      
yc<-scan("yri.ceu.ml")
yj<-scan("yri.jpt.ml")
jc<-scan("jpt.ceu.ml")
    source("plot2dSFS.R")
plot2<-function(s,...){
    dim(s)<-c(21,21)
    s[1]<-NA
    s[21,21]<-NA
s<-s/sum(s,na.rm=T)

    pal <- color.palette(c("darkgreen","#00A600FF","yellow","#E9BD3AFF","orange","red4","darkred","black"), space="rgb")
    pplot(s/sum(s,na.rm=T),pal=pal,...)
}

plot2(yc,ylab="YRI",xlab="CEU")
x11()
plot2(yj,ylab="YRI",xlab="JPT")
x11()
plot2(jc,ylab="JPT",xlab="CEU")
```
If you had problems running the above commands the plots can be found [here](results/2dsfs.pdf)

Due to the very limited amount of data the plots are very noisy. However
they are still informative.The colors indicate the density. High density
means many sites will look like this and low density (green) means that
few sites looks like this.

Based on the plots try to guess

  - Which populations has most private SNPs (sites that are only
    polymorphic in this population)

  - Which two populatons are most closely related?

close R

In order to get a measure of this populations are most closely related
we willl estimate the pairwise Fst

    #first will will index the sample so the same sites are analysed for each population
    $REALSFS fst index jpt.saf.idx ceu.saf.idx -sfs jpt.ceu.ml -fstout jpt.ceu
    $REALSFS fst index yri.saf.idx ceu.saf.idx -sfs yri.ceu.ml -fstout yri.ceu
    $REALSFS fst index yri.saf.idx jpt.saf.idx -sfs yri.jpt.ml -fstout yri.jpt
    
    #get the global estimate
    $REALSFS fst stats jpt.ceu.fst.idx
    $REALSFS fst stats yri.jpt.fst.idx
    $REALSFS fst stats yri.ceu.fst.idx 

look at the weigthed Fst (Fst.Weight).

  - which two populations are most closely related?

  - which two populations are most distantly related?

Lets see how the Fst and PBS varies between different regions of the
genome my using a sliding windows approach (windows site of 50kb)

    $REALSFS fst index yri.saf.idx jpt.saf.idx ceu.saf.idx -fstout yri.jpt.ceu -sfs yri.jpt.ml -sfs yri.ceu.ml -sfs jpt.ceu.ml
    $REALSFS fst stats2 yri.jpt.ceu.fst.idx -win 50000 -step 10000 >slidingwindowBackground

read the data into R

``` 
 ##run in R                      
r<-read.delim("slidingwindowBackground",as.is=T,head=T)
names(r)[-c(1:4)] <- c("wFst_YRI_JPT","wFst_YRI_CEU","wFst_JPT_CEU","PBS_YRI","PBS_JPT","PBS_CEU")


head(r) #print the results to the screen

#plot the distribution of Fst
mmax<-max(c(r$wFst_YRI_JPT,r$wFst_YRI_CEU,r$wFst_JPT_CEU),na.rm=T)
par(mfcol=c(3,2))
hist(r$wFst_YRI_JPT,col="lavender",xlim=c(0,mmax),br=20)
hist(r$wFst_YRI_CEU,col="mistyrose",xlim=c(0,mmax),br=20)
hist(r$wFst_JPT_CEU,col="hotpink",xlim=c(0,mmax),br=20)

mmax<-max(c(r$PBS_CEU,r$PBS_YRI,r$PBS_JPT),na.rm=T)

#plot the distribution of PBS
mmax<-max(c(r$PBS_CEU,r$PBS_YRI,r$PBS_JPT),na.rm=T)
hist(r$PBS_YRI,col="lavender",xlim=c(0,mmax),br=20)
hist(r$PBS_CEU,col="mistyrose",xlim=c(0,mmax),br=20)
hist(r$PBS_JPT,col="hotpink",xlim=c(0,mmax),br=20)

```
If you had problems running the above commands you can find the result [here](results/fst_pbs.pdf)

note the maximum observed values for both the pairwise fst and the PBS

Lets do the same for not so randomly selection 1Mb region of on chr 5.
Remember to close R

``` 
                                                                                                                       
#a African population for a region on chr 5                                                 
find $BAMFOLDERchr5 | grep bam$ | grep YRI > YRIchr5.filelist
#a Asian population for a region on chr 5                                                                                       
find $BAMFOLDERchr5 | grep bam$ | grep JPT > JPTchr5.filelist
#a European population for a region on chr 5                                                                                        
find $BAMFOLDERchr5 |  grep bam$ | grep CEU > CEUchr5.filelist

#use the same filters and options as before
FILTERS="-minMapQ 30 -minQ 20 -baq 1 -C 50 -minInd 8"
OPT=" -dosaf 1 -gl 2"

#get site frequency likelihoods
$ANGSD -b  YRIchr5.filelist  -anc $ANC -out yriChr5 $FILTERS $OPT -ref $REF
$ANGSD -b  JPTchr5.filelist  -anc $ANC -out jptChr5 $FILTERS $OPT -ref $REF
$ANGSD -b  CEUchr5.filelist  -anc $ANC -out ceuChr5 $FILTERS $OPT -ref $REF

#estimate the 1D SFS
$REALSFS yriChr5.saf.idx ceuChr5.saf.idx >yri.ceuChr5.ml
$REALSFS yriChr5.saf.idx jptChr5.saf.idx >yri.jptChr5.ml
$REALSFS jptChr5.saf.idx ceuChr5.saf.idx >jpt.ceuChr5.ml

#get FST and PBS in sliding window
$REALSFS fst index yriChr5.saf.idx jptChr5.saf.idx ceuChr5.saf.idx -fstout yri.jpt.ceuChr5 -sfs yri.jptChr5.ml -sfs yri.ceuChr5.ml -sfs jpt.ceuChr5.ml
$REALSFS fst stats2 yri.jpt.ceuChr5.fst.idx -win 50000 -step 10000 >slidingwindowChr5
```

Lets view how it looks in this region

    #run in R
    r<-read.delim("slidingwindowChr5",as.is=T,head=T)
    names(r)[-c(1:4)] <- c("wFst_YRI_JPT","wFst_YRI_CEU","wFst_JPT_CEU","PBS_YRI","PBS_JPT","PBS_CEU")
    
    
    par(mfrow=1:2)
    plot(r$midPos,r$wFst_YRI_CEU,ylim=c(0,max(r$wFst_YRI_CEU)),type="b",pch=18,ylab="Fst",xlab="position on Chr 5")
    points(r$midPos,r$wFst_YRI_JPT,col=2,type="b",pch=18)
    points(r$midPos,r$wFst_JPT_CEU,col=3,type="b",pch=18)
    legend("topleft",fill=1:3,c("YRI vs. CEU","YRI vs. JPT","JPT vs CEU"))
    
    plot(r$midPos,r$PBS_YRI,ylim=c(0,max(r$PBS_CEU)),type="b",pch=18,ylab="PBS",xlab="position on Chr 5")
    points(r$midPos,r$PBS_JPT,col=2,type="b",pch=18)
    points(r$midPos,r$PBS_CEU,col=3,type="b",pch=18)
    legend("topleft",fill=1:3,c("YRI","JPT","CEU"))

If you problems running the above commands the results can be found [here](results/fst_pbs_chr5.pdf)

  - Compare the values you observed on this part of the genome with the
    random pars of the genome you looked at [earlier](results/fst_pbs.pdf)). Is this region extreme?

  - Why is there two peak for the Fst and only one for the PBS?

  - In which of the populations are this loci under selection?

Find out what genes is in this region by going to the [UCSC
browser](https://genome.ucsc.edu/index.html). Choose Genome browser.
Choose human GRCh37/hg19 and find the region. Read about this gene on
wikipedia and see if this fits PBS results.

### What happens if we try to call genotypes?

We can compare with what happens if we try to call genotypes by calling
SNPs and genotypes like GATK. If you are running out of time then skip
this part

``` 
FILTERS2="-minMapQ 30 -minQ 20 -minInd 10"                           

OPT2="-gl 2 -doGeno 2 -doPost 2 -doMajorMinor 4 -doMaf 1 -SNP_pval 1e-6 -postCutoff 0.95"
$ANGSD -b  YRI.filelist  -out yri $FILTERS2 $OPT2 -ref $REF &  
$ANGSD -b  JPT.filelist  -out jpt $FILTERS2 $OPT2 -ref $REF &
$ANGSD -b  CEU.filelist  -out ceu $FILTERS2 $OPT2 -ref $REF  
```

While it runs you can look at the options we choose:

  - *minInd 10*: minimum individuals with data (in this case it means
    with called genotypers). Why do we need this?

  - *-doGeno 2* Print only the count (0,1,2) and not the based e.g.
    AA,AT,TT

  - *-doPost 2* Use uniform prior for genotype i.e. call the genotype
    with the highest likelihood

  - *-DoMajorMinor 4* Use the Ancestral allele from the chimp

  - *-doMaf 1 -SNP\_pval 1e-6* Use this p-value cutoff to call SNPs.
    What would happend to the SFS if you change this threshold?

  - *-PostCutoff 0.95* only call genotype with a propability above 0.95

Plot the results in R

``` 
 ##run in R                      
#plot the results
nnorm <- function(x) x/sum(x)
getSFS<-function(x) table(factor(rowSums(read.table(x)[,-c(1:2)]),levels=1:20))

res <- rbind(
  YRI=getSFS("yri.geno.gz"),
  JPI=getSFS("jpt.geno.gz"),
  CEU=getSFS("ceu.geno.gz")
)
colnames(res) <- 1:20

# density instead of expected counts
res <- t(apply(res,1,nnorm))

#plot the none ancestral sites
barplot(res,beside=T,legend=c("YRI","JPT","CEU"),names=1:20,main="SFS from called genotypes")

#plot the polymorphic sites. 
resPoly <- t(apply(res[,-20],1,nnorm))
barplot(resPoly,beside=T,legend=c("YRI","JPT","CEU"),names=1:19,main="SFS from call\
ed genotypes")


#down sample to 5 individuals (10 chromosome) and exclude fixed derived
downsampleSFS <- function(x,chr){ #x 1:2n , chr < 2n
    n<-length(x)
    mat <- sapply(1:chr,function(i) choose(1:n,i)*choose(n- (1:n),chr-i)/choose(n,chr))
    nnorm( as.vector(t(mat) %*% x)[-chr] )
}
resDown <- t(apply(res,1,downsampleSFS,chr=10))
barplot(resDown,beside=T,legend=c("YRI","JPT","CEU"),names=1:9)

```
If you had problems running the above commands the result can be found [here](results/yri.jpt.ceu.gc.1dsfs.pdf)
  - How does this compare to the likelhood based estimates
    ([pdf](results/yri.jpt.ceu.1dsfs.pdf))
