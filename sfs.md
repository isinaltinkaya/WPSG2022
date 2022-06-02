# Sitefrequency spectrum and VCF files

###
First setup some paths and environment variables

```
DATA="BCF"
ANGSD="angsd"
REALSFS="realsfs"
```

Validate that we have setup our variables correctly
```
ls ${DATA}
```

### Site allele frequencies

```
${ANGSD} -vcf-gl ${DATA}/POP1.bcf -doSaf 1 -anc ${DATA}/chr20.fa.gz -out POP1
${ANGSD} -vcf-gl ${DATA}/POP2.bcf -doSaf 1 -anc ${DATA}/chr20.fa.gz -out POP2
##if the above runs takes forever, (it took 5minutes on my desktop),
##we can limit the analyses to 30megabases in the central part of chromosome20
${ANGSD} -vcf-gl ${DATA}/POP1.bcf -doSaf 1 -anc ${DATA}/chr20.fa.gz -out POP1 -r chr20:20000000-50000000
${ANGSD} -vcf-gl ${DATA}/POP2.bcf -doSaf 1 -anc ${DATA}/chr20.fa.gz -out POP2 -r chr20:20000000-50000000
```

Which files was generated?

|Filetype     | Explanation                                           |
| --- | ------------------------------------------ |
| arg | arguments used for the analysis   |
| saf.gz | containing the sample allele frequencies for all sites   |
| saf.pos.gz | containing the position         |
| saf.idx | index file containing the binary offset |



### Site frequency spectrum
The data are the sample allele frequency loglikelihoods these can be viewed with:
```
realSFS print FILE.saf.idx|head
```
The first two columns are the chromosome and position followed by the saf for each bin. We can obtain an estimate of the global site frequency spectrum for each popoulation using the following commands
 
```
${REALSFS} POP1.saf.idx > POP1.sfs
${REALSFS} POP2.saf.idx > POP2.sfs
```

Have a look at the .sfs files :

```
cat POP1.sfs
cat POP2.sfs
```
We need to plot these, we will use R

```
p1 <- scan(POP1.sfs)
p2 <- scan(POP2.sfs)
barplot(rbind(p1,p2))
barplot(rbind(p1,p2)[,-1])
```


### Allele frequency posterior probabilities and associated statistics (`-doThetas`)



```
realSFS saf2theta POP1.saf.idx -sfs POP1.sfs -outname POP1
```



- OUTFILE.thetas.idx

We can view the theta statistics using `./thetaStat print thetas.idx`. This file contains log scaled per site estimates of the thetas.


```
thetaStat print POP1.thetas.idx
```
```
$thetaStat print testout.thetas.idx 2>/dev/null |head                        
#Chromo Pos     Watterson       Pairwise        thetaSingleton  thetaH  thetaL                   
chr20   1       -13.837903      -15.382814      -12.393384      -19.039749      -16.050478
chr20   2       -14.297906      -15.843701      -12.852455      -19.502541      -16.511412
chr20   3       -13.446123      -14.991596      -12.001015      -18.649746      -15.659290
chr20   4       -12.615373      -14.158298      -11.172954      -17.810963      -14.825852
chr20   5       -14.952734      -16.499620      -13.506134      -20.160820      -17.167391
chr20   6       -11.360343      -12.901918      -9.919370       -16.551733      -13.569401
chr20   7       -14.651113      -16.197880      -13.204640      -19.858820      -16.865644
chr20   8       -14.741365      -16.288082      -13.294944      -19.948916      -16.955843
chr20   9       -8.865955       -10.400315      -7.432686       -14.034883      -11.067410
```


| Column index | 1           | 2        | 3                                                           | 4                                                           | 5              | 6                                                         | 7                                          |
|--------------|-------------|----------|-------------------------------------------------------------|-------------------------------------------------------------|----------------|-----------------------------------------------------------|--------------------------------------------|
| Column ID    | #Chromo     | Pos      | Watterson                                                   | Pairwise                                                    | thetaSingleton | thetaH                                                    | thetaL                                     |
| Example data | chr20       | 1        | -13.837903 | -15.382814 | -12.393384     | -19.039749 | -16.050478 |
|              | Contig name | Position | Watterson's theta                                           | ThetaD Nucleotide diversity                                 |                | FayH  | L |
|              |             |          | $$\sum _{i=1}^{n-1} \eta _i/a^{-1}, a= \sum _{i=1}^{n-1}i $$ | $${{n} \choose {2}}^{-1} \sum _{i=1}^{n-1}i(n-i) \eta _i $$ | $$\eta _ 1$$ | $${{n} \choose {2}}^{-1} \sum _{i=1}^{n-1}i^ 2 \eta _i$$ | $${n-1}^{-1} \sum _{i=1}^{n-1}i \eta _i $$ |



### Sliding window

We can do a sliding window analysis using a window size of 50kb and a step size of 10kb:

```
thetaStat do_stat out.thetas.idx -win 50000 -step 10000  -outnames theta.thetasWindow.gz
```

`pestPG` contains the sum of the per site estimates for a region

```
#(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)   Chr     WinCenter       tW      tP      tF      tH      tL    Tajima   fuf     fud     fayh    zeng    nSites
(0,63025519)(1,63025520)(0,63025520)    chr20   31512760        29084.489811    29094.351398    29120.408460    34251.072423    31672.711913  0.001278 -0.001687       -0.003197       -0.142269       0.072371        63025519
```



```
./misc/thetaStat do_stat out.thetas.idx
```




```
d<-read.table("POP1.pestPG",header=F)
colnames(d)<-c("Index","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
plot(t$Tajima)
plot(t$tW)
```
