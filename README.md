# draft_teaching


## Site frequency spectrum


We can now use the SFS as prior information using `-pest`

Compute the allele frequency posterior probabilities and associated statistics (`-doThetas`)

```
realSFS saf2theta sim_rep0_d10.saf.idx -sfs simulations/sfs/sim_rep0_d10.sfs -outname o
```


```
thetaStat print out.thetas.idx
```
```
$ /maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/angsd/misc/thetaStat print testout.thetas.idx 2>/dev/null |head                        
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


```
./misc/thetaStat do_stat out.thetas.idx
```

### Sliding window

We can do a sliding window analysis using a window size of 50kb and a step size of 10kb:

```
thetaStat do_stat out.thetas.idx -win 50000 -step 10000  -outnames theta.thetasWindow.gz
```

`./thetaStat print thetas.idx`  log scaled per site estimates of the thetas
`pestPG` sum of the per site estimates for a region

```
#(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)   Chr     WinCenter       tW      tP      tF      tH      tL    Tajima   fuf     fud     fayh    zeng    nSites
(0,63025519)(1,63025520)(0,63025520)    chr20   31512760        29084.489811    29094.351398    29120.408460    34251.072423    31672.711913  0.001278 -0.001687       -0.003197       -0.142269       0.072371        63025519
```


