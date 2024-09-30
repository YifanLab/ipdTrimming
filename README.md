# 6mA_footprint
##Here is the comparison of pipelines between ipdSummary and ipdTrimming.

![image](https://github.com/user-attachments/assets/610a723a-c0bf-4203-a37e-1d0d260944c3)

For ipdTrimming, a Local trimming of IPD outliers at the subread level was introduced into IPD conversion: for each site in a DNA molecule, the top 10% of subread IPD values were removed; the remaining subread IPD values were then averaged to generate the CCS IPD value.  This was implemented by a modified CCS module. 

```console
usc@hpc:~$ ccs -j 20 --hifi-kinetics  subread.bam hifi.bam
    
```

For this modified CCS module, we get the CCS bam file with trimmed IPD value. Following a movie-time normalization was applied to compensate for variations in the polymerase elongation rate on individual DNA molecules, each CCS IPD value was normalized against the CCS IPD value averaged across all sites in a DNA molecule.  This CCS IPD value was then compared to that of an unmodified base with the same sequence context, provided by a kmer model pretrained on the Sequel data. The last two steps were implemented by the script ipdRatiocalculator_FromCCS.py 

```console
usc@hpc:~$ python ipdRatiocalculator_FromCCS.py hifi.bam hifi.withIPDr.bam
#the cpu number can be set in this script.
```

