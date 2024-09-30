# 6mA_footprint
##Here is the comparison of pipelines between ipdSummary and ipdTrimming.

![image](https://github.com/user-attachments/assets/610a723a-c0bf-4203-a37e-1d0d260944c3)

For ipdTrimming, a Local trimming of IPD outliers at the subread level was introduced into IPD conversion: for each site in a DNA molecule, the top 10% of subread IPD values were removed; the remaining subread IPD values were then averaged to generate the CCS IPD value.  This was implemented by a modified CCS module. 
