

cat ccs.np_ec.xls |perl -ne 'chomp;if($_=~/np:i:(\d+)/){$cov=$1+1;if($cov>=20){print "$_\n"}}' |cut -f 1 > ccs.20xall.id

cat ccs2mm10.blastnout cc2mm10.remained.bnout |perl -ne 'chomp;@ar=split(/\t/,$_);if($ar[2]/$ar[10]>=0.98){print "$_\n"}'|sort -k1,1 -k14,14gr|sort -k1,1 -u --merge|cut -f 1-2,7-10|perl -ne 'chomp;@ar=split(/\t/,$_);$len=$ar[3]-$ar[2];if($ar[4]>$ar[5]){$end=$ar[5]+$len;print "$ar[0]\t$ar[1]\t$ar[2]\t$ar[3]\t$ar[5]\t$end\t-\n"}else{$end=$ar[4]+$len;print "$ar[0]\t$ar[1]\t$ar[2]\t$ar[3]\t$ar[4]\t$end\t+\n"}'|perl -ne 'chomp;@ar=split(/\t/,$_);if($_!~/:/){print "$_\n"}else{if($ar[0]=~/:(\d+)\-(\d+)/){$start=$1;$tmpend=$2;$end=$tmpend-$start;print "$ar[0]\t$ar[1]\t$start\t$end\t$ar[4]\t$ar[5]\t$ar[6]\n"}}' >ccs2mm10.gt0.98.transwithnew.xls

nohup cat sdcal_sm/*.sd.xls|grep A|grep -v all|perl -ne 'chomp;@ar=split(/\t/,$_);if($ar[6]> 0.4){print "$_\n"}'|cut -f 1|sort|uniq -c |grep '      2 '|sed 's/      2 //'|cat ccs.20xall.id -|sed 's/"//g'|sort|uniq -c|grep '      1 '|sed 's/      1 //'|cat ccs2mm10.gt0.98.transwithnew.xls - | cut -f 1 | sort | uniq -c|grep '      2 '|sed 's/      2 //' |perl /project/yalidou_405/Pacbio_pipeline/YifanLab-main/Brdu_CCS/brdu_syn2h/gettargetsm.pl - ipd20x.allraw.A.csv >ipd20x.mapped2genome.sdlt0.4.allA.csv 2>err.log&

python3.8 /project/yalidou_405/Pacbio_pipeline/YifanLab-main/Brdu_CCS/6ma_gaussianAxfor6mAwithPeak.py ipd20x.sdg.ax_ipd.csv

nohup cat ipd20x.mapped2genome.sdlt0.4.allA.csv |perl -ne 'chomp;@ar=split(/,/,$_);if($ar[-3]>=2.144){print "$_\n"}' > ipd20x.mapped2genome.sdlt0.4.methyA.csv 2>err.log&


 nohup cat ipd20x.mapped2genome.sdlt0.4.methyA.Percentage.xls |perl -ne 'chomp;@ar=split(/\t/,$_);if($ar[-1]>=10){print "$_\n"}' |cut -f 1| perl /project/yalidou_405/Pacbio_pipeline/YifanLab-main/PG3683HMMmouse_batch2/outputs/demultiplexing_files/PstI_minusTam/getargetblast.pl - ccs2mm10.gt0.98.transwithnew.xls |perl -ne 'chomp;@ar=split(/\t/,$_);print "$ar[1]\t$ar[4]\t$ar[5]\t$ar[0]\n"' |windowBed -a /project/yalidou_405/Pacbio_pipeline/YifanLab-main/PBG3447_mouseHMM/PBG8302_Ref/mm10.enhancerfmtar500.bed -b - -w 0|perl -ne 'chomp;@ar=split(/\t/,$_);if($ar[1]>=$ar[5] and $ar[2]<=$ar[6]){print "$_\n"}' >enhanceroverlapwithccsar500.xls 2>err.log&
 
 nohup cat ipd20x.mapped2genome.sdlt0.4.methyA.Percentage.xls |perl -ne 'chomp;@ar=split(/\t/,$_);if($ar[-1]>=10){print "$_\n"}'|cut -f 1|perl /project/yalidou_405/Pacbio_pipeline/YifanLab-main/PG3683HMMmouse_batch2/outputs/demultiplexing_files/PstI_minusTam/getargetblast.pl - ccs2mm10.gt0.98.transwithnew.xls|perl -ne 'chomp;@ar=split(/\t/,$_);print "$ar[1]\t$ar[4]\t$ar[5]\t$ar[0]\n"' |sort -k1,1 -k2,2n|windowBed -b /project/yalidou_405/Pacbio_pipeline/YifanLab-main/PBG3447_mouseHMM/PBG8302_Ref/mm10.enhancerfmtar500.bed -a - -w 0 |perl -ne 'chomp;@ar=split(/\t/,$_);if($ar[1] <=$ar[5] and $ar[2]>=$ar[6]){print "$_\n"}' >enhanceroverlapwithccsar500.xls 2>err.log&
 
#ccs to genome##
nohup perl /project/yalidou_405/Pacbio_pipeline/YifanLab-main/sb210_seq2_rep2/5e318406-a881-11ec-a98e-b07b25d42266/blast_mapping/ccsanalysi_Perl/ccs2refpostion_blastn.pl ccs2mm10.gt0.98.transwithnew.xls ipd20x.mapped2genome.sdlt0.4.methyA.xls >ipd20x.mapped2genome.sdlt0.4.methyAongenome.xls 2>err.log&

#motif around 6nAcount#
nohup perl  /project/yalidou_405/Pacbio_pipeline/YifanLab-main/PBG3447_mouseHMM/HindIII_plusTam/gettargetregionmethypos_onlongmolfull.pl ipd20x.mapped2genome.sdlt0.4.methyAongenome.xls enhanceroverlapwithccsar500.xls >enhancer.ar500m6Acount.xls 2>err.log&
nohup cat mouseminusoverlapwithctcf.xls |cut -f 1-4,7-12|perl /project/yalidou_405/Pacbio_pipeline/YifanLab-main/PB000486_20220706/NF_EcoGII/gettargetregionmethypos_onlongmolfull.pl mouse.mergedminus.methyongenomesort.xls - >ctcf.mouse6mAminus.xls 2>err.log&
nohup less mouseplusoverlapwithctcf.xls|cut -f 1-4,7-12|perl /project/yalidou_405/Pacbio_pipeline/YifanLab-main/PB000486_20220706/NF_EcoGII/gettargetregionmethypos_onlongmolfull.pl mouse.mergedplus.methyongenomesort.xls - >ctcf.mouse6mAplus.xls 2>err.log&



nohup computeMatrix scale-regions -S bowtie2ref.HMM_H3K4me1_1.bw  -R mm10.PEar500.bed -o mm10.PEh3k4m1.tag.gz 2>err.log&


##Centremeric region analysi#

nohup less ecog2.centrall.xls|perl -ne 'chomp;@ar=split(/\t/,$_);@ele=split(/,/,$ar[3]);@pos=split(/,/,$ar[2]);if($#ele>=4){for($i=0;$i<=$#ele-4;$i++){$subele=join(",", @ele[$i..$i+4]);$subpos=join(",",@pos[$i..$i+5])  ;print "$ar[0]\t$subele\t$subpos\n"}}'|grep 'D2,D1,D2,D1,D1'|perl -ne 'chomp;@ar=split(/\t/,$_);@pos=split(/,/,$ar[-1]);$len=$pos[-1]-$pos[0]+1;print "$_\t$pos[0],$pos[-1]\n"'|perl formatcentrosm.pl - |perl encodeseq.withipdATall.pl /project/xwang787_650/PB_000352_human/WGA_EcoGII/ccs.hifi.fa /project/xwang787_650/PB_000352_human/WGA_EcoGII/ccsallA.methyAgt2.0.csv - >ecog2.D21211encode.xls 2>err.log&

