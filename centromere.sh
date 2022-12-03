
##enrichment of centromere reads#
#using blastn to map CCS reads to the T2T reference genome sequences##

nohup perl get20xreads.pl ccs2genomew0.fa ccs.20xall.id >ccs20xall.fa 2>err.log&
nohup blastn  -query ccs20xall.fa -num_threads 10 -max_target_seqs 2  -max_hsps 2   -db /project/xwang787_650/PB_000352_human/WGA_EcoGII/T2T_refgenome/chm13.draft_v1.1.fasta  -word_size 50 -outfmt "6 qacc sacc length nident  mismatch gaps qstart qend sstart send qlen slen bitscore score evalue" >ccs20xtoT2Tgenome.bnout 2>err.log&

