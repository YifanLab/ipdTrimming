
##Pheatmap for autocorr#
smph = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/PB000486_20220706/nfhia5.sm.allidencodeforauto.sub1000.xls", row.names = 1, sep="\t")
smphmad <- smph[rev(order(apply(smph[, c(50,999)], 1,mad))),]
breaks <- 0.0005*(-50:50)
pheatmap(smphmad, show_rownames = F, show_colnames = F, cluster_cols = F, cluster_rows = F, breaks = breaks,color=colorRampPalette(c('blue4', 'white', 'red'))(100))

##lineplot for autocorr##
smlin = read.csv("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/PB000486_20220706/metHia5.smautocorr.mean.xls", sep=",", header = T)
ggplot(smlin, aes(x=X, y= X0))+geom_point(size=1)+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+ theme_bw()+theme(axis.title.x = element_blank(), axis.title.y = element_blank())


##

ctcfnf = read.csv("/var/www/visul_jb/sb210_tetrahymena/methia5.ccsbedoverlapwithctcfaround50.count.xls", sep="\t", header = F, row.names = 1)
ctcfnffilter  = ctcfnf[rowSums(ctcfnf) >= 5,]
ctcfnfinnercount = as.data.frame(rowSums(ctcfnffilter[,c(52:71)]))
ctcfnfoutcount = as.data.frame(rowSums(ctcfnffilter[,c(1:51,72:120)]))
ctcfcountsum = cbind(ctcfnfinnercount, ctcfnfoutcount)
colnames(ctcfcountsum) <- c('Inner', 'Outer')
ctcfcountsum$InnerPer = ctcfcountsum$Inner/19
ctcfcountsum$OuterPer = ctcfcountsum$Outer/111
p3 = ggplot(ctcfcountsum, aes(x= InnerPer, y= OuterPer))+stat_density_2d( aes(fill = ..density..), geom='raster', contour = F)+geom_density_2d(color='red', alpha=0.3) +scale_fill_distiller(palette = 4, direction = 1, limits=c(0,100))+scale_x_continuous(expand = c(0,0),limits = c(0,0.5))+scale_y_continuous(expand = c(0,0), limits = c(0,0.5))



#p2 = ggplot(ctcfcountsum, aes(x=OuterPer))+geom_density(binwidth = 0.01, position = 'identity', fill='#66c2a4', alpha=0.4)+theme_void()+theme(plot.margin = margin())+scale_x_continuous(expand=c(0,0), limits = c(0,0.5))+scale_y_continuous(expand=c(0,0))+coord_flip()

p1 = ggplot(ctcfcountsum, aes(x=InnerPer))+geom_histogram(binwidth = 0.01, bins=50 , fill='#66c2a4', alpha=0.4)+scale_x_continuous(expand=c(0,0), limits = c(-0.02, 0.5))+scale_y_continuous(expand=c(0,0))+theme_void()
p2 = ggplot(ctcfcountsum, aes(x=OuterPer))+geom_histogram(binwidth = 0.01, bins=50 , fill='#66c2a4', alpha=0.4)+theme_void()+theme(plot.margin = margin())+scale_x_continuous(expand=c(0,0), limits = c(-0.02,0.5))+scale_y_continuous(expand=c(0,0))+coord_flip()
alignxhist <- align_plots(p1, p3, align = "v",  axis="tblr")[[1]]
alignyhist <- align_plots(p2, p3, align = "h", axis="tblr")[[1]]
plot_grid(alignxhist, NULL, p3, alignyhist, ncol=2, nrow=2, rel_heights = c(0.2,1), rel_widths = c(1, 0.2))


#innercountper = as.data.frame(ggplot_build(p1)$data[[1]])[,c('count', 'x')]

#View(as.data.frame(ggplot_build(p1)$data[[1]])[, c(2,3,4,5)])
innercountper = as.data.frame(ggplot_build(p1)$data[[1]])[, c(2,3,4,5)]
innercountper[innercountper$count>0,]

ctctgroup1 = ctcfcountsum[ctcfcountsum$InnerPer>=0 & ctcfcountsum$InnerPer<=0.005,]
ctctgroup2 = ctcfcountsum[ctcfcountsum$InnerPer>=0.045 & ctcfcountsum$InnerPer<=0.055,]
ctctgroup3 = ctcfcountsum[ctcfcountsum$InnerPer>=0.105 & ctcfcountsum$InnerPer<=0.115,]
ctctgroup4 = ctcfcountsum[ctcfcountsum$InnerPer>=0.155 & ctcfcountsum$InnerPer<=0.165,]
ctctgroup5 = ctcfcountsum[ctcfcountsum$InnerPer>=0.205 & ctcfcountsum$InnerPer<=0.215,]
ctctgroup1$Group <- 'Group1'
ctctgroup2$Group <- 'Group2'
ctctgroup3$Group <- 'Group3'
ctctgroup4$Group <- 'Group4'
ctctgroup5$Group <- 'Group5'
ggplot(rbind(ctctgroup1, ctctgroup2, ctctgroup3, ctctgroup4, ctctgroup5), aes(x= OuterPer, color=Group))+geom_density(alpha=0.3)+scale_fill_brewer()+theme_bw()





##compared mean and plot##
obs = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/Motif_Jasper/hg19_MA0139.1.ar50.s45.methyobs.mean.xls", sep="\t", header = T, row.names = 1)
exp = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/Motif_Jasper/hg19_MA0139.1.ar50.s45.methyexpect.mean.xls", sep="\t", header = T, row.names = 1)
ovse = cbind(obs, exp)
colnames(ovse) = c('Obs', 'Exp')
ovse$OvE_D = ovse$Obs - ovse$Exp
p1 = ggplot()+geom_line(data=ovse, aes(x=seq(1:nrow(ovse)), y=Obs, color='Obs'))+geom_line(data=ovse, aes(x=seq(1:nrow(ovse)), y=Exp, color='Exp'))+scale_color_manual(values=c('red', 'blue'))
p2 = ggplot()+geom_line(data=ovse, aes(x=seq(1:nrow(ovse)), y=OvE_D, color='OvE'))+scale_color_manual(values = c('black'))


#Figure7s for G&H#
library(cowplot)
smrawdata <- read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/Genome_region/chr_144.34800_35106.fmt.xls", sep="\t", header=F)
smsortid <- read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/Genome_region/chr_144.34800_35106.sort.id", sep="\t", header = F)
orderid <-  factor(1:length(smsortid$V1), labels =  rev(smsortid$V1))
smrawdata$V1 <-  factor(smrawdata$V1, levels = levels(orderid))

smmethyplot <- ggplot(smrawdata) + geom_rect(aes(xmin=V3, xmax=V4,ymin=-1, ymax=1),fill='#EBEBEB')+ geom_point(aes(x=V5, y=0, color=V7),size=1)  +theme_bw()+theme(panel.grid =element_blank(), panel.border = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_color_manual(values=c('black', 'blue',  'red'))    +xlab("chr_144") + facet_wrap(~V1, ncol=1)+theme(strip.background = element_blank(), strip.text.x = element_blank())+theme(legend.position = 'none', axis.title.y = element_blank())


smsigdata <- read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/Genome_region/chr_144.34800_35106.sigmap.xls", sep="\t", header=F)
smsigplot <-  ggplot(smsigdata, aes(x=V2, y= V5))+geom_bar(stat='identity', width = 0.5)+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand = c(0,0), breaks = c(0,1))+theme_classic()+theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank())
plot_grid(smsigplot, smmethyplot, ncol=1, rel_heights = c(0.1, 1))

##figure 7c##
dyadpoint <- read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/Genome_region/chr_057.253670_253973.dyad.xls", sep="\t", header = F)
dyadpoint$V4 <- rev(0.01*(1:nrow(dyadpoint)))

dyadplot <- ggplot(dyadpoint, aes(x= V3, y=V4))+geom_point(size=0.5) + theme_nothing()

smrawdata <- read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/Genome_region/chr_057.253670_253973.fmt.xls", sep="\t", header=F)
smsortid <- read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/Genome_region/chr_057.253670_253973.sort.id", sep="\t", header = F)
orderid <-  factor(1:length(smsortid$V1), labels =  rev(smsortid$V1))
smrawdata$V1 <-  factor(smrawdata$V1, levels = levels(orderid))
smmethyplot <- ggplot(smrawdata) + geom_rect(aes(xmin=V3, xmax=V4,ymin=-1, ymax=1),fill='#EBEBEB')+ geom_point(aes(x=V5, y=0),color='red',size=1)  +theme_bw()+theme(panel.grid =element_blank(), panel.border = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())  +xlab("chr_144") + facet_wrap(~V1, ncol=1)+theme(strip.background = element_blank(), strip.text.x = element_blank())+theme(legend.position = 'none', axis.title.y = element_blank())

smsigdata <- read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/Genome_region/chr_057.253670_253973.sigmap.xls", sep="\t", header=F)
smsigplot <- ggplot(smsigdata, aes(x=V2, y= V5))+geom_bar(stat='identity')+theme_nothing()
smsigplot <-  ggplot(smsigdata, aes(x=V2, y= V5))+geom_bar(stat='identity', width = 0.5)+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand = c(0,0), breaks = c(0,1))+theme_classic()+theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank())

plot_grid(smsigplot, smmethyplot, dyadplot, ncol=1, rel_heights = c(0.1, 1,0.1))


##plot random rank auto##
randomrank1 = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/CTCF_NF_EcoG2/Auto_corr/rank1.randeom1000.autocorr.xls", sep="\t", row.names = 1)
randomrank2 = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/CTCF_NF_EcoG2/Auto_corr/rank2.randeom1000.autocorr.xls", sep="\t", row.names = 1)
randomrank3 = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/CTCF_NF_EcoG2/Auto_corr/rank3.randeom1000.autocorr.xls", sep="\t", row.names = 1)
randomrank4 = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/CTCF_NF_EcoG2/Auto_corr/rank4.randeom1000.autocorr.xls", sep="\t", row.names = 1)
randomrank5 = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/CTCF_NF_EcoG2/Auto_corr/rank5.randeom1000.autocorr.xls", sep="\t", row.names = 1)
randomrank6 = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/CTCF_NF_EcoG2/Auto_corr/rank6.randeom1000.autocorr.xls", sep="\t", row.names = 1)
randomrank7 = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/CTCF_NF_EcoG2/Auto_corr/rank7.randeom1000.autocorr.xls", sep="\t", row.names = 1)
randomrank8 = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/CTCF_NF_EcoG2/Auto_corr/rank8.randeom1000.autocorr.xls", sep="\t", row.names = 1)
randomrank9 = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/CTCF_NF_EcoG2/Auto_corr/rank9.randeom1000.autocorr.xls", sep="\t", row.names = 1)
randomrank10 = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/CTCF_NF_EcoG2/Auto_corr/rank10.randeom1000.autocorr.xls", sep="\t", row.names = 1)


randomrank1 = randomrank1[rev(order(apply(randomrank1[, c(1,999)], 1,mad))),]
randomrank2 = randomrank2[rev(order(apply(randomrank2[, c(1,999)], 1,mad))),]
randomrank3 = randomrank3[rev(order(apply(randomrank3[, c(1,999)], 1,mad))),]
randomrank4 = randomrank4[rev(order(apply(randomrank4[, c(1,999)], 1,mad))),]
randomrank5 = randomrank5[rev(order(apply(randomrank5[, c(1,999)], 1,mad))),]
randomrank6 = randomrank6[rev(order(apply(randomrank6[, c(1,999)], 1,mad))),]
randomrank7 = randomrank7[rev(order(apply(randomrank7[, c(1,999)], 1,mad))),]
randomrank8 = randomrank8[rev(order(apply(randomrank8[, c(1,999)], 1,mad))),]
randomrank9 = randomrank9[rev(order(apply(randomrank9[, c(1,999)], 1,mad))),]
randomrank10 = randomrank10[rev(order(apply(randomrank10[, c(1,999)], 1,mad))),]

pheatmap(rbind(randomrank1, randomrank2, randomrank3, randomrank4, randomrank5, randomrank6, randomrank7, randomrank8, randomrank9, randomrank10), show_rownames = F, show_colnames = F, cluster_cols = F, cluster_rows = F, breaks = breaks,color=colorRampPalette(c('blue4', 'white', 'red'))(100), gaps_row = c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000))

pheatmap(rank1full,cluster_rows = F,  cluster_cols = F, show_rownames = F, show_colnames = F, annotation_row = rank1bindall, annotation_colors = ann_colors, color=colorRampPalette(c('blue4', 'white', 'red'))(50))



##clustering matirx##
rank1bind =cbind(rank1fullid, cluster= sort(cutree(rank1fulpheat$tree_row, h=15.5)))
rank1bindall = cbind(rank1fullid, rank1bind)
rank1bindall = rank1bindall[,c(1,5)]
colnames(rank1bindall) = c('Rank', 'Cluster')
rank1bindall$Cluster = gsub('^', 'cluster',rank1bindall$Cluster)

##Assign color to each cluster##
clucolorpan = brewer.pal(n=7, name='Dark2')
cluid = c(paste0('cluster', 1:7))
names(clucolorpan) = cluid

ann_colors = list(
  Rank = c(Rank1 = "#1B9E77", Rank10 = "#D95F02"),
  Cluster = clucolorpan
)

#pheatmap(rank1full,cluster_rows = F,  cluster_cols = F, show_rownames = F, show_colnames = F, annotation_row = rank1bindall, annotation_colors = ann_colors)


pheatmap(rank1full,cluster_rows = F,  cluster_cols = F, show_rownames = F, show_colnames = F, annotation_row = rank1bind, annotation_colors = ann_colors, color=colorRampPalette(c('blue4', 'white', 'red'))(50))


#rank1fulpheat = pheatmap(rank1fullcor ,cluster_rows = T,  cluster_cols = T, show_rownames = F, show_colnames = F, annotation_row = rank1bind ,  color=colorRampPalette(c('blue4', 'white', 'red'))(120), breaks  = breaks, annotation_col =rank1bind)


###chip vs m6A##
#plus DE#
DEpaclusort = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/DE_ar500.plusTonlyclu.sortbysrr1373.xls", sep = "\t", row.names = 1)
pheatmap(DEpaclusort , show_colnames = F, show_rownames = F, cluster_rows = F, cluster_cols = F,color=colorRampPalette(c('blue4', 'white', 'red'))(50), clustering_distance_rows = "correlation")
DEpaclum = as.data.frame(colMeans(DEpaclusort))
colnames(DEpaclum) = 'DE_P'
ggplot()+geom_line(data = DEpaclum, aes(x=seq(1:1001), y = DE_P), color='red')+theme_bw()

DEpacluchip = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/DE_ar500.plusTonlyclu.srrchip.xls", sep = "\t")
DEpacluchipm = as.data.frame(colMeans(DEpacluchip[,c(407:806)] ))
pheatmap(DEpacluchip[,c(407:806)] , show_colnames = F, show_rownames = F, cluster_rows = F, cluster_cols = F,color=colorRampPalette(c('blue4', 'white', 'red'))(100), clustering_distance_rows = "correlation", breaks = breaks)
colnames(DEpacluchipm) = 'DE_Pchip'
ggplot()+geom_line(data = DEpacluchipm, aes(x=seq(1:400), y = DE_Pchip), color='red')+theme_bw()


#minus DE#
DEmaclusort = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB0003l
                         52/DE_vs_PETssgene/DE_ar500.minusTonlyclu.sortbysrr1373.xls", sep = "\t", row.names = 1)
pheatmap(DEmaclusort , show_colnames = F, show_rownames = F, cluster_rows = F, cluster_cols = F,color=colorRampPalette(c('blue4', 'white', 'red'))(50), clustering_distance_rows = "correlation")
DEmaclum = as.data.frame(colMeans(DEmaclusort))
colnames(DEmaclum) = 'DE_M'
ggplot()+geom_line(data = DEmaclum, aes(x=seq(1:1001), y = DE_M), color='blue')+theme_bw()

DEmacluchip = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/DE_ar500.minusTonlyclu.srrchip.xls", sep = "\t")
DEmacluchipm = as.data.frame(colMeans(DEmacluchip[,c(407:806)] ))
pheatmap(DEmacluchip[,c(407:806)] , show_colnames = F, show_rownames = F, cluster_rows = F, cluster_cols = F,color=colorRampPalette(c('blue4', 'white', 'red'))(100), clustering_distance_rows = "correlation", breaks = breaks)
colnames(DEmacluchipm) = 'DE_Mchip'
ggplot()+geom_line(data = DEmacluchipm, aes(x=seq(1:400), y = DE_Mchip), color='blue')+theme_bw()


#plus PE#
PEpaclusort = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/PE_ar500.plusTonlyclu.sortbysrr1373.xls", sep = "\t", row.names = 1)
pheatmap(PEpaclusort , show_colnames = F, show_rownames = F, cluster_rows = F, cluster_cols = F,color=colorRampPalette(c('blue4', 'white', 'red'))(50), clustering_distance_rows = "correlation")
PEpaclum = as.data.frame(colMeans(PEpaclusort))
colnames(PEpaclum) = 'PE_P'
ggplot()+geom_line(data = PEpaclum, aes(x=seq(1:1001), y = PE_P), color='red')+theme_bw()

PEpacluchip = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/PE_ar500.plusTonlyclu.srrchip.xls", sep = "\t")
PEpacluchipm = as.data.frame(colMeans(PEpacluchip[,c(407:806)] ))
pheatmap(PEpacluchip[,c(407:806)] , show_colnames = F, show_rownames = F, cluster_rows = F, cluster_cols = F,color=colorRampPalette(c('blue4', 'white', 'red'))(100), clustering_distance_rows = "correlation", breaks = breaks)
colnames(PEpacluchipm) = 'PE_Pchip'
ggplot()+geom_line(data = PEpacluchipm, aes(x=seq(1:400), y = PE_Pchip), color='red')+theme_bw()


#minus PE#
PEmaclusort = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/PE_ar500.minusTonlyclu.sortbysrr1373.xls", sep = "\t", row.names = 1)
pheatmap(PEmaclusort , show_colnames = F, show_rownames = F, cluster_rows = F, cluster_cols = F,color=colorRampPalette(c('blue4', 'white', 'red'))(50), clustering_distance_rows = "correlation")
PEmaclum = as.data.frame(colMeans(PEmaclusort))
colnames(PEmaclum) = 'PE_M'
ggplot()+geom_line(data = PEmaclum, aes(x=seq(1:1001), y = PE_M), color='blue')+theme_bw()

PEmacluchip = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/PE_ar500.minusTonlyclu.srrchip.xls", sep = "\t")
PEmacluchipm = as.data.frame(colMeans(PEmacluchip[,c(407:806)] ))
pheatmap(PEmacluchip[,c(407:806)] , show_colnames = F, show_rownames = F, cluster_rows = F, cluster_cols = F,color=colorRampPalette(c('blue4', 'white', 'red'))(100), clustering_distance_rows = "correlation", breaks = breaks)
colnames(PEmacluchipm) = 'PE_Mchip'
ggplot()+geom_line(data = PEmacluchipm, aes(x=seq(1:400), y = PE_Mchip), color='blue')+theme_bw()


##

fmtfile <- function(input, output, grp){
  cludemgt = read.table(input, row.names = 1, sep = "\t")
  output =  as.data.frame(colMeans(cludemgt))
  output$Group = grp
  output$Pos = seq(1:1001)
  colnames(output) = c('NormC', 'Group', 'Pos')
  return(output)
}

DEm50m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt50.xls", "DE_clu50", 'DEm50m')
DEm60m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt60.xls", "DE_clu60", 'DEm60m')
DEm70m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt70.xls", "DE_clu70", 'DEm70m')
DEm80m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt80.xls", "DE_clu80", 'DEm80m')
DEm90m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt90.xls", "DE_clu90", 'DEm90m')
DEm100m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt100.xls", "DE_clu100", 'DEm100m')
DEm110m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt110.xls", "DE_clu110", 'DEm110m')
DEm120m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt120.xls", "DE_clu120", 'DEm120m')
DEm130m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt130.xls", "DE_clu130", 'DEm130m')
DEm140m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt140.xls", "DE_clu140", 'DEm140m')
DEm150m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt150.xls", "DE_clu150", 'DEm150m')
DEm160m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt160.xls", "DE_clu160", 'DEm160m')
DEm170m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt170.xls", "DE_clu170", 'DEm170m')
DEm180m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt180.xls", "DE_clu180", 'DEm180m')
DEm190m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt190.xls", "DE_clu190", 'DEm190m')
DEm200m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.minusTonlyclu.srrchip.clugt200.xls", "DE_clu200", 'DEm200m')

DEmsitecbind = rbind(DEm50m, DEm60m, DEm70m, DEm80m, DEm90m, DEm100m, DEm110m, DEm120m, DEm130m, DEm140m, DEm150m, DEm160m, DEm170m, DEm180m, DEm190m, DEm200m)

DEp50m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt50.xls", "DE_clu50", 'DEp50m')
DEp60m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt60.xls", "DE_clu60", 'DEp60m')
DEp70m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt70.xls", "DE_clu70", 'DEp70m')
DEp80m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt80.xls", "DE_clu80", 'DEp80m')
DEp90m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt90.xls", "DE_clu90", 'DEp90m')
DEp100m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt100.xls", "DE_clu100", 'DEp100m')
DEp110m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt110.xls", "DE_clu110", 'DEp110m')
DEp120m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt120.xls", "DE_clu120", 'DEp120m')
DEp130m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt130.xls", "DE_clu130", 'DEp130m')
DEp140m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt140.xls", "DE_clu140", 'DEp140m')
DEp150m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt150.xls", "DE_clu150", 'DEp150m')
DEp160m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt160.xls", "DE_clu160", 'DEp160m')
DEp170m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt170.xls", "DE_clu170", 'DEp170m')
DEp180m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt180.xls", "DE_clu180", 'DEp180m')
DEp190m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt190.xls", "DE_clu190", 'DEp190m')
DEp200m = fmtfile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/m6AGroup_compara/clu_size/DE_ar500.plusTonlyclu.srrchip.clugt200.xls", "DE_clu200", 'DEp200m')

DEpsitecbind = cbind(DEp50m, DEp60m, DEp70m, DEp80m, DEp90m, DEp100m, DEp110m, DEp120m, DEp130m, DEp140m, DEp150m, DEp160m, DEp170m, DEp180m, DEp190m, DEp200m)


rollrowinmatrix <- function(x){
  df = list()
  for(i in 1:nrow(x)){
    #print(i)
    df[[i]] <- rollmean(as.numeric(x[i,]), 10)
  }
  df = as.data.frame(do.call(rbind, df))
  row.names(df) = row.names(x)
  df = as.matrix(df)
  return(df)
}



##one DE##
mvspplot <- function(x){
DE145m = DEmsm[grep(paste0(x,'$'), row.names(DEmsm)),]
DE145p = DEpsm[grep(paste0(x,'$'), row.names(DEpsm)),]
if(nrow(DE145m) >= 2 && nrow(DE145p) >= 2 ){
DE145msmooth = rollrowinmatrix(DEmsm[rownames(DE145m),])
DE145psmooth = rollrowinmatrix(DEpsm[rownames(DE145p),])



DE145msmplot = pheatmap(DE145msmooth[ names(sort(rowSums(DE145msmooth[, c(400:600)]))),], show_rownames = T, show_colnames = F, cluster_rows = F, cluster_cols = F, color=colorRampPalette(c('blue4', 'white', 'red'))(100))

DE145psmplot = pheatmap(DE145psmooth[ names(sort(rowSums(DE145psmooth[, c(400:600)]))),], show_rownames = T, show_colnames = F, cluster_rows = F, cluster_cols = F, color=colorRampPalette(c('blue4', 'white', 'red'))(100))


ggsave2(paste0(x,'.pvsm.png'), plot_grid(DE145msmplot[[4]], DE145psmplot[[4]], ncol=1, rel_heights = c(nrow(DE145msmooth), nrow(DE145psmooth))))
}
}

##culmulative sum##
#culall = list()
mvspcul <- function(xtag){
  culall = c()
  DE145m = DEmsm[grep(paste0(xtag,'$'), row.names(DEmsm)),]
  DE145p = DEpsm[grep(paste0(xtag,'$'), row.names(DEpsm)),]
  if(nrow(DE145m) >= 6 && nrow(DE145p) >= 6 ){
    DE145mcm = as.data.frame(colMeans(rollrowinmatrix((DE145m))))
    DE145pcm = as.data.frame(colMeans(rollrowinmatrix((DE145p))))
    colnames(DE145mcm) = 'NormC'
    colnames(DE145pcm) = 'NormC'
    DE145mcm$Pos = seq(1:nrow(DE145mcm))
    DE145pcm$Pos = seq(1:nrow(DE145pcm))
    DE145mcm$Group = 'Minus_T'
    DE145pcm$Group = 'Plus_T'
    DE145rb = rbind(DE145mcm, DE145pcm)
    lineplot = ggplot(DE145rb, aes(x= Pos,  y = NormC, color=Group))+geom_line()+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0), limits = c(0,0.6))+scale_color_manual(values = c('blue', 'red'))+theme_cowplot()
    #ggsave2(paste0(x,'.pvsm.line.png'),lineplot)
    
    DE145msmooth = rollrowinmatrix(DEmsm[rownames(DE145m),])
    DE145psmooth = rollrowinmatrix(DEpsm[rownames(DE145p),])
    
    DE145msmplot = pheatmap(DE145msmooth[ names(sort(rowSums(DE145msmooth[, c(400:600)]))),], show_rownames = F, show_colnames = F, cluster_rows = F, cluster_cols = F, breaks = breakcolor, color=colorRampPalette(c('blue4', 'white', 'red'))(100))
    DE145psmplot = pheatmap(DE145psmooth[ names(sort(rowSums(DE145psmooth[, c(400:600)]))),], show_rownames = F, show_colnames = F, cluster_rows = F, cluster_cols = F, breaks = breakcolor, color=colorRampPalette(c('blue4', 'white', 'red'))(100))
    ##epimarkplot##
    h3k27acDEt = as.data.frame(t(h3k27acDE[xtag,7:106]))
    h3k4m1DEt = as.data.frame(t(h3k4m1DE[xtag,7:106]))
    cenpaDEt = as.data.frame(t(cenpaDE[xtag,7:106]))
    hoxa9DEt = as.data.frame(t(hoxa9DE[xtag,7:106]*2))
    hoxa9DEt$Group = 'Hoxa9'
    cenpaDEt$Group = 'CEPBA'
    h3k27acDEt$Group = 'H3K27ac'
    h3k4m1DEt$Group = 'H3K4m1'
    
    hoxa9DEt$Pos = seq(1:nrow(hoxa9DEt))
    cenpaDEt$Pos = seq(1:nrow(cenpaDEt))
    h3k27acDEt$Pos = seq(1:nrow(h3k27acDEt))
    h3k4m1DEt$Pos = seq(1:nrow(h3k4m1DEt))
    
    #demarkall = rbind(hoxa9DEt, cenpaDEt, h3k27acDEt, h3k4m1DEt)
    demarkall = rbind(hoxa9DEt, cenpaDEt, h3k27acDEt)
    colnames(demarkall) = c('Markvalue', 'Group', 'Pos')
    demarkplot = ggplot(demarkall, aes(x= Pos, y = Markvalue, color=Group))+geom_line()+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0) , sec.axis = sec_axis( trans=~./2, name="For Hoax9"))+theme_cowplot()
    
    ##save figure##
    ggsave2(paste0(xtag,'.pvsm.lineaddheatmap.png'), plot_grid(demarkplot, lineplot, DE145msmplot[[4]], DE145psmplot[[4]], ncol=1, rel_heights = c(round((nrow(DE145msmooth)+nrow(DE145psmooth))/2), round((nrow(DE145msmooth)+nrow(DE145psmooth))/2),  nrow(DE145msmooth), nrow(DE145psmooth))), bg='white')
    culm = sum(seq(1:length(colMeans(DE145p))) * colMeans(DE145m))
    culp = sum(seq(1:length(colMeans(DE145p))) * colMeans(DE145p))
    culall = c(culm, culp)
  }
  return(culall)
}

culall = list()
for(i in 1:length(DEallname)){
  culall[[DEallname[i]]] = mvspcul(DEallname[i])
}

culalldf = as.data.frame(do.call(rbind, culall))



##load single molecular#
DEmsm = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/DE_ar500.minusT.xls", sep= "\t", row.names = 1)
DEpsm = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/DE_ar500.plusT.xls", sep= "\t", row.names = 1)

DEmsm = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/PE_ar500.minusT.xls", sep= "\t", row.names = 1http://10.185.31.113:8787/graphics/a1842a93-08c9-41f1-9cbc-a363e65422f5.png)
DEpsm = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/PE_ar500.plusT.xls", sep= "\t", row.names = 1)

##load mark##

cenpaDE = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/Epimark_enhar500/mm10.DEar500.cenpa.xls", row.names = 1, sep = "\t")
hoxa9DE = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/Epimark_enhar500/mm10.DEar500.hoxa9.xls", row.names = 1, sep = "\t")
h3k27acDE = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/Epimark_enhar500/mm10.DEar500.H3K27ac.xls", row.names = 1, sep = "\t")
h3k4m1DE = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/Epimark_enhar500/mm10.DEar500.H3K4m1.xls", row.names = 1, sep = "\t")

cenpaDE = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/Epimark_enhar500/mm10.PEar500.cenpa.xls", row.names = 1, sep = "\t")
hoxa9DE = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/Epimark_enhar500/mm10.PEar500.hoxa9.xls", row.names = 1, sep = "\t")
h3k27acDE = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/Epimark_enhar500/mm10.PEar500.H3K27ac.xls", row.names = 1, sep = "\t")
h3k4m1DE = read.table("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/DE_vs_PETssgene/Epimark_enhar500/mm10.PEar500.H3K4m1.xls", row.names = 1, sep = "\t")



h3k27acDEt = as.data.frame(t(h3k27acDE['DE733',7:106]))
h3k4m1DEt = as.data.frame(t(h3k4m1DE['DE733',7:106]))
cenpaDEt = as.data.frame(t(cenpaDE['DE733',7:106]))
hoxa9DEt = as.data.frame(t(hoxa9DE['DE733',7:106]))
hoxa9DEt$Group = 'Hoxa9'
cenpaDEt$Group = 'CEPBA'
h3k27acDEt$Group = 'H3K27ac'
h3k4m1DEt$Group = 'H3K4m1'

hoxa9DEt$Pos = seq(1:nrow(hoxa9DEt))
cenpaDEt$Pos = seq(1:nrow(cenpaDEt))
h3k27acDEt$Pos = seq(1:nrow(h3k27acDEt))
h3k4m1DEt$Pos = seq(1:nrow(h3k4m1DEt))

demarkall = rbind(hoxa9DEt, cenpaDEt, h3k27acDEt, h3k4m1DEt)

demarkplot = ggplot(demarkall, aes(x= Pos, y = DE733, color=Group))+geom_line()+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+theme_cowplot()








