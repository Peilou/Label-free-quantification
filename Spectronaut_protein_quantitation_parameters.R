# This program is for Spectronaut DIA data quantitation analysis
# DDA acquisition of HighpH 8 fractions for library built-up
# library and dia search was completed by Spectronaut
# export the protein quant result report containing stripped pep and PG.quant
# PG.proteinGroups, PG.Genes,PG.Organisms, PG.ProteinDescripitions

library(tidyverse)
library(RColorBrewer)
library(reshape2)

# set the working directory to export file
raw <- read_csv("DIA_spectronaut_Report.csv")
meta<-read_csv("metadata.csv") # metadata file for the project
 
annot<-raw[,1:4]
#colnames(annot)[1]<-'access'

pep_raw<-raw[,c(1,5:22)]
int_raw<-raw[,c(1,23:40)]
pep_raw<-column_to_rownames(pep_raw,var = 'PG.ProteinGroups')
int_raw<-column_to_rownames(int_raw,var = 'PG.ProteinGroups')
colnames(int_raw)<-meta$names
colnames(pep_raw)<-meta$names

int_df<-int_raw

# determine the intensity cut-off threshold

int_df[int_df == "Filtered" ] <- NA  
int_df<-as.data.frame(sapply(int_df,as.numeric))
row.names(int_df)<-row.names(int_raw)

hist(log2(int_df[, "LP2_R4_462"]), 100, main = "Histogram of log2 intensities", 
     col = "steelblue", border = "steelblue", freq = TRUE)

int_q001<-stack(lapply(int_df[1:18], quantile, na.rm=TRUE, prob = 0.01, names = FALSE))
int_q001_m<-mean(int_q001$values) ## get the mean is 23

# export int_q001 data
write_csv(int_q001,"intensity_quantile.csv")

int_df<-int_df[,order(colnames(int_df))]
raw<-int_df


##--------set parameters for the data analysis------------------------------
groups<-as.list(unique(meta$group))
groups<-groups[order(unlist(groups))]

name.group<-groups #group names, the names must be also shown in the csv file
reps<-6 # how many replicates in each group
min.quant<-23  # the mean determined by quantile analysis above
n.rep<-3 # define the minimum replicate number of each group that identified the protein
na.imput<-0.1 # define the missing value imputation number
#quant<-'Spec'

#----------1.1 prepare other data----------------------------------------

tissue <- gl(length(name.group),reps,labels = name.group) # define group
color.j<-color.factor(name.group,reps) # define groups color

# manually difine the color
#mytissue<-factor(rep(c('notrap_bb','notrap_hela','trap_bb','trap_hela'),times=c(5,5,6,5)))
#mycolor<-factor(rep(c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072'),times=c(5,5,6,5)))

#-----------2. data filtering, cleaning-------------------------------------

raw.clean.all<-cleanup.data(raw,name.group,min.quant,n.rep,na.imput) # data filter
raw.clean<-raw.clean.all[,!grepl("count",names(raw.clean.all))]

write_csv(rownames_to_column(raw.clean, var = 'access') ,"raw.clean.csv")

plot.pro.num(raw.clean,fname = '23int',color.j, min.quant)

plotbox("23int",dat = raw.clean,color = color.j)
plotbar("23int",dat = raw.clean,color = color.j)

#-----------histogram distribution---------------------------------------------

plot_missingval(raw.clean,min.quant,"23int")

# histogram
plot_hist(raw.clean,"23int",bin=100)

#----------3. normalization -------------------------------------------------
Norm.m<-'SUM' # choice "SUM", "MEDIAN", "MEAN"
df_sl <- Norm_SUM(raw.clean, tissue) # data sum normalization

write_csv(rownames_to_column(df_sl,var = 'access'),'df_sl_sum_int23.csv')

plot_PCA(df_sl,"norm_all_samples_int30")

library(pheatmap)
plot_cor_heatmap(df_sl,'df_slc')   

#-----------5. df_ave, df_ratio calculation----------------------------------

df_ave<-cal_ave(df_sl,name.group)
df_ratio<-cal_ratio(df_ave,control = 'LP1',name.group)
ratio.dis(df_ratio,'df_sl','LP1')


#----------------------6. ANOVA analysis----------------

annot<-column_to_rownames(annot,var = 'PG.ProteinGroups')

df_stat<-single.ANOVA(df_sl,meta,group = 'group',annot=annot,df_ratio)


#-------------7.significant filtering and heatmap------------------------------------------

df_sig<-df_stat%>%
  filter(pval<=0.05 & (LP2>1.5|LP2<0.67|LP3>1.5|LP3<0.67))
write_csv(df_sig,"LP_int23_sig_p005_FC1pt5.csv")

#---------------7.2 post hoc test-------------------------------------------
df_avo<-read_csv("df_anova_result.csv")
df_sig<-read_csv("LP_int23_sig_p005_FC1pt5.csv")
df_sig_sl<-df_sig[,c(1,7:24)]
df_sig_sl<-column_to_rownames(df_sig_sl,var = 'Accession')
meta<-read_csv("metadata.csv")
meta<-meta[order(meta$names),]


HSD<-single.HSD(df_sig_sl,meta,"group",df_sig)

HSD_test<-function(x,meta,group){
  m<-cbind(intensity=x,meta)
  lm<-aov(intensity ~group ,data = m)
  out <- HSD.test(lm, "group",alpha = 0.2, group = TRUE, console = TRUE)
  df_hop<-as.data.frame(t(out$groups[2]))
  rownames(df_hop)<-row.names(x)
  df_hop<-df_hop[ , order(names(df_hop))]
  return(df_hop)
}

single.HSD <- function(df_sl,meta,group,df_sig){
  
  group.res<-apply(df_sl,1,HSD_test,meta,group)
  group.res <- do.call("rbind", group.res)
  
  df_stat<- df_sig %>%
    left_join(rownames_to_column(group.res,var='Accession'),by='Accession')
  
  write_csv(df_stat,"df_HSD_result_2.csv")
  
  return(group.res)
}

#-----------7.3 t-test analysis----------------------------------------
df_sl<-df_avo[,c(1,7:24)]
df_sl<-column_to_rownames(df_sl,var = 'Accession')

df_stat<-cal_stat(df_sl,meta,df_avo)
write_csv(df_stat,"df_ttest.csv")

sig_G2<-df_stat%>%filter(pval.y<=0.05 & (LP2>1.5|LP2<0.67))
df_stat<-df_stat%>% 
  left_join(sig_G2 %>% transmute(Accession,sig.G2='yes'))%>%
  replace_na(list(sig.G2 = 'no'))

sig_G3<-df_stat%>%filter(pval.x.x<=0.05 & (LP3>1.5|LP3<0.67))
df_stat<-df_stat%>% 
  left_join(sig_G3 %>% transmute(Accession,sig.G3='yes'),by='Accession')%>%
  replace_na(list(sig.G3 = 'no'))

df_stat$G2_G3<-df_stat$LP3/df_stat$LP2
sig_G2_G3<-df_stat%>%filter(pval.y.y<=0.05 & (G2_G3>1.5|G2_G3<0.67))
df_stat<-df_stat%>% 
  left_join(sig_G2_G3 %>% transmute(Accession,sig.G2.G3='yes'),by='Accession')%>%
  replace_na(list(sig.G2.G3 = 'no'))

cal_stat<-function(df_sl,meta,df_avo){
  df_sl_G2<-df_sl[,c(1:12)]
  df_sl_G3<-df_sl[,c(1:6,13:18)]
  df_sl_G2_G3<-df_sl[,c(7:18)]
  
  meta_G2<-meta %>% filter(group %in% c("LP1","LP2"))
  meta_G3<-meta %>% filter(group %in% c("LP1","LP3"))
  meta_G2_G3<-meta %>% filter(group %in% c("LP2","LP3"))
  
  df_stat_G2<-single_ttest(df_sl_G2,meta_G2,group = 'group')
  df_stat_G3<-single_ttest(df_sl_G3,meta_G3,group = 'group')
  df_stat_G2_G3<-single_ttest(df_sl_G2_G3,meta_G2_G3,group = 'group')
  
  df_stat<-df_avo %>% left_join(rownames_to_column(df_stat_G2,var='Accession'),by='Accession')%>%
    left_join(rownames_to_column(df_stat_G3,var='Accession'),by='Accession')%>%
    left_join(rownames_to_column(df_stat_G2_G3,var='Accession'),by='Accession')
  return(df_stat)
  
}

single_ttest <- function(df_sl,metadata,group){
  
  group.res<-data.frame(pval=apply(df_sl,1,test.ttest,metadata,group))
  group.res$padj<-p.adjust(group.res$pval,method='BH',n=nrow(group.res))
  
  df_stat<- group.res
  
  return(df_stat)
}


#------------7. volcano plot-------------------------------------------------------
df<-read_csv("t-test.csv")
df_stat_b<-df

df_stat_b$log2fc<-log2(df_stat_b$`LP2/LP1`)
df_stat_b$pval<-df_stat_b$pval.LP2_LP1
df_stat_b$sig<-df_stat_b$`?sig.LP2/LP1`
myVolcano(df_stat_b,"g2")

df_stat_b$log2fc<-log2(df_stat_b$`LP3/LP1`)
df_stat_b$pval<-df_stat_b$pval.LP3_LP1
df_stat_b$sig<-df_stat_b$`?sig.LP3/LP1`
myVolcano(df_stat_b,"g3")

df_stat_b$log2fc<-log2(df_stat_b$`LP3/LP2`)
df_stat_b$pval<-df_stat_b$pval.LP3_LP2
df_stat_b$sig<-df_stat_b$`?sig.LP3/LP2`
myVolcano(df_stat_b,"g2_g3")


myVolcano<-function(data,fname){
  tiff(str_c(fname,"_volvano.tiff"),width = 6, height = 6, units = 'in', res = 300, compression = 'lzw')
  
  plot(data$log2fc,(-10)*log10(data$pval),col="#00000033",pch=19,las=1,
       xlab="Log2 (Fold change)",
       ylab="-10log10 (P value)"
  )  
  #up<-subset(data,data$pval<p_threshold & data$log2fc>log2(fc_threshold))
  up<-data%>%filter(sig=='yes' & (log2fc>0))
  
  points(up$log2fc,(-10)*log10(up$pval),col=1,bg=brewer.pal(9,"YlOrRd")[6],pch=21,cex=1.5)
  #text(up[,1],-log10(up[,2]),rownames(up),adj = -0.1)
  #down<-subset(data,data$pval<p_threshold & data$log2fc<(-log2(fc_threshold)))
  down<-data%>%filter(sig=='yes' & (log2fc<0))
  
  points(down$log2fc,(-10)*log10(down$pval),col=1,bg=brewer.pal(11,"RdBu")[9],pch=21,cex=1.5)
  
  abline(h=(-10)*log10(0.05),v=c(-log2(1.5),log2(1.5)),lty=2,lwd=1)
  
  dev.off()
  
}

#df_vol<-myVolcanoData(df_stat_b,p_threshold=0.05,fc_threshold=2)
#myVolcano2(df_vol,p_threshold=0.05,fc_threshold=2,"brain_int")
#----------8. venn diagram------------------------------
library(VennDiagram)
myCol <- c('red','blue','black')
myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(sig_G2$Accession, sig_G3$Accession,sig_G2_G3$Accession),
  category.names = c("LP2/LP1" , "LP2/LP1","LP3/LP2"),
  filename = 'venn_diagramm.png',
  output=TRUE,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
)
##----------------9. heatmap--------------------------------------------
sig<-df %>% filter(`?sig.LP2/LP1`=='yes' | `?sig.LP3/LP1`=='yes'| `?sig.LP3/LP2`=='yes')

sig_m<-sig[,c(1,21:38)]
sig_m<-column_to_rownames(sig_m,var = 'Accession')
sig_m$LP1<-rowMeans(sig_m[,1:6])
sig_m$LP2<-rowMeans(sig_m[,7:12])
sig_m$LP3<-rowMeans(sig_m[,13:8])

sig_heat<-sig_m[,19:21]

library(gplots)

plot_heatmap(sig_heat,myCol)

#heatmap.2(data.matrix(sig_heat))



while (!is.null(dev.list())) dev.off()
dev.off()
