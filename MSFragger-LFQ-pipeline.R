# LFQ pipeline after MSFragger database search.
# Protein quantitation is based on spectrum counts reported from MSFragger.
# The generated <combined_protein.tsv> is used. 

# load library
library(tidyverse)
library(RColorBrewer)
library(reshape2)

# set working directory to MSFragger tsv file location.
raw <- read_tsv("combined_protein.tsv")

# filter protein probability 0.99, peptide probability 0.99, no duplicate protein group
raw1<- raw %>% filter(`Protein Probability`>=0.99 & `Top Peptide Probability`>=0.99) %>%
  distinct(`Protein Group`,.keep_all = TRUE)

# To keep protein group information 
raw1<- raw %>% filter(`Protein Probability`>=0.99 & `Top Peptide Probability`>=0.99)


#
raw0<-raw
raw1<-raw

# extract SPC information
raw_spc<- raw1[,grepl("Total Spectral Count",names(raw1))]
raw_spc<-raw_spc[,-1]

# If using protein intensity
raw_int<- raw1[,grepl("Total Intensity",names(raw1))]

# extract protein description information
annot<- raw1[,c(1:14)]
annot['Indistinguishable Proteins']<-raw1['Indistinguishable Proteins']

# protein unique id: `Protein ID` or `Protein`
#raw_spc$access<-annot$`Protein ID`
#raw_int$access<-annot$`Protein ID`

annot<-column_to_rownames(annot,var = 'Protein')

raw_spc$access<-annot$Protein
raw_spc<-column_to_rownames(raw_spc,var = 'access')

# read metadata of the project 
metadata<-read_csv("metadata.csv")

names(raw_spc)<-metadata$name

raw_spc_all<-rownames_to_column(annot,var = 'access') %>%
  left_join(rownames_to_column(raw_spc,var = 'access'))

#raw_spc_all<-raw_spc %>% left_join(annot,by = c("access" = "Protein ID"))

# export cleaning dataset
write_csv(raw_spc_all,"raw_spc_all.csv")


#-----------1.0set parameters for the data analysis------------------------------

name.group<-list('A','B','C') #group names, the names must be also shown in the csv file
reps<-3 # how many replicates in each group
min.quant<-2  # minimum spectrum counts cut off for protein identification/quantification
n.rep<-2 # define the minimum replicate number of each group that identified the protein
na.imput<-0.1 # define the missing value imputation number
#quant<-'Spec'

#----------1.1 prepare other data----------------------------------------

tissue <- gl(length(name.group),reps,labels = name.group) # define group
color.j<-color.factor(name.group,reps) # define groups color

# Manually define color 
#mytissue<-factor(rep(c('notrap_bb','notrap_hela','trap_bb','trap_hela'),times=c(5,5,6,5)))
#mycolor<-factor(rep(c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072'),times=c(5,5,6,5)))

raw<-raw_spc
#raw<-raw_int

#-----------2. data filtering, cleaning-------------------------------------
raw.clean.all<-cleanup.data(raw,name.group,min.quant,n.rep,na.imput) # data filter
raw.clean<-raw.clean.all[,!grepl("count",names(raw.clean.all))]

plot.pro.num(raw_spc,fname = 'spc',color.j,min.quant=2)

plotbox("raw_spc",dat = raw.clean,color = color.j)
plotbar("raw_spc",dat = raw.clean,color = color.j)

#-----------histogram distribution---------------------------------------------

plot_missingval(raw.clean,min.quant,"Msfragger_spc")

# histogram
plot_hist(raw.clean,"Msfragger_spc",bin=20)

#----------3. normalization -------------------------------------------------
Norm.m<-'SUM' # choice "SUM", "MEDIAN", "MEAN"
df_sl <- Norm_SUM(raw.clean, tissue) # data sum normalization

write_csv(rownames_to_column(df_sl_m,var = 'access'),'df_sl_MSFragger_SPC.csv')


library(pheatmap)

plot_cor_heatmap(df_sl,'df_sl_spc')           # samples correlation plot
plot_PCA(df_sl,"norm_all_samples_spc")

#----------4. remove certain samples based on PCA quality control----------

#toMatch <- c("A","B")
#df_sl2<-df_sl[,!grepl(paste(toMatch, collapse="|"),names(df_sl))]

#-----------5. df_ave, df_ratio calculation----------------------------------
df_sl2<-df_sl

df_ave<-cal_ave(df_sl2,name.group)

#metadata_2<-metadata %>% filter(!name %in% toMatch) 
df_ratio<-cal_ratio(df_ave,control = 'C',list('A','B'))

ratio.dis(df_ratio,'SPC','C')

#------------6. statistical analysis------------------------------------------

#metadata_2<-metadata %>% filter(!name %in% toMatch)

df_stat<-single.ANOVA (df_sl2,metadata,group = 'group',annot=annot,df_ratio,raw_spc)

df_stat_b<-df_stat

#------------7. volcano plot-------------------------------------------------------

df_stat_b$log2fc<-log2(df_stat_b$B)
myVolcano(df_stat_b,p_threshold=0.01,fc_threshold=2,"B_C")

df_stat_b$log2fc<-log2(df_stat_b$A)
myVolcano(df_stat_b,p_threshold=0.01,fc_threshold=2,"A_C")

#-------------7.significant filtering and heatmap------------------------------------------

df_sig<-df_stat%>%
  filter(padj<=0.01 & (B>2|B<0.5|A>2|A<0.5))

df_stat <- df_stat %>%
  left_join(df_sig %>% transmute(Accession, sig= 'yes')) %>%
  replace_na(list(sig = 'no'))


df_sig_B<-df_stat%>%
  filter(padj<=0.01 & B>2|B<0.5)
df_sig_A<-df_stat%>%
  filter(padj<=0.01 & A>2|A<0.5)

write_csv(df_stat,"sum_norm_anova_results.csv")



while (!is.null(dev.list())) dev.off()
dev.off()
