#------define color function---------------------------------------------
color.factor<-function(name.group,reps){
  ## this function generate color factor based on groups
  if(length(name.group)==2){
    color.j<-gl(2,reps,labels = c('black','red'))
  }else{
    color.j <- I(brewer.pal(length(name.group), name = 'Set3'))
    color.j <- gl(length(color.j),reps,labels = color.j)
  }
  return(color.j)
}


##--------percentage missing values------------------------------
plot_missingval<-function(df,min.quant,fname){
  
  df[df < min.quant] <- NA
  
  missing.values <- df %>%
    gather(key = "key", value = "val") %>%
    mutate(isna = is.na(val)) %>%
    group_by(key) %>%
    mutate(total = n()) %>%
    group_by(key, total, isna) %>%
    summarise(num.isna = n()) %>%
    mutate(pct = num.isna / total * 100)
  
  levels <-(missing.values  %>% filter(isna == T) )$key
  
  percentage.plot <- missing.values %>%
    ggplot() +
    geom_bar(aes(x = key, 
                 y = pct, fill=isna), 
             stat = 'identity', alpha=0.8) +
    scale_x_discrete(limits = levels) +
    coord_flip() +
    labs(title = "Percentage of missing values", x =
           'Variable', y = "% of missing values")
  
  tiff(str_c(fname,"_missingVal.tiff"),width = 6, height = 5, units = 'in', res = 300, compression = 'lzw')
  print(percentage.plot)
  dev.off()
}

##------------------histogram distribution-----------------------------
plot_hist<-function(df,fname,bin){
  p<-ggplot(gather(df), aes(log2(value))) + 
    geom_histogram(bins = bin) + 
    facet_wrap(~key, scales = 'free_x')
  tiff(str_c(fname,"_histogram_distribution.tiff"),width = 6, height = 6, units = 'in', res = 300, compression = 'lzw')
  print(p)
  dev.off()
}


#----Data clean up for normalization analysis-------------------
cleanup.data<-function(raw,name.group,min.quant,n.rep,na.imput){
  
  raw.uniq<- raw
  raw.uniq[is.na(raw.uniq)]<-na.imput
  raw.uniq[raw.uniq==0]<-na.imput
  
  #colnames(raw.uniq) <- gsub('#Spec ', '', colnames(raw.uniq))
  i<-''
  for (i in name.group) {
    temp=data.frame(rowSums(raw.uniq[,grepl(i,names(raw.uniq))]>=min.quant))
    colnames(temp)<-paste(i,'count',sep = '.')
    raw.uniq=cbind(raw.uniq,temp)
  }
  raw.uniq.spc.all<-raw.uniq[apply(raw.uniq[,(ncol(raw.uniq)-length(name.group)+1):ncol(raw.uniq)]>=n.rep,1,any),]
  
  #raw.uniq.spc.all<-raw.uniq.spc.all[ , order(names(raw.uniq.spc.all))]
  #write_csv(raw.uniq.spc.all, "raw.uniq.spc.all.cleaning.csv")
  
  return(raw.uniq.spc.all)
}

#---------------------protein number----------------------------------------

plot.pro.num<-function(raw.clean,fname,color.j,min.quant){
  df<-raw.clean
  pro.num<-data.frame(protein=colSums(df>= min.quant))
  write_csv(rownames_to_column(pro.num,var = 'sample'),str_c(fname,"_protein.number.csv"))
  
  tiff(str_c(fname,"proteinNumber_barplot.tiff"),width = 10, height = 7, units = 'in', res = 300, compression = 'lzw')
  barplot(pro.num$protein,names=rownames(pro.num),col=color.j,las=2,ylab='proteinNumber',main = "protein number",cex.names=0.5)
  abline(h=4000,col='red')
  dev.off()
  
}
#----------------------boxplot, barplot---------------------------------------------
plotbox<-function(fname,dat,color){
  tiff(str_c(fname,"_boxplot.tiff"),width = 10, height = 7, units = 'in', res = 300, compression = 'lzw')
  boxplot(log2(dat),col=color,las=2, notch = TRUE, main = str_c(fname," data"),ylab='log2(SpC)',cex.axis=0.5)
  dev.off()
}

plotbar<-function(fname,dat,color){
  tiff(str_c(fname,"_barplot.tiff"),width = 10, height = 7, units = 'in', res = 300, compression = 'lzw')
  barplot(colSums(dat),col=color,las=2,ylab='SUM SpC',main = str_c(fname," data"),cex.names=0.5)
  #abline(h=sum(dat$WT_00h_R1_TMT1),col='red')
  #abline(h=2*sum(dat$WT_00h_R1_TMT1),col='red')
  dev.off()
}

#---check sample similarity, correlation analysis---------------

plot_cor_heatmap<-function(df,fname){
  cormat <- round(cor(df),2)
  #head(cormat)
  
  tiff(str_c(fname," heatmap.tiff"),width = 10, height = 7, units = 'in', res = 300, compression = 'lzw')
  pheatmap(cormat)
  dev.off()
}   

#-----clustering PCA plot---------------------------------------------
## need to modify the sample name 
plot_PCA<-function(df,fname){
  
  p <- prcomp(t(df), scale=TRUE)
  #s <- summary(p)
  #screeplot(p)
  #screeplot(p, type="lines")
  
  df_out <- as.data.frame(p$x)
  df_out<- rownames_to_column(df_out,var='rowname')
  df_out<-separate(df_out,col='rowname',into=c('group','rep'),sep = '_')
  percentVar <- round(p$sdev^2 / sum( p$sdev^2 ),digits = 3)
  
  tiff(str_c(fname,"_PCA.tiff"),width = 7, height = 5, units = 'in', res = 300, compression = 'lzw')
  p<-ggplot(df_out, aes(x = PC1, y = PC2)) + 
    geom_point(aes(color = factor(group)),size=4) + 
    geom_text(aes(label=rep),size=2) +
    scale_shape_manual(values=c(15,16,17,18,3,7,8)) +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    ggtitle(str_c(fname,"_PCA"))
  print(p)
  dev.off()
  
}
#----function for SUM simple normalization------------------
library(robustbase)
Norm_Method <- function(df, Norm.m, color = NULL, plot = TRUE) {
  
  
  # compute scaling factors to make colsums match the average sum
  if(Norm.m=='SUM'){
    norm_facs <-  colSums(df,na.rm = TRUE)/mean(colSums(df,na.rm = TRUE))
  }else if(Norm.m=='MEDIAN'){
    df.mat<-as.matrix(df)
    norm_facs <- colMedians(df.mat,na.rm = TRUE)/mean(colMedians(df.mat,na.rm = TRUE))
  }else if(Norm.m=='MEAN'){
    norm_facs <- colMeans(df,na.rm = TRUE)/mean(colMeans(df,na.rm = TRUE))
  }
  
  # cat("SL Factors:\n", sprintf("%-5s -> %f\n", colnames(df), norm_facs))
  norm.factor<-data.frame('sample'=colnames(df),'norm.factor'=norm_facs)
  df_sl  <- sweep(df, 2, norm_facs, FUN = "/")
  
  # visualize results and return data frame
  if(plot == TRUE) {
    tiff("Normalization.boxplot.barplot.tiff",width = 10, height = 10, units = 'in', res = 300, compression = 'lzw')
    par(mfrow=c(2,2))
    boxplot(log10(df), col = color, notch = TRUE, main = "Original data",ylab='log10(SPC)')
    boxplot(log10(df_sl), col = color, notch = TRUE, main = "Normalized data",ylab='log10(SPC)')
    
    if(Norm.m=='SUM'){
      barplot(colSums(df), col = color, main = "Original data",ylab='Sum(SPC)')
      barplot(colSums(df_sl), col = color, main = "Normalized data",ylab='Sum(SPC)')
    } else if(Norm.m=='MEDIAN'){
      df.mat<-as.matrix(df)
      df_sl.mat<-as.matrix(df_sl)
      barplot(colSums(df.mat), col = color, main = "Original data",ylab='sum(SPC)')
      barplot(colSums(df_sl.mat), col = color, main = "Normalized data",ylab='sum(SPC)')
    } else {
      barplot(colMeans(df), col = color, main = "Original data",ylab='Mean(SPC)')
      barplot(colMeans(df_sl), col = color, main = "Normalized data",ylab='Mean(SPC)')
    }
    
    dev.off()
    par(mfrow=c(1,1))
  }
  write_csv(norm.factor,'normalization.factor.csv')
  return(df_sl)
  
}

Norm_SUM <- function(df, color = NULL, plot = TRUE) {
  # This makes each channel sum to the average grand total
  # df - data frame contain spc
  # returns a new data frame with normalized values
  
  # compute scaling factors to make colsums match the average sum
  norm_facs <-  colSums(df)/mean(colSums(df))
  # cat("SL Factors:\n", sprintf("%-5s -> %f\n", colnames(df), norm_facs))
  norm.factor<-data.frame('sample'=colnames(df),'norm.factor'=norm_facs)
  df_sl  <- sweep(df, 2, norm_facs, FUN = "/")
  
  # visualize results and return data frame
  if(plot == TRUE) {
    tiff("Normalization.boxplot.barplot.tiff",width = 10, height = 10, units = 'in', res = 300, compression = 'lzw')
    par(mfrow=c(2,2))
    boxplot(log10(df), col = color, notch = TRUE, main = "Original data",ylab='log10(SPC)',las=2, cex.names=.5)
    boxplot(log10(df_sl), col = color, notch = TRUE, main = "Normalized data",ylab='log10(SPC)',las=2, cex.names=.5)
    barplot(colSums(df), col = color, main = "Original data",ylab='Sum(SPC)',las=2, cex.names=.5)
    barplot(colSums(df_sl), col = color, main = "Normalized data",ylab='Sum(SPC)',las=2, cex.names=.5)
    dev.off()
    par(mfrow=c(1,1))
  }
  write_csv(norm.factor,'normalization.factor.csv')
  return(df_sl)
  
}

#------------ave,ratio analysis----------------------------

cal_ave<-function(df,name.group){
  g_list <- list()
  j<-1
  
  for(i in name.group){
    tmp<- rowMeans(df[,grepl(i,names(df))])
    g_list[[j]] <- tmp
    j=j+1
  }
  
  df_ave<-do.call(cbind.data.frame, g_list)
  colnames(df_ave) <- name.group
  return(df_ave)
}

cal_ratio<-function(df,control,name.group){
  g_list <- list()
  ave<-df[control]
  j<-1
  
  for(i in name.group){
    tmp<- df[,grepl(i,names(df))]/ave
    g_list[[j]] <- tmp
    j=j+1
  }
  
  df_ratio<-do.call(cbind.data.frame, g_list) 
  
  colnames(df_ratio) <- name.group
  
  return(df_ratio)
}

ratio.dis<-function(df,fname,control){
  df_long<-melt(df,id.vars = NULL)
  df_long<-subset(df_long,variable != control)
  
  tiff(str_c(fname,"_ratio_distribution.tiff"),width = 7, height = 5, units = 'in', res = 300, compression = 'lzw')
  p<-ggplot(aes(x=log2(value), colour=variable), data=df_long)+ geom_density()
  
  print(p)
  dev.off()
}

#one-way anova x~group

oneway.ANOVA<-function(x,meta,group){
  m<-cbind(intensity=x,meta)
  c(rownames(x),anova(aov(intensity~group,data = m))$`Pr(>F)`[1])
}

test.ttest<-function(x,meta,group){
  m<-cbind(intensity=x,meta)
  c(rownames(x),t.test(intensity~group,data = m)$p.value)
}

single.ANOVA <- function(df_sl,metadata,group,annot,df_ratio,raw_spc){
  
  group.res<-data.frame(pval=apply(df_sl,1,oneway.ANOVA,metadata,group))
  group.res$padj<-p.adjust(group.res$pval,method='BH',n=nrow(group.res))
  
  df_stat<- rownames_to_column(annot,var='Accession') %>%
    right_join(rownames_to_column(group.res,var='Accession')) %>%
    left_join(rownames_to_column(df_sl,var='Accession'))%>%
    left_join(rownames_to_column(df_ratio,var='Accession'))%>%
    left_join(rownames_to_column(raw_spc,var='Accession'),by='Accession')
  
  write_csv(df_stat,"df_anova_result.csv")
  
  return(df_stat)
}

single_ttest <- function(df_sl,metadata,group,annot,df_ratio,fname){
  
  group.res<-data.frame(pval=apply(df_sl,1,test.ttest,metadata,group))
  group.res$padj<-p.adjust(group.res$pval,method='BH',n=nrow(group.res))
  
  df_stat<- rownames_to_column(annot,var='Accession') %>%
    right_join(rownames_to_column(group.res,var='Accession')) %>%
    left_join(rownames_to_column(df_sl,var='Accession'))%>%
    left_join(rownames_to_column(df_ratio,var='Accession'))
  
  write_csv(df_stat,str_c(fname,"_df_anova_result.csv"))
  
  return(df_stat)
}

myVolcano<-function(data,p_threshold=0.05,fc_threshold=2,fname){
  tiff(str_c(fname,"_volvano.tiff"),width = 6, height = 6, units = 'in', res = 300, compression = 'lzw')
  
  plot(data$log2fc,(-10)*log10(data$padj),col="#00000033",pch=19,las=1,
       xlab="Log2 (Fold change)",
       ylab="-10log10 (adjust P value)"
  )  
  #up<-subset(data,data$pval<p_threshold & data$log2fc>log2(fc_threshold))
  up<-data%>%filter(padj<=p_threshold & (log2fc>log2(fc_threshold)))
    #filter_at(vars(18:20), any_vars(.>=2))%>%
    #filter_at(vars(21:23), any_vars(.>=2))
  
  points(up$log2fc,(-10)*log10(up$padj),col=1,bg=brewer.pal(9,"YlOrRd")[6],pch=21,cex=1.5)
  #text(up[,1],-log10(up[,2]),rownames(up),adj = -0.1)
  #down<-subset(data,data$pval<p_threshold & data$log2fc<(-log2(fc_threshold)))
  down<-data%>%filter(padj<=p_threshold & (log2fc<(-log2(fc_threshold))))
    #filter_at(vars(18:20), any_vars(.>=2))%>%
    #filter_at(vars(21:23), any_vars(.>=2))
  
  points(down$log2fc,(-10)*log10(down$padj),col=1,bg=brewer.pal(11,"RdBu")[9],pch=21,cex=1.5)
  
  abline(h=(-10)*log10(p_threshold),v=c(-log2(fc_threshold),log2(fc_threshold)),lty=2,lwd=1)
  
  dev.off()
  
}


