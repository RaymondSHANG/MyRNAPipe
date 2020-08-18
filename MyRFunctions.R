myPCA <-function(mat_set_f,title_f ="My PCA", condition_f,log_f=FALSE, center_f=TRUE, scale_f=FALSE,pc_x_f=1L,pc_y_f=2L,color_set_f=NA){
  #Notice, the colnames of mat_set_f should be the sample names, the rownames of mat_set_f should be gene names
  #pc_x_f=1L
  #pc_y_f=2L
  if(log_f){
    mat_set_f <- log2(mat_set_f+1)
  }
  mat_set_f.pca <- prcomp(t(mat_set_f),center = center_f,scale. = scale_f) 
  #computation of variances
  eigenvalues <- (mat_set_f.pca$sdev) ^ 2
  var_explained <- eigenvalues * 100 / sum(eigenvalues)
  
  #set label names
  x_lab <- paste0('PC',pc_x_f,' (', round(var_explained[pc_x_f],digits=1))
  x_lab <- paste0(x_lab, '%)')
  y_lab <- paste0('PC',pc_y_f,' (', round(var_explained[pc_y_f],digits=1))
  y_lab <- paste0(y_lab, '%)')
  
  #Extract PC1 and PC2 to pcs, you can also change this PC1 and PC2 to others
  pcs_set_f <- data.frame(mat_set_f.pca$x[, c(pc_x_f, pc_y_f)])
  pcs_set_f$Sample <- rownames(pcs_set_f)
  rownames(pcs_set_f) <- NULL
  #add 'Group' information from experimental design
  pcs_set_f <- dplyr::left_join(pcs_set_f, condition_f,by = 'Sample')
  #ggplot
  pc_x_f <- paste0('PC', pc_x_f)
  pc_y_f <- paste0('PC', pc_y_f)
  p1 <-ggplot(pcs_set_f, aes_string(pc_x_f, pc_y_f, colour = 'Group'))+  
    geom_point(size = 5, alpha = 0.8)  +
    ggtitle(title_f)+
    xlab(x_lab)+
    ylab(y_lab)+
    theme(text = element_text(size=10),axis.text = element_text(size=9),plot.title = element_text(size = 10)) +
    theme(plot.title = element_text(hjust = 0.5))
  if(! is.na(color_set_f[1])){
    p1 <- p1 + scale_color_manual(values=c(color_set_f))
  }
  result <- list(plot = p1,data = mat_set_f.pca)
  return (result)
}

myPCAerrorBar <- function(myPCA_f,groupvars_f="Group",title_f="PCA with error bar",bar_f="se",color_set_f=NA){
  library("Rmisc")
  if(FALSE){
    pcs_f <-test$plot$data#p2$data
    groupvars_f <- "condition"
    title_f <-"PCA with error bar"
    bar_f <-"se"
  }
  pcs_f <- myPCA_f$plot$data
  label_x_f <- colnames(pcs_f)[1]
  label_y_f <- colnames(pcs_f)[2]
  se_x_f <- summarySE(pcs_f,measurevar = label_x_f,groupvars=groupvars_f,conf.interval = 0.95)
  se_y_f <- summarySE(pcs_f,measurevar = label_y_f,groupvars=groupvars_f,conf.interval = 0.95)
  temp_x <- dplyr::select(se_x_f,Group = groupvars_f,PCx=label_x_f,ci_PCx='ci',se_PCx='se',sd_PCx='sd')
  #colnames(temp_x)[2:5] <- gsub("PC1",label_x_f, colnames(temp_x)[2:5])
  colnames(temp_x)[1] <- groupvars_f
  temp_y <- dplyr::select(se_y_f,Group=groupvars_f,PCy=label_y_f,ci_PCy='ci',se_PCy='se',sd_PCy='sd')
  #colnames(temp_y)[2:5] <- gsub("PC1",label_y_f, colnames(temp_y)[2:5])
  colnames(temp_y)[1] <- groupvars_f 
  
  df_f_group <- dplyr::left_join(temp_x, temp_y,by=groupvars_f)
  
  title_f <- paste0(title_f," (",bar_f,")")
  ybar <- 0
  xbar <-0
  
  ybar <- as.name(paste0(bar_f,"_PCy"))
  xbar <- as.name(paste0(bar_f,"_PCx"))
  
  eigenvalues <- (myPCA_f$data$sdev) ^ 2
  var_explained <- eigenvalues * 100 / sum(eigenvalues)
  x_lab <- paste0(label_x_f,' (', round(var_explained[as.numeric(gsub("PC","",label_x_f))],digits=1))
  x_lab <- paste0(x_lab, '%)')
  y_lab <- paste0(label_y_f,' (', round(var_explained[as.numeric(gsub("PC","",label_y_f))],digits=1))
  y_lab <- paste0(y_lab, '%)')
  
  p<-ggplot(data = df_f_group,aes_string( x="PCx", y="PCy",colour=groupvars_f)) + 
    geom_point(size=7,alpha=0.5) + 
    geom_errorbar(aes(ymin = PCy-eval(ybar),ymax = PCy+eval(ybar))) + 
    geom_errorbarh(aes(xmin = PCx-eval(xbar),xmax = PCx+eval(xbar)))+
    ggtitle(title_f)+
    xlab(x_lab)+
    ylab(y_lab)+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(text = element_text(size=10),axis.text = element_text(size=9),plot.title = element_text(size = 10))
  
  if(! is.na(color_set_f[1])){
    p <- p + scale_color_manual(values=c(color_set_f))
    p <- p + scale_fill_manual(values=c(color_set_f))
    
  }
  return (p)
}


myPCAplot <- function(geneset_f,name_f,shortname_f,condition_f,vsd_data_f,directory_f,log_f=FALSE,center_f=TRUE,scale_f=FALSE,pc_x_f=1L,pc_y_f=2L,color_set_f=NA){
  #
  if(FALSE){
    geneset_f<- geneset2
    name_f <- geneset_name2
    shortname_f <- "test"
    condition_f<-conditions
    vsd_data_f<-vsd_data
    directory_f<-directory_save
  }
  
  mat_set_vsd_f <- vsd_data_f[geneset_f,condition_f$Sample]
  title_f <- paste0("PCA based on normalized expressions\n",name_f,"")
  p_1 <- myPCA(mat_set_f = mat_set_vsd_f ,title_f=title_f,condition_f=condition_f, log_f=log_f,center_f=center_f,scale_f = scale_f,pc_x_f=pc_x_f,pc_y_f=pc_y_f,color_set_f=color_set_f)
  p_2<-myPCAerrorBar(p_1,title_f = title_f ,bar_f="se",color_set_f=color_set_f)
  filename1 <- paste0(directory_f,"/PCA_",shortname_f,".jpg")
  filename2 <- paste0(directory_f,"/PCA_",shortname_f,"_ErrorBar",".jpg")
  ggsave(filename1,plot=p_1$plot, device = "jpeg",width = 12, height = 8, units = "cm",dpi = 300)
  ggsave(filename2,plot=p_2, device = "jpeg",width = 12, height = 8, units = "cm",dpi = 300)
}




###
library(genefilter)
library(ggplot2)
library(ggrepel)
#install.packages("ggrepel")
#mycolors
plotPCA.Raymond <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE,pc_1=1,pc_2=2,label=FALSE) 
{
  #
  #object=vsd
  #intgroup="summary"
  #ntop=500
  #pc_1=2
  #pc_2=3
  if(FALSE){
    object=vsd
    intgroup="Summary"
    ntop=500
    pc_1=2
    pc_2=3
    returnData=FALSE
  }
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  if (length(intgroup) > 1) {
    group <-factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }  else {
    group <-colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, pc_1], PC2 = pca$x[, pc_2], group = group, 
                  intgroup.df, name = colData(object)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(pc_1,pc_2)]
    return(d)
  }
  # ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0(paste0("PC",pc_1,": "), round(percentVar[pc_1] * 100), "% variance")) + ylab(paste0(paste0("PC",pc_2,": "), round(percentVar[pc_2] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
  if(label){
    label_tag="name"
  }else{
    label_tag=NA
  }
  p <- ggplot(data = d, aes_string("PC1", "PC2", color=intgroup)) +
    geom_point(size=3) +
    xlab(paste0(paste0("PC",pc_1,": "),round(percentVar[pc_1] * 100),"% variance")) +
    ylab(paste0(paste0("PC",pc_2,": "),round(percentVar[pc_2] * 100),"% variance")) + 
    coord_fixed()+geom_text(aes_string(label = label_tag, color=intgroup),hjust=0, vjust=0)
  #scale_color_manual(values = mycolors)
  return (p)
  
}