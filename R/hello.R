# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
library(GSVA)
library(ggplot2)
library(ggpubr)
library(e1071)
library(survival)
library(survminer)
library(KMsurv)
library(ggpubr)
library(gridExtra)
library(easyGgplot2)

T_FRE<<-c()
sample_TIL=c()
group_fre=c()
TAG=0
ICB_response=NULL
group_content=c()

#' @title self-build reference file
#' @description
#' @details Please input the expression matrix separated by tab, rownames of the file must be gene symbol.
#' @param sig_file self-build cell type gene signature file.
#' @return A list of cell type gene signatures.

self_sig=function(sig_file){
  ref_sig=read.csv(sig_file,row.names=1,fill = TRUE,check.names = F,header=F)
  paper_marker<-list()
  #compensation_matrix<<-compensation_matrix[,-6]
  for(cell in row.names(ref_sig)){paper_marker[[cell]]<-unique(as.vector(unlist(ref_sig[cell,])))}
  return(paper_marker)
}


#' Title
#'
#' @param sample
#' @param data_type
#' @param group_tag
#' @param response_tag
#' @param customer
#' @param sig_file
#' @param exp_file
#' @param save_abun
#' @param save_group
#' @param save_icb
#' @param save_plot
#' @param save_bar
#'
#' @return

getResult = function(sample,data_type,group_tag,response_tag,customer,sig_file=NULL,exp_file=NULL,save_abun="Sample_abundance",save_group="Group_abundance",save_icb="ICB_response",save_plot="Group_comparison.png",save_bar="Abundance_barplot.png"){
  print(save_bar)
  if (customer==TRUE){
    paper_marker<-self_sig(sig_file)
    marker_exp=read.csv(exp_file,row.names=1,fill = TRUE,check.names = F,header=T)
  }
  sample=read.table(sample,header=T,sep="\t",check.names = F)
  sample=as.matrix(sample)
  if(ncol(sample)==2){
    sam_name=colnames(sample)[2]
    row.names(sample)=sample[,1]
    sample=as.matrix(sample[,-1])
    colnames(sample)=sam_name
  }else{
    row.names(sample)=sample[,1]
    sample=sample[,-1]
  }
  #  sample=as.matrix(sample)
  #  row.names(sample)=sample[,1]
  #  sample=sample[,-1]
  group_index=0
  if (group_tag){
    #group_index=as.numeric(as.vector(unlist(grep("group",row.names(sample)))))
    group_column<-sample[1,]
    group_content<<-sample[1,]
    sample=sample[-1,]
  }
  sam = apply(sample,2,as.numeric)
  row.names(sam) = row.names(sample)
  tt = intersect(row.names(sam),as.vector(unlist(paper_marker)))
  genes = intersect(tt,row.names(marker_exp))
  sam_exp = as.matrix(sam[genes,])
  colnames(sam_exp) = colnames(sample)
  marker_exp = marker_exp[genes,]
  marker_tag_mat = c()
  for(cell in names(paper_marker)){
    tag = marker_tag(genes,as.vector(unlist(paper_marker[cell])))
    marker_tag_mat = cbind(marker_tag_mat,tag)
  }
  row.names(marker_tag_mat) = row.names(marker_exp)
  colnames(marker_tag_mat) = names(paper_marker)
  # progress$set(value = 20,detail = "Reconstructing sample ratio profile")
  # Sys.sleep(0.1)
  #  shinyWidgets::updateProgressBar(title = 'Reconstructing sample ratio profile',id = "pb2",value=20,session=getDefaultReactiveDomain())
  exp_new = apply(sam_exp,2,function(x) sample_ratio(x,marker_exp,marker_tag_mat,data_type))
  # progress$set(value = 30,detail = "ssgsea for abundance prediction")
  # Sys.sleep(0.1)
  # shinyWidgets::updateProgressBar(title = 'ssgsea for abundance prediction',id = "pb2",value=40,session=getDefaultReactiveDomain())
  result = gsva(exp_new,paper_marker,method="ssgsea",ssgsea.norm=TRUE,parallel.sz=20)
  if(ncol(result)<3){
    result[which(result<0)]=0
  }else{
    result = result - apply(result,1,min)
    #result[which(result<0)]=0
  }
  compensation_matrix_num = apply(compensation_matrix,2,as.numeric)
  # progress$set(value = 20,detail = "Adjusting result by Compensation matrix")
  # incProgress(0.2, detail = "Immune infiltration calculating")
  # Sys.sleep(0.5)
  # shinyWidgets::updateProgressBar(title = 'Immune infiltration calculating',id = "pb2",value=70,session=getDefaultReactiveDomain())
  if(customer==0){
    row.names(compensation_matrix_num) = row.names(compensation_matrix)
    result_norm = compensation(result,compensation_matrix_num)
    #result_norm = norm_func(result_adjust)
    # incProgress(0.2, detail = "Immune infiltration calculating")
    #Sys.sleep(0.5)
    #shinyWidgets::updateProgressBar(title = 'Adjusting result by Compensation matrix',id = "pb2",value=90,session=getDefaultReactiveDomain())
    if(ncol(result_norm)==1){
      InfiltrationScore=sum(result_norm)
    }else{
      # print("error test")
      if(nrow(result_norm)==24){
        InfiltrationScore = apply(result_norm[c('Bcell','CD4_T','CD8_T','DC','Macrophage','Monocyte','Neutrophil','NK'),],2,sum)
      }else{
        InfiltrationScore = apply(result_norm,2,sum)
      }
    }
    InfiltrationScore = (InfiltrationScore/max(InfiltrationScore))*0.9
    result_mat = rbind(result_norm,InfiltrationScore)
  }else{
    result_mat=result
    # print(head(result_norm))
  }
  if(ncol(result_mat)>1){
    result_mat=apply(result_mat,1,function(x) round(x,3))
  }else{
    result_mat=t(round(result_mat,3))
  }

  if(group_tag){
    group_name<-sort(unique(as.vector(unlist(group_column))))
    p_va=c()
    group_column<-as.numeric(as.factor(as.vector(unlist(group_column))))
    result_group=cbind(result_mat,group_column)
    result_tt=apply(result_group,2,as.numeric)
    if (length(group_name)>2){
      for (cell in colnames(result_group)){
        result_group_new<-result_group[,c(cell,"group_column")]
        t=aov(group_column~.,data.frame(result_group_new))
        p_va=c(p_va,round(summary(t)[[1]][["Pr(>F)"]],2))
      }
    }else{
      g1_index=grep(1,group_column)
      g2_index=grep(2,group_column)
      for (cell in colnames(result_mat)){
        c_=wilcox.test(result_mat[g1_index,cell],result_mat[g2_index,cell])
        p_va=c(p_va,round(c_$p.value,2))
      }
    }
    p_va=p_va[-26]
    row.names(result_tt)=row.names(result_group)
    result_tt=data.frame(result_tt)
    exp_median=aggregate(.~group_column,data=result_tt,median)
    exp_median=rbind(exp_median[,-1],p_va)
    #row.names(exp_median)=c(group_name,"p value")
    row.names(exp_median)=c(group_name,"p value")
    group_fre<<-exp_median
    write.table(group_fre,save_group,sep="\t",quote=F,col.names = NA)
    T_FRE<<-result_mat
    plot_fun("boxplot",1,save_plot,customer)
  }
  if(response_tag){
    feature=names(paper_marker[-c(23,24)])
    df=data.frame(cbind(train_data[,feature],train_tag))
    rlt1<- svm(as.numeric(as.factor(train_tag)) ~ ., data=df, kernel="radial",cross=5,type="eps-regression")

    rlt<- svm(train_tag ~ ., data=df, kernel="radial",cross=5,type="C-classification")
    pred_data=result_mat[,feature]
    pred_result1=predict(rlt1,pred_data)
    #print(pred_result1)
    pred_result=predict(rlt,pred_data)
    Response=as.vector(unlist( pred_result))
    Score=round(as.vector(unlist( pred_result1)),3)
    colnames(result_mat)=c("CD4_naive","CD8_naive","Cytotoxic","Exhausted","Tr1","nTreg","iTreg","Th1","Th2","Th17",
                           "Tfh","Central_memory","Effector_memory","NKT","MAIT","DC","B_cell","Monocyte","Macrophage","NK","Neutrophil","Gamma_delta","CD4_T","CD8_T","InfiltrationScore")
    ICB_response<<-cbind(Response,Score,result_mat)
    colnames(result_mat)=c("CD4_naive","CD8_naive","Cytotoxic","Exhausted","Tr1","nTreg","iTreg","Th1","Th2","Th17",
                           "Tfh","Central_memory","Effector_memory","NKT","MAIT","DC","B_cell","Monocyte","Macrophage","NK","Neutrophil","Gamma_delta","CD4_T","CD8_T","InfiltrationScore")
    write.table(ICB_response,save_icb,sep="\t",quote=F,col.names = NA)
  }
  result_mat=as.matrix(result_mat)
  #print(ncol(result_mat))
  if(ncol(result_mat)==1){
    colnames(result_mat)=colnames(sample)
    T_FRE<<-t(result_mat)
  }else{
    T_FRE<<-result_mat
  }

  #  shinyWidgets::updateProgressBar(title = 'Done',id = "pb2",value=100,session=getDefaultReactiveDomain())
  TAG<<-1
  # incProgress(20, detail = "Done"
  #print(head(result_mat))
 # print(save_bar)
  plot_fun("barplot",1,save_bar,customer)
  write.table(result_mat,save_abun,sep="\t",quote=F,col.names = NA)
  # return(result_mat)
  # })
}



#' Title
#'
#' @param exp
#' @param gene_name
#' @param data_type
#' @param pre_result
#'
#' @return sample immune cell infilatration score

immune_infiltate_calculate=function(exp,gene_name,data_type,pre_result){
  inf = 0
  names(exp) = gene_name
  for (cell in names(immune_infiltate_marker)){
    abun = 0
    markers = as.vector(unlist(immune_infiltate_marker[cell]))
    for (gene in markers){
      if(data_type == "microarray"){
        abun = abun + as.numeric(exp[gene])/marker_exp_T[gene,cell]
      }else{
        abun = abun + as.numeric(log2(exp[gene]+1))/marker_exp_T[gene,cell]
      }
    }
    inf = inf+abun/length(as.vector(unlist(immune_infiltate_marker[cell])))
  }
  return(inf)
}


#' Title
#'
#' @param comgenes
#' @param tag_gene
#'
#' @return cell signature gene tag matrix
marker_tag = function(comgenes,tag_gene){
  a = comgenes
  a[which(comgenes%in%tag_gene)] = 1
  a[which(a!=1)] = 0
  a = as.numeric(a)
  return(a)
}


#' Title
#'
#' @param data
#' @param marker_exp
#' @param marker_tag_mat
#' @param data_type
#'
#' @return sample expression deviation matrix


sample_ratio = function(data,marker_exp,marker_tag_mat,data_type){
  exp = 0
  if(data_type == "microarray"){
    for (cell in colnames(marker_exp)){
      exp = exp+data/(marker_exp[,cell])*marker_tag_mat[,cell]
    }
  }else{
    for (cell in colnames(marker_exp)){
      exp = exp+log2(data+1)/marker_exp[,cell]*marker_tag_mat[,cell]
    }
  }
  return(exp)
}

#' Title
#'
#' @param raw_score
#' @param compensation_matrix
#'
#' @return corrected NES

compensation = function(raw_score,compensation_matrix){
  raw_score=as.matrix(raw_score)
  compensation_matrix = compensation_matrix
  diag(compensation_matrix) = 1
  rows <- rownames(raw_score)[rownames(raw_score) %in%  rownames(compensation_matrix)]
  if(ncol(raw_score)==1){
    scores <- as.matrix(pracma::lsqlincon(compensation_matrix[rows,rows], raw_score, lb = 0))
  }else{
    scores <- apply(raw_score[rows,], 2, function(x) pracma::lsqlincon(compensation_matrix[rows,rows], x, lb = 0))

  }
  scores[scores < 0] = 0
  rownames(scores) <- rows
  return(scores)
}


#' Title
#'
#' @param figure_type
#' @param group_tag
#' @param save_plot
#' @param customer
#'
#' @return
#' @export
#'
#' @examples

plots<-list()
plot_fun=function(figure_type,group_tag,save_plot,customer){
 # print(ncol(T_FRE))
  if(customer==0){
    data=T_FRE[,seq(1,ncol(T_FRE)-1)]
    colnames(data)<-c("CD4 naive","CD8 naive","Tc","Tex","Tr1","nTreg","iTreg","Th1","Th2","Th17","Tfh","Tcm","Tem","NKT","MAIT","DC","B cell","Monocyte","Macrophage","NK","Neutrophil","Tgd","CD4 T","CD8 T")
    p=NULL
    N=nrow(data)
    count=ncol(data)
    groups=as.factor(as.vector(unlist(group_content)))
  }else{
    data=T_FRE
    count=ncol(data)
    p=NULL
    N=nrow(data)
  }
  if ((group_tag)&&(figure_type=="boxplot")){
    plot_index<-0
    for (cell in (colnames(data))){
      #  Sys.sleep(0.55)
      plot_index=plot_index+1
      cell_plot(data[,plot_index],plot_index,groups,cell)
      #incProgress(0.04, detail = cell)
    }
    #incProgress(0.02, detail = "Multiplot generating")
    # pdf(paste(save_plot,".pdf",sep=""),16,4)
    # ggplot2.multiplot(plotlist = plots,cols = 12)
    # dev.off()
    png(save_plot,bg="transparent",width =1500 ,height=400)
    ggplot2.multiplot(plotlist = plots,cols = 12)
    dev.off()
    #group_box_figure<<-ggarrange(plotlist = plots,ncol=12)
    group_tag<<-1
    # withProgress(message = 'Creating plot', value = 0.1, {

    #incProgress(0.02, detail = "Done")
  }
  #print(group_box_figure)

  if(figure_type=="barplot"){
    #  print("Draw barplot")
    x=rep(row.names(data),each=count)
    y=c()
    p=NULL
    for(i in seq(1,N)){
      y=c(y,as.numeric(data[i,]))
    }
    df=data.frame(x=x,y=y)
    cell_type=rep(colnames(data),times=N)
    #print(head(df))
    #p1=ggplot(df, mapping = aes(x = x, y = y, fill = cell_type)) + xlab("Sample")+ylab("Abundance")+geom_bar(stat = 'identity',position = 'fill',width=0.3)+theme(legend.text=element_text(size=12),axis.text.y=element_text(size=14),axis.text.x=element_text(angle=90,hjust = 0.5,vjust=0.5,size=14))+scale_fill_discrete(name="Cell type")+theme(panel.background = element_rect(fill='#EDEDED', colour='#EDEDED'))+
    # theme(plot.margin=unit(x=c(0,0,0,0),units="mm"),plot.background = element_rect(fill="#EDEDED"))+theme(axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16))

    png(save_plot,bg="transparent",width =1500 ,height=400)
    p<-ggplot(df, mapping = aes(x = x, y = y, fill = cell_type)) + xlab("Sample")+ylab("Abundance")+geom_bar(stat = 'identity',position = 'fill',width=0.3)+theme(legend.text=element_text(size=12),axis.text.y=element_text(size=14),axis.text.x=element_text(angle=90,hjust = 0.5,vjust=0.5,size=14))+scale_fill_discrete(name="Cell type")+theme(panel.background = element_rect(fill='#EDEDED', colour='#EDEDED'))+theme(plot.margin=unit(x=c(0,0,0,0),units="mm"),plot.background = element_rect(fill="#EDEDED"))+theme(axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16))
    print(p)
    dev.off()

    #return(p1)
  }

}

#' Title
#'
#' @param fra
#' @param plot_index
#' @param group
#' @param celltype
#'
#' @return group comparison result

cell_plot=function(fra,plot_index,group,celltype){
  #df$fra=as.numeric(c((pre_fra[,"iTreg"]+pre_fra[,"nTreg"]+pre_fra[,"Tr1"]),(on_fra[,"iTreg"]+on_fra[,"nTreg"]+on_fra[,"Tr1"])))
  df=data.frame(x=as.factor(as.numeric(group)),y=fra)
  if(length(unique(group))>2){
    p=ggplot(df,mapping= aes(x=x, y=y,fill=as.factor(group)))+theme(axis.text.x=element_text(angle=50,hjust = 0.5,vjust=0.5))+
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black",size = 0.2),plot.title = element_text(hjust = 0.5,size = 8))+
      labs(title =celltype,x=NULL,y=NULL)+theme(legend.position='none',axis.text.x =element_text(size=16), axis.text.y=element_text(size=16),plot.title = element_text(size=16))+geom_boxplot(width=0.5)+stat_compare_means(method = "anova",label = "p.signif")

  }else{
    p=ggplot(df,mapping= aes(x=x, y=y,fill=as.factor(group)))+theme(axis.text.x=element_text(angle=50,hjust = 0.5,vjust=0.5))+scale_fill_manual(values=c("#0077c8","#e84118"))+
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black",size = 0.2),plot.title = element_text(hjust = 0.5,size = 8))+
      labs(title =celltype,x=NULL,y=NULL)+theme(legend.position='none',axis.text.x =element_text(size=16), axis.text.y=element_text(size=16),plot.title = element_text(size=16))+geom_boxplot(width=0.5)+stat_compare_means(method = "wilcox",label = "p.format")
  }
  plots[[plot_index]]<<-p
}
#print(args)
#getResult(args$sample,args$dataType,args$groupTag,args$customer,args$signatureGene,args$signatureExpression,args$abundance,args$groupDiff,
#                  args$response,args$abundancePlot,args$groupPlot)
