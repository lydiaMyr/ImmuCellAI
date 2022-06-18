

#' Title
#'
#' @param sample
#' @param data_type
#' @param group_tag
#' @param response_tag
#' @param customer
#' @param sig_file
#' @param exp_file
#'
#' @return
#' @export
#'
#' @examples
#'

ImmuCellAI = function(sample,data_type,group_tag,response_tag,customer,sig_file=NULL,exp_file=NULL){
  data("marker_exp")
  data("marker_exp_T")
  data("paper_marker")
  data("train_data")
  data("train_tag")
  data("compensation_matrix")
  data("immune_infiltate_marker")
  if (customer==TRUE){
    paper_marker<-sig_file
    marker_exp=exp_file
  }
  group_index=0
  if (group_tag){
    #group_index=as.numeric(as.vector(unlist(grep("group",row.names(sample)))))
    group_column<-sample[1,]
    group_content<-sample[1,]
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
  print("work_done")
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
  #print("test")
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
   # write.table(group_fre,save_group,sep="\t",quote=F,col.names = NA)
    T_FRE<<-result_mat
   # plot_fun("boxplot",1,save_plot,customer)
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
    ICB_response<<-cbind(Response,Score)
    colnames(result_mat)=c("CD4_naive","CD8_naive","Cytotoxic","Exhausted","Tr1","nTreg","iTreg","Th1","Th2","Th17",
                           "Tfh","Central_memory","Effector_memory","NKT","MAIT","DC","B_cell","Monocyte","Macrophage","NK","Neutrophil","Gamma_delta","CD4_T","CD8_T","InfiltrationScore")
    #write.table(ICB_response,save_icb,sep="\t",quote=F,col.names = NA)
  }
  result_mat=as.matrix(result_mat)
  #print(ncol(result_mat))
  if(ncol(result_mat)==1){
    colnames(result_mat)=colnames(sample)
    T_FRE<<-t(result_mat)
  }else{
    T_FRE<<-result_mat
  }
  return(list(Sample_abundance=T_FRE, Group_result=group_fre, Response=ICB_response))

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


