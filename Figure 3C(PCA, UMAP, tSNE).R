library(optparse)
option_list <- list(
  make_option(c("-i", "--infile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/tool_runing/a7f2e171cab47364d8d8ab6c7fc99169/input.json',
              action = "store", help = "Input a exp file path!"
  ),
  make_option(c("-o", "--outfile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/web_file_catche/runing/dae42f7681eca28355c6ba9105a3b5ff',
              action = "store", help = "Input a outfolder path!"
  )
)
logs=c()

tryCatch({
  Args <- commandArgs()
  opt = parse_args(OptionParser(option_list = option_list, usage = "Data press"))
  #logs=c(logs,paste0('geting data:',paste0(paste0(names(opt),'=',opt),collapse = ',')))
  logs=c(logs,paste0('run pcg.R-',basename(opt$outfile)))
  
  library(jsonlite)
  data<-jsonlite::stream_in(file(opt$infile),pagesize = 1000)
  method=unlist(data$method)
  exp_path=unlist(data$exp_path)
  groups=unlist(data$groups)
  samples=unlist(data$samples)
  outfolder=opt$outfile
  logs=c(logs,paste0('reading data:',basename(exp_path)))
  #exp_path='/pub1/data/mg_projects/projects/web_script/tool_runing/test_data/Merge_RNAseqCount.txt'
  dat=data.table::fread(exp_path, sep = "\t",header = T,stringsAsFactors = F,check.names = F
                        ,na.strings="NA",data.table = F)  
  dat=dat[match(unique(dat[,1]),dat[,1]),]
  row.names(dat)=as.character(dat[,1])
  dat=dat[,-1]
  #which((dat[,1]==''))
  logs=c(logs,paste0('readed data,ncol=',ncol(dat),'nrow=',nrow(dat)))
  exp=apply(dat, 2, as.numeric)
  row.names(exp)=row.names(dat)
  logs=c(logs,'converted numeric')
  dat=(exp[apply(t(exp), 2, function(x){
    return(sd(x,na.rm=T))
  })>0,])
  nnames=row.names(dat)
  if(sum(is.na(dat))>0){
    dat=t(impute::impute.knn(as.matrix(t(dat)))$data)
  }
  colnames(dat)=colnames(exp)
  row.names(dat)=nnames
  
  if(method=='MVT'){
    lmfit=limma::lmFit(dat)
    logs=c(logs,'outping')
    write.table(cbind(lmfit$Amean,sqrt(lmfit$sigma)),file = paste0(outfolder,'/stat.mtx'),quote = F,row.names = F,col.names = F,sep = '\t')
    logs=c(logs,'outputed')
  }else if(method=='boxplot'){
    bx=apply(dat, 2, function(x){
      q1=quantile(x, c(0,0.25, 0.5, 0.75,1), na.rm=T)
      otl=1.5*(q1[4]-q1[2])
      otl.u=ifelse(q1[4]+otl>q1[5],q1[5],q1[4]+otl)
      otl.d=ifelse(q1[2]-otl<q1[1],q1[1],q1[2]-otl)
      return(c(otl.d,q1[2:4],otl.u,sd(x,na.rm = T),mean(x,na.rm = T)))
    })
    logs=c(logs,'outping')
    write.table(cbind(Tag=colnames(bx),t(bx)),file = paste0(outfolder,'/stat.mtx'),quote = F,row.names = F,col.names = F,sep = '\t')
    logs=c(logs,'outputed')
  }else if(method=='PCAplot'){
    nn=ceiling(ncol(dat)/3)
    if(nn>15) nn=15
    if(nn<2) nn=2
    library(umap)
    p_input=t(dat)
    p_input=as.data.frame(scale(p_input))
    #colnames(p_input)=row.names(dat)
    #colnames(p_input)
    #class(p_input)
    #head(p_input[,1:10])
    logs=c(logs,paste0('runing umap,N=',nn))
    iris.umap = umap(scale(p_input), n_neighbors = nn, random_state = 123)
    umaply=iris.umap$layout
    logs=c(logs,'outping umap')
    write.table(cbind(Tag=row.names(umaply),umaply),file = paste0(outfolder,'/stat_UMAP.mtx'),quote = F,row.names = F,col.names = F,sep = '\t')
    logs=c(logs,'outputed umap')
    
    logs=c(logs,paste0('runing PCA',''))
    pca1=stats::prcomp(p_input)
    logs=c(logs,'pca runed')
    as.matrix(p_input)
    head(p_input)
    unscale=function (data, center = NULL, scale = NULL) 
    {
      if (is.null(scale)) {
        scale <- attr(data, "scaled:scale")
      }
      if (is.null(center)) {
        center <- attr(data, "scaled:center")
      }
      if (!is.null(scale) && !is.logical(scale)) {
        data <- base::scale(data, center = FALSE, scale = 1/scale)
      }
      if (!is.null(center) && !is.logical(center)) {
        data <- base::scale(data, center = -center, scale = FALSE)
      }
      as.data.frame(data)
    }
    
    fortify=function (model) 
    {
      if (is(model, "prcomp")) {
        d <- as.data.frame(model$x)
        values <- model$x %*% t(model$rotation)
      }
      else if (is(model, "princomp")) {
        d <- as.data.frame(model$scores)
        values <- model$scores %*% t(model$loadings[, ])
      }
      else {
        stop(paste0("Unsupported class for fortify.pca_common: ", 
                    class(model)))
      }
      values <- unscale(values, center = model$center, 
                                   scale = model$scale)
      #values <- cbind_wraps(data, values)
      #cbind(data,values[match(row.names(data),row.names(values)),])
      return(d);
      #d <- cbind_wraps(values, d)
      #post_fortify(d)
    }
    
    plot.data <- fortify(pca1)
    #dim(plot.data)
    #dim(plot.data2)
    
    logs=c(logs,'fortify PCA')
    plot.data$rownames <- rownames(plot.data)
    ve <- pca1$sdev^2/sum(pca1$sdev^2)
    PC <- paste0("PC", 1:length(pca1$sdev))
    lam <- pca1$sdev
    lam <- lam * sqrt(nrow(plot.data))
    plot.data.pc <- t(t(plot.data[,PC])/lam)
    otd=rbind(summary(pca1)[[6]][2,],plot.data.pc)
    logs=c(logs,'outping PCA')
    write.table(cbind(Tag=row.names(otd),otd)
                ,file = paste0(outfolder,'/stat_PCA.mtx')
                ,quote = F,row.names = F,col.names = F,sep = '\t')
    logs=c(logs,'outputed PCA')
    logs=c(logs,paste0('runing tSNE',''))
    floor((nrow(p_input)-1 )/3)->tpk
    if(tpk>0){
      if(tpk>50) tpk=50
      tsne <- Rtsne::Rtsne(as.matrix(p_input), check_duplicates = FALSE, pca = T,
                         perplexity=tpk, theta=0.0, dims=2)
      logs=c(logs,'outping tSNE')
      write.table(cbind(Tag=row.names(p_input), tsne$Y)
                  ,file = paste0(outfolder,'/stat_tSNE.mtx')
                  ,quote = F,row.names = F,col.names = F,sep = '\t')
      logs=c(logs,'outputed tSNE')
    }
    
  }
  
},error = function(e) {
  print(conditionMessage(e))
  logs=c(logs,paste0('error:',conditionMessage(e)))
}, finally = {
  write.table(logs,file = paste0(opt$outfile,'/run.log'),quote = F,row.names = T,col.names = T,sep = '\t')
})
  