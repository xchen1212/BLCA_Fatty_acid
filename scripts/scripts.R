dir.create('BLCA_Fatty_acid')
setwd('BLCA_Fatty_acid')
if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("PDFs")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
}
library(reshape2)
library(ggpubr)
library(ggsci)
library(maftools)
library(tidyr)
library(pheatmap)
library(clusterProfiler)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(corrplot)
library(survminer)
library(survival)
options(stringsAsFactors = F)
source('base.R')
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}
my_boxplot=function(dat,group,group_cols,test_method='kruskal.test',fill= "Group",label=c("p.format",'p.signif')[1],
                    xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat=ARGs.score[tcga.subtype$Samples,'score']
  # group=tcga.subtype$Cluste
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  p=ggplot(data, aes(x=Group, y=value,fill=Group)) +
    geom_boxplot()+
    #scale_color_manual(values = sub.col) +  #箱线图颜色
    scale_fill_manual(values = group_cols)+   #箱线图填充颜色
    ggpubr::stat_compare_means(aes(group=Group), label = label, method = test_method)+
    labs(x=xlab, y = ylab, fill = fill) +
    theme_bw()+
    theme(legend.position =legend.position,                 #图例位置
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = x.size),
          axis.text.y = element_text(size = y.size)) # 将图表标题居中
  return(p)
}


my_mutiboxplot=function(dat,group,group_cols,test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],bw=T,
                        xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,legend.position='top',fill='group'){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  colnames(dat.melt)=c('Group','type','value')
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot()+
    #scale_color_manual(values = sub.col) +  #箱线图颜色
    scale_fill_manual(values =group_cols)+   #箱线图填充颜色
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    theme_bw()+
    #theme_classic()
    theme(legend.position = legend.position,                 #图例位置
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust)) # 将图表标题居中
  return(p)
}

my_mutiboxplot_seg=function(dat,group,group_cols,test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],fill='group',
                            xlab='',ylab='score',title='',xsize=10,xangle = 45,xhjust = 1,ysize=10,
                            legend.position='top',nrow,ncol){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  colnames(dat.melt)=c('Group','type','value')
  p=dat.melt %>%
    ggplot(aes(x=Group, y=value,fill=Group)) +
    geom_boxplot()+facet_wrap(~type,scales = 'free',nrow = nrow,ncol = ncol)+
    #scale_color_manual(values = sub.col) +  #箱线图颜色
    scale_fill_manual(values =group_cols)+   #箱线图填充颜色
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill = fill,title =title) +
    #theme_light()+
    #theme_classic()+
    theme(legend.position = legend.position,                 #图例位置
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = xsize,angle = xangle, hjust = xhjust),
          axis.text.y = element_text(size= ysize)) # 将图表标题居中
  return(p)
}


mg_violin_1=function(data,xangle=0,ylab='value',xlab='',leg.title='Group',test_method='anova',cmp_test_method='t.test',
                     legend.pos='r',melt=F,jitter=T,ylim=NULL,show_compare=NULL,point_size=NULL,group_col){
  library(ggplot2)
  if(is.null(ylim)){
    
  }
  if(melt){
    data_m=data
    colnames(data_m)=c('Group','value')
  }else{
    data_m=reshape2::melt(data)
    colnames(data_m)=c('Group','value')
  }
  if(!is.null(ylim)){
    data_m$value[data_m$value>ylim[2]]<-NA
  }
  data_m=data_m[which(!is.na(data_m[,1])),]
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos='right'
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='t'){
    pos='top'
  }else if(legend.pos=='r'){
    pos='right'
  }else if(legend.pos=='b'){
    pos='bottom'
  }
  uni.group=unique(data_m[,1])
  ct=length(uni.group)
  
  p1<-ggplot(data_m,aes(x=Group,y=value))+geom_violin(alpha=0.7)
  # if(ct<=4){
  #   p1=p1+ggsci::scale_fill_lancet()
  # }else if(ct<=10){
  #   p1=p1+ggsci::scale_fill_npg(name=leg.title)
  # }else if(ct<=20){
  #   p1=p1+ggsci::scale_fill_d3(palette = "category20",name=leg.title)
  # }else if(ct<=30){
  #   cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10),ggsci::pal_d3("category20", alpha = 0.6)(20))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }else if(ct<=38){
  #   cbPalette=c(ggsci::pal_lancet()(10)
  #               ,ggsci::pal_npg("nrc", alpha = 0.6)(10)
  #               ,ggsci::pal_d3("category20", alpha = 0.6)(20)
  #               ,ggsci::pal_nejm("default", alpha = 0.6)(8))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }
  p1=p1+scale_fill_manual(values=group_col)
  if(jitter){
    if(is.null(point_size)){
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2)
    }else{
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2,size=point_size)
    }
  }
  
  p1=p1+theme_bw()+geom_boxplot(width=0.2,aes(fill=Group),outlier.shape = NA)
  p1=p1+theme(axis.text.x=tx, 
              axis.text.y=element_text(family="Times",face="plain"), 
              axis.title.y=element_text(family="Times",face="plain"), 
              legend.text=element_text(face="plain", family="Times", colour="black" 
              ),
              legend.title=element_text(face="plain", family="Times", colour="black"
              ),
              legend.justification=pos, legend.position=pos
              ,legend.background = element_rect(fill = NA, colour = NA)
  )+ylab(ylab)+xlab(xlab)
  til=''
  if(test_method=='anova'){
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=t.test(x1,x2)$p.value 
      til=paste0('t-tests p=',signif(pv,2))
    }else{
      fit <- aov(value~Group, data = data_m)
      pv=summary(fit)[[1]][5][[1]]
      fv=summary(fit)[[1]][4][[1]]
      til=paste0('ANOVA tests p=',signif(pv,2))
    }
  }else{
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=wilcox.test(x1,x2)$p.value 
      til=paste0('wilcox.tests p=',signif(pv,2))
    }else{
      fit=kruskal.test(value~Group, data = data_m)
      pv=fit$p.value
      til=paste0('Kruskal-Wallis test p=',signif(pv,2))
    }
  }
  p1=p1+ggtitle(til) 
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(is.null(show_compare)){
    if(length(uni.group)>5){
      show_compare=F
    }else{
      show_compare=T
    }
  }
  if(show_compare){
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,method = cmp_test_method,label= "p.signif")
  }
  return(p1)
}
coxFun <- function(dat){
  library(survival)
  colnames(dat)=c('time','status','gene')
  fmla=as.formula("Surv(time,status)~gene")
  cox=coxph(fmla,data=dat)
  p=summary(cox)[[7]][5]
  result=c(p,summary(cox)[[8]][1],summary(cox)[[8]][3],summary(cox)[[8]][4])
  return(result)
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  set.seed(2021)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Parial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
get_riskscore<-function(dat,os,os.time,step=T,direction=c("both", "backward", "forward")[1]){
  tcga_dat1 <- cbind(time=os.time,
                     status=os,
                     dat)
  tcga_dat1=crbind2DataFrame(tcga_dat1)
  colnames(tcga_dat1)=gsub('-','__',colnames(tcga_dat1))
  gene111=gsub('-','__',colnames(dat))
  fmla <- as.formula(paste0("Surv(time, status) ~"
                            ,paste0(gene111,collapse = '+')))
  cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
  if(step==T){
    cox1 <- step(cox,direction =direction)
  }else{
    cox1=cox
  }
  lan <- coef(cox1)
  #round(lan, 3)
  genes <- names(cox1$coefficients)
  mult_results=paste0(round(lan, 3), '*', names(lan),collapse = '+')
  risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))
  
  data_gene_score_final<-tcga_dat1
  data_gene_score_final$Samples<-rownames(data_gene_score_final)
  data_gene_score_final$riskscore=risk.tcga
  data_gene_score_final$riskscorez=mosaic::zscore(risk.tcga)
  #最佳截断
  optimalCutoff <- survminer::surv_cutpoint(data.frame(time=data_gene_score_final$time/365,
                                                       event=data_gene_score_final$status,
                                                       risk=data_gene_score_final$riskscore), 
                                            time = "time", event = "event",variables = c("risk"))
  optimalCutoff=optimalCutoff$cutpoint$cutpoint[1]
  #print(optimalCutoff)
  #optimalCutoff=median(data_gene_score_final$riskscore)
  #optimalCutoff=0
  data_gene_score_final$Risk=ifelse(data_gene_score_final$riskscore>optimalCutoff,'High','Low')
  table(data_gene_score_final$Risk)
  data_gene_score_final$cutoff=optimalCutoff
  return(list(result=data_gene_score_final,module.gene=genes,model=mult_results))
}
bioForest=function(rt=null,col){
  #读取输入文件
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #输出图形
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #绘制右边森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  ## 未过滤任何基因
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 降序排列
  return(geneList)
}
mg_nomogram=function(clinical_riskscore,
                     os,
                     status,
                     title='Nomogram',
                     quick=T,
                     mks = c(1,3,5)){
  #clinical_riskscore=dat1[,3:5]
  #os=dat1[,1]
  #status=dat1[,2]
  #sum(is.na(norm.stat.al[,3]))
  norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
  norm.stat.al=as.data.frame(norm.stat.al)
  library(rms)
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
  options(datadist='MG_Grobal_DDSet')
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
  cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
  #summary(cox2)
  #surv=Survival(cox2)
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  #cindex.orig=1-rcorr.cens(fp,Surv(norm.stat.al$time,norm.stat.al$status))
  cut.time=c()
  if(quantile(os[!is.na(os)])['75%']<12){
    cut.time=mks
  }else if(quantile(os[!is.na(os)])['75%']<365){
    cut.time=c(12*mks[1],12*mks[2],12*mks[3])
  }else{
    cut.time=c(365*mks[1],365*mks[2],365*mks[3])
  }
  cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['100%'])]
  print(cut.time)
  #regplot(cox2)
  #  print(regplot(cox3#对观测2的六个指标在列线图上进行计分展示
  #,observation=pbc[2,] #也可以不展示
  #预测3年和5年的死亡风险，此处单位是day
  #              ,title=title
  #              ,failtime = cut.time
  #              ,prfail = TRUE #cox回归中需要TRUE
  #              ,showP = T #是否展示统计学差异
  #              ,droplines = F#观测2示例计分是否画线
  #,colors = mg_colors[1:3] #用前面自己定义的颜色
  #,rank="decreasing") #根据统计学差异的显著性进行变量的排序
  #,interval="confidence"
  #,rank="decreasing"
  #,clickable=T
  #              ,points=TRUE)) #展示观测的可信区间
  
  #  plot(nom)
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
    #    surv0 <- function(x)surv(cut.time[i],lp=x) 
    #    survs=c(survs,list(surv0))
  }
  #f2<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[2]) 
  #cal3<-calibrate(f2, cmethod="KM", method="boot",u=cut.time[2],m=100,B=200) 
  #f3<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[3]) 
  #cal5<-calibrate(f3, cmethod="KM", method="boot",u=cut.time[3],m=100,B=200) 
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c(paste0(mks[1], '-Year Survival'),
                           paste0(mks[2], '-Year Survival'),
                           paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
    #cal3=get_best_calibrate(cox2,cut.time[2])
    #cal5=get_best_calibrate(cox2,cut.time[3])
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  #customLayout::lay_show(cl)
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  #lines(cal_all[[1]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[1], pch = 16)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
      #lines(cal_all[[i]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[i], pch = 16)
    }
    #plot(cal3,lwd = 2,lty = 0,errbar.col = mg_colors[3],xlim = c(0,1),ylim= c(0,1),col = mg_colors[3],add = T)
    #lines(cal3[,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[3], pch = 16)
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c(paste0(mks[1], '-Year'),
                               paste0(mks[2], '-Year'),
                               paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=crbind2DataFrame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  #plot(cal1,xlim=c(0,1),ylim=c(0,1))
  #plot(cal3,xlim=c(0,1),ylim=c(0,1))
  #plot(cal5,xlim=c(0,1),ylim=c(0,1))
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}
mg_RunGSEA_use=function(mod=c('exp_group','exp_gene','rank')[1],exp_Path=NULL, sample_group_path=NULL, outFolder=NULL,gene=NULL, column=NULL,lower=50,upper=50, gmt_Path=c("KEGG",'GO_BP','GO_CC','GO_MF','reactome','HALLMARK','TF')[1],plot_svg=FALSE,top=10,min=5,max=5000,outLog=T){
  
  if(is.null(exp_Path)|is.null(mod)|is.null(outFolder)|is.null(gmt_Path)){
    return(NULL)
  }
  if(plot_svg){
    svg='true'
  }else{
    svg='false'
  }
  if(gmt_Path=='KEGG'){
    #gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c2.cp.kegg.v7.0.symbols.gmt')
    gmt_Path="/pub1/data/mg_projects/users/yanj/project/public/gmt/c2.cp.kegg.v7.5.1.symbols.gmt"
  }
  else if(gmt_Path=='GO_BP'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c5.bp.v7.0.symbols.gmt')
  }
  else if(gmt_Path=='GO_CC'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c5.cc.v7.0.symbols.gmt')
  }
  else if(gmt_Path=='GO_MF'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c5.mf.v7.0.symbols.gmt')
  }
  else if(gmt_Path=='reactome'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c2.cp.reactome.v7.0.symbols.gmt')
  }
  else if(gmt_Path=='HALLMARK'){
    #gmt_Path=paste0(MG_Grobal_baseFolder,'/source/h.all.v7.0.symbols.gmt')
    gmt_Path="/pub1/data/mg_projects/users/yanj/project/public/gmt/h.all.v7.5.1.symbols.gmt"
  }
  else if(gmt_Path=='TF'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c3.tft.v7.0.symbols.gmt')
  }
  if(file.exists(paste0(getwd(),'/',outFolder))){
    outFolder=paste0(getwd(),'/',outFolder)
  }else if(!file.exists(outFolder)){
    dir.create(outFolder)
    if(file.exists(paste0(getwd(),'/',outFolder))){
      outFolder=paste0(getwd(),'/',outFolder)
    }
  }
  
  if(file.exists(paste0(getwd(),'/',gmt_Path))){
    gmt_Path=paste0(getwd(),'/',gmt_Path)
  }
  if(file.exists(paste0(getwd(),'/',exp_Path))){
    exp_Path=paste0(getwd(),'/',exp_Path)
  }
  
  command=NULL
  if(mod=='exp_group'){
    if(!is.null(exp_Path)&!is.null(sample_group_path)&!is.null(outFolder)){
      if(file.exists(paste0(getwd(),'/',sample_group_path))){
        sample_group_path=paste0(getwd(),'/',sample_group_path)
      }
      command=paste0(MG_Grobal_baseFolder,'/jre/GSEA/MG_GSEA.jar exp_group '
                     ,exp_Path,' ',sample_group_path,' ',outFolder,' ',gmt_Path,' ',svg,' ',top,' ',min,' ',max) 
    }
  }else if(mod=='exp_gene'){
    if(!is.null(exp_Path)&!is.null(gene)&!is.null(outFolder)){
      command=paste0(MG_Grobal_baseFolder,'/jre/GSEA/MG_GSEA.jar exp_gene '
                     ,exp_Path,' ',gene,' ',outFolder,' ',gmt_Path,' ',svg,' ',top,' ',min,' ',max,' ',lower,' ',upper) 
    }
  }else if(mod=='rank'){
    if(!is.null(exp_Path)&!is.null(column)&!is.null(outFolder)){
      command=paste0(MG_Grobal_baseFolder,'/jre/GSEA/MG_GSEA.jar rank '
                     ,exp_Path,' ',column-1,' ',outFolder,' ',gmt_Path,' ',svg,' ',top,' ',min,' ',max) 
    }
  }
  
  if(!is.null(command)){
    if(MG_Grobal_System=='win'){
      command=paste0(MG_Grobal_baseFolder,'/jre/bin/java -jar ',command)
    }else{
      command=paste0('java -jar ',command)
    }
    print(paste0('RunGSEA CMD:',command))
    logs=system(command, intern = !outLog, 
                ignore.stdout = FALSE, ignore.stderr = FALSE, 
                wait = TRUE, input = NULL, show.output.on.console = TRUE, 
                minimized = FALSE, invisible = TRUE)
    if(outLog){
      if(logs==0){
        print('Run GSEA succ')
      }else{
        print('Run GSEA error')
      }
    }else{
      print('Runed GSEA')
      print(logs)
      logs=logs[grep('######/',logs)]
      if(length(logs)==1){
        logs=unlist(strsplit(logs[1],'/'))
        if(length(logs)>1){
          return(logs[2:length(logs)])
        }
      }
    }
  }
  return(NULL)
}
#####data.pre#########
#TCGA####
#临床数据
tcga_cli<-read.delim('origin_datas/TCGA/Merge_BLCA_clinical.txt',sep='\t',header = T)
colnames(tcga_cli)[1:30]
table(tcga_cli$A22_histological_type)
table(tcga_cli$A23_Tumor_tissue_site)

tcga_cli=data.frame(Samples=tcga_cli$A0_Samples,
                    Age=tcga_cli$A17_Age,
                    Gender=tcga_cli$A18_Sex,
                    T.stage=tcga_cli$A3_T,
                    N.stage=tcga_cli$A4_N,
                    M.stage=tcga_cli$A5_M,
                    Stage=tcga_cli$A6_Stage,Grade=tcga_cli$A7_Grade,
                    OS.time=tcga_cli$A1_OS,Status=tcga_cli$A2_Event)
head(tcga_cli)
table(tcga_cli$T.stage)
tcga_cli$T.stage=gsub('[ab]','',tcga_cli$T.stage)
tcga_cli$T.stage[tcga_cli$T.stage==''|tcga_cli$T.stage=='TX']<-NA

table(tcga_cli$N.stage)
tcga_cli$N.stage[tcga_cli$N.stage==''|tcga_cli$N.stage=='NX']<-NA

table(tcga_cli$M.stage)
tcga_cli$M.stage[tcga_cli$M.stage==''|tcga_cli$M.stage=='MX']<-NA

table(tcga_cli$Stage)
#tcga_cli$Stage=gsub('[ABC]','',tcga_cli$Stage)
tcga_cli$Stage[tcga_cli$Stage=='']<-NA
tcga_cli$Stage=gsub('Stage ','',tcga_cli$Stage)

table(tcga_cli$Grade)
tcga_cli$Grade[tcga_cli$Grade=='Unknown']<-NA

fivenum(as.numeric(tcga_cli$Age))
tcga_cli$Age1=ifelse(tcga_cli$Age>69,'>69','<=69')

tcga_cli$OS.time
tcga_cli=tcga_cli[which(tcga_cli$OS.time>0),]
table(tcga_cli$Status)
tcga_cli$OS[tcga_cli$Status=='Alive']=0
tcga_cli$OS[tcga_cli$Status=='Dead']=1
tcga_cli$Samples=paste0(tcga_cli$Samples,'-01')
rownames(tcga_cli)=tcga_cli$Samples
head(tcga_cli)
dim(tcga_cli)

#筛选编码蛋白基因
genecode=read.delim('GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
head(mrna_genecode)


#表达谱
tcga_data<-read.delim('origin_datas/TCGA/BLCA_TPM.txt',sep='\t',header = T,row.names = 1,check.names = F)
tcga_data[1:4,1:4]
table(substr(colnames(tcga_data),14,15))

sample_T=colnames(tcga_data)[which(as.numeric(substr(colnames(tcga_data),14,15))==1)]#肿瘤样本
sample_N=colnames(tcga_data)[which(as.numeric(substr(colnames(tcga_data),14,15))==11)]#正常样本
tcga_type=data.frame(Samples=c(sample_T,sample_N),Type=rep(c('Tumor','Normal'),c(length(sample_T),length(sample_N))))
rownames(tcga_type)=tcga_type$Samples
table(tcga_type$Type)

range(tcga_data)
tcga_tpm=tcga_data[intersect(mrna_genecode$SYMBOL,rownames(tcga_data)),c(sample_N,sample_T)]
tcga_tpm_log=log2(tcga_tpm+1)
dim(tcga_tpm_log)

tcga_tpm_log_T=tcga_tpm_log[,intersect(tcga_cli$Samples,sample_T)]
dim(tcga_tpm_log_T)
tcga_cli=tcga_cli[intersect(tcga_cli$Samples,colnames(tcga_tpm_log_T)),]
dim(tcga_cli)


###########GSE32894##################
load('origin_datas/GEO/GSE32894.RData')
GSE32894_cli=GSE32894$Sample
table(GSE32894_cli$Platform)
GSE32894_cli=data.frame(Samples=substr(GSE32894_cli$Title,1,9),
                        Gender=GSE32894_cli$gender,
                        Age=GSE32894_cli$age,
                        Stage=GSE32894_cli$tumor_stage,
                        Grade=GSE32894_cli$tumor_grade,
                        Status=GSE32894_cli$`dod_event_(yes/no)`,
                        OS.time=GSE32894_cli$`time_to_dod_(months)`)
rownames(GSE32894_cli)=GSE32894_cli$Samples

table(GSE32894_cli$Stage)
GSE32894_cli$Stage[GSE32894_cli$Stage=='Tx'|GSE32894_cli$Stage=='Ta']<-NA
GSE32894_cli$Stage=gsub('[ab]','',GSE32894_cli$Stage)

table(GSE32894_cli$Grade)
GSE32894_cli$Grade[GSE32894_cli$Grade=='Gx']<-NA

table(GSE32894_cli$Status)
GSE32894_cli=GSE32894_cli[which(GSE32894_cli$Status!='NULL'),]
GSE32894_cli$OS=ifelse(GSE32894_cli$Status=='no',0,1)
GSE32894_cli=crbind2DataFrame(GSE32894_cli)
GSE32894_cli$OS.time=GSE32894_cli$OS.time*30
head(GSE32894_cli)


load('origin_datas/GEO/GSE32894_exp.RData')
range(GSE32894_exp)
GSE32894_exp[1:5,1:5]
dim(GSE32894_exp)
#20201   216
GSE32894_cli=GSE32894_cli[intersect(colnames(GSE32894_exp_lumi),GSE32894_cli$Samples),]
dim(GSE32894_cli)


######HALLMARK GSVA###########
h.gmt <- getGmt("h.all.v7.5.1.symbols.gmt",
                collectionType=BroadCollection(category="h"),
                geneIdType=SymbolIdentifier())
tcga.hall.ssGSEA <- gsva(as.matrix(tcga_tpm_log_T),
                          h.gmt,
                          method='ssgsea',
                          min.sz=10,
                          max.sz=500,
                          verbose=TRUE)
rownames(tcga.hall.ssGSEA)=gsub('HALLMARK_','',rownames(tcga.hall.ssGSEA))
tcga.hall.ssGSEA[1:5,1:5]


####脂肪酸代谢相关基因############
fatty.acid.gene=read.gmt('origin_datas/HALLMARK_FATTY_ACID_METABOLISM.v2022.1.Hs.gmt')
fatty.acid.gene=fatty.acid.gene$gene

##########单因素
fatty.acid.gene.cox=cox_batch(dat = tcga_tpm_log_T[intersect(fatty.acid.gene,rownames(tcga_tpm_log_T)),
                                                   tcga_cli$Samples],
                              time = tcga_cli$OS.time,
                              event = tcga_cli$OS)
table(fatty.acid.gene.cox$p.value<0.05)
fatty.acid.gene.cox.fit=fatty.acid.gene.cox[fatty.acid.gene.cox$p.value<0.05,]
table(fatty.acid.gene.cox.fit$HR>1)
length(rownames(fatty.acid.gene.cox.fit))

fatty.acid.fit.gene=rownames(fatty.acid.gene.cox.fit)
length(fatty.acid.fit.gene)
#30


#######01.基因的基因组景观（CNV,MAF）############
dir.create('results/01.Genome Landscape')
###########MAF
tcga.maf=getTCGAMAFByCode('BLCA')
pdf('results/01.Genome Landscape//Fig1a.pdf',height = 7,width = 7)
oncoplot(maf = tcga.maf,genes =fatty.acid.fit.gene)
dev.off()


############CMV
library("ComplexHeatmap")
library("circlize")
library(org.Hs.eg.db)
library(clusterProfiler)
get_CNV_Preprocess=function(df_cnv){
  df_cnv$`Gene Symbol`=gsub("\\..*","",df_cnv$`Gene Symbol`)
  rownames(df_cnv)=df_cnv$`Gene Symbol`
  
  library(TCGAutils)
  aliquot_id_to_submitter_id=UUIDtoBarcode(colnames(df_cnv)[-c(1:3)]
                                           ,from_type = 'aliquot_ids')
  colnames(df_cnv)[-c(1:3)]=aliquot_id_to_submitter_id[,2]
  colnames(df_cnv)=substr(colnames(df_cnv),1,15)
  df_cnv=df_cnv[,-which(duplicated(colnames(df_cnv)))]
  
  df_cnv=dplyr::distinct(df_cnv,`Gene Symbol`,.keep_all=TRUE)
  
  ensg=df_cnv$`Gene Symbol`
  
  idmap=bitr(ensg,fromType="ENSEMBL",toType="SYMBOL",OrgDb="org.Hs.eg.db")
  idmap=dplyr::distinct(idmap,ENSEMBL,.keep_all=TRUE)## 一个基因匹配到多个geneid,随机取一个
  
  cnv.inds=which(!is.na(idmap$SYMBOL)) ## 去掉没有注释到gene symbol的行
  idmap=idmap[cnv.inds,]
  df_cnv=df_cnv[idmap$ENSEMBL,]
  df_cnv$`Gene Symbol`=idmap$SYMBOL
  
  df_cnv=dplyr::distinct(df_cnv,`Gene Symbol`,.keep_all=TRUE)
  dim(df_cnv)
  
  # rownames(df_cnv)=df_cnv$`Gene Symbol`
  df_cnv=df_cnv[,-c(2:3)]
  return(df_cnv)
}
get_CNV_Freq=function(df_cnv,genes_custom,genes_type=NULL){
  df_cnv=reshape2::melt(df_cnv,id.vars='Gene Symbol',measure.vars=colnames(cnv.all)[-1]
                        ,variable.name='Sample',value.name='CopyNum')
  head(df_cnv)
  cnv_frq=as.data.frame(matrix(data=0,nrow=1,ncol=5))
  colnames(cnv_frq)<-c('type','gene','CNV_gain','CNV_loss','none_CNV')
  r<-1
  for(i in 1:length(genes_custom))
  {
    cnv_gene<-df_cnv[which(df_cnv$`Gene Symbol` == genes_custom[i]),]
    total_counts<-dim(cnv_gene)[1]
    if(is.null(genes_type)){
      cnv_frq[r,1]<-NA
    }else{
      cnv_frq[r,1]<-genes_type[i]
    }
    cnv_frq[r,2]<-genes_custom[i]
    
    cnv_frq[r,3]<-(dim(cnv_gene[which(cnv_gene$CopyNum >0),])[1])/total_counts ## 在整个拷贝数数据集中
    cnv_frq[r,4]<-(dim(cnv_gene[which(cnv_gene$CopyNum <0),])[1])/total_counts ## CNV_loss
    cnv_frq[r,5]<-(total_counts/total_counts)-cnv_frq[r,3]-cnv_frq[r,4] ## none_CNV
    r<-r+1
  }
  cnv_frq<-cnv_frq[order(cnv_frq$CNV_loss),]
}
###########不同拷贝数的比较
get_CNV_Custom=function(df_cnv,df_exp,genes_custom,abbr){
  df_cnv=reshape2::melt(df_cnv,id.vars='Gene Symbol',measure.vars=colnames(df_cnv)[-1]
                        ,variable.name='Sample',value.name='CopyNum')
  head(df_cnv)
  
  df_cnv$Sample=as.character(df_cnv$Sample)
  df_cnv=df_cnv[df_cnv$`Gene Symbol` %in% genes_custom,]
  dim(df_cnv)
  
  exp.samples=colnames(df_exp)
  print(table(substr(exp.samples,14,15)))
  
  sampl.num=length(unique(df_cnv$Sample))
  
  all.mt=rbind()
  dat.al=rbind()
  
  for(g in genes_custom){
    # g=genes[1]
    inds=which(df_cnv$`Gene Symbol`==g)
    sub_matrix<-df_cnv[inds,]
    dim(sub_matrix)
    
    all.mt=rbind(all.mt,
                 c(sum(sub_matrix$CopyNum>0)/sampl.num,
                   sum(sub_matrix$CopyNum<0)/sampl.num))
    ##### 分别统计每个样本的突变类型和基因的表达情况
    amplification.samp=sub_matrix[which(sub_matrix$CopyNum>0),2]
    deletion.samp=sub_matrix[which(sub_matrix$CopyNum<0),2]
    diploid.samp=setdiff(unique(df_cnv$Sample),c(amplification.samp,deletion.samp))
    
    length(c(amplification.samp,deletion.samp,diploid.samp))
    
    ######## CNV+RNA-seq均有数据的样本名称
    amplification.samp=intersect(amplification.samp,exp.samples)
    deletion.samp=intersect(deletion.samp,exp.samples)
    diploid.samp=intersect(diploid.samp,exp.samples)
    
    ###
    ind1=match(deletion.samp,exp.samples)## del
    ind2=match(diploid.samp,exp.samples) ## amp
    ind3=match(amplification.samp,exp.samples) ## diploid
    ind4=grep('-11$',exp.samples) ## 原发组织样本
    
    sp=c(rep('deletion',length(ind1)),rep('diploid',length(ind2))
         ,rep('amplification',length(ind3)),rep('Normal',length(ind4)))
    # print(length(sp))
    samples=c(deletion.samp,diploid.samp,amplification.samp,exp.samples[grep('-11$',exp.samples)])
    ####
    dat=data.frame(Expression=as.numeric(df_exp[which(row.names(df_exp)==g), c(ind1,ind2,ind3,ind4)])
                   ,category=sp,type=rep(g,length(sp)),sample=samples)
    dim(dat)
    dat.al=rbind(dat.al,dat)  ## 一行代表一个样本
  }
  dat.al=cbind(dat.al,fa=rep(abbr,nrow(dat.al)))
  rownames(all.mt)=genes_custom
  colnames(all.mt)=c('Amp','Del')
  return(dat.al)
}

##########cnv
cnv.all=read.delim('origin_datas/TCGA/Merge_GeneLevelCopyNumber.txt',sep = '\t',stringsAsFactors = F,header = T,check.names = F)
cnv.all=get_CNV_Preprocess(cnv.all)
dim(cnv.all)#19013   412
cnv.all[1:5,1:5]
###########肿瘤样本中的基因cnv频率
tcga.cnv.sp.selected=colnames(cnv.all)[-1][(which(substr(colnames(cnv.all)[-1],13,15)=="-01"))]
length(tcga.cnv.sp.selected)#410
cnv.all=cnv.all[,c("Gene Symbol",tcga.cnv.sp.selected)]

cnv_freq=get_CNV_Freq(df_cnv=cnv.all,genes_custom = fatty.acid.fit.gene)

cnv_freq_1<-melt(cnv_freq[,-1])
colnames(cnv_freq_1)<-c("gene","CNV","value")
head(cnv_freq_1)

cnv_freq[1:5,1:5]
dat=cnv_freq[,c(2:4)]
dat1=melt(dat)
dat1[1:5,]
colnames(dat1)=c('Genes','CNV type','value')
pdf('results/01.Genome Landscape/Fig1b.pdf',height = 5,width = 10)
ggplot(dat1,aes(x= value, y= Genes)) +  
  geom_line(aes(group = Genes)) +
  xlab("Frequency(%)")+ylab('')+
  geom_point(aes(color = `CNV type`), size = 3)+coord_flip()+
  theme_classic()+
  theme(legend.position = 'top',
    axis.text.x = element_text(size = 10,angle = 45,hjust = 1))
dev.off()  




###############基因在肿瘤样本中不同CNV比较
tcga.exp.cmp.for_CNV=get_CNV_Custom(df_cnv=cnv.all,df_exp=tcga_tpm_log_T
                                    ,genes_custom = fatty.acid.fit.gene
                                    ,abbr = 'TCGA-BLCA')
table(tcga.exp.cmp.for_CNV$category)
tcga.exp.cmp.for_CNV$type=factor(tcga.exp.cmp.for_CNV$type,levels = fatty.acid.fit.gene)
head(tcga.exp.cmp.for_CNV)

pdf('results/01.Genome Landscape/Fig1c.pdf',height = 6,width = 18)
ggboxplot(tcga.exp.cmp.for_CNV,
          x="type", y="Expression", color = "category",
          scales = "free_x",
          facet.by = "fa",ncol=3, short.panel.labs = F)+ylab('Gene Expression Level \nlog2(TPM+1)')+
  scale_color_manual(values = ggsci::pal_lancet('lanonc')(9)[c(1,2,3)])+
  stat_compare_means(aes(group=category), label = "p.signif", method = "kruskal.test",label.y = 10)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()


#######02.一致性聚类#############
dir.create('results/02.Cluster')
library(ConsensusClusterPlus)
clusterAlg_name=c('hc','pam','km','kmdist')[2]
distance_name=c('pearson','spearman','euclidean','binary','maximum','canberra','minkowski')[2]
consen_gene=fatty.acid.fit.gene
length(consen_gene)
##############TCGA#######################
tcga_consen_data=as.matrix(tcga_tpm_log_T[intersect(consen_gene,rownames(tcga_tpm_log_T)),])
tcga_consen_data=sweep(tcga_consen_data,1,apply(tcga_consen_data, 1, median))  #22,3
tcga_consen_data=as.matrix(tcga_consen_data)
dim(tcga_consen_data)
tcga_clust_subtype <- ConsensusClusterPlus(tcga_consen_data
                                           , maxK = 10, reps = 500, pItem = 0.8
                                           , pFeature = 1
                                           , title = "TCGA_subtype"
                                           , clusterAlg = clusterAlg_name
                                           , distance = distance_name
                                           , plot = "pdf"
                                           , writeTable = T
                                           , seed = 123456)
k=3
cluster.color=pal_jco()(10)[c(1,2,4)]
tcga.subtype <- data.frame(Samples = names(tcga_clust_subtype[[k]]$consensusClass),
                           Cluster=tcga_clust_subtype[[k]]$consensusClass)
tcga.subtype$Cluster=paste0('C',tcga.subtype$Cluster)
table(tcga.subtype$Cluster)
writeMatrix(tcga.subtype,'results/02.Cluster/tcga.subtype.txt',row=F)
tcga.cluster.km=ggsurvplot(fit=survfit( Surv(OS.time, OS) ~ Cluster,
                                   data = data.frame(OS.time = tcga_cli[rownames(tcga.subtype),]$OS.time/365
                                                     , OS = tcga_cli[rownames(tcga.subtype),]$OS
                                                     , Cluster=tcga.subtype$Cluster)),
                      data=data.frame(OS.time = tcga_cli[rownames(tcga.subtype),]$OS.time/365
                                      , OS = tcga_cli[rownames(tcga.subtype),]$OS
                                      , Cluster=tcga.subtype$Cluster),
                      conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,surv.median.line = 'hv',
                      title='TCGA-BLCA',ggtheme=custom_theme(),
                      linetype = c("solid", "dashed","strata")[1],
                      palette = cluster.color,
                      #legend = c('top', 'bottom', 'left', 'right', 'none')[4],
                      legend = c(0.8,0.75), # 指定图例位置
                      legend.title = "",
                      legend.labs = c("C1","C2","C3"))
fig2d=mg_merge_plot(tcga.cluster.km$plot,tcga.cluster.km$table,ncol = 1,nrow = 2,heights = c(3.5,1),align = 'v')

fam.sub.dat=data.frame(FATTY_ACID_METABOLISM= as.numeric(t(tcga.hall.ssGSEA['FATTY_ACID_METABOLISM',tcga.subtype$Samples])),
                       Cluster=tcga.subtype$Cluster)
fig2e=fam.sub.dat %>%
  ggplot(aes(x=Cluster, y=FATTY_ACID_METABOLISM,fill=Cluster)) +
  geom_violin()+  
  scale_fill_manual(values =cluster.color)+
  geom_boxplot(width=0.2,position=position_dodge(1),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+
  ggpubr::stat_compare_means(aes(group=Cluster), label = "p.signif", method = 'kruskal.test')+
  theme_classic()+labs(x='',y='FATTY ACID METABOLISM score')+
  theme(legend.position = 'none',axis.text = element_text(color = 'black'),
        axis.text.x = element_text(size = 10),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),axis.title.y = element_text(size = 12))
fig2de=mg_merge_plot(fig2d,fig2e,labels = c('D','E'))
savePDF('results/02.Cluster/Fig2de.pdf',fig2de,height = 6,width = 10)



########03.临床特征#########################
dir.create('results/03.subtype.cli')
plotMutiBar_use=function(dat,ist=F,margin=T,xlb='',ylb='',lineCol='black',lineW=0.5,legTitle='Group',showValue=F,showLine=T,xangle=0,isAuto=T){
  library(ggplot2)
  #library(tidyverse)
  #library(reshape2)
  #library(optparse)
  if(ist){
    dat=t(dat)
  }
  lbc=colnames(dat)
  lbr=row.names(dat)
  bk_dat=dat
  if(margin){
    dat=dat%*%diag(1/c(apply(t(dat), 1, sum)))
  }
  row.names(dat)=paste0('R',1:(nrow(dat)))
  colnames(dat)=paste0('C',1:(ncol(dat)))
  row.names(bk_dat)=paste0('R',1:(nrow(bk_dat)))
  colnames(bk_dat)=paste0('C',1:(ncol(bk_dat)))
  #df=cbind(bg=paste0('R',1:nrow(dat)),dat)
  #colnames(df)=c('bg',paste0('C',1:(ncol(dat))))
  tp.dat=as.data.frame(cbind(bg=row.names(dat),dat))
  tp.dat[,1]=as.character(tp.dat[,1])
  for(i in 2:ncol(tp.dat)){
    tp.dat[,i]=as.numeric(as.character(tp.dat[,i]))
  }
  mt.df=reshape2::melt(tp.dat)
  colnames(mt.df)=c('bg','variable','value')
  
  pg=ggplot(mt.df, aes(x=variable, y=value, fill=bg))+geom_bar(stat = "identity", width=lineW, col=lineCol) 
  if(showLine){
    for (i in 2:(ncol(tp.dat)-1)) {
      tmp=tp.dat[order(tp.dat[,1],decreasing = T),]
      tmp[,i]=base::cumsum(tmp[,i])
      tmp[,i+1]=base::cumsum(tmp[,i+1])
      colnames(tmp)[c(i,i+1)]=c('STY','ED')
      tmp1=cbind(tmp,STX=rep(i-1+lineW/2,nrow(tmp))
                 ,EDX=rep(i-lineW/2,nrow(tmp)))
      pg=pg+geom_segment(data=tmp1,aes(x=STX, xend=EDX, y=STY, yend=ED))
    }
  }
  
  if(showValue){
    pg=pg+geom_text(data=mt.df,aes(label=sprintf("%0.2f", round(value, digits = 2))),position=position_stack(vjust=0.5))
  }
  pg=pg+scale_x_discrete(breaks = paste0('C',1:(ncol(dat))),label = lbc)
  pg=pg+labs(x=xlb, y=ylb)+ggsci::scale_fill_npg()+theme(legend.position = "bottom")
  pg=pg+ggsci::scale_fill_npg()+scale_fill_discrete(breaks = paste0('R',1:nrow(dat)),label = lbr,name=legTitle)
  if(xangle>0){
    pg=pg+theme(axis.text.x = element_text(angle = xangle, hjust = 1),legend.position = "bottom")
  }
  
  g.tb=matrix(0,nrow=ncol(dat),ncol=ncol(dat))
  for(i in 1:(ncol(dat))){
    for(j in 1:ncol(dat)){
      if(i!=j){
        g.tb[i,j]=round(-log10((chisq.test(bk_dat[,c(i,j)])$p.value)),2)
      }
    }
  }
  colnames(g.tb)=lbc
  row.names(g.tb)=lbc
  g.tb=reshape2::melt(g.tb) 
  colnames(g.tb)=c('A1','A2','A3')
  g.tb$A4=paste0(g.tb[,3],ifelse(g.tb[,3]>-log10(0.05),'(*)',''))
  stable.p=ggplot(g.tb, aes(A1, A2)) + geom_tile(aes(fill = A3),colour = "white") +xlab('')+ylab('')+ scale_fill_gradient(low = "white",high = "steelblue")+geom_text(aes(x=A1,y=A2,label=A4))+theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank())
  stable.p=stable.p+ggtitle('-log10(Chi-Square test p.value)')
  if(isAuto){
    g1=ggpubr::ggarrange(stable.p,pg, ncol = 1, nrow = 2,heights = c(0.5,1),align = "hv")
    return(g1)
  }else{
    return(list(Bar=pg,Table=stable.p))
  }
}

tcga.subtype.cli=merge(tcga_cli,tcga.subtype,bu='Samples')
rownames(tcga.subtype.cli)=tcga.subtype.cli$Samples
head(tcga.subtype.cli)

chisq.test(table(tcga.subtype.cli$Age1,tcga.subtype.cli$Cluster)   )$p.value #ns
chisq.test(table(tcga.subtype.cli$Gender,tcga.subtype.cli$Cluster) )$p.value #ns
chisq.test(table(tcga.subtype.cli$T.stage,tcga.subtype.cli$Cluster))$p.value #*
chisq.test(table(tcga.subtype.cli$N.stage,tcga.subtype.cli$Cluster))$p.value #ns
chisq.test(table(tcga.subtype.cli$M.stage,tcga.subtype.cli$Cluster))$p.value #*
chisq.test(table(tcga.subtype.cli$Stage,tcga.subtype.cli$Cluster)  )$p.value #***
chisq.test(table(tcga.subtype.cli$Grade,tcga.subtype.cli$Cluster)  )$p.value #***
chisq.test(table(tcga.subtype.cli$Status,tcga.subtype.cli$Cluster) )$p.value #***

table(fatty.acid.gene.cox.fit$HR>1)
fatty.acid.gene.cox.fit$type=ifelse(fatty.acid.gene.cox.fit$HR>1,'Risk','Protective')

gene_anno=data.frame(type=fatty.acid.gene.cox.fit$type[order(fatty.acid.gene.cox.fit$type)])
rownames(gene_anno)=rownames(fatty.acid.gene.cox.fit[order(fatty.acid.gene.cox.fit$type),])

type.col=c( "#FBB4AE" ,"#B3CDE3" )
names(type.col)=c('Risk','Protective')

cluster.color.use=cluster.color
names(cluster.color.use)=c('C1','C2','C3')

age.col=c("#B3E2CD", "#FDCDAC" )
names(age.col)=c('>69','<=69')

sex.col=c("#8DD3C7" ,"#FFFFB3" )
names(sex.col)=c('FEMALE','MALE')

tstage.col=brewer.pal(4,"Set1")
names(tstage.col)=c('T1','T2','T3','T4')

nstage.col=brewer.pal(4,"Set1")
names(nstage.col)=c('N0','N1','N2','N3')

mstage.col=pal_aaas()(10)[1:2]
names(mstage.col)=c('M0','M1')

stage.col=brewer.pal(4,"Set1")
names(stage.col)=c('I','II','III','IV')

grade.col=c("#E41A1C", "#377EB8")
names(grade.col)=c('High Grade','Low Grade')

status.col=c("blue", "red")
names(status.col)=c('Alive','Dead')

colnames(tcga.subtype.cli)
cli_anno=tcga.subtype.cli[order(tcga.subtype.cli$Cluster),c(13,10,11,3:8)]
color_anno=list(type=type.col,Cluster=cluster.color.use,Age1=age.col,Gender=sex.col,
                T.stage=tstage.col,N.stage=nstage.col,M.stage=mstage.col,
                Stage=stage.col,Grade=grade.col,Status=status.col)

pdf('results/03.subtype.cli/Fig3.pdf',height = 9,width =12)
pheatmap(tcga_tpm_log_T[rownames(gene_anno),rownames(cli_anno)],
         scale = 'row',
         color =  circlize::colorRamp2(c(-4, 0, 4), c('#3B4992FF', 'white', '#EE0000FF')),
         main="", # 设置图形标题
         annotation_col = cli_anno,
         annotation_row = gene_anno,
         annotation_colors = color_anno,
		 display_numbers = F, # 热图上显示数值
         cluster_cols = F, # 去掉横向、纵向聚类
         cluster_rows = F,
         show_rownames = T, #去掉横、纵坐标id
         show_colnames = F,
         fontsize_row = 10) 
dev.off()

########04.亚型间免疫#########
dir.create('results/04.subtype.immu')
load('tcga.ciber.RData')
pdf('results/04.subtype.immu/Fig4A.pdf',height =6,width = 15)
my_mutiboxplot(tcga.ciber[tcga.subtype$Samples,1:22]
                   , group = tcga.subtype$Cluster
                   , legend.pos = 'top'
                   , group_cols = cluster.color
                   , test_method = c('kruskal.test','wilcox.test')[1]
                   , fill = 'Group'
                   , ylab = 'Fraction',size = 10)	
dev.off()

load('tcga.est.RData')
pdf('results/04.subtype.immu/Fig4C.pdf',height =4,width = 6)
my_mutiboxplot(tcga.est[tcga.subtype$Samples,]
                   , group = tcga.subtype$Cluster
                   , legend.pos = 'top'
                   , group_cols = cluster.color
                   , test_method =  c('kruskal.test','wilcox.test')[1]
                   ,  fill = 'Group'
                   , ylab = 'Score',size = 10,angle = 0,hjust = 0.5)
dev.off()


load('tcga.immu.ssgsea.RData')
tcga.mcp=read.delim('tcga.mcp.txt')


tcga.immu.ssgsea.use=tcga.immu.ssgsea[tcga.subtype$Samples,]
colnames(tcga.immu.ssgsea.use)=paste0(colnames(tcga.immu.ssgsea.use),'_ssGSEA')
tcga.mcp.use=tcga.mcp[tcga.subtype$Samples,]
colnames(tcga.mcp.use)=paste0(colnames(tcga.mcp.use),'_MCP-COUNT')

immu.dat.sub=cbind(tcga.immu.ssgsea.use[tcga.subtype$Samples,],tcga.mcp.use[tcga.subtype$Samples,])
cli_anno1=data.frame(Cluster=tcga.subtype.cli[order(tcga.subtype.cli$Cluster),'Cluster'])
rownames(cli_anno1)=tcga.subtype.cli[order(tcga.subtype.cli$Cluster),'Samples']

immu.anno=data.frame(Method=rep(c('ssGSEA','MCP-Count'),c(28,10)))
rownames(immu.anno)=colnames(immu.dat.sub)
immu.col=ggsci::pal_jama()(7)[1:2]
names(immu.col)=c('ssGSEA','MCP-Count')

pdf('results/04.subtype.immu/Fig4B.pdf',height = 10,width = 10)
pheatmap(t(immu.dat.sub[rownames(cli_anno1),rownames(immu.anno)]),
         scale = 'row',
         #border="white", # 设置边框为白色
         color =  colorRampPalette(c("navy", "white", "red"))(200),
         main="Differences in immunity scores among cluster", # 设置图形标题
         annotation_col = cli_anno1,
         annotation_row = immu.anno,
         annotation_colors = list(Cluster=cluster.color.use,Method=immu.col),
         display_numbers = F, # 热图上显示数值
         cluster_cols = F, # 去掉横向、纵向聚类
         cluster_rows = F,
         show_rownames = T, #去掉横、纵坐标id
         show_colnames = F,
         #legend_breaks=c(-5,0,5),
         gaps_row = 28,
         gaps_col = c(180,180+124),
         fontsize_row = 10, # 分别设置横向和纵向字体大小
         fontsize_col = 12)
dev.off()


library(tinyarray)
pdf('results/04.subtype.immu/Fig4E.pdf',height =4,width = 7)
draw_heatmap(t(tcga_tpm_log_T[c('CTLA4','PDCD1','PDCD1LG2','CD80','CD28','HAVCR2'),order(tcga.subtype$Cluster)]),
             group_list = factor(tcga.subtype$Cluster[order(tcga.subtype$Cluster)]),
             cluster_cols = F,
             show_rownames = T,
             legend = T,
             annotation_legend = T,
             #n_cutoff = 3,
             color_an = cluster.color)
dev.off()

tcga.Angiogenesis=read.delim('tcga.Angiogenesis.txt')
pdf('results/04.subtype.immu/Fig4D.pdf',height =4,width = 5)
my_boxplot(dat = tcga.Angiogenesis[tcga.subtype$Samples,'Angiogenesis_score'],
         group = tcga.subtype$Cluster,
         group_cols = cluster.color,#label = 'p.signif',
         ylab = 'Angiogenesis score',fill = 'Cluster',
         test_method =  c('kruskal.test','wilcox.test')[1])
dev.off()


#####TIDE########
tcga_tide_res<-read.csv('tcga.tide.res.csv',row.names = 1,stringsAsFactors = F)
head(tcga_tide_res)

tide_sel=c('TIDE','IFNG','Exclusion','Dysfunction','MDSC')
fig4f=my_mutiboxplot_seg(dat = tcga_tide_res[tcga.subtype$Samples,tide_sel[1:4]],
                   group = tcga.subtype$Cluster,
                   test_method = c('t.test','wilcox.test','paired_t.test','paired_wilcox.test','anova','kruskal.test')[6],
                   ylab = 'Score',group_cols =cluster.color,
                   legend.pos = 'none',nrow = 1,ncol = 4,xangle = 0,xhjust = 0.5)
fig4f

tcga_tide_res.sub=cbind(tcga_tide_res[tcga.subtype$Samples,],tcga.subtype)
tide.per=prop.table(table(tcga_tide_res.sub$Cluster,tcga_tide_res.sub$Responder),margin=2)
tide.per=reshape2::melt(tide.per)
colnames(tide.per)<-c("Cluster","Responder","Percentage")
tide.per$Percentage<-round(tide.per$Percentage,digits=2)

fig4g=ggplot(tide.per,aes(x=Cluster,y=Percentage,fill=Responder))+
  geom_bar(position = "fill",stat="identity")+
  theme(legend.position = 'top',axis.text.x = element_text(size = 10))

fig4fg=mg_merge_plot(fig4f,fig4g,labels = c('F','G'),widths = c(4,1))
savePDF('results/04.subtype.immu/Fig4FG.pdf',fig4fg,height = 5,width = 15)

############05.亚型间通路################
dir.create('results/05.subtype.pathway')
############亚型间差异激活的通路
tcga.subtype.use1=tcga.subtype
tcga.subtype.use1$Cluster1=ifelse(tcga.subtype.use1$Cluster=='C1','C1','Other')
tcga.subtype.use1$Cluster2=ifelse(tcga.subtype.use1$Cluster=='C2','C2','Other')
tcga.subtype.use1$Cluster3=ifelse(tcga.subtype.use1$Cluster=='C3','C3','Other')

tcga.geneList1=getGeneFC(gene.exp=tcga_tpm_log_T[,tcga.subtype.use1$Samples],
                         group=tcga.subtype.use1$Cluster1,
                         ulab='C1',
                         dlab='Other')
tcga.geneList2=getGeneFC(gene.exp=tcga_tpm_log_T[,tcga.subtype.use1$Samples],
                         group=tcga.subtype.use1$Cluster2,
                         ulab='C2',
                         dlab='Other')
tcga.geneList3=getGeneFC(gene.exp=tcga_tpm_log_T[,tcga.subtype.use1$Samples],
                         group=tcga.subtype.use1$Cluster3,
                         ulab='C3',
                         dlab='Other')

h.all.gmt<-read.gmt("/pub1/data/mg_projects/users/yanj/project/public/gmt/h.all.v7.5.1.entrez.gmt")
tcga.hallmark.gsea1<-GSEA(tcga.geneList1,TERM2GENE = h.all.gmt,seed=T)
tcga.hallmark.gsea2<-GSEA(tcga.geneList2,TERM2GENE = h.all.gmt,seed=T)
tcga.hallmark.gsea3<-GSEA(tcga.geneList3,TERM2GENE = h.all.gmt,seed=T)

library(enrichplot)
library(ggplot2)

fig5a=clusterProfiler::dotplot(tcga.hallmark.gsea1,split=".sign",showCategory=nrow(tcga.hallmark.gsea1@result),
                         title='TCGA C1 vs Other',font.size=10)+facet_grid(~.sign)+
  scale_y_discrete(labels=function(x) str_remove(x,"HALLMARK_"))
fig5b=clusterProfiler::dotplot(tcga.hallmark.gsea2,split=".sign",showCategory=nrow(tcga.hallmark.gsea2@result),
                         title='TCGA C2 vs Other',font.size=10)+facet_grid(~.sign)+
  scale_y_discrete(labels=function(x) str_remove(x,"HALLMARK_"))
fig5c=clusterProfiler::dotplot(tcga.hallmark.gsea3,split=".sign",showCategory=nrow(tcga.hallmark.gsea3@result),
                         title='TCGA C3 vs Other',font.size=10)+facet_grid(~.sign)+
  scale_y_discrete(labels=function(x) str_remove(x,"HALLMARK_"))
fig5abc=mg_merge_plot(fig5a,fig5b,fig5c,labels = LETTERS[1:3],ncol=3,nrow=1,common.legend = T)
savePDF('results/05.subtype.pathway/Fig5abc.pdf',fig5abc,height = 6,width = 16)

tcga.hallmark.gsea.res1=tcga.hallmark.gsea1@result
tcga.hallmark.gsea.res2=tcga.hallmark.gsea2@result
tcga.hallmark.gsea.res3=tcga.hallmark.gsea3@result
writeMatrix(tcga.hallmark.gsea.res1,'results/05.subtype.pathway/tcga.hallmark.gsea.res1.txt')
writeMatrix(tcga.hallmark.gsea.res2,'results/05.subtype.pathway/tcga.hallmark.gsea.res2.txt')
writeMatrix(tcga.hallmark.gsea.res3,'results/05.subtype.pathway/tcga.hallmark.gsea.res3.txt')

rownames(tcga.hallmark.gsea.res1)=gsub("HALLMARK_","",rownames(tcga.hallmark.gsea.res1))
rownames(tcga.hallmark.gsea.res2)=gsub("HALLMARK_","",rownames(tcga.hallmark.gsea.res2))
rownames(tcga.hallmark.gsea.res3)=gsub("HALLMARK_","",rownames(tcga.hallmark.gsea.res3))

tcga.hallmark.union=Reduce(union,list(rownames(tcga.hallmark.gsea.res1),
                                      rownames(tcga.hallmark.gsea.res2),
                                      rownames(tcga.hallmark.gsea.res3)
))
length(tcga.hallmark.union)#35

tcga.hallmark.heatmap.dat=matrix(0,nrow = 3,ncol = length(tcga.hallmark.union))
rownames(tcga.hallmark.heatmap.dat)=c('C1 vs OTHER', 'C2 vs OTHER','C3 vs OTHER')
colnames(tcga.hallmark.heatmap.dat)=tcga.hallmark.union

tcga.hallmark.heatmap.dat[1,match(rownames(tcga.hallmark.gsea.res1),colnames(tcga.hallmark.heatmap.dat))]=tcga.hallmark.gsea.res1$NES
tcga.hallmark.heatmap.dat[2,match(rownames(tcga.hallmark.gsea.res2),colnames(tcga.hallmark.heatmap.dat))]=tcga.hallmark.gsea.res2$NES
tcga.hallmark.heatmap.dat[3,match(rownames(tcga.hallmark.gsea.res3),colnames(tcga.hallmark.heatmap.dat))]=tcga.hallmark.gsea.res3$NES
range(tcga.hallmark.heatmap.dat)

pdf('results/05.subtype.pathway/Fig5d.pdf',height = 6,width = 16)
Heatmap(as.matrix(tcga.hallmark.heatmap.dat)
        , name = "NES"
        , col = circlize::colorRamp2(c(-3, 0, 3), c('blue', 'white', 'red'))
        , border = T
        , show_column_names = T
        , show_column_dend = F
        , show_row_dend =F
        , cluster_columns=T
        , cluster_rows=F
        , rect_gp = gpar(col = "white", lwd =0.1)
        , row_names_gp = gpar(fontsize = 8)
        ,column_names_gp = gpar(fontsize=10)
        )
dev.off()

#############06.差异基因##########
dir.create('results/06.DEGs')
p_fit=0.05
fc_fit=log2(2)
tcga.degs.c1=mg_limma_DEG(exp = tcga_tpm_log_T,group = as.character(tcga.subtype.use1$Cluster1),ulab = 'C1',dlab = 'Other')
tcga.degs.c1$Summary
writeMatrix(round(tcga.degs.c1$DEG[which(tcga.degs.c1$DEG$adj.P.Val<p_fit & abs(tcga.degs.c1$DEG$logFC)>fc_fit),],3),
           'results/06.DEGs/tcga.c1.degs.txt')
tcga.degs.c1.res=rownames(tcga.degs.c1$DEG[which(tcga.degs.c1$DEG$adj.P.Val<p_fit & abs(tcga.degs.c1$DEG$logFC)>fc_fit),])
length(tcga.degs.c1.res)

tcga.degs.c2=mg_limma_DEG(exp = tcga_tpm_log_T,group = as.character(tcga.subtype.use1$Cluster2),ulab = 'C2',dlab = 'Other')
tcga.degs.c2$Summary
writeMatrix(round(tcga.degs.c2$DEG[which(tcga.degs.c2$DEG$adj.P.Val<p_fit & abs(tcga.degs.c2$DEG$logFC)>fc_fit),],3),
           'results/06.DEGs/tcga.c2.degs.txt')
tcga.degs.c2.res=rownames(tcga.degs.c2$DEG[which(tcga.degs.c2$DEG$adj.P.Val<p_fit & abs(tcga.degs.c2$DEG$logFC)>fc_fit),])
length(tcga.degs.c2.res)

tcga.degs.c3=mg_limma_DEG(exp = tcga_tpm_log_T,group = as.character(tcga.subtype.use1$Cluster3),ulab = 'C3',dlab = 'Other')
tcga.degs.c3$Summary
writeMatrix(round(tcga.degs.c3$DEG[which(tcga.degs.c2$DEG$adj.P.Val<p_fit & abs(tcga.degs.c3$DEG$logFC)>fc_fit),],3),
           'results/06.DEGs/tcga.c3.degs.txt')
tcga.degs.c3.res=rownames(tcga.degs.c3$DEG[which(tcga.degs.c3$DEG$adj.P.Val<p_fit & abs(tcga.degs.c3$DEG$logFC)>fc_fit),])
length(tcga.degs.c3.res)

tcga.sub.degs=unique(c(tcga.degs.c1.res,tcga.degs.c2.res,tcga.degs.c3.res))
length(tcga.sub.degs)
writeMatrix(tcga.sub.degs,'results/06.DEGs/tcga.subtype.degs.txt',header = F)


#######07.PPI#########
dir.create('results/07.PPI')
###将tcga.sub.degs基因输入STRING,置信度得分设置为>0.7
mcode_clust1=read.csv('results/07.PPI/MCODE_Cluster1.csv')
mcode_clust1=mcode_clust1[,c('Degree','name')]
mcode_clust1=mcode_clust1[order(mcode_clust1$Degree,decreasing = T),]

mcode_clust1_genes=mcode_clust1$name
length(mcode_clust1_genes)
#172

gene.enrich=enrichmentORA(mcode_clust1_genes,
                                 mp_dbs=c('pathway_KEGG',
                                          'geneontology_Biological_Process',
                                          'geneontology_Cellular_Component',
                                          'geneontology_Molecular_Function'))
gene.enrich.fit=gene.enrich[gene.enrich$FDR<0.05,]
writeMatrix(gene.enrich.fit,'results/07.PPI/mcode1.gene.enrich.txt')
dim(gene.enrich.fit)
table(gene.enrich.fit$DB)
pdf('results/07.PPI/enrich_dotplot.pdf',height = 10,width = 14)
dotplot_batch(gene.enrich.fit,dbs = c('pathway_KEGG',
                                      'geneontology_Biological_Process',
                                      'geneontology_Cellular_Component',
                                      'geneontology_Molecular_Function'),top = 10,FDR = T)

dev.off()

#######08.风险模型###############
dir.create('results/08.Module')
pre.genes=mcode_clust1_genes

tcga_model_data=t(tcga_tpm_log_T[pre.genes,tcga_cli$Samples])
colnames(tcga_model_data)=gsub('-','__',colnames(tcga_model_data))
tcga_model_data=merge(data.frame(Samples=tcga_cli$Samples,OS=tcga_cli$OS,OS.time=tcga_cli$OS.time),
                      data.frame(Samples=rownames(tcga_model_data),tcga_model_data),by='Samples')
rownames(tcga_model_data)=tcga_model_data$Samples
tcga_model_data=tcga_model_data[,-1]
tcga_model_data=crbind2DataFrame(tcga_model_data)
dim(tcga_model_data)

tcga.cox=cox_batch(dat = tcga_tpm_log_T[pre.genes,tcga_cli$Samples],time = tcga_cli$OS.time,event = tcga_cli$OS)
tcga.cox=na.omit(tcga.cox)
tcga.cox$coef=log(tcga.cox$HR)
head(tcga.cox)
rownames(tcga.cox)=gsub('-','__',rownames(tcga.cox))
p_cutoff=0.05
table(tcga.cox$p.value<p_cutoff)
tcga.cox$Type=rep('None',nrow(tcga.cox))
tcga.cox$Type[which(tcga.cox$p.value<p_cutoff & tcga.cox$HR>1)]='Risk'
tcga.cox$Type[which(tcga.cox$p.value<p_cutoff & tcga.cox$HR<1)]='Protective'
tcga.cox.fit=tcga.cox[tcga.cox$p.value<p_cutoff,]
table(tcga.cox.fit$Type)
writeMatrix(round(tcga.cox[rownames(tcga.cox.fit),-6],3),'results/08.Module/tcga.cox.txt')
tcga.lasso=get_riskscore.lasso(dat = tcga_model_data[tcga_cli$Samples,rownames(tcga.cox.fit)],
                               os = tcga_model_data[tcga_cli$Samples,]$OS,
                               os.time = tcga_model_data[tcga_cli$Samples,]$OS.time)
length(tcga.lasso$lasso.gene)#13
fig8ab=tcga.lasso$plot

tcga.module.risk=get_riskscore(dat = tcga_model_data[tcga_cli$Samples,tcga.lasso$lasso.gene],
                               #dat=tcga_model_data[tcga_cli$Samples,rownames(tcga.cox.fit)],
                               os = tcga_model_data[tcga_cli$Samples,]$OS,
                               os.time = tcga_model_data[tcga_cli$Samples,]$OS.time,
                               step = T,direction = c("both", "backward", "forward")[1])
length(tcga.module.risk$module.gene)
#5
tcga.module.risk$model
# "-0.171*SRC+-0.303*CTLA4+-0.092*CTSE+0.29*LAMA2+0.12*ADAMTSL4"

fmla <- as.formula(paste0("Surv(OS.time, OS) ~"
                          ,paste0(tcga.module.risk$module.gene,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga_model_data))
fig8c=ggforest(cox, data =tcga_model_data, 
         main = "Hazardratio", 
         fontsize =1.0, 
         noDigits = 2)#画森林图
fig8c

tcga.risktype.cli=data.frame(tcga_cli,
                             Cluster=tcga.subtype[tcga_cli$Samples,]$Cluster,
                             Riskscore=tcga.module.risk$result[tcga_cli$Samples,]$riskscore,
                             Risktype=ifelse(tcga.module.risk$result[tcga_cli$Samples,]$riskscore>
                                               median(tcga.module.risk$result$riskscore),'High','Low'))
writeMatrix(tcga.risktype.cli,'results/08.Module/tcga.risktype.cli.txt',row = F)

fig8d=ggplotTimeROC(time = tcga.module.risk$result$time/365
             ,status = tcga.module.risk$result$status
             ,score = tcga.module.risk$result$riskscore
             ,mks = c(1,3,5))


tcga.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga.risktype.cli),
                   data=tcga.risktype.cli,
                   conf.int = T,pval = T,fun = "pct",risk.table = T, size = 0.7,
                   title='TCGA-BLCA',ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   palette = ggsci::pal_lancet()(8),
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   legend.title = "",
                   legend.labs = c("High","Low"))
fig8e=mg_merge_plot(tcga.km$plot,tcga.km$table,ncol = 1,nrow = 2,heights = c(3.5,1),align = 'v')


tcga.os.per=prop.table(table(tcga.risktype.cli$Status,tcga.risktype.cli$Risktype),margin=2)
tcga.os.per=reshape2::melt(tcga.os.per)
colnames(tcga.os.per)<-c("Status","Risktype","Percentage")
tcga.os.per$Percentage<-round(tcga.os.per$Percentage,digits=2)
fig8f=ggplot(tcga.os.per,aes(x=Risktype,y=Percentage,fill=Status))+
  geom_bar(position = "fill",stat="identity")+
  scale_fill_manual(values = c('red','blue'))+
  labs(x="",y = "Percentage",fill = "Status",title= paste0('-log10(Chi-Squared Test pvalue) = ',round(-log10(chisq.test(table(tcga.risktype.cli$Status,tcga.risktype.cli$Risktype))$p.value),3)))+ theme_bw()+coord_flip()



##############GSE32894########
GSE32894_model_data <- t(GSE32894_exp[intersect(rownames(GSE32894_exp),pre.genes),])
colnames(GSE32894_model_data)=gsub('-','__',colnames(GSE32894_model_data))
GSE32894_model_data=merge(data.frame(Samples=GSE32894_cli$Samples,OS.time=GSE32894_cli$OS.time,OS=GSE32894_cli$OS),
                          data.frame(Samples=rownames(GSE32894_model_data),GSE32894_model_data),
                          by='Samples')
#colnames(GSE32894_model_data)
rownames(GSE32894_model_data)=GSE32894_model_data$Samples
GSE32894_model_data[1:5,1:5]
GSE32894_model_data=GSE32894_model_data[,-1]
GSE32894_model_data=crbind2DataFrame(GSE32894_model_data)
dim(GSE32894_model_data)

GSE32894.module.risk=get_riskscore(dat = GSE32894_model_data[,intersect(tcga.module.risk$module.gene,colnames(GSE32894_model_data))],
                                   os=GSE32894_model_data$OS,
                                   os.time = GSE32894_model_data$OS.time,
                                   step=F,
                                   direction=c("both", "backward", "forward")[1])
length(GSE32894.module.risk$module.gene)

GSE32894.risktype.cli=data.frame(Samples=rownames(GSE32894.module.risk$result),
                                 OS=GSE32894.module.risk$result$status,
                                 OS.time=GSE32894.module.risk$result$time,
                                 Riskscore=GSE32894.module.risk$result$riskscore,
                                 Risktype=ifelse(GSE32894.module.risk$result$riskscore>
                                                   median(GSE32894.module.risk$result$riskscore),'High','Low')
                                 )
fig8g=ggplotTimeROC(GSE32894_model_data$OS.time,
              GSE32894_model_data$OS,
              GSE32894.module.risk$result$riskscorez,mks = c(1,3,5))
fig8g
GSE32894.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = GSE32894.risktype.cli),
           data=GSE32894.risktype.cli,
           conf.int = T,pval = T,fun = "pct",risk.table = T, size = 0.7,
           title='GSE32894',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           palette = ggsci::pal_lancet()(8),
           legend = c('top', 'bottom', 'left', 'right', 'none')[1],
           legend.title = "",
           legend.labs = c("High","Low"))
fig8h=mg_merge_plot(GSE32894.km$plot,GSE32894.km$table,ncol = 1,nrow = 2,heights = c(3.5,1),align = 'v')
fig8h

GSE32894.os.per=prop.table(table(GSE32894.risktype.cli$OS,GSE32894.risktype.cli$Risktype),margin=2)
GSE32894.os.per=reshape2::melt(GSE32894.os.per)
colnames(GSE32894.os.per)<-c("Status","Risktype","Percentage")
GSE32894.os.per$Status=ifelse(GSE32894.os.per$Status==0,'Alive','Dead')
GSE32894.os.per$Percentage<-round(GSE32894.os.per$Percentage,digits=2)
fig8i=ggplot(GSE32894.os.per,aes(x=Risktype,y=Percentage,fill=Status))+
  geom_bar(position = "fill",stat="identity")+
  scale_fill_manual(values = c('red','blue'))+
  labs(x="",y = "Percentage",fill = "Status",title= paste0('-log10(Chi-Squared Test pvalue) = ',round(-log10(chisq.test(table(GSE32894.risktype.cli$OS,GSE32894.risktype.cli$Risktype))$p.value),3)))+ theme_bw()+coord_flip()


fig8=mg_merge_plot(mg_merge_plot(fig8ab,fig8c,widths = c(1.5,1)),
                   mg_merge_plot(fig8e,fig8d,fig8f,fig8h,fig8g,fig8i,ncol=3,nrow=2,labels = LETTERS[4:9]),
                   nrow=2,ncol=1,heights = c(1,2))
savePDF('results/08.Module/Fig8.pdf',fig8,height = 15,width = 16)



###########09.风险分组间临床特征###############
dir.create('results/09.risktype.cli')
chisq.test(table(tcga.risktype.cli$Age1,tcga.risktype.cli$Risktype)   )$p.value 
chisq.test(table(tcga.risktype.cli$Gender,tcga.risktype.cli$Risktype) )$p.value 

chisq.test(table(tcga.risktype.cli$T.stage,tcga.risktype.cli$Risktype))$p.value 
chisq.test(table(tcga.risktype.cli$N.stage,tcga.risktype.cli$Risktype))$p.value 
chisq.test(table(tcga.risktype.cli$M.stage,tcga.risktype.cli$Risktype))$p.value 
chisq.test(table(tcga.risktype.cli$Stage,tcga.risktype.cli$Risktype)  )$p.value 
chisq.test(table(tcga.risktype.cli$Grade,tcga.risktype.cli$Risktype)  )$p.value 
chisq.test(table(tcga.risktype.cli$Status,tcga.risktype.cli$Risktype) )$p.value 


fig9a=list()
cli.per1=prop.table(table(tcga.risktype.cli$T.stage,tcga.risktype.cli$Risktype),margin=2)
cli.per1=reshape2::melt(cli.per1)
colnames(cli.per1)<-c("T.stage","Risktype","Percentage")
cli.per1$Percentage<-round(cli.per1$Percentage,digits=2)
fig9a[[1]]=ggplot(cli.per1,aes(x=Risktype,y=Percentage,fill=T.stage))+
  geom_bar(position = "fill",stat="identity")+
  labs(x="",y = "Percentage",fill = "T.stage",title= paste0('-log10(Chi-Squared Test pvalue) = ',round(-log10(chisq.test(table(tcga.risktype.cli$T.stage,tcga.risktype.cli$Risktype))$p.value),3)))+
  theme_bw()#+ 


cli.per2=prop.table(table(tcga.risktype.cli$N.stage,tcga.risktype.cli$Risktype),margin=2)
cli.per2=reshape2::melt(cli.per2)
colnames(cli.per2)<-c("N.stage","Risktype","Percentage")
cli.per2$Percentage<-round(cli.per2$Percentage,digits=2)
fig9a[[2]]=ggplot(cli.per2,aes(x=Risktype,y=Percentage,fill=N.stage))+
  geom_bar(position = "fill",stat="identity")+
  labs(x="",y = "Percentage",fill = "N.stage",title= paste0('-log10(Chi-Squared Test pvalue) = ',round(-log10(chisq.test(table(tcga.risktype.cli$N.stage,tcga.risktype.cli$Risktype))$p.value),3)))+
  theme_bw()#+ geom_text(aes(label = Percentage),position=position_stack(vjust =0.5),size = 5)

cli.per3=prop.table(table(tcga.risktype.cli$M.stage,tcga.risktype.cli$Risktype),margin=2)
cli.per3=reshape2::melt(cli.per3)
colnames(cli.per3)<-c("M.stage","Risktype","Percentage")
cli.per3$Percentage<-round(cli.per3$Percentage,digits=2)
fig9a[[3]]=ggplot(cli.per3,aes(x=Risktype,y=Percentage,fill=M.stage))+
  geom_bar(position = "fill",stat="identity")+
  labs(x="",y = "Percentage",fill = "N.stage",title= paste0('-log10(Chi-Squared Test pvalue) = ',round(-log10(chisq.test(table(tcga.risktype.cli$M.stage,tcga.risktype.cli$Risktype))$p.value),3)))+
  theme_bw()#+ geom_text(aes(label = Percentage),position=position_stack(vjust =0.5),size = 5)


cli.per4=prop.table(table(tcga.risktype.cli$Stage,tcga.risktype.cli$Risktype),margin=2)
cli.per4=reshape2::melt(cli.per4)
colnames(cli.per4)<-c("Stage","Risktype","Percentage")
cli.per4$Percentage<-round(cli.per4$Percentage,digits=2)
fig9a[[4]]=ggplot(cli.per4,aes(x=Risktype,y=Percentage,fill=Stage))+
  geom_bar(position = "fill",stat="identity")+
  labs(x="",y = "Percentage",fill = "Stage",title= paste0('-log10(Chi-Squared Test pvalue) = ',round(-log10(chisq.test(table(tcga.risktype.cli$Stage,tcga.risktype.cli$Risktype))$p.value),3)))+
  theme_bw()#+ geom_text(aes(label = Percentage),position=position_stack(vjust =0.5),size = 5)


cli.per5=prop.table(table(tcga.risktype.cli$Grade,tcga.risktype.cli$Risktype),margin=2)
cli.per5=reshape2::melt(cli.per5)
colnames(cli.per5)<-c("Grade","Risktype","Percentage")
cli.per5$Percentage<-round(cli.per5$Percentage,digits=2)
fig9a[[5]]=ggplot(cli.per5,aes(x=Risktype,y=Percentage,fill=Grade))+
  geom_bar(position = "fill",stat="identity")+
  labs(x="",y = "Percentage",fill = "Grade",title= paste0('-log10(Chi-Squared Test pvalue) = ',round(-log10(chisq.test(table(tcga.risktype.cli$Grade,tcga.risktype.cli$Risktype))$p.value),3)))+ theme_bw()
  
fig9a=mg_merge_plot(fig9a,ncol = 5,nrow = 1,labels = LETTERS[1:5])
fig9a



tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)
table(tcga_cox_datas$Age1)
table(tcga_cox_datas$T.stage)
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T1'|tcga_cox_datas$T.stage=='T2']<-'T1+T2'
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T3'|tcga_cox_datas$T.stage=='T4']<-'T3+T4'

table(tcga_cox_datas$N.stage)
tcga_cox_datas$N.stage[tcga_cox_datas$N.stage=='N1'|tcga_cox_datas$N.stage=='N2'|tcga_cox_datas$N.stage=='N3']<-'N1+N2+N3'

table(tcga_cox_datas$M.stage)

table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='I'|tcga_cox_datas$Stage=='II']<-'I+II'
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='III'|tcga_cox_datas$Stage=='IV']<-'III+IV'

table(tcga_cox_datas$Grade)


fig9b=list()
cli.dat1=data.frame(Group=tcga_cox_datas$T.stage,value=tcga_cox_datas$Riskscore)
cli.dat1=crbind2DataFrame(cli.dat1[which(cli.dat1$Group!='NA'),])
cli.dat1=melt(cli.dat1)

fig9b[[1]]=
  cli.dat1 %>%
  ggplot(aes(x=Group, y=value,fill=Group)) +
  geom_violin()+  
  scale_fill_manual(values = pal_nejm()(10))+
  geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+labs(x='',y='Risk Score',title = 'T.stage')+
  ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method =c('wilcox.test','anova')[1])+
  theme(legend.position = 'none',axis.text = element_text(color = 'black'),
        title = element_text(size = 12),
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


cli.dat2=data.frame(Group=tcga_cox_datas$N.stage,value=tcga_cox_datas$Riskscore)
cli.dat2=crbind2DataFrame(cli.dat2[which(cli.dat2$Group!='NA'),])
cli.dat2=melt(cli.dat2)

fig9b[[2]]=cli.dat2 %>%
  ggplot(aes(x=Group, y=value,fill=Group)) +
  geom_violin()+  
  scale_fill_manual(values = pal_nejm()(10))+
  geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+labs(x='',y='Risk Score',title = 'N.stage')+
  ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method =c('wilcox.test','anova')[1])+
  theme(legend.position = 'none',axis.text = element_text(color = 'black'),
        title = element_text(size = 12),
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


cli.dat3=data.frame(Group=tcga_cox_datas$M.stage,value=tcga_cox_datas$Riskscore)
cli.dat3=crbind2DataFrame(cli.dat3[which(cli.dat3$Group!='NA'),])
cli.dat3=melt(cli.dat3)

fig9b[[3]]=cli.dat3 %>%
  ggplot(aes(x=Group, y=value,fill=Group)) +
  geom_violin()+  
  scale_fill_manual(values = pal_nejm()(10))+
  geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+labs(x='',y='Risk Score',title = 'M.stage')+
  ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method =c('wilcox.test','anova')[1])+
  theme(legend.position = 'none',axis.text = element_text(color = 'black'),
        title = element_text(size = 12),
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


cli.dat4=data.frame(Group=tcga_cox_datas$Stage,value=tcga_cox_datas$Riskscore)
cli.dat4=crbind2DataFrame(cli.dat4[which(cli.dat4$Group!='NA'),])
cli.dat4=melt(cli.dat4)

fig9b[[4]]=cli.dat4 %>%
  ggplot(aes(x=Group, y=value,fill=Group)) +
  geom_violin()+  
  scale_fill_manual(values = pal_nejm()(10))+
  geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+labs(x='',y='Risk Score',title = 'Stage')+
  ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method =c('wilcox.test','anova')[2])+
  theme(legend.position = 'none',axis.text = element_text(color = 'black'),
        title = element_text(size = 12),
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))


cli.dat5=data.frame(Group=tcga_cox_datas$Grade,value=tcga_cox_datas$Riskscore)
cli.dat5=crbind2DataFrame(cli.dat5[which(cli.dat5$Group!='NA'),])
cli.dat5=melt(cli.dat5)

fig9b[[5]]=cli.dat5 %>%
  ggplot(aes(x=Group, y=value,fill=Group)) +
  geom_violin()+  
  scale_fill_manual(values = pal_nejm()(10))+
  geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+labs(x='',y='Risk Score',title = 'Grade')+
  ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method =c('wilcox.test','anova')[1])+
  theme(legend.position = 'none',axis.text = element_text(color = 'black'),
        title = element_text(size = 12),
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))

fig9b=mg_merge_plot(fig9b,ncol = 5,nrow = 1,labels = LETTERS[6:10])
fig9b

fig9c=list()
table(tcga_cox_datas$Age1)
cli.km1=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$Age1=='<=69'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$Age1=='<=69'),],
           conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
           title='Age <=69',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           #palette =risktype.color,
           #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "",
           legend.labs = c("High","Low"))
cli.km2=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$Age1=='>69'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$Age1=='>69'),],
           conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
           title='Age >69',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           #palette =risktype.color,
           #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "",
           legend.labs = c("High","Low"))

table(tcga_cox_datas$Gender)
cli.km3=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$Gender=='FEMALE'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$Gender=='FEMALE'),],
           conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
           title='FEMALE',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           #palette =risktype.color,
           #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "",
           legend.labs = c("High","Low"))
cli.km4=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$Gender=='MALE'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$Gender=='MALE'),],
           conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
           title='MALE',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           #palette =risktype.color,
           #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "",
           legend.labs = c("High","Low"))

table(tcga_cox_datas$Stage)
cli.km5=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$Stage=='I+II'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$Stage=='I+II'),],
           conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
           title='Stage I+II',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           #palette =risktype.color,
           #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "",
           legend.labs = c("High","Low"))
cli.km6=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$Stage=='III+IV'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$Stage=='I+II'),],
           conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
           title='Stage III+IV',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           #palette =risktype.color,
           #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "",
           legend.labs = c("High","Low"))

table(tcga_cox_datas$Grade)
cli.km7=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$Grade=='High Grade'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$Grade=='High Grade'),],
           conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
           title='High Grade',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           #palette =risktype.color,
           #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "",
           legend.labs = c("High","Low"))


table(tcga_cox_datas$T.stage)
cli.km8=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$T.stage=='T1+T2'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$T.stage=='T1+T2'),],
           conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
           title='T1+T2',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           #palette =risktype.color,
           #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "",
           legend.labs = c("High","Low"))
cli.km9=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$T.stage=='T3+T4'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$T.stage=='T3+T4'),],
           conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
           title='T3+T4',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           #palette =risktype.color,
           #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "",
           legend.labs = c("High","Low"))

table(tcga_cox_datas$N.stage)
cli.km10=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$N.stage=='N0'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$N.stage=='N0'),],
           conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
           title='N0',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           #palette =risktype.color,
           #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "",
           legend.labs = c("High","Low"))
cli.km11=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                        data = tcga_cox_datas[which(tcga_cox_datas$N.stage=='N1+N2+N3'),]),
           data=tcga_cox_datas[which(tcga_cox_datas$N.stage=='N1+N2+N3'),],
           conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
           title='N1+N2+N3',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           #palette =risktype.color,
           #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "",
           legend.labs = c("High","Low"))

table(tcga_cox_datas$M.stage)
cli.km12=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                 data = tcga_cox_datas[which(tcga_cox_datas$M.stage=='M0'),]),
                    data=tcga_cox_datas[which(tcga_cox_datas$M.stage=='M0'),],
                    conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
                    title='M0',ggtheme=custom_theme(),
                    linetype = c("solid", "dashed","strata")[1],
                    #palette =risktype.color,
                    #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
                    legend = c(0.8,0.75), # 指定图例位置
                    legend.title = "",
                    legend.labs = c("High","Low"))

fig9c=mg_merge_plot(cli.km1$plot,cli.km2$plot,cli.km3$plot,cli.km4$plot,
                    cli.km5$plot,cli.km6$plot,cli.km7$plot,cli.km12$plot,
                    cli.km8$plot,cli.km9$plot,cli.km10$plot,cli.km11$plot,
                    ncol=4,nrow=3)
fig9c


fig9=mg_merge_plot(fig9a,fig9b,fig9c,nrow=3,ncol=1,heights = c(1,1,3),labels = c('','','K'))
savePDF('results/09.risktype.cli/Fig9.pdf',fig9,height = 16,width = 14)


##########10.风险分组间免疫特征#####################
risktype.col= ggsci::pal_aaas()(10)[1:2]
dir.create('results/10.risktype.immu')
pdf('results/10.risktype.immu/Fig10a.pdf',height = 5,width = 15)
mg_PlotMutiBoxplot(tcga.ciber[tcga.risktype.cli$Samples,1:22]
           , group = tcga.risktype.cli$Risktype
           , legend.pos = 'top'
           ,group_cols = risktype.col
           , test_method = c('kruskal.test','wilcox.test')[2]
           , add = 'boxplot'
           , ylab = 'Fraction')	
dev.off()


############29TMEsignature
tme.type=readxl::read_excel('origin_datas/TME.geneSets.classification.PMID34019806.xlsx',sheet = "Sheet2")
tme.type=data.frame(tme.type,check.names = F,stringsAsFactors = F)
tme.type$Group=factor(tme.type$Group,levels=unique(tme.type$Group))
sort(tme.type$`Process/Signature`)
tme.anno=data.frame(group=tme.type$Group)
rownames(tme.anno)=tme.type$`Process/Signature`
####################TME geneset
tme.genesets=readxl::read_excel('origin_datas/29signatrue.PMID34019806.xlsx')
tme.genesets=data.frame(tme.genesets)
head(tme.genesets)
names(table(tme.genesets$Gene.signature))
sort(names(table(tme.genesets$Gene.signature)))
tme.genesets.list=split(x=tme.genesets,f=tme.genesets$Gene.signature)
tme.genesets.list=sapply(tme.genesets.list, function(x){subset(x,select='Gene',drop=TRUE)})

#tcga.tme.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = tcga_tpm_log_T
#                                                ,genelist = tme.genesets.list)
#save(tcga.tme.ssgsea,file='results/tcga.tme.ssgsea.RData')
load('results/tcga.tme.ssgsea.RData')
pdf('results/10.risktype.immu/Fig10b.pdf',height = 10,width = 6)
mg_PlotMutiBoxplot(data = t(tcga.tme.ssgsea[rownames(tme.anno),tcga.risktype.cli$Samples])
                   ,group = tcga.risktype.cli$Risktype
                   ,test_method = 'wilcox.test'
                   ,legend.pos = 'top'
                   ,ylab = 'score',xangle = 0
                   ,group_cols = risktype.col
                   ,add = 'boxplot')+coord_flip()
dev.off()

immu.dat.RS=cbind(tcga.risktype.cli$Riskscore,
                  tcga_model_data[tcga.risktype.cli$Samples,tcga.module.risk$module.gene],
                  t(tcga.tme.ssgsea[rownames(tme.anno),tcga.risktype.cli$Samples]))
colnames(immu.dat.RS)[1]='Riskcsore'

cor_res <- Hmisc::rcorr(as.matrix(immu.dat.RS),type = 'spearman')
cor_res$P[is.na(cor_res$P)] <- 0
p.dat=cor_res$P
p.dat[1:5,1:5]
p.dat<-ifelse(p.dat<0.0001,'****',ifelse(p.dat<0.001,'***', ifelse(p.dat<0.01,'**',ifelse(p.dat<0.05,'*',''))))

tme.col=pal_jama()(9)[1:4]
names(tme.col)=names(table(tme.anno$group))

pdf('results/10.risktype.immu/Fig10c.pdf',height = 7,width = 10,onefile = F)
pheatmap(cor_res$r[-c(1:6),c(tcga.module.risk$module.gene,'Riskcsore')],
         color = circlize::colorRamp2(c(-1, 0, 1), c('#3B4992FF', 'white', '#EE0000FF')),
         annotation_row = tme.anno,
         annotation_colors = list(group=tme.col),
         display_numbers = p.dat[-c(1:6),c(tcga.module.risk$module.gene,'Riskcsore')], # 热图上显示数值
         cluster_cols = F, # 去掉横向、纵向聚类
         cluster_rows = F,
         show_rownames = T, #去掉横、纵坐标id
         show_colnames = T,
         gaps_row = c(8,8+7,8+7+12),
         fontsize_row =10, # 分别设置横向和纵向字体大小
         fontsize_col = 10)
dev.off()


##########PMID:28650338
T_cell_GEP=readMatrix('origin_datas/T_cell-inflamed.txt')
T_cell_GEP=T_cell_GEP$Gene
T_cell_GEP_score=ssGSEAScore_by_genes(gene.exp = tcga_tpm_log_T,genes = T_cell_GEP)
T_cell_GEP_score=data.frame(Samples=colnames(T_cell_GEP_score),score=as.numeric(T_cell_GEP_score))
rownames(T_cell_GEP_score)=T_cell_GEP_score$Samples
writeMatrix(T_cell_GEP_score,'results/10.risktype.immu/T_cell_GEP_score.txt',row = F)
fig10d=my_boxplot(dat = T_cell_GEP_score[tcga.risktype.cli$Samples,'score'],
           group = tcga.risktype.cli$Risktype,
           group_cols =  pal_nejm()(8),x.size = 10,
           ylab = 'T-cell–inflamed GEP score',fill = 'Risktype',
           test_method =  c('kruskal.test','wilcox.test')[2])
fig10d

###########血管生成
fig10e=my_boxplot(dat = tcga.Angiogenesis[tcga.risktype.cli$Samples,'Angiogenesis_score'],
           group = tcga.risktype.cli$Risktype,
           group_cols =  pal_nejm()(8),x.size = 10,
           ylab = 'Angiogenesis score',fill = 'Risktype',
           test_method =  c('kruskal.test','wilcox.test')[2])
fig10e
fig10de=mg_merge_plot(fig10e,fig10d,labels = c('E','F'),common.legend = T)
savePDF('results/10.risktype.immu/Fig10de.pdf',fig10de,height = 4,width = 9)


############TIDE
fig10f=my_mutiboxplot_seg(dat = tcga_tide_res[tcga.risktype.cli$Samples,tide_sel],
                   group = tcga.risktype.cli$Risktype,
                   test_method = c('t.test','wilcox.test','paired_t.test','paired_wilcox.test','anova','kruskal.test')[2],
                   ylab = 'Score',group_cols =risktype.col,
                   legend.pos = 'none',nrow = 1,ncol = 5,xangle = 0,xhjust = 0.5)
fig10f
library(pRRophetic)
library(ggplot2)
############### Cisplatin,顺铂
# set.seed(12345)
# predictedPtype_Cisplatin <- pRRopheticPredict(as.matrix(tcga_tpm_log_T)
#                                               , "Cisplatin"
#                                               , selection=1
#                                               ,dataset = "cgp2016")
# predictedPtype_Cisplatin <- data.frame(predictedPtype_Cisplatin)
# 
# tcga_durg_ic50_res <- predictedPtype_Cisplatin

drugs <- c("Cisplatin","Erlotinib","Rapamycin","Sunitinib","PHA-665752","MG-132","Paclitaxel","Cyclopamine","AZ628","Sorafenib","VX-680","Imatinib","TAE684","Crizotinib","Saracatinib","S-Trityl-L-cysteine","Z-LLNle-CHO","Dasatinib","GNF-2","CGP-60474","CGP-082996","A-770041","WH-4-023","WZ-1-84","BI-2536","BMS-509744","CMK","Pyrimethamine","JW-7-52-1","A-443654","GW843682X","MS-275","Parthenolide","KIN001-135","TGX221","Bortezomib","XMD8-85","Roscovitine","Salubrinal","Lapatinib","Vinorelbine","NSC-87877","QS11","CP466722","Midostaurin","Shikonin","AKT inhibitor VIII","Embelin","Bexarotene","Bleomycin","Phenformin")
length(drugs)
# for (drug in drugs) {
#  print(drug)
#  set.seed(12345)
#  tmpic50 <- pRRopheticPredict(as.matrix(tcga_tpm_log_T)
#                               , drug
#                               , selection=1
#                               , dataset = "cgp2016")
#  tmpic50 <- data.frame(tmpic50)
#  colnames(tmpic50) <- drug
#  tcga_durg_ic50_res <- cbind(tcga_durg_ic50_res, tmpic50)
# }
# save(tcga_durg_ic50_res,file='results/tcga_durg_ic50.RData')
load('results/tcga_durg_ic50.RData')
head(tcga_durg_ic50_res)
tcga_durg_ic50_res=tcga_durg_ic50_res[,-1]
writeMatrix(tcga_durg_ic50_res,'results/10.risktype.immu/tcga_durg_ic50_res.txt')

tcga.ic50.dat=cbind(tcga.risktype.cli[,'Riskscore'],tcga_durg_ic50_res[tcga.risktype.cli$Samples,])
head(tcga.ic50.dat)
colnames(tcga.ic50.dat)[1]='Riskscore'
ic50.cor.RS=Hmisc::rcorr(as.matrix(tcga.ic50.dat),type = 'spearman')
ic50.cor.RS.res=data.frame(Names=names(ic50.cor.RS$r['Riskscore',]),
                           cor=as.numeric(ic50.cor.RS$r['Riskscore',]),
                           p.val=as.numeric(ic50.cor.RS$P['Riskscore',]))
ic50.cor.RS.res=ic50.cor.RS.res[-1,]
head(ic50.cor.RS.res)
colnames(ic50.cor.RS.res)=c('Drugs','cor','pvalue')
ic50.cor.RS.res$Drugs=factor(ic50.cor.RS.res$Drugs,
                             levels = ic50.cor.RS.res$Drugs[order(ic50.cor.RS.res$cor,decreasing = T)], ordered=TRUE)
ic50.cor.RS.res$pvalue=ifelse(ic50.cor.RS.res$pvalue==0,1e-16,ic50.cor.RS.res$pvalue)
head(ic50.cor.RS.res)

high.suit.drug=as.character(ic50.cor.RS.res$Drugs[which(ic50.cor.RS.res$pvalue<0.05 & ic50.cor.RS.res$cor < -0.195)])
length(high.suit.drug)
low.suit.drug=as.character(ic50.cor.RS.res$Drugs[which(ic50.cor.RS.res$pvalue<0.05 & ic50.cor.RS.res$cor  > 0.23)])
length(low.suit.drug)

fig10g=my_mutiboxplot_seg(dat = tcga_durg_ic50_res[tcga.risktype.cli$Samples,c(high.suit.drug,low.suit.drug)],
                   group = tcga.risktype.cli$Risktype,legend.position = 'none',
                   group_cols = ggsci::pal_aaas()(10)[1:2],ylab = 'IC50',xsize = 12,ysize = 12,
                   test_method = 'wilcox.test',xangle = 0,xhjust = 0.5,ncol = 4,nrow=2)+theme_bw()
fig10g

fig10fg=mg_merge_plot(fig10f,fig10g,ncol=1,nrow=2,labels = c('F','G'),heights = c(1,2),common.legend = T)
savePDF('results/10.risktype.immu/Fig10fg.pdf',fig10fg,height = 8,width = 9)


################Riskscore与免疫的相关性
library(tidyverse)
library(ggcor)
library(vegan)
cr=psych::corr.test(x=tcga.risktype.cli$Riskscore,y=tcga.ciber[tcga.risktype.cli$Samples,1:22]
                    ,method = 'spearman')
df=t(rbind(cr$r,cr$p))
colnames(df)=c('r','p.value')
df=data.frame(RiskScore='RiskScore',gene.signatures=rownames(df),df)
df
df <- df %>%
  mutate(lty = cut(r, breaks = c(-1, 0, 1),
                   labels = c("r <= 0", "r > 0")),
         col = cut(p.value, breaks = c(0, 0.01, 0.05, 1),
                   labels = c("< 0.01", "< 0.05", ">= 0.05"),
                   right = FALSE, include.lowest = TRUE))
head(df)
writeMatrix(df[,c(1:4)],'results/10.risktype.immu/RS_cor_cibersort.txt',row = F)


corrmat.color=colorRampPalette(c('blue', 'white','red'))(100)
pdf('results/10.risktype.immu/Fig10h.pdf',height = 10,width =6)
quickcor(tcga.ciber[tcga.risktype.cli$Samples,1:22], cor.test = TRUE,type = "lower") + #upper
  geom_square(data = get_data(p.value < 0.05, type = "lower")) + 
  anno_link(df, mapping = aes(colour = col,
                              size = abs(r),
                              linetype = lty),diag.label = TRUE) +
  scale_fill_gradient2n(colours = corrmat.color) +
  remove_x_axis()+
  scale_size_area(max_size = 2) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  guides(
    fill = guide_colourbar(title = "corr", order = 1),
    colour = guide_legend(title = "spearman's p", order = 2),
    size = guide_legend(title = "spearman's r", order = 3),
    linetype = "none")
dev.off()

######11.风险分组间通路特征###########
dir.create('results/11.risktype.pathway')
###########GSEA###########
dir.create('results/11.risktype.pathway/GSEA')
writeMatrix(dat = tcga_tpm_log_T[,tcga.risktype.cli$Samples],outpath = 'results/11.risktype.pathway/tcga_exp.txt',row = T,header = T)
writeMatrix(dat = tcga.risktype.cli[,c('Samples','Risktype')],outpath = 'results/11.risktype.pathway/tcga.risktype.txt',row = F,header = T)
mg_RunGSEA_use(mod = 'exp_group',exp_Path = 'results/11.risktype.pathway/tcga_exp.txt'
               ,sample_group_path = 'results/11.risktype.pathway/tcga.risktype.txt'
               ,outFolder = 'results/11.risktype.pathway/GSEA'
               ,gmt_Path = 'HALLMARK',outLog=F)
tcga_GSEA_PTEN=parseGSEAResult('results/11.risktype.pathway/GSEA/my_analysis.Gsea.1670223236201/')
head(tcga_GSEA_PTEN$EnrichTable)

table(tcga_GSEA_PTEN$EnrichTable$FDR<0.25 & tcga_GSEA_PTEN$EnrichTable$NES<0 )
table(tcga_GSEA_PTEN$EnrichTable$FDR<0.25 & tcga_GSEA_PTEN$EnrichTable$NES>0 )
######high富集
inds1=which(tcga_GSEA_PTEN$EnrichTable$FDR<0.25 & tcga_GSEA_PTEN$EnrichTable$NES>0 & tcga_GSEA_PTEN$EnrichTable$NP<0.05)
tcga_GSEA_PTEN$EnrichTable[inds1,]
writeMatrix(tcga_GSEA_PTEN$EnrichTable[inds1,],'results/11.risktype.pathway/Hig_GSEA.txt')
tcag.plot1<-plot_GSEA_By_nodes(tcga_GSEA_PTEN,indexs=inds1)
tcag.plot1
savePDF('results/11.risktype.pathway/Fig11b.pdf',tcag.plot1,height = 6,width = 8)




load('results/tcga.onco.RData')
writeMatrix(tcga.onco,'results/11.risktype.pathway/tcga.onco.txt')
mg_PlotMutiBoxplot(data = tcga.onco[tcga.risktype.cli$Samples,]
                   ,group = tcga.risktype.cli$Risktype
                   ,test_method = 'wilcox.test'
                   ,legend.pos = 'top'
                   ,ylab = 'score',xangle = 0
                   ,group_cols = risktype.col
                   ,add = 'boxplot')#+coord_flip()
RS_cor_path=cbind(data.frame(Riskscore=tcga.risktype.cli$Riskscore),
                  tcga.onco[tcga.risktype.cli$Samples,])
library(ggcor)
pdf('results/11.risktype.pathway/Fig11a.pdf',height = 6,width = 6)
corrplot(cor(RS_cor_path),
         tl.col = 'black',
         # cl.align.text = 'l',
         tl.srt=45, 
         col=colorRampPalette(c('blue', 'white','red'))(100),
         order = c("original", "AOE", "FPC", "hclust", "alphabet")[1],
         type="upper",)

dev.off()

#######预后模型#################
dir.create('results/13.Nomogram')

#单因素和多因素####
#Age
tcga_cox_datas=crbind2DataFrame(tcga_cox_datas)
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat

#Gender
Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                                data=tcga_cox_datas))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_sig_cox[[8]]),
                                 HR = round(Gender_sig_cox[[7]][,2],3),
                                 lower.95 = round(Gender_sig_cox[[8]][,3],3),
                                 upper.95 = round(Gender_sig_cox[[8]][,4],3),
                                 p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat

#T.stage
T.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~T.stage,
                                 data=tcga_cox_datas))
T.stage_sig_cox_dat <- data.frame(Names=rownames(T.stage_sig_cox[[8]]),
                                  HR = round(T.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(T.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(T.stage_sig_cox[[8]][,4],3),
                                  p.value=round(T.stage_sig_cox[[7]][,5],3))
T.stage_sig_cox_dat

#N.stage
N.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~N.stage,
                                 data=tcga_cox_datas))
N.stage_sig_cox_dat <- data.frame(Names=rownames(N.stage_sig_cox[[8]]),
                                  HR = round(N.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(N.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(N.stage_sig_cox[[8]][,4],3),
                                  p.value=round(N.stage_sig_cox[[7]][,5],3))
N.stage_sig_cox_dat

#M.stage
M.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~M.stage,
                                 data=tcga_cox_datas))
M.stage_sig_cox_dat <- data.frame(Names=rownames(M.stage_sig_cox[[8]]),
                                  HR = round(M.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(M.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(M.stage_sig_cox[[8]][,4],3),
                                  p.value=round(M.stage_sig_cox[[7]][,5],3))
M.stage_sig_cox_dat

#Stage
Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat

#Grade
Grade_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Grade,
                               data=tcga_cox_datas))
Grade_sig_cox_dat <- data.frame(Names=rownames(Grade_sig_cox[[8]]),
                                HR = round(Grade_sig_cox[[7]][,2],3),
                                lower.95 = round(Grade_sig_cox[[8]][,3],3),
                                upper.95 = round(Grade_sig_cox[[8]][,4],3),
                                p.value=round(Grade_sig_cox[[7]][,5],3))
Grade_sig_cox_dat

#riskscore
riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                 data=tcga_cox_datas))
riskscore_sig_cox_dat <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                    HR = round(riskscore_sig_cox[[7]][,2],3),
                                    lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                    upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                    p.value=round(riskscore_sig_cox[[7]][,5],3))
riskscore_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     T.stage_sig_cox_dat,
                     N.stage_sig_cox_dat,
                     M.stage_sig_cox_dat,
                     Stage_sig_cox_dat,
                     Grade_sig_cox_dat,
                     riskscore_sig_cox_dat)
data.sig <- data.frame(Names=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
data.sig$Names
rownames(data.sig) <- c("Age",
                        "Gender",
                        "T.stage",
                        "N.stage",
                        "M.stage",
                        "Stage",
                        "Grade",
                        "RiskScore")
data.sig$Names <- rownames(data.sig)
data.sig
data.sig$p.value=ifelse(data.sig$p.value<0.001,'<0.001',data.sig$p.value)
data.sig.use=data.sig
data.sig.use$HR=paste0(data.sig.use$sig_cox_dat.HR,'(',data.sig.use$sig_cox_dat.lower.95,',',data.sig.use$sig_cox_dat.upper.95,')')
data.sig.use=data.sig.use[,c(1,2,6,3:5)]
data.sig.use
colnames(data.sig.use)[3]=c('Hazard Ratio(95%CI)')
data.sig.use

sig.dat<- read.csv("results/13.Nomogram/data.sig.csv",header = F)
sig.dat
library(forestplot)
pdf('results/13.Nomogram/Fig13a.pdf',height = 5,width = 8,onefile = F)
forestplot(labeltext=as.matrix(sig.dat[,1:3]), #直接显示在图中的文本设置
           mean=sig.dat$V4, #森林图中中间点的位置数据设置
           lower=sig.dat$V5, #森林图的最小值设置
           upper=sig.dat$V6, #森林图的最大值设置
           is.summary=c(T,F,F,F,F,F,F,F,F), #设置重点突出显示的文本（标题等）
           hrzl_lines=list("2"=gpar(col="black",lwd="2"), "10"=gpar(col="black",lwd="2")), #在森林图中添加线条并设置线条的颜色,粗细
           xlab="Hazard Ratio", #设置森林图底部的标签名称
           zero=1, #设置森林图中的标准线,可设置多条
           lwd.zero=2,#设置标准线的粗细
           #clip = c(0.8, 3), #设置森林图的误差线范围,超出部分显示为箭头
           #xticks = c(1,1.5,2,2.5,3),#设置森林图的下标刻度
           txt_gp = fpTxtGp(label = gpar(fontfamily = "serif"),ticks = gpar(cex = 1), xlab = gpar(cex = 1), cex = 1),#设置文本的字体,字体大小
           boxsize = 0.3, #设置平均数点的大小
           #graph.pos = 4, #设置森林图的位置,第四列
           fn.ci_norm = fpDrawDiamondCI, #设置平均数点的类型,有标准、圆等可选
           col=fpColors(line = "black", box="coral1",zero = "cadetblue",), #设置线,点的颜色
           lty.ci = 7, #设置误差线的类型
           lwd.ci = 2, #设置误差线的粗细
           #colgap=unit(0.05,"npc"), #设置列之间的距离
           #graphwidth=unit(0.2,"npc")
) #设置森林图在整个图中的大小
dev.off()




#多因素
#T.stage +  N.stage +M.stage
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age+Gender+Stage+Grade+Riskscore, 
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Names=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
data.muti
rownames(data.muti) <- c("Age",
                         "Gender",
                         #"T.stage",
                         "Stage",
                         "Grade",
                         "RiskScore")
data.muti$Names <- rownames(data.muti)
data.muti
data.muti$p.value=ifelse(data.muti$p.value<0.001,'<0.001',data.muti$p.value)
#mg_forestplot_v2(data.muti,xlog = T,colgap = 15,lineheight = 10)

data.muti.use=data.muti
data.muti.use$HR=paste0(data.muti.use$muti_cox_dat.HR,'(',data.muti.use$muti_cox_dat.lower.95,',',data.muti.use$muti_cox_dat.upper.95,')')
data.muti.use=data.muti.use[,c(1,2,6,3:5)]
data.muti.use
colnames(data.muti.use)[3]=c('Hazard Ratio(95%CI)')
data.muti.use

muti.dat<- read.csv("results/13.Nomogram/data_muti.csv",header = F)
muti.dat
library(forestplot)
pdf('results/13.Nomogram/Fig13b.pdf',height = 5,width = 8,onefile = F)
forestplot(labeltext=as.matrix(muti.dat[,1:3]), #直接显示在图中的文本设置
           mean=muti.dat$V4, #森林图中中间点的位置数据设置
           lower=muti.dat$V5, #森林图的最小值设置
           upper=muti.dat$V6, #森林图的最大值设置
           is.summary=c(T,F,F,F,F,F,F,F,F), #设置重点突出显示的文本（标题等）
           hrzl_lines=list("2"=gpar(col="black",lwd="2"), "7"=gpar(col="black",lwd="2")), #在森林图中添加线条并设置线条的颜色,粗细
           xlab="Hazard Ratio", #设置森林图底部的标签名称
           zero=1, #设置森林图中的标准线,可设置多条
           lwd.zero=2,#设置标准线的粗细
           #clip = c(0.8, 3), #设置森林图的误差线范围,超出部分显示为箭头
           #xticks = c(1,1.5,2,2.5,3),#设置森林图的下标刻度
           txt_gp = fpTxtGp(label = gpar(fontfamily = "serif"),ticks = gpar(cex = 1), xlab = gpar(cex = 1), cex = 1),#设置文本的字体,字体大小
           boxsize = 0.3, #设置平均数点的大小
           #graph.pos = 4, #设置森林图的位置,第四列
           fn.ci_norm = fpDrawDiamondCI, #设置平均数点的类型,有标准、圆等可选
           col=fpColors(line = "black", box="coral1",zero = "cadetblue",), #设置线,点的颜色
           lty.ci = 7, #设置误差线的类型
           lwd.ci = 2, #设置误差线的粗细
           #colgap=unit(0.05,"npc"), #设置列之间的距离
           #graphwidth=unit(0.2,"npc")
) #设置森林图在整个图中的大小
dev.off()

#列线图####
pdf('results/13.Nomogram/Fig13c.pdf', width = 12, height = 10)
nom.plot=mg_nomogram(data.frame(RiskScore=tcga_cox_datas$Riskscore,
                                Age=tcga_cox_datas$Age,
                                Stage=tcga_cox_datas$Stage),
                     os = tcga_cox_datas$OS.time,
                     status = tcga_cox_datas$OS,
                     mks = c(1,3,5)
)
dev.off()
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))

################AUC
library("timeROC")
tcga_cox_auc=tcga_cox_datas
tcga_cox_auc$Riskscore=as.numeric(tcga_cox_auc$Riskscore)
tcga_cox_auc$T.stage=as.numeric(as.factor(tcga_cox_auc$T.stage))
tcga_cox_auc$N.stage=as.numeric(as.factor(tcga_cox_auc$N.stage))
tcga_cox_auc$M.stage=as.numeric(as.factor(tcga_cox_auc$M.stage))
tcga_cox_auc$Stage=as.numeric(as.factor(tcga_cox_auc$Stage))
tcga_cox_auc$Grade=as.numeric(as.factor(tcga_cox_auc$Grade))
tcga_cox_auc$Age=as.numeric(as.factor(tcga_cox_auc$Age))
tcga_cox_auc$Gender=as.numeric(as.factor(tcga_cox_auc$Gender))

head(tcga_cox_auc)
ROC.DSST.Age=timeROC(T=tcga_cox_auc$OS.time/365,
                     delta=tcga_cox_auc$OS,
                     marker=tcga_cox_auc$Age,
                     cause=1,weighting="marginal",
                     times=c(1,2,3,4,5),
                     iid=TRUE)
ROC.DSST.Age$AUC

ROC.DSST.Gender=timeROC(T=tcga_cox_auc$OS.time/365,
                        delta=tcga_cox_auc$OS,
                        marker=tcga_cox_auc$Gender,
                        cause=1,weighting="marginal",
                        times=c(1,2,3,4,5),
                        iid=TRUE)
ROC.DSST.Gender$AUC

ROC.DSST.T=timeROC(T=tcga_cox_auc$OS.time/365,
                   delta=tcga_cox_auc$OS,
                   marker=tcga_cox_auc$T.stage,
                   cause=1,weighting="marginal",
                   times=c(1,2,3,4,5),
                   iid=TRUE)
ROC.DSST.T$AUC

ROC.DSST.N=timeROC(T=tcga_cox_auc$OS.time/365,
                   delta=tcga_cox_auc$OS,
                   marker=tcga_cox_auc$N.stage,
                   cause=1,weighting="marginal",
                   times=c(1,2,3,4,5),
                   iid=TRUE)
ROC.DSST.N$AUC

ROC.DSST.M=timeROC(T=tcga_cox_auc$OS.time/365,
                   delta=tcga_cox_auc$OS,
                   marker=tcga_cox_auc$M.stage,
                   cause=1,weighting="marginal",
                   times=c(1,2,3,4,5),
                   iid=TRUE)
ROC.DSST.M$AUC

ROC.DSST.stage=timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$Stage,
                       cause=1,weighting="marginal",
                       times=c(1,2,3,4,5),
                       iid=TRUE)
ROC.DSST.stage$AUC

ROC.DSST.Grade=timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$Grade,
                       cause=1,weighting="marginal",
                       times=c(1,2,3,4,5),
                       iid=TRUE)
ROC.DSST.Grade$AUC

ROC.DSST.Risk=timeROC(T=tcga_cox_auc$OS.time/365,
                      delta=tcga_cox_auc$OS,
                      marker=tcga_cox_auc$Riskscore,
                      cause=1,weighting="marginal",
                      times=c(1,2,3,4,5),
                      iid=TRUE)
ROC.DSST.Risk$AUC

ROC.DSST.Nomo<-timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$Riskscore,
                       other_markers=as.matrix(tcga_cox_auc[,c("Stage")]),
                       cause=1,
                       weighting="cox",
                       times=c(1,2,3,4,5),
                       iid=F)
ROC.DSST.Nomo$AUC

pdf('results/13.Nomogram/Fig13e.pdf',height = 5,width = 12)
par(mfrow=c(1,3))
plot(ROC.DSST.Age,time=1, col=pal_d3()(10)[6], lwd=2, title = "",)   #time是时间点，col是线条颜色
plot(ROC.DSST.Gender,time=1, col=pal_d3()(10)[5], add=TRUE, lwd=2)
plot(ROC.DSST.stage,time=1, col=pal_d3()(10)[4], add=TRUE, lwd=2) 
plot(ROC.DSST.Grade,time=1, col=pal_d3()(10)[2], add=TRUE, lwd=2) 
plot(ROC.DSST.Risk,time=1, col=pal_d3()(10)[2], add=TRUE, lwd=2) 
plot(ROC.DSST.Nomo,time=1, col=pal_d3()(10)[1], add=TRUE, lwd=2) 
legend("bottomright",
       c(paste0("Age AUC at 1 year: ",round(ROC.DSST.Age[["AUC"]][1],2)), 
         paste0("Gender AUC at 1 year: ",round(ROC.DSST.Gender[["AUC"]][1],2)),
         paste0("Stage AUC at 1 year: ",round(ROC.DSST.stage[["AUC"]][1],2)), 
         paste0("Grade AUC at 1 year: ",round(ROC.DSST.Grade[["AUC"]][1],2)), 
         paste0("RiskScore AUC at 1 year: ",round(ROC.DSST.Risk[["AUC"]][1],2)), 
         paste0("Nomogram AUC at 1 year: ",round(ROC.DSST.Nomo[["AUC"]][1],2))
       ),
       col=pal_d3()(10)[c(6,5,4,3,2,1)],
       lty=1, lwd=2,bty = "n",cex = 1)   


#pdf('results/08.nomogram/Fig8f.pdf',height = 5,width = 5)
plot(ROC.DSST.Age,time=3, col=pal_d3()(10)[6], lwd=2, title = "",)   #time是时间点，col是线条颜色
plot(ROC.DSST.Gender,time=3, col=pal_d3()(10)[5], add=TRUE, lwd=2)
plot(ROC.DSST.stage,time=3, col=pal_d3()(10)[4], add=TRUE, lwd=2) 
plot(ROC.DSST.Grade,time=3, col=pal_d3()(10)[3], add=TRUE, lwd=2) 
plot(ROC.DSST.Risk,time=3, col=pal_d3()(10)[2], add=TRUE, lwd=2) 
plot(ROC.DSST.Nomo,time=3, col=pal_d3()(10)[1], add=TRUE, lwd=2) 
legend("bottomright",
       c(paste0("Age AUC at 3 year: ",round(ROC.DSST.Age[["AUC"]][3],2)), 
         paste0("Gender AUC at 3 year: ",round(ROC.DSST.Gender[["AUC"]][3],2)),
         paste0("Stage AUC at 3 year: ",round(ROC.DSST.stage[["AUC"]][3],2)), 
         paste0("Grade AUC at 3 year: ",round(ROC.DSST.Grade[["AUC"]][3],2)), 
         paste0("RiskScore AUC at 3 year: ",round(ROC.DSST.Risk[["AUC"]][3],2)), 
         paste0("Nomogram AUC at 3 year: ",round(ROC.DSST.Nomo[["AUC"]][3],2))
       ),
       col=pal_d3()(10)[c(6,5,4,3,2,1)],
       lty=1, lwd=2,bty = "n",cex = 1)    


#pdf('results/08.nomogram/Fig8g.pdf',height = 5,width = 5)
plot(ROC.DSST.Age,time=5, col=pal_d3()(10)[6], lwd=2, title = "",)   #time是时间点，col是线条颜色
plot(ROC.DSST.Gender,time=5, col=pal_d3()(10)[5], add=TRUE, lwd=2) 
plot(ROC.DSST.stage,time=5, col=pal_d3()(10)[4], add=TRUE, lwd=2) 
plot(ROC.DSST.Grade,time=5, col=pal_d3()(10)[3], add=TRUE, lwd=2) 
plot(ROC.DSST.Risk,time=5, col=pal_d3()(10)[2], add=TRUE, lwd=2) 
plot(ROC.DSST.Nomo,time=5, col=pal_d3()(10)[1], add=TRUE, lwd=2) 
legend("bottomright",
       c(paste0("Age AUC at 5 year: ",round(ROC.DSST.Age[["AUC"]][5],2)), 
         paste0("Gender AUC at 5 year: ",round(ROC.DSST.Gender[["AUC"]][5],2)), 
         paste0("Stage AUC at 5 year: ",round(ROC.DSST.stage[["AUC"]][5],2)), 
         paste0("Grade AUC at 5 year: ",round(ROC.DSST.Grade[["AUC"]][5],2)), 
         paste0("RiskScore AUC at 5 year: ",round(ROC.DSST.Risk[["AUC"]][5],2)), 
         paste0("Nomogram AUC at 5 year: ",round(ROC.DSST.Nomo[["AUC"]][5],2))),
       col=pal_d3()(10)[c(6,5,4,3,2,1)],
       lty=1, lwd=2,bty = "n",cex = 1)   
dev.off()

save.image(file = 'BLCA_Fatty_Acid_metabolic.RData')

