ssGSEAScore_by_genes=function(gene.exp,genes){
  #library('GSVA')
  #library(GSEABase)
  #all.list=list()
  gs=GSEABase::GeneSet(setName='GeneSet', setIdentifier=paste0("101"),geneIds=unique(genes),GSEABase::SymbolIdentifier()) 
  
  gsc <- GSEABase::GeneSetCollection(list(gs))
  fl <- tempfile()
  GSEABase::toGmt(gsc, fl)
  cgeneset=GSEABase::getGmt(fl)
  ssGSEA.geneset <- GSVA::gsva(as.matrix(gene.exp), cgeneset,method='ssgsea',
                               min.sz=1, max.sz=5000, verbose=TRUE)
  #detach('package:GSVA')
  #detach('package:GSEABase')
  #row.names(ssGSEA.geneset)
  return(ssGSEA.geneset)
}
ssGSEAScore_by_muti_group_genes=function(gene.exp,genelist){
  #library('GSVA')
  #library(GSEABase)
  all.list=list()
  for(i in 1:length(genelist)){
    gs=GSEABase::GeneSet(setName=names(genelist)[i], setIdentifier=paste0("101"),geneIds=unique(genelist[[i]]),GSEABase::SymbolIdentifier()) 
    #gs
    all.list=c(all.list,list(gs))
  }
  #length(all.list)
  gsc <- GSEABase::GeneSetCollection(all.list)
  fl <- tempfile()
  GSEABase::toGmt(gsc, fl)
  cgeneset=GSEABase::getGmt(fl)
  ssGSEA.geneset <- GSVA::gsva(as.matrix(gene.exp), cgeneset,method='ssgsea',
                               min.sz=1, max.sz=5000, verbose=TRUE)
  #detach('package:GSVA')
  #detach('package:GSEABase')
  #row.names(ssGSEA.geneset)
  return(ssGSEA.geneset)
}
cox_batch=function(dat,time,event){
  t.inds=which(!is.na(time)&!is.na(event))
  dat1=dat[,t.inds]
  os=time[t.inds]
  ev=event[t.inds]
  
  ct=sum(ev%in%c(0,1))
  if(ct!=length(ev)){
    print('event must be 0(alive) or 1(dead)')
    return(NULL)
  }
  
  res=t(apply(dat1, 1, function(x){
    ep=as.numeric(as.character(x))
    ind2=which(!is.na(ep))
    print(length(ind2))
    if(length(ind2)>1){
      os1=os[ind2]
      ev1=ev[ind2]
      ep1=ep[ind2]
      return(coxRun(data.frame(os1,ev1,ep1)))
    }else{
      return(c(NA,NA,NA,NA))
    }
  }))
  colnames(res)=c('p.value','HR','Low 95%CI','High 95%CI')
  row.names(res)=row.names(dat1)
  return(as.data.frame(res))
}

ggplotKMCox=function(dat,title='Groups',labs=NULL,add_text=NULL,palette='npg',show_confint=T,show_median_text=T){
  library(ggplot2)
  library(survival)
  colnames(dat)=c('time','status','groups')
  #sdf<-survdiff(Surv(time,status) ~ groups,data=dat)
  #print((sdf))
  #summary(sdf)
  #p<-pchisq(sdf$chisq,length(sdf$n)-1,lower.tail=FALSE)
  sf<-survival::survfit(Surv(time,status) ~ groups,data=dat)
  surv=survminer::ggsurvplot(sf, data = dat, palette = palette, #jco palette 
                             pval = TRUE,surv.median.line='hv'
                             #,conf.int = T
                             ,conf.int.style ='step'
                             , pval.coord=c(0, 0.2), #Add p-value 
                             risk.table = TRUE, 
                             legend.title = title
                             ,legend.labs = labs
                             ,conf.int=show_confint
  )
  p1=surv$plot+theme_bw()+theme(axis.text.y=element_text(family="Times",face="plain")
                                ,axis.text.x=element_blank()
                                ,axis.title.x=element_blank()
                                ,plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches")
                                #,axis.title.y=element_blank()
                                ,legend.position=c(1,1), legend.justification=c(1,1)
                                ,legend.background = element_rect(fill = NA, colour = NA)
                                ,legend.title = element_text(family="Times",face="plain")
                                ,legend.text = element_text(family="Times",face="plain"))
  if(show_median_text){
    median_labels=c()
    for(st in unique(surv$data.survplot$strata)){
      st1=surv$data.survplot[which(surv$data.survplot$strata==st),]
      x_m=-1
      if(min(st1$surv)<0.5){
        inds=which(st1$surv==0.5)
        if(length(inds)>0){
          x_m=st1$time[inds[1]]
        }else{
          x_m=max(st1$time[st1$surv>=0.5])
        }
      }
      if(x_m>0){
        median_labels=c(median_labels,round(x_m,1))
      }
    }
    if(length(median_labels)>0){
      txt_median=surv$data.survplot[1:length(median_labels),]
      txt_median[,5]=rep(0.5,length(median_labels))
      txt_median[,1]=median_labels
      txt_median$Text=median_labels
      p1=p1+geom_text(data=txt_median,aes(x=time, y=surv, label=Text),color="black",hjust =1,angle=90,alpha=0.5,vjust=0,nudge_y = -0.01)
    }
  }
  #p1=p1+text()
  #tms=data.frame(Group=tms.gp,value=tms.tps,Attribute=rep(data_m[1,1],length(tms.gp))
  #               ,ymax=rep(max(ylim),length(tms.gp)))
  #p4=p4+geom_text(data=tms,aes(x=Group, y=ymax, label=value),color="black")
  if(is.null(add_text)){
    gp=unique(dat[,3])
    vls=1:length(gp)
    gvls=vls[match(dat[,3],gp)]
    g.cox=coxRun(data.frame(dat[,1],dat[,2],gvls))
    add_text=paste0('HR=',round(g.cox[2],2)
                    ,', 95%CI(',round(g.cox[3],2)
                    ,', ',round(g.cox[4],2),')'
                    ,'\nHR By:',paste0(gp,collapse = '<'))
    
  }
  text.tb=surv$data.survplot[1,]
  text.tb[1,1]=0
  text.tb[1,5]=0
  text.tb$Text=add_text
  p1=p1+geom_text(data=text.tb,aes(x=time, y=surv, label=Text),color="black",hjust =0)
  
  p2=surv$table+theme_bw()+theme(axis.text.y=element_text(family="Times",face="plain")
                                 #,axis.text.x=element_blank()
                                 #,axis.title.x=element_blank()
                                 #,axis.title.y=element_blank()
                                 ,plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches")
                                 ,plot.title=element_blank()
                                 ,legend.position=c(1,1), legend.justification=c(1,1)
                                 #,legend.background = element_rect(fill = NA, colour = NA)
                                 ,legend.title = element_text(family="Times",face="plain")
                                 ,legend.text = element_text(family="Times",face="plain"))
  
  g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(1,0.3),align = "v")
  return(g2)
}

mg_violin_1=function(data,xangle=0,ylab='value',xlab='',leg.title='Group',test_method='anova',cmp_test_method='t.test',legend.pos='r',melt=F,jitter=T,ylim=NULL,show_compare=NULL,point_size=NULL,group_col){
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

mg_PlotMutiBoxplot=function(data,group,group_cols='jco'
                            ,test_method=c('t.test','wilcox.test','paired_t.test','paired_wilcox.test','anova','kruskal.test')[1]
                            ,order=NULL,size=1,fill=F,outlier.shape=NA,yscale=c('none','log2','log10')[1]
                            ,xangle=45,ylab='Value',xlab='',box_span=0.7
                            ,orientation = c("vertical", "horizontal", "reverse")[1]
                            ,legend.pos=NULL,melt=F,ylim=NULL,binwidth=0.05
                            ,add=c("none", "dotplot", "jitter", "boxplot", "point", "mean"
                                   , "mean_se", "mean_sd", "mean_ci", "mean_range", "median"
                                   , "median_iqr", "median_mad", "median_range")[3]){
  paired=FALSE
  if(test_method=='paired_t.test'|test_method=='paired_wilcox.test'){
    test_method=gsub('paired_','',test_method)
    paired=TRUE
  }
  print(class(data))
  if(add=='jitter'){
    fill=F
  }
  
  library(ggplot2)
  if(!melt){
    #print(class(data))
    if(class(data)=='numeric'|class(data)=='integer'){
      data=as.numeric(data)
      vd1.sbs=data.frame(group,rep('Tag',length(group)),data)
      #print(vd1.sbs)
    }else{
      data=as.data.frame(data)
      data$ID=group
      vd1.sbs <- reshape2::melt(data, id.vars=c("ID"))
    }
    colnames(vd1.sbs)=c('category','type','Score')
    Data=vd1.sbs
  }else{
    vd1.sbs=data
    colnames(vd1.sbs)=c('category','type','Score')
    Data=vd1.sbs
  }
  
  #vd1.sbs[,2]=paste0('C',as.numeric(as.character(vd1.sbs[,2])))
  if(is.null(order)){
    order=unique(vd1.sbs[,2])
  }
  
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos=c(0,0)
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
  }else if(legend.pos=='top'){
    pos='top'
  }else if(legend.pos=='buttom'){
    pos='buttom'
  }else{
    pos='right'
  }
  print(pos)
  if(fill){
    p <- ggpubr::ggboxplot(vd1.sbs, x="type", y="Score", fill = "category", yscale = yscale
                           ,palette = group_cols,width = box_span,size = size,order = order,outlier.shape=outlier.shape
                           ,orientation=orientation,add=add,add.params = list(binwidth=binwidth)
                           ,short.panel.labs = T)#按dose进
  }else{
    p <- ggpubr::ggboxplot(vd1.sbs, x="type", y="Score", color = "category", yscale = yscale
                           ,palette = group_cols,width = box_span,size = size,order = order,outlier.shape=outlier.shape
                           ,orientation=orientation,add=add,add.params = list(binwidth=binwidth)
                           ,short.panel.labs = T)#按dose进
  }
  
  p=p+ggpubr::stat_compare_means(aes(group=category), label = "p.signif", method = test_method,paired=paired
                                 #,label.y = max(vd1.sbs[,1])
  )
  #p=p+ylab(ylab)+xlab(xlab)
  #p=p+theme(axis.text.x = element_text(angle = xangle, hjust = 1))
  p=p+theme(axis.text.x=tx, #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
            axis.text.y=element_text(family="Times",face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
            axis.title.y=element_text(family="Times",face="plain"), #设置y轴标题的字体属性
            #panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
            legend.text=element_text(face="plain", family="Times", colour="black"  #设置图例的子标题的字体属性
            ),
            legend.title=element_text(face="plain", family="Times", colour="black" #设置图例的总标题的字体属性
            ),
            legend.justification=pos, legend.position=pos
            ,legend.background = element_rect(fill = NA, colour = NA)
            #,panel.grid.major = element_blank(),   #不显示网格线
            #panel.grid.minor = element_blank()
  )+ylab(ylab)+xlab(xlab) #设置x轴和y轴的标题
  if(!is.null(ylim)){
    p=p+ylim(ylim)
  }
  return(p)
}

mg_limma_DEG=function(exp,group,ulab,dlab){
  library(limma)
  ind1=which(group==ulab)
  ind2=which(group==dlab)
  
  sml <- c(rep('G1',length(ind1)),rep('G0',length(ind2)))    # set group names
  eset=exp[,c(ind1,ind2)]
  fl <- as.factor(sml)
  
  design <- model.matrix(~fl+0)
  colnames(design) <- levels(fl)
  cont.matrix<-makeContrasts(contrasts='G1-G0',levels=design)
  #print(head(eset))
  fit<-lmFit (eset,design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  #print(sml)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(eset))
  regulated=ifelse(tT$logFC>0,'Up','Down')
  lfcs=c(log2(1.2),log2(1.3),log2(1.5),1)
  all.deg.cnt=cbind()
  for(lfc in lfcs){
    deg1=regulated[which(abs(tT$logFC)>lfc&tT$P.Value<0.05)]
    deg2=regulated[which(abs(tT$logFC)>lfc&tT$P.Value<0.01)]
    deg3=regulated[which(abs(tT$logFC)>lfc&tT$adj.P.Val<0.05)]
    deg4=regulated[which(abs(tT$logFC)>lfc&tT$adj.P.Val<0.01)]
    all.deg.cnt=cbind(all.deg.cnt,c(paste0(sum(deg1=='Up'),'|',sum(deg1=='Down'))
                                    ,paste0(sum(deg2=='Up'),'|',sum(deg2=='Down'))
                                    ,paste0(sum(deg3=='Up'),'|',sum(deg3=='Down'))
                                    ,paste0(sum(deg4=='Up'),'|',sum(deg4=='Down'))))
  }
  row.names(all.deg.cnt)=c('p<0.05','p<0.01','FDR<0.05','FDR<0.01')
  colnames(all.deg.cnt)=paste0(c('1.2','1.3','1.5','2'),'-fold')
  return(list(Exp=eset,Group=group[c(ind1,ind2)],DEG=tT,Summary=all.deg.cnt))
}

mg_volcano=function(logfc,pvalue,symbol=NULL,cutFC=1,cutPvalue=0.05
                    ,showText=NULL
                    ,colors=c(mg_colors[2],'grey',mg_colors[1])
                    ,xlim=NULL,ylim=NULL
                    ,legend.pos='right'
                    ,ylab='-log10(FDR)',leg='State',xlab='log2(FoldChange)'){
  library(ggplot2)
  pos=c(0,0)
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
  }else{
    pos='right'
  }
  cange=rep('None',length(logfc))
  cange[which(logfc>cutFC&pvalue<cutPvalue)]='Up'
  cange[which(logfc< -cutFC&pvalue<cutPvalue)]='Down'
  if(is.null(symbol)){
    symbol=rep('',length(logfc))
    showText=NULL
  }
  vo.input=data.frame(logFC=logfc,FDR=pvalue,change=cange,SYMBOL=symbol)
  #print(head(vo.input))
  p1 <- ggplot(data = vo.input, 
               aes(x = logFC, 
                   y = -log10(FDR)))
  p1=p1+geom_point(alpha=0.4, size=3.5, aes(color=change))  
  p1=p1+scale_color_manual(values=colors,limits = c("Down",'None', "Up"),name=leg) 
  p1=p1+geom_vline(xintercept=c(-cutFC,cutFC),lty=4,col="black",lwd=0.8)  
  p1=p1+geom_hline(yintercept = -log10(cutPvalue),lty=4,col="black",lwd=0.8)  
  p1=p1+ylab(ylab)+xlab(xlab)
  p1=p1+theme_bw()
  p1=p1+theme(
    axis.text.y=element_text(family="Times",face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
    axis.title.y=element_text(family="Times",face="plain"), #设置y轴标题的字体属性
    legend.text=element_text(face="plain", family="Times", colour="black"  #设置图例的子标题的字体属性
    ),
    legend.title=element_text(face="plain", family="Times", colour="black" #设置图例的总标题的字体属性
    ),
    legend.justification=pos, legend.position=pos
    ,legend.background = element_rect(fill = NA, colour = NA)
  )
  if(is.null(showText)|is.null(symbol)){
    showText=c()
  }
  
  if(length(showText)>0){
    for_label <-vo.input[match(intersect(showText,vo.input$SYMBOL),vo.input$SYMBOL),]
    p1=p1+geom_point(size = 3, shape = 1, data = for_label)+ggrepel::geom_label_repel(
      aes(label = SYMBOL),
      data = for_label,
      color="black"
    )
  }
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(!is.null(xlim)){
    p1=p1+xlim(xlim)
  }
  
  return(p1)
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



getC_index=function(riskscore,os,status){
  inds=which(!is.na(riskscore)&!is.na(os)&!is.na(status))
  riskscore=riskscore[inds]
  os=os[inds]
  status=status[inds]
  c1=survcomp::concordance.index(x=riskscore, surv.time=os, surv.event=status,
                                 method="noether")
  #c2=concordance.index(x=riskscore[order(rnorm(length(riskscore)))], surv.time=os, surv.event=status,
  #                     method="noether")
  #p=min(cindex.comp(c1, c2)$p.value,cindex.comp(c2, c1)$p.value)  
  return(c1)
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



mg_nomogram_buti=function(cox_model,cut.time,title='Nomogram'){
  library(regplot)
  regplot(cox_model#对观测2的六个指标在列线图上进行计分展示
          ,observation=pbc[2,] #也可以不展示
          #预测3年和5年的死亡风险，此处单位是day
          ,title=title
          ,failtime = cut.time
          ,prfail = TRUE #cox回归中需要TRUE
          ,showP = T #是否展示统计学差异
          ,droplines = F#观测2示例计分是否画线
          #,colors = mg_colors[1:3] #用前面自己定义的颜色
          #,rank="decreasing") #根据统计学差异的显著性进行变量的排序
          #,interval="confidence"
          #,rank="decreasing"
          #,clickable=T
          ,points=TRUE)
  
}

crbind2DataFrame=function(dat,full=F){
  print(class(dat))
  if(class(dat)=='table'){
    if(!is.na(ncol(dat))){
      dat=apply(dat,2,function(x){
        return(x)
      })
    }
  }
  if(class(dat)!='data.frame'){
    dat1=as.data.frame(as.matrix(dat))
  }else{
    dat1=dat
  }
  dat1.class=apply(dat1, 2, class)
  #which(dat1.class!='numeric')
  #print(head(dat1))
  for(i in which(dat1.class!='numeric')){
    dat1[,i]=as.character(dat1[,i])
    if(full){
      dat1[,i]=as.numeric(dat1[,i])
    }else{
      dt=dat1[which(gsub(' ','',dat1[,i])!=''&!is.na(dat1[,i])),i]
      dt=dt[which(dt!='Inf'&dt!='NaN'&dt!='NA')]
      #dt[which(is.na(as.numeric(dt)))]
      if(sum(is.na(as.numeric(dt)))<length(dt)*0.1){
        #print(dat1[,i])
        dat1[,i]=as.numeric(dat1[,i])
      }
    }
  }
  return(dat1)  
}
plotMutiBar=function(dat,ist=F,margin=T,xlb='',ylb='',lineCol='black',lineW=0.5,legTitle='Group',showValue=F,showLine=T,xangle=0,isAuto=T){
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
  stable.p=stable.p+ggtitle('-log10(anova p value)')
  if(isAuto){
    g1=ggpubr::ggarrange(stable.p,pg, ncol = 1, nrow = 2,heights = c(0.5,1),align = "hv")
    return(g1)
  }else{
    return(list(Bar=pg,Table=stable.p))
  }
}
mg_merge_plot=function(...,ncol = NULL, nrow = NULL,
                       labels = NULL, label.x = 0, label.y = 1, hjust = -0.5,
                       vjust = 1.5, font.label = list(size = 14, color = "black", face ="bold"
                                                      , family = NULL), align = c("none", "h", "v", "hv")
                       ,widths = 1, heights = 1, legend = NULL, common.legend = FALSE){
  ml=list(...)
  if(length(ml)==1){
    if(is.list(ml[[1]])){
      ml=ml[[1]]
    }
  }
  gal=ggpubr::ggarrange(plotlist = ml, ncol = ncol, nrow = nrow
                        ,labels = labels,label.x = label.x, label.y = label.y
                        , hjust = hjust,
                        vjust = vjust, font.label = font.label, align = align
                        ,widths = widths, heights = heights, legend = legend
                        , common.legend = common.legend
  )
  return(gal)
}


coxRun=function(dat){
  library(survival)
  colnames(dat)=c('time','status','AS')  
  dat=dat[which(!is.na(dat[,1])&!is.na(dat[,3])&!is.na(dat[,2])),]
  #print(nrow(dat))
  if(nrow(dat)<10){
    print(paste0('Sample Num is small:',nrow(dat)))
    return(c(NA,NA,NA,NA))
  }
  #if(quantile(dat[,3])['25%']==quantile(dat[,3])['50%']) return(c(NA,NA,NA,NA))
  fmla <- as.formula("Surv(time, status) ~AS")
  if(table(dat[,2])[1]>1&table(dat[,2])[2]>1){
    cox <- survival::coxph(fmla, data = dat)
    re=c(summary(cox)[[7]][5],summary(cox)[[7]][2],summary(cox)[[8]][3],summary(cox)[[8]][4])
    return(re)
  }else{
    return(c(NA,NA,NA,NA))
  }
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
coxFun <- function(dat){
  library(survival)
  colnames(dat)=c('time','status','gene')
  fmla=as.formula("Surv(time,status)~gene")
  cox=coxph(fmla,data=dat)
  p=summary(cox)[[7]][5]
  result=c(p,summary(cox)[[8]][1],summary(cox)[[8]][3],summary(cox)[[8]][4])
  return(result)
}

enrichmentORA=function(genes,minNum=10,maxNum = 500, sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.5,
                       topThr = 10, reportNum = 20, perNum = 1000,mp_dbs=c('pathway_KEGG','pathway_Wikipathway','pathway_Reactome','geneontology_Biological_Process','geneontology_Cellular_Component','geneontology_Molecular_Function'
                                                                           ,'disease_Disgenet','disease_OMIM','drug_DrugBank','phenotype_Human_Phenotype_Ontology')){
  #library('WebGestaltR')
  refFile <- system.file("extdata", "referenceGenes.txt", package="WebGestaltR")
  refGenes=read.csv(refFile,sep = '\t',stringsAsFactors = F,header = F)[,1]
  cm.genes=intersect(genes,refGenes)
  if(length(cm.genes)<length(genes)*0.1){
    library(hgu133plus2.db)
    columns(hgu133plus2.db)
    tryCatch({
      eid=AnnotationDbi::select(hgu133plus2.db,keys = genes,keytype = 'ENTREZID',columns = c('SYMBOL'))
      genes1=unique(eid$SYMBOL)
      cm.genes1=intersect(genes1,refGenes)
    },error=function(e){cm.genes1=c();
      genes1=c()
    })
    tryCatch({
      eid=AnnotationDbi::select(hgu133plus2.db,keys = genes,keytype = 'ENSEMBL',columns = c('SYMBOL'))
      genes2=unique(eid$SYMBOL)
      cm.genes2=intersect(genes2,refGenes)
    },error=function(e){
      cm.genes2=c();
      genes2=c()
    })
    if(length(genes1)>length(genes2)){
      genes3=genes1
    }else{
      genes3=genes2      
    }
    cm.genes1=intersect(genes3,refGenes)
    if(length(cm.genes)<length(cm.genes1)){
      genes=genes3
      cm.genes=cm.genes1
    }
  }
  
  #mp_dbs=c('pathway_KEGG','pathway_Wikipathway','pathway_Reactome','geneontology_Biological_Process','geneontology_Cellular_Component','geneontology_Molecular_Function'
  #,'disease_Disgenet','disease_OMIM','drug_DrugBank','phenotype_Human_Phenotype_Ontology')
  if(length(cm.genes)>0){  
  
  fld=paste0(tempfile())
  dir.create(fld) 
  tmp=paste0(fld,'/gene.txt')
  write.table(data.frame(genes),file=tmp,row.names = F,col.names = F,sep = '\t')
  result=rbind()
  for(db in mp_dbs){
    print(paste0('starting ',db))
    tryCatch({
    enrichResult <- WebGestaltR::WebGestaltR(enrichMethod="ORA", organism="hsapiens"
                                ,enrichDatabase=db,minNum=minNum,maxNum = maxNum, sigMethod = sigMethod, fdrMethod = fdrMethod, fdrThr = fdrThr,
                                topThr = topThr, reportNum = reportNum, perNum = perNum
                                , interestGeneFile=tmp,
                                interestGeneType="genesymbol", 
                                referenceGeneFile=refFile,
                                referenceGeneType="genesymbol", isOutput=F,
                                outputDirectory=fld, projectName=NULL)
    if(!is.null(enrichResult)){
      enrichResult$DB=rep(db,nrow(enrichResult))
      result=rbind(result,enrichResult)
    }
    },error=function(e){cat("Error:",'enrichment fail:',db,"\n");})
  }
  unlink(fld,recursive = T)
  #detach('package:WebGestaltR')
  return(result)
  }else{
    #detach('package:WebGestaltR')
    return(data.frame())
  }
  
}

