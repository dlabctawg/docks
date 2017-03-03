#utilities
lintran<-function(x,s1=c(0,1),s2=c(0,1)) {a=diff(s2)/diff(s1);b=s2[1]-a*s1[1];return(a*x+b)}
ec<-function(x) strsplit(x,',')[[1]]
tscore<-function(x) (x-mean(x))/sd(x)

#' google books ngram 2 token time series
#'
#' @param pipe hot stuff
#' @param cfso
#' @param query A character vector of terms or a long format database where the last column is the term and leading columns are factors. Factors will be returned on the result.
#' @param ys
#' @param slug
#' @param sm
#' @param ye
#'
#' @return
#' @export
#'
#' @examples
gbng2tts.f<-function(
  query=c('chicken','egg')
  ,dump=paste0('d',.Platform$file.sep,'d')
  ,pipe=paste0('d',.Platform$file.sep,'p')
  ,ys=1900
  ,ye=2000
  ,corp='eng_us_2012'
  ,cfso=T
  ,slug='auto'
  ,sm=0
){
  require(magrittr)
  require(data.table)
  require(ngramr)
  if(!file.exists(dump)) stop(paste('Data dump directory',dump,'does not exist.'))
  if(!file.exists(pipe)) stop(paste('Data pipe directory',pipe,'does not exist.'))
  query<-data.table(query)
  pfs<-.Platform$file.sep
  n<-dim(query) %>% as.list()
  names(n)<-c('r','c')
  if(slug=='auto') {
    if(n$c==1) {
      slug<-query[,n$c,with=F] %>% gsub(pattern='[^A-Za-z ]',replacement='') %>% strsplit(split=' +')
    } else {
      slug<-query[,!n$c,with=F]
    }
    slug<-unlist(slug) %>% substr(0,3) %>% unique() %>% paste0(collapse='-')
  }
  fso<-sub(paste0(pfs,'+'),pfs,paste(pipe,pfs,'gbng',ys,ye,'-',slug,'.RData',sep=''))
  if(cfso) {
    if(!file.exists(fso)) {cat('\nNo saved ngramr output. Set cfso to FALSE to submit new query, and back to TRUE to use saved output.');return(NULL)}
    cat('Loading \"',fso,'\"\n',sep='')
    load(fso)
  } else {
    qu<-unique(query[[n$c]])
    lq<-length(qu)
    b<-ceiling(lq/12) # 12 term search limit
    if(b>1) {p<-split(qu,f=cut(1:lq,b,include.lowest=T))} else {p<-list(qu)}
    cat('Querying Google Ngrams in',length(p),ifelse(length(p)==1,'batch.\n','batches.\n'))
    tts<-lapply(p,function(x) ngram(phrases=x,corpus = corp,year_start=ys,year_end=ye,smoothing=sm,count=T))
    if(file.exists(dump)) {
      cat('Saving:\n')
      lapply(tts,function(x) {
        dfso<-sub(paste0(pfs,'+'),pfs,paste0(dump,pfs,'gbng',paste(range(x$Year),collapse=''),'-',paste0(unique(x$Phrase),collapse='-'),'.txt',sep='',collapse=''))
        cat('\"',dfso,'\"\n',sep='')
        write.table(x,file = dfso,sep='\t',quote = F,na = '',row.names = F,col.names = T)
      })
    }
    tts<-rbindlist(lapply(tts,data.table))
    tts<-merge(x=tts,y=query,by.x='Phrase',by.y=names(query)[n$c])
    setattr(tts,'date.queried',Sys.Date())
    if(!file.exists(dump)) cat('Saving:\n')
    cat('\"',fso,'\"\n',sep='')
    save(tts,file=fso)
  }
  tts
}

tts2grgr.f<-function(
  gbng2tts
  ,order=unique(gbng2tts$Phrase)
  ,mxlg=5
){
  require(data.table)
  require(forecast)
  require(lmtest)
  require(ggplot2)
  if(missing(order)) {d<-copy(gbng2tts)} else {d<-gbng2tts[,Phrase:=factor(Phrase,levels=order)]}
  setkey(d,Phrase,Year)
  lv<-c('Identity','First Difference')
  d<-rbindlist(list(
    d[,list(
      Year=Year[-1]
      ,Frequency=Frequency[-1] %>% lintran(s1=range(Frequency[-1]))
      ,Diffs=factor('Identity',levels = lv)
    ),by=Phrase]
    ,d[,list(
      Year=Year[-1]
      ,Frequency=diff(Frequency)  %>% lintran(s1=c(0,max(Frequency[-1]))) %>% BoxCox(BoxCox.lambda(Frequency))
      ,Diffs=factor('First Difference',levels = lv)
    ),by=Phrase]
  ))
  
  p <- ggplot(d, aes(Year,Frequency)) + geom_line() +
    geom_text(
      aes(x, y, label=Phrase)
      ,data=data.frame(x=-Inf,y=Inf,Phrase=unique(d$Phrase),Diffs='Identity')
      , hjust=0,vjust=1,size=3) +
    facet_wrap( ~ Phrase + Diffs,scales='free_y',ncol=2,strip.position='left') +
    theme(
      axis.text.x = element_text(angle = 90,vjust=.5,debug=F)
      ,legend.position="bottom"
      ,panel.grid.minor.x = element_blank()
      ,panel.grid.minor.y = element_blank()
      ,strip.background = element_blank()
      #,strip.placement = "outside"
      ,strip.text.y = element_blank()
      ,axis.text.y = element_blank()
      ,axis.ticks = element_blank()
    ) +
    scale_x_continuous(breaks=seq(round(min(d$Year),-1),round(max(d$Year),-1),10))  +
    # scale_y_continuous(breaks=function(x) ifelse(max(x)==1,list(c(0,.5,1)),list(c(-.2,-.1,0,.1,.2)))[[1]]) +
    ylab('Scaled Frequency')
  
  if(length(unique(d$Phrase))>2) {
    warning('Granger test not implemented for more than two phrases.')
    return(p)
  } else {
    setkey(d,Diffs,Phrase,Year)
    rl<-select.lags(d[list('First Difference',order[2]),Frequency],d[list('First Difference',order[2]),Frequency],mxlg)
    lr<-select.lags(d[list('First Difference',order[1]),Frequency],d[list('First Difference',order[2]),Frequency],mxlg)
  }
  ret<-list(
    test=data.table(
      lag=c(-(mxlg:1),0,1:mxlg)
      ,pval=1-c(rev(rl$pvals),rep(NA,3),lr$pvals)
      ,aic=c(rev(rl$ic[,'aic']),NA,lr$ic[,'aic'])
      ,bic=c(rev(rl$ic[,'bic']),NA,lr$ic[,'bic'])
    )
    ,plot=list(d=p)
  )
  r<-dim(ret$test)[1]
  ret$plot$t<-ggplot(data=melt(ret$test,id.vars = 'lag',variable.name = 'criterion')
                     ,mapping = aes(x=lag,y=value)) +
    geom_line() + geom_point(size=2,shape=21,fill=c(
      replace(rep('black',r),which.min(ret$test$pval),'white')
      ,replace(rep('black',r),which.min(ret$test$aic),'white')
      ,replace(rep('black',r),which.min(ret$test$bic),'white')
    )) +  facet_grid(criterion~.,scales = 'free_y') +
    labs(title='Model Selection Criteria',x=paste('(',paste(paste(rev(order),collapse=' > '),')          Lags          (',paste(order,collapse=' > '),')'))) +
    scale_x_continuous(breaks=ret$test$lag,minor_breaks = NULL) + scale_y_continuous(minor_breaks = NULL)
  ret
}

tts2arima.f<-function(
  gbng2tts
  ,by=c('batch','Phrase')
)
{
  require(forecast)
  tts2arima<-gbng2tts[,list(aa=list(
    auto.arima(
      ts(Frequency,start = min(Year),frequency = 1)
      ,lambda = BoxCox.lambda(Frequency)
      ,stepwise = F
      ,seasonal = F
      ,trace=F)
  )),keyby=by]
  tts2arima[,Predicted:=list(lapply(aa,fitted))]
  tts2arima[,`Predicted t-Score`:=list(lapply(fit,tscore))]
  tts2arima
}

# from http://stats.stackexchange.com/questions/160671/estimate-lag-for-granger-causality-test
select.lags<-function(x,y,max.lag=8) {
  y<-as.numeric(y)
  y.lag<-embed(y,max.lag+1)[,-1,drop=FALSE]
  x.lag<-embed(x,max.lag+1)[,-1,drop=FALSE]
  
  t<-tail(seq_along(y),nrow(y.lag))
  
  ms=lapply(1:max.lag,function(i) lm(y[t]~y.lag[,1:i]+x.lag[,1:i]))
  
  pvals<-mapply(function(i) anova(ms[[i]],ms[[i-1]])[2,"Pr(>F)"],max.lag:2)
  ind<-which(pvals<0.05)[1]
  ftest<-ifelse(is.na(ind),1,max.lag-ind+1)
  
  aic<-as.numeric(lapply(ms,AIC))
  bic<-as.numeric(lapply(ms,BIC))
  structure(list(ic=cbind(aic=aic,bic=bic),pvals=pvals,
                 selection=list(aic=which.min(aic),bic=which.min(bic),ftest=ftest)))
}
