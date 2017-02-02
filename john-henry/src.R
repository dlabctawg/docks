# utilities
lintran<-function(x,s1=c(0,1),s2=c(0,1)) {a=diff(s2)/diff(s1);b=s2[1]-a*s1[1];return(a*x+b)}
ec<-function(x) strsplit(x,',')[[1]]

#'Converts full text data to stm style document term matrix. Redundant with stm function textProcessor()
#'
#' @param source.dir A path to a folder containing a plain text file for each document.
#' @param string A character string where each element is a full document.
#' @param lt
#' @param save.to.disk
#' @param check.for.saved.output
#' @param out.dir
#' @param sample.docs
#' @param names A vector of document names. Coerced to character. If unspecified no names will be given. If unspecified and source given, files names will be assigned instead.
ftxt2stmbow.f<-function(
	source.dir=NULL
	,string=NULL
	,names=NULL
	,lt=1
	,save.to.disk=F
	,check.for.saved.output=F
	,out.dir
	,sample.docs
){
	require(stm,quietly = T)
	require(data.table,quietly = T)
	require(SnowballC,quietly = T)
	require(tm,quietly = T)
	pfs<-.Platform$file.sep
	sfn<-"ftxt2stmbow.RData"
	if(check.for.saved.output) if(any(grepl(sfn,dir(recursive=T,full.names=T,ignore.case=T)))) {
		warning(paste('Loading and returning first saved',sfn),call.=F)
		l<-ls()
		load(dir(pattern=sfn,full.names=T,recursive=T,ignore.case=F)[1])
		return(get(setdiff(ls(),c(l,'l'))))
	}

	if(all(sapply(list(source.dir,string),is.null))) stop('Specify either source.dir or string.')
	if(all(!sapply(list(source.dir,string),is.null))) stop('Specify either source.dir or string but not both.')


	# 1. Preprocessing functions in base R
#something stupid is happening between docs and string
	if(!is.null(source.dir)) {
		docs<-list() # container for docs
		files<-list.files(source.dir,full.names=T,recursive = T)
		if(!missing(sample.docs)) files<-sample(files,sample.docs)
		files<-files[order(sub(paste('.*','(.+)',sep=pfs),'',files))] # helps later to have files in alpha order by document name
		for(i in files) docs[[i]]<-readLines(i,warn=F)
		docs<-lapply(docs,FUN=paste,collapse=' ') # each doc is one long character string
		txt<-unlist(docs)
		string<-sapply(docs,FUN=strsplit,split='\\s') # split docs into words "\\s+" is a regex. "[[:space:]]+" also works but is R specific
		if(is.null(names)) names<-sub(paste0('^.+',pfs,'(.+)\\..+'),'\\1',files)
		names(string)<-names
		docs<-string
	} else {
		docs<-string
		txt<-unlist(docs)
		string<-sapply(docs,FUN=strsplit,split='\\s') # split docs into words "\\s+" is a regex. "[[:space:]]+" also works but is R specific
		if(is.null(names)) {names(string)<-names} else {names<-names(string)}
	}

	docs<-lapply(docs,FUN=tolower) # transform to lower case
	docs<-lapply(docs,FUN=removePunctuation) # ...
	docs<-lapply(docs,FUN=removeNumbers)
	docs<-lapply(docs,FUN=removeWords,stopwords('english'))
	docs<-lapply(docs,FUN=stemDocument,language='english')

	## here we end with a list of untabulated tokenized character vectors in original document order, with blanks preserved.

	# 2. Match tokens to originals for applying stm results for qualitative cross validation.
	ftxt2stmbow<-list()
	ftxt2stmbow$map<-mapply(function(o,t) {t<-c(t,rep('',length(o)-length(t)));data.table(o,t,ord=1:length(o))},o=string,t=docs,SIMPLIFY = F) #original and tokens
	for(i in 1:length(ftxt2stmbow$map)) ftxt2stmbow$map[[i]]$o[!sapply(ftxt2stmbow$map[[i]]$o,length)]<-'' # replace character(0) with ""
	docs<-lapply(docs,FUN=function(x) x[!!nchar(x)]) #remove blanks

	# 3. Sparse matrix format expected by stm.

	ftxt2stmbow$vocab<-sort(unique(unlist(docs))) # stm expects as input a list of matrices where each element of the list is a document and where the first row of each matrix is the index of a word and the second row is the count of the word. This is a more memory efficient form since zeros are not stored.
	for(i in 1:length(docs)){
		tb<-table(docs[[i]])
		ftxt2stmbow$documents[[i]]<-rbind(vocab.index=which(ftxt2stmbow$vocab%in%names(tb)),frequency=tb)
	}

	names(ftxt2stmbow$documents)<-names
	p<-prepDocuments(
		ftxt2stmbow$documents
		,ftxt2stmbow$vocab
		,meta = data.table(names,sapply(ftxt2stmbow$map,function(x) sum(x$t!='')/sum(x$o!='')),1:length(ftxt2stmbow$documents))
		,lower.thresh=lt
	)
setnames(p$meta,c('names','count.prop','ord'))
setkey(p$meta,ord)
p$meta[,ord:=1:length(ord)]
	if(length(p$docs.removed)) {
		cat('\nRemoved:',sub(paste0('.+',pfs),'',files[p$docs.removed]))
		ftxt2stmbow$map[p$docs.removed]<-NULL
	}
  if('docs.removed'%in%names(p)){
    drt<-T
    if(is.null(p$docs.removed)) drt<-F
  }
	ftxt2stmbow<-do.call(list,c(txt=ifelse(drt,list(txt[-p$docs.removed]),list(txt)),ftxt2stmbow['map'],p))

	if(save.to.disk) try(save(ftxt2stmbow,file=gsub(paste0(pfs,'+'),pfs,paste(out.dir,sfn,sep=pfs))))
	ftxt2stmbow
}

#'Imports from ICPSR Congressional Record files saved as a CTAWG list, and performs conventional preprocessing to yield an stm formatted term-document "sparse matrix".
#' @param icpsr.cong A CTAWG list of 3 data frames named according to ICPSR plain text tables and containing the "GPOspeech", "SpeakerID", and "GPOdescr" substrings.
icpsr2stmbow.f<-function(icpsr,sample.size=100,wrdmin=0,wrdmax=0,lt=1,save.to.disk=F,check.for.saved.output=F,out.dir=NULL){

	sfn<-paste("icpsr2stmbow-samp",ifelse(sample.size,sample.size,'all'),".RData",sep="")
	if(check.for.saved.output) if(any(grepl(sfn,dir(recursive=T,full.names=T,ignore.case=T)))) {
		warning(paste('Loading and returning first saved',sfn),call.=F)
		l<-ls()
		load(dir(pattern=sfn,full.names=T,recursive=T,ignore.case=F)[1])
		return(get(setdiff(ls(),c(l,'l'))))
	}

	require(stm,quietly = TRUE)
	s<-1:nrow(icpsr[[grep('GPOspeech',names(icpsr))]])
	if(wrdmin) {wn<-icpsr[[grep('GPOdescr',names(icpsr))]]$word.count>=wrdmin} else {wn<-rep(T,length(s))}
	if(wrdmax) {wx<-icpsr[[grep('GPOdescr',names(icpsr))]]$word.count<=wrdmax} else {wx<-rep(T,length(s))}
	s<-s[wn&wx]
	if(sample.size) s<-sample(s,sample.size)
	GPOspeech<-grep('GPOspeech',names(icpsr))
	icpsr2stmbow<-textProcessor(
		documents=icpsr[[GPOspeech]]$speech[s]
		,lowercase=TRUE
		,removestopwords=TRUE
		,removenumbers=TRUE
		,removepunctuation=TRUE
		,stem=TRUE
		,wordLengths=c(3,Inf)
		,sparselevel=1
		,language="en"
		,verbose=TRUE
		,onlycharacter= FALSE
		,striphtml=FALSE
		,customstopwords=NULL)
	nm<-icpsr[[GPOspeech]]$speechID[s]
	if(length(icpsr2stmbow$docs.removed)) nm<-nm[-icpsr2stmbow$docs.removed]
	p<-prepDocuments(icpsr2stmbow$documents,icpsr2stmbow$vocab,meta = nm,lower.thresh = lt)
	icpsr2stmbow$documents<-p$documents
	icpsr2stmbow$vocab<-p$vocab
	names(icpsr2stmbow$documents)<-p$meta
	if(save.to.disk) try(save(icpsr2stmbow,file=paste(grep(paste(out.dir,'$',sep=''),dir(include.dirs = T,full.names=T,recursive=F,ignore.case=F),value = T)[1],sfn,sep=pfs)))
	icpsr2stmbow
}

#' Browse ICPSR records using either speech or speaker IDs.
#' @param icpsr A CTAWG list of 3 data frames named according to ICPSR plain text tables and containing the "GPOspeech", "SpeakerID", and "GPOdescr" substrings.
#' @param speechID The name of the field containing the speech ID. First three digits indicate the nth Congress, next 4 indicate the year, following 7 are the numerical index in order of appearance.
#' @param speakerID Speaker ID, which the function links via the "GPOdesc" table to merge speaker and speech data.
browse.icpsr<-function(icpsr,speechID=NULL,speakerID=NULL,print.max=50){
	require(data.table,quietly = T)
	if(is.null(speechID)&is.null(speakerID)) stop('Must specify either speechID or speakerID.')

	icpsr<-icpsr.cong107
	for(i in grep('GPOspeech',names(icpsr))) icpsr[[i]]$speechID<-as.character(icpsr[[i]]$speechID)
	GPOspeech<-rbindlist(icpsr[grep('GPOspeech',names(icpsr))])
	icpsr[grep('GPOspeech',names(icpsr))]<-NULL
	setkey(GPOspeech,speechID)

	for(i in grep('SpeakerID',names(icpsr))) icpsr[[i]]$id<-as.character(icpsr[[i]]$id)
	SpeakerID<-rbindlist(icpsr[grep('SpeakerID',names(icpsr))])
	icpsr[grep('SpeakerID',names(icpsr))]<-NULL
	setnames(SpeakerID,'id','speakerID')
	setkey(SpeakerID,speakerID)

	for(i in grep('GPOdesc',names(icpsr))) {
		icpsr[[i]]$speechID<-as.character(icpsr[[i]]$speechID)
		icpsr[[i]]$speakerID<-as.character(icpsr[[i]]$speakerID)
	}
	GPOdesc<-rbindlist(icpsr[grep('GPOdesc',names(icpsr))])
	icpsr[grep('GPOdesc',names(icpsr))]<-NULL

	if(!is.null(speechID)){
		speechID<-as.character(speechID)
		setkey(GPOdesc,speechID)
		ret<-merge(GPOdesc,GPOspeech[speechID])
		setkey(ret,speakerID)
		ret<-merge(SpeakerID,ret)
		if(nrow(ret)<print.max) {cat(paste('Printing first',floor(print.max),'documents.\n'));for(i in 1:ifelse(nrow(ret)>print.max,print.max,nrow(ret))) {print(t(ret[i]));cat('\n')}}
		return(ret)
	}
	if(!is.null(speakerID)){
		speakerID<-as.character(speakerID)
		setkey(GPOdesc,speakerID)
		ret<-merge(GPOdesc,SpeakerID[speakerID])
		setkey(GPOdesc,speechID)
		setkey(ret,speechID)
		ret<-merge(ret,GPOspeech)
		cat(paste('Printing first',floor(print.max),'documents.'))
		if(nrow(ret)<print.max) {cat(paste('Printing first',floor(print.max),'documents.\n'));for(i in 1:ifelse(nrow(ret)>print.max,print.max,nrow(ret))) {print(t(ret[i]));cat('\n')}}
		return(ret)
	}
}

#' Takes stm style bag of words as input and gives LDA.
#'
#' @param stmbow stm formatted bag of words, as given by textProcessor command or stmbow importer.
#' @param k Number of topics.
#' @param alpha The alpha parameter must be greater than 0. Alpha < 1 assumes that each document is constructed from relatively few topics. Alpha > 1 assumes that each is constructed from many topics. If you aren't sure choose the convention: 50/k, which will also be the default if you specify nothing.
#' @param visualize.results Logical, if TRUE provides visualization tools to help interpret LDA output.
#' @param verbose Print iteration progress (very verbose!).
#' @param out.dir
#' @param sig.pri
#' @param it
#' @param check.for.saved.output
#' @param ... Other arguments to stm.
#' @param save.to.disk Save an .RData file of the object; helpful for lengthy estimations. May be buggy across environments.
stmbow2topmod.f<-function(
	stmbow,out.dir
	,k=0,alpha=NULL,sig.pri=.5,it='Spectral'
	,visualize.results=F
	,verbose=F,save.to.disk=F,check.for.saved.output=F
	,...
)
{
  pfs<-.Platform$file.sep
  if(is.null(alpha)) alpha<-50/k
	if(check.for.saved.output) if(any(grepl(paste("stm-model-k",k,"-alpha",round(alpha,3),".RData",sep=""),dir(recursive=T,full.names=T,ignore.case=T)))) {
		warning(paste('Loading and returning first saved',paste("stm-model-k",k,"-alpha",round(alpha,3),".RData.",sep="")),call.=F)
		load(dir(pattern=paste("stm-model-k",k,"-alpha",round(alpha,3),".RData",sep=""),full.names=T,recursive=T,ignore.case=F)[1])
		pdf(paste(out.dir,'topicquality.pdf',sep=pfs))
		topicQuality(stmbow2topmod$mod,stmbow$documents)
		dev.off()
		return(stmbow2topmod)
	}
	# Check package requirements and arguments
	require(stm,quietly = T)
	require(data.table)
	
	stmbow2topmod<-selectModel(
		documents=stmbow$documents
		,vocab=stmbow$vocab
		,K=k
		,init.type = it
		,control=list(alpha=alpha)
		,sigma.prior=sig.pri
		,verbose=verbose
		,runs = 10
		,...
		)
	stmbow2topmod$model<-stmbow2topmod$runout[[1]]
	stmbow2topmod$top.word.phi.beta<-sapply(data.frame(stmbow2topmod$model$beta$logbeta),function(x) sapply(x,function(y) ifelse(is.infinite(y),.Machine$double.eps,exp(y)))) # called beta by stm, epsilon closest thing to zero the machine can represent, necessary to prevent error
	colnames(stmbow2topmod$top.word.phi.beta)<-stmbow2topmod$model$vocab
	stmbow2topmod$doc.top.theta<-stmbow2topmod$model$theta
	rownames(stmbow2topmod$doc.top.theta)<-names(stmbow$documents)
	stmbow2topmod$doc.length<-sapply(stmbow$documents,ncol)
	stmbow2topmod$vocab<-stmbow2topmod$model$vocab
	tn<-data.table(do.call(rbind,sapply(stmbow$documents, t)))
	setnames(tn,c('ix','freq'))
	setkey(tn,ix)
	tn<-tn[,list('freq'=sum(freq)),by=ix]
	stmbow2topmod$term.frequency<-tn$freq
	names(stmbow2topmod$term.frequency)<-stmbow$vocab[tn$ix]
	stmbow2topmod$documents<-stmbow$documents
	pdf(paste(out.dir,'topicquality.pdf',sep=pfs))
	stmbow2topmod$tq<-topicQuality(stmbow2topmod$mod,stmbow$documents)
	dev.off()
	if(save.to.disk){
		f<-paste(out.dir,paste("stm-model-k",k,"-alpha",round(alpha,3),".RData",sep=""),sep=pfs)
		save(stmbow2topmod,file=f)
	}
	stmbow2topmod
}

#' Different approaches to reduce and visualize the standard topic-document and topic-word word matrices output by LDA and other topic modelling estimations.
lda2viz.f<-function(stmbow2topmod,out.dir,rt=F,ob=F,launch=T){
	# from http://cpsievert.github.io/LDAvis/newsgroup/newsgroup.html
	# http://nlp.stanford.edu/events/illvi2014/papers/sievert-illvi2014.pdf
	# http://glimmer.rstudio.com/cpsievert/xkcd/
	require(LDAvis,quietly = T)
	require(servr,quietly = T)
	library(RJSONIO,quietly = T)

	json <- createJSON(
		phi = stmbow2topmod$top.word.phi.beta
		,theta = stmbow2topmod$doc.top.theta
		,vocab = stmbow2topmod$vocab
		,doc.length = stmbow2topmod$doc.length
		,term.frequency = stmbow2topmod$term.frequency
		,reorder.topics = rt
	)
	save(json,file=paste0(out.dir,.Platform$file.sep,'viz.RData'))
	cat('Topic index maps to propability index like so:\n')
	to<-fromJSON(json)$topic.order
	print(rbind(viz=to))
	od<-paste(out.dir,'viz',sep=.Platform$file.sep)
	serVis(json,open.browser = ob,out.dir=od)
	call<-paste('python -m SimpleHTTPServer && open http://localhost:8000',sep='')
	if(launch) {setwd(od);cat('You may need to manually kill the Python process');system(call,wait=F)}
	list(call=call,top.ord=to)
	}

#'Network visualizations and clustering.
lda2netviz.f<-function(stmbow2topmod,thresh="choose"){
	require(igraph,quietly=T)
	require(network,quietly = T)

	bam<-stmbow2topmod$doc.top.theta
	b<-quantile(bam,seq(0,1,.05))
	h<-hist(bam,breaks=b,col="black")
	abline(v=b,col=gray(0,.5))
	text(
		x=rev(h$breaks)
		,y=seq(0,max(h$density),length.out=21)
		,labels=rev(paste('|',names(b),' <= ',round(b,3),sep=''))
		,pos=4
		,offset=0
		,cex=.5
	)
	if(thresh=="choose"){
		cat("\nPlease choose an edge weight threshold by clicking on the histogram near the x-axis where you would like to cut the distribution of probabilities that a document draws words from a particular topic (i.e. theta or the document-topic probability matrix). Relationships between documents and topics that fall below this threshold will be ignored.\n")
		thresh<-locator(n=1,type="p",col="red")
		abline(v=thresh$x,col="red")
		text(
			x=thresh$x
			,y=thresh$y
			,labels=rev(paste('|',round(mean(bam<thresh$x)*100,2),'% <= ',round(thresh$x,3),sep=''))
			,pos=4
			,offset=0
			,cex=1
			,col="red"
		)
	}
	bam[bam<thresh$x]<-0
	m1am<-bam%*%t(bam)
	m1net<-network(m1am,directed=F,loops=F)
	browser()
	network.vertex.names(m1net)<-sub(paste('.*','(.+)',sep=.Platform$file.sep),'\\1',rownames(m1am))
	#	pdf('doc-by-top-net.pdf')
	plot(m1net
			 ,displaylabels=T
			 ,label=paste('T',1:nrow(m1am),sep='')
			 ,label.pos=5
			 ,label.col="white"
			 ,label.cex=.25
			 ,vertex.col="black"
			 ,vertex.cex=2
	)
	#	dev.off()

	m2am<-t(bam)%*%bam
	m2net<-network(m2am,directed=F,loops=F)

	#	pdf('top-by-doc-net.pdf')
	plot(m2net
			 ,displaylabels=T
			 ,label=paste('D',1:nrow(m2am),sep='')
			 ,label.pos=5
			 ,label.col="black"
			 ,label.cex=.75
			 ,vertex.col="white"
			 ,vertex.cex=2
	)
	#	dev.off()

	# 	b<-quantile(m1am,seq(0,1,.1))
	# 	h<-hist(m1am,breaks=b)
	# 	abline(v=b,col'pink')
	# 	h<-hist(m2am,breaks=quantile(m2am,seq(0,1,.1)))

	bel<-which(bam>0,arr.ind=T)
	w<-bam[bel]
	bel<-cbind(
		sub(paste(".+","(.+)",sep=.Platform$file.sep),"\\1",rownames(bel))
		,paste("t",bel[,2],sep="")
		)
	browser()
	o<-order(bel[,1],bel[,2])
	bel<-data.frame(bel[o,])
	w<-w[o]
	colnames(bel)<-c("document","topic")
	nm1<-length(unique(bel$document))
	nm2<-length(unique(bel$topic))
	bnet<-network(bel,bipartite=nm1,matrix.type="edgelist")
	bnet%e%"w"<-w
	pdf('bimodal-net.pdf')
	plot(
		bnet
		,displaylabels=T
		,label=c(
			paste(1:nm1)
			,network.vertex.names(bnet)[-(1:nm1)]
		)
		,label.pos=5
		,label.col=c(rep("white",nm1),rep("black",nm2))
		,label.cex=.75
		,vertex.col=c(rep("black",nm1),rep("white",nm2))
		,vertex.cex=2
		#,vertex.sides=c(rep(3,nm1),rep(20,nm2))
	)
	dev.off()


}

#'Community detection.
lda2netcom.f<-function(stmbow2topmod,out.dir,freq.weight,reps=10,mx.grp=30){
	require(data.table)
	require(igraph)
	require(magrittr)

	if(!missing(freq.weight)) stmbow2topmod$doc.top.theta<-freq.weight * stmbow2topmod$doc.top.theta
	m<-stmbow2topmod$doc.top.theta%*%t(stmbow2topmod$doc.top.theta)
	g<-graph_from_adjacency_matrix(
		adjmatrix = m
		,mode = 'undirected'
		,weighted = 'ew'
		,diag = F
	)

	md<-list()
	bmod=-Inf
	for(i in 1:reps){
		cat('.')
		sg<-spinglass.community(g,spins=mx.grp,weights=E(g)$ew)
		md[[i]]<-modularity(sg,weights=E(g)$ew)
		if(md[[i]]>bmod) {
			bmod<-md[[i]]
			sgb<-sg
		}
	}
	save(sgb,file=paste(out.dir,'sgb.RData',sep=.Platform$file.sep))

	# diagnostics
	md<-unlist(md)
	pdf(paste(out.dir,'com.pdf',sep=.Platform$file.sep),w=8.5,h=11)
	par(mfrow=c(2,1),mar=rep(2,4))
	plot(md,type='l',main=paste('Modularity of ',reps,' trial',ifelse(reps==1,'','s'),sep=''))
	abline(h=mean(md),col='blue')
	points(x=which.max(md),y=max(md),col='red')
	plot(density(md),main='Modularity density')
	abline(v=max(md),col='red')
	abline(v=mean(md),col='blue')
	dev.off()

	# sort by topic association
	rk<-data.table(
		group=membership(sgb)
		,stmbow2topmod$doc.top.theta
	)
	setnames(rk,sub('^V','T',names(rk)))
	cols<-2:ncol(rk)
	rm<-apply(rk[,-1],1,function(x) list(names(sort(-x)[1:2])))
	rm<-rk[,list(group,rm)]
	rm<-rm[,list(m=list(unique(unlist(rm)))),keyby=group]
	rm[,group:=as.character(group)]
	setkey(rm,group)

	rk<-rk[ ,lapply(.SD, mean), .SDcols = cols,keyby=group]
	rs<- rank(-apply(rk[,-1],1,max))
	rk[,group:=factor(group,levels=order(rs))]
	rk<-melt(data = rk,id.vars = 'group',variable.name = 'name',value.name = 'p')
	rk<-rk[,rank:=rank(-p),by=group]
	setkey(rk,group)

	for(i in as.character(1:nrow(rm))) rk[i,lb:=ifelse(name%in%rm[list(i),m][[1]],as.character(name),'')]

	setkey(rk,lb)
	ix<-rk[!'',do.call(seq,as.list(range(rank)))]
	setkey(rk,rank)
	xpl<-rk[,range(p)]
	xpl<-xpl+c(-.06*xpl[2],.06*xpl[2])
	tp<-ggplot(rk[list(ix)],aes(x=rank,y=p,label=lb)) +
		geom_line(color='gray') +
		geom_text(size=3,angle=90,color='red') +
		ggtitle('Scree plot for average topic probabilities within groups') +
		facet_wrap(~group) + theme(
			legend.position="bottom"
			#,panel.grid.minor.x = element_blank()
			#,panel.grid.major.x = element_blank()
			#,panel.grid.minor.y = element_blank()
			#,panel.grid.major.y = element_blank()
			#,strip.background = element_blank()
			#,strip.placement = "outside"
			#,strip.text.y = element_blank()
			#,axis.text.x = element_blank()
			#,axis.text.y = element_blank()
			#,axis.ticks.x = element_blank()
			#,axis.ticks.y = element_blank()
		) + expand_limits(y=xpl)
	ggsave(device = 'pdf',filename = 'screecom.pdf',path=out.dir)

	com<-list(mem=membership(sgb),c=sgb,g=g,mod=md,r=rk,p=tp)
}

#'Color coding original documents in topic highlighting.
lda2ftxt.f<-function(
	map,doc.top,top.word,lda2rel
	,intensify=T,intensity=.3
	,sample=10
	,out.dir=NULL
	,index=0,spacing=.4,ptsize=10,axes=F
	,pdf=F,mxdoc.word.prop=0
	,fname
	,top.ord
)
{
	# d<-dir(include.dirs = T,recursive=T,full.names = T)
	# out.dir<-grep(paste(out.dir,'$',sep=''),d,value=T)
	# if(!length(out.dir)) {warning('out.dir not found, saving to current working directory.');out.dir<-getwd()}

	require(colorspace)
	library(RColorBrewer)
	require(data.table)

	mx<-max(sapply(1:nrow(doc.top), function(x) max(doc.top[x,]*top.word)))
	s<-1:length(map)
	if(sample) s<-sort(sample(s,sample))
	if(index[1]) s<-index
	ret<-list()
	for(i in s){

		cat(c('\nRendering document',i))
		w<-order(doc.top[i,],decreasing=T)[1:2]
		tn<-doc.top[i,w]
		names(tn)<-paste('T',w,sep='')
		p<-round(sum(tn)*100,1)

		cat('\nCalculating document\'s topic by word probability matrix =\n(document by topic prob vector) * (global topic by word prob matrix)')
		m<-t(doc.top[i,w]*top.word[w,])

		cat('\nOriginal range of predicted document\'s topic by word probabilities:')
		print(range(m))

		if(intensify) {
			m<-lintran(m,c(0,mx),c(intensity,1))
			cat('Range of document\'s topic by word probabilities after linear amplification:')
			print(range(m))
		}

		colnames(m)<-names(tn)
		m<-data.table(
			t=rownames(m),m
			,r=apply(m,1,function(x) x[2]/sum(x))
			,w=apply(m,1,function(x) sum(x))
		)
		setkey(m,t)
		setkey(map[[i]],t)
		m<-merge(map[[i]],m,all.x=T,all.y=F)
		if(all(is.na(m[,4:7,with=F]))) {cat('Too short, skipping.');next}

		pal<-rev(brewer.pal(n = 7,name = 'RdYlGn'))
		col<-round(range(na.omit(m$r))*1000)-500
		col<-max(abs(col))+1
		col<-(500-col):(500+col)
		col<-cut(col,breaks = length(pal),label=pal,include.lowest = T) #hex(HSV(H=240+120*m$r,S=1,V=m$w))
		rd<-round(m$r*1000)
		rd[!is.na(rd)]<-rd[!is.na(rd)]-min(rd[!is.na(rd)])+2
		rd[is.na(rd)]<-1
		col<-c('gray80',as.character(col))[rd]
		m[,'col':=col]
		setkey(m,ord)
		mr<-.6

		par(mar=rep(0,4))
		plot.new()
		pw<-8
		ph<-1
		par(family='Times',ps=ptsize,fin=c(pw,ph))
		plot.window(xlim=c(0,pw),ylim=c(0,1))

		o<-as.vector(rbind(m$o,' '))
		t<-as.vector(rbind(m$t,' '))
		c<-as.vector(rbind(m$col,'pink'))
		shu<-strheight(LETTERS,units = 'inches')
		swu<-strwidth(o,units = 'inches')*1.25

		fl<-list()
		x<-list()
		sww<-swu
		ct<-0
		f<-cumsum(sww)%/%pw
		while(max(f)>0){
			ft<-f==0
			x[[length(x)+1]]<-cumsum(sww[ft])
			sww<-sww[!ft]
			fl[[length(fl)+1]]<-rep(ct,sum(ft))
			f<-cumsum(sww)%/%pw
			ct<-ct+1
		}
		ft<-f==0
		x[[length(x)+1]]<-cumsum(sww[f==0])
		fl[[length(fl)+1]]<-rep(ct,sum(ft))
		f<-unlist(fl)
		x<-unlist(sapply(x,function(x) c(0,x[-length(x)])))
		l<-max(f)*2+1.5

		lda2relc<-copy(lda2rel)
		setkey(lda2relc,Category)
		lda2relc<-lda2relc[names(tn)]
		lda2relc[,inc:=Term%in%t]

		leg<-sum(c(
			head=4
			,rows=nrow(lda2relc)/length(unique(lda2relc$Category))/length(unique(lda2relc$lambda))
			,foot=0
		))

		ph<-sum(c(lines=l,legend=leg,divider=1))

		par(mar=c(b=2.25,l=2.25,t=3.25,r=2.25),fin=c(pw,ph*spacing),yaxs='i')
		plot.window(xlim=c(0,pw),ylim=c(ph,1))
		box()
		if(axes){
			axis(1)
			axis(2)
		}
		abline(h=leg)

		#Topic Super Column headings
		text(x=seq(0,pw,length.out = 9)[c(3,7)],y=2,labels=names(tn),col=pal[c(1,length(pal))],font=2,pos=1,offset=0)
		#Rank listed
		text(x=seq(0,pw,length.out = 9)[-c(1,5,9)],y=3,labels=c('Anchor','~','Common'),font=3,pos=1,offset=0)

		loc<-data.table(expand.grid(Category=names(tn),lambda=unique(lda2relc$lambda)),col=c(2,6,3,7,4,8))
		setkey(lda2relc,Category,lambda)
		setkey(loc,Category,lambda)
		lda2relc<-merge(lda2relc,loc)
		lda2relc[,clr:=ifelse(inc,'black','gray')]

		# List of terms
		text(x=seq(0,pw,length.out = 9)[lda2relc$col],y=lda2relc$ord+3,labels=lda2relc$Term,col = lda2relc$clr,pos=1,offset=0)
		#Original text
		text(x=x,y=f*2+leg+1,labels=o,pos=4,offset=0,cex=1,col='black',ps=11)
		#Tokenized text
		text(x=x,y=f*2+leg+1.5,labels=t,pos=4,offset=0,cex=1,font=3,col=c)
		#Line numbers
		text(x=0,y=unique(f)*2+leg+1.25,labels=unique(f)+1,pos=4,offset=-0.5,cex=.75,col='lightgray',ps=11,srt=90)

		if(missing(top.ord)) {to<-''} else {to<-paste0('-',top.ord[as.numeric(sub('T','',names(tn)))])}

		title(main=paste0('Green ',names(tn)[1],to[1],' ~ ',ifelse(missing(fname),'Document',fname[i]),'-',sprintf('%03d',i),' ~ ','Red ',names(tn)[2],to[2],'\n',round(tn*100,2)[1],' + ',round(tn*100,2)[2],' = ',round(sum(tn)*100,2),' %',sep=''))

		if(pdf) {
			try(dev.copy2pdf(
				file=paste(out.dir,paste0(sprintf('%02d',round(tn[1]*100)),'-',ifelse(missing(fname),'doc',fname[i]),'-',sprintf('%03d',i),'.pdf'),sep=.Platform$file.sep)
				, width=pw, height=ph*spacing),silent = T)
			dev.off(dev.prev())
		}
		try(dev.off())
		cat('Range of ',ifelse(intensify,'amplified ',''),'probabilites for ',names(m)[[4]],sep='')
		print(range(na.omit(m[[4]])))
		cat('Range of ',ifelse(intensify,'amplified ',''),'probabilites for ',names(m)[[5]],sep='')
		print(range(na.omit(m[[5]])))

	}

	m
}

#' Original authors code to compute relevance. 'Forked' from https://github.com/cpsievert/LDAvis/blob/6f93aa85499b705c9ae6c56e5985df637f9f5132/R/createJSON.R
lda2rel.f<- function(stmbow2topmod,R = 10,lambda.step = 0.5,reorder.topics = FALSE,save.to.disk=F,check.for.saved.output=F,out.dir,...) {

	sfn<-'lda2rel.RData'
	if(check.for.saved.output) if(any(grepl(sfn,dir(recursive=T,full.names=T,ignore.case=T)))) {
		warning(paste('Loading and returning first saved',sfn),call.=F)
		l<-ls()
		load(dir(pattern=sfn,full.names=T,recursive=T,ignore.case=F)[1])
		return(get(setdiff(ls(),c(l,'l'))))
	}

	require(data.table,quietly = T)
	phi = stmbow2topmod$top.word.phi.beta
	theta = stmbow2topmod$doc.top.theta
	vocab = stmbow2topmod$vocab
	doc.length = stmbow2topmod$doc.length
	term.frequency = stmbow2topmod$term.frequency
	#rm(stmbow2topmod)
	# Set the values of a few summary statistics of the corpus and model:
	dp <- dim(phi)  # should be K x W
	dt <- dim(theta)  # should be D x K

	N <- sum(doc.length)  # number of tokens in the data
	W <- length(vocab)  # number of terms in the vocab
	D <- length(doc.length)  # number of documents in the data
	K <- dt[2]  # number of topics in the model

	# check that certain input dimensions match
	if (dp[1] != K) stop("Number of rows of phi does not match
											 number of columns of theta; both should be equal to the number of topics
											 in the model.")
	if (D != dt[1]) stop("Length of doc.length not equal
											 to the number of rows in theta; both should be equal to the number of
											 documents in the data.")
	if (dp[2] != W) stop("Number of terms in vocabulary does
											 not match the number of columns of phi (where each row of phi is a
											 probability distribution of terms for a given topic).")
	if (length(term.frequency) != W) stop("Length of term.frequency
																				not equal to the number of terms in the vocabulary.")
	if (any(nchar(vocab) == 0)) stop("One or more terms in the vocabulary
																	 has zero characters -- all terms must have at least one character.")

	# check that conditional distributions are normalized:
	phi.test <- all.equal(rowSums(phi), rep(1, K), check.attributes = FALSE)
	theta.test <- all.equal(rowSums(theta), rep(1, dt[1]),
													check.attributes = FALSE)
	if (!isTRUE(phi.test)) stop("Rows of phi don't all sum to 1.")
	if (!isTRUE(theta.test)) stop("Rows of theta don't all sum to 1.")

	# compute counts of tokens across K topics (length-K vector):
	# (this determines the areas of the default topic circles when no term is
	# highlighted)
	topic.frequency <- colSums(theta * doc.length)
	topic.proportion <- topic.frequency/sum(topic.frequency)

	# re-order the K topics in order of decreasing proportion:
	if(reorder.topics) {o <- order(topic.proportion, decreasing = TRUE)} else {o <- seq_along(topic.proportion)}

	phi <- phi[o, ]
	theta <- theta[, o]
	topic.frequency <- topic.frequency[o]
	topic.proportion <- topic.proportion[o]

	# compute intertopic distances using the specified multidimensional
	# scaling method:
	# 	mds.res <- mds.method(phi)
	# 	if (is.matrix(mds.res)) {
	# 		colnames(mds.res) <- c("x", "y")
	# 	} else if (is.data.frame(mds.res)) {
	# 		names(mds.res) <- c("x", "y")
	# 	} else {
	# 		warning("Result of mds.method should be a matrix or data.frame.")
	# 	}
	# 	mds.df <- data.frame(mds.res, topics = seq_len(K), Freq = topic.proportion*100,
	# 											 cluster = 1, stringsAsFactors = FALSE)
	# note: cluster (should?) be deprecated soon.

	# token counts for each term-topic combination (widths of red bars)
	term.topic.frequency <- phi * topic.frequency

	# compute term frequencies as column sums of term.topic.frequency
	# we actually won't use the user-supplied term.frequency vector.
	# the term frequencies won't match the user-supplied frequencies exactly
	# this is a work-around to solve the bug described in Issue #32 on github:
	# https://github.com/cpsievert/LDAvis/issues/32
	term.frequency <- colSums(term.topic.frequency)
	stopifnot(all(term.frequency > 0))

	# marginal distribution over terms (width of blue bars)
	term.proportion <- term.frequency/sum(term.frequency)

	# Old code to adjust term frequencies. Deprecated for now
	# adjust to match term frequencies exactly (get rid of rounding error)
	#err <- as.numeric(term.frequency/colSums(term.topic.frequency))
	# http://stackoverflow.com/questions/3643555/multiply-rows-of-matrix-by-vector
	#term.topic.frequency <- sweep(term.topic.frequency, MARGIN=2, err, `*`)

	# Most operations on phi after this point are across topics
	# R has better facilities for column-wise operations
	phi <- t(phi)

	# compute the distinctiveness and saliency of the terms:
	# this determines the R terms that are displayed when no topic is selected
	topic.given.term <- phi/rowSums(phi)  # (W x K)
	kernel <- topic.given.term * log(sweep(topic.given.term, MARGIN=2,
																				 topic.proportion, `/`))
	distinctiveness <- rowSums(kernel)
	saliency <- term.proportion * distinctiveness

	# Order the terms for the "default" view by decreasing saliency:
	default.terms <- vocab[order(saliency, decreasing = TRUE)][1:R]
	counts <- as.integer(term.frequency[match(default.terms, vocab)])
	Rs <- rev(seq_len(R))
	default <- data.frame(Term = default.terms, logprob = Rs, loglift = Rs,
												Freq = counts, Total = counts, Category = "Default",
												stringsAsFactors = FALSE,lambda=NA)
	topic_seq <- rep(seq_len(K), each = R)
	category <- paste0("Topic", topic_seq)
	lift <- phi/term.proportion

	# Collect R most relevant terms for each topic/lambda combination
	# Note that relevance is re-computed in the browser, so we only need
	# to send each possible term/topic combination to the browser
	find_relevance <- function(i) {
		relevance <- i*log(phi) + (1 - i)*log(lift)
		idx <- apply(relevance, 2,
								 function(x) order(x, decreasing = TRUE)[seq_len(R)])
		# for matrices, we pick out elements by their row/column index
		indices <- cbind(c(idx), topic_seq)
		data.frame(Term = vocab[idx], Category = category,
							 logprob = round(log(phi[indices]), 4),
							 loglift = round(log(lift[indices]), 4),
							 stringsAsFactors = FALSE)
	}
	lambda.seq <- seq(0, 1, by=lambda.step)
	#if (missing(cluster)) {
	tinfo <- lapply(as.list(lambda.seq), function(x) {x<-data.frame(find_relevance(x),lambda=as.character(x));data.frame(x,ord=1:R)})
	#} else {
	#	tinfo <- parallel::parLapply(cluster, as.list(lambda.seq), find_relevance)
	#}
	tinfo <- unique(do.call("rbind", tinfo))
	tinfo$Total <- term.frequency[match(tinfo$Term, vocab)]
	rownames(term.topic.frequency) <- paste0("Topic", seq_len(K))
	colnames(term.topic.frequency) <- vocab
	tinfo$Freq <- term.topic.frequency[as.matrix(tinfo[c("Category", "Term")])]
	#tinfo <- rbind(default, tinfo)
	tinfo$Category<-sub('opic','',tinfo$Category)

	lda2rel<-data.table(tinfo)

	if(save.to.disk) try(save(lda2rel,file=paste(grep(paste(out.dir,'$',sep=''),dir(include.dirs = T,full.names=T,recursive=F,ignore.case=F),value = T)[1],sfn,sep=.Platform$file.sep)))

	lda2rel
}

source2drive2R.f <- function(
	title.search = stop("Specify something to search for")
	,filter.extension = ""
	,makedt=T
	,ignore.trash=T
)
{
	library(RGoogleDrive)
	library(httr)
	library(data.table)
	raw.list <- GET("https://www.googleapis.com/drive/v2/files?maxResults=100000",
									config(token = getOption("drive.auth")))
	parsed.list <- httr::content(raw.list, as = "parsed")
	all.files <- sapply(parsed.list$items, function(x) x$title, simplify = F)
	if(ignore.trash) all.files<-all.files[!sapply(parsed.list$items, function(x) x$labels$trashed, simplify = T)]
	query<-paste(title.search, ".*",filter.extension,"$", sep = "")
	target.files <- sapply(
		query
		, function(x) grep(x, all.files, value = T)
		, simplify = F
	)
	target.w <- sapply(
		query
		, function(x) grep(x,all.files, value = F)
		, simplify = F)
	ret <- list()
	for (i in names(target.w)) for (j in 1:length(target.files[[i]])) {
		filesize <- paste(round(as.integer(parsed.list$items[[target.w[[i]][j]]]$fileSize)/(2^20),
														3), "MB")
		cat(c("\nDownloading ", filesize, " \"", target.files[[i]][j], "\""),
				sep = "")

				t0 <- proc.time()
		ret[[i]][target.files[[i]][j]] <- list(httr::content(GET(parsed.list$items[[target.w[[i]][j]]]$exportLinks$`text/csv`,
																											 config(token = getOption("drive.auth"))), as = "text"))
		t1 <- proc.time()
		cat("\nDownloaded in", round((t1 - t0)[3]/60, 2), "minutes.")
		if (makedt) {
			cat("\nConverting target text file to data.table.")
			ret[[i]][target.files[[i]][j]] <- list(fread(ret[[i]][[target.files[[i]][j]]],verbose = T))
			cat(" Done.")
		}
		attr(ret[[i]][[target.files[[i]][j]]], "source.file.size") <- filesize
	}
	ret
}

azlyrics2ftxt.f<-function(
  query='John Henry' # alphanumeric and spaces only, no funny business
  ,pgmx=10 # maximum number of pages
  ,pause='rpois'
  ,dq='~/john-henry/d/q'
  ,dd='~/john-henry/d/d'
  ,cfso=T
){
  library(data.table)
  library(rvest)
  library(magrittr)
  library(tools)
  
  slug<-paste0('azlyrics2ftxt-',gsub(' ','_',query))
  if(cfso){
    fso<-dir(pattern=slug,recursive = T,full.names = T,include.dirs = T)
    lf<-length(fso)
    if(!lf) stop('No saved output. Set cfso=F to query www.azlyrics.com.')
    if(length(lf)>1) {
      cat('Read which file?')
      print(data.frame(file=fso))
      fc<-''
      while(!fc%in%as.character(1:lf)) {cat('\nType index number and press -return-');fc<-readLines(n=1)}
      fso<-fso[as.integer(fc)]
      cat('\nLoading',fso)
      load(fso)
    } else {cat('\nLoading',fso);load(fso)}
  } else {  
    #list for return
    ret<-list()
    
    #get pages
    for(i in 1:pgmx){
      h<-paste0('http://search.azlyrics.com/search.php?q=',gsub(' +','+',query),'&p=',i,'&w=songs')
      if(i==1) cat('\nScraping',h,'\n')
      s<-read_html(h) %>%
        html_nodes(xpath="//a") %>% html_attr('href') %>% grep(pattern='/lyrics/',value=T)
      if(!length(s)) break
      for(j in s){
        if(pause=='rpois') Sys.sleep(1+rpois(1,50)/100)
        t<-read_html(j)
        cat('.')
        ret[[length(ret)+1]]<-data.table(
          url=j
          ,title=t %>% html_nodes(xpath="//b") %>% html_text %>% extract(2) %>% sub('^.(.+).$','\\1',.) %>% ifelse(length(.),.,NA_character_)
          ,artist=t %>% html_nodes(xpath="//b") %>% html_text %>% extract(1) %>% tolower %>% toTitleCase %>% sub(' [^ ]+$','',.) %>% ifelse(length(.),.,NA_character_)
          ,featuring=t %>% html_nodes(xpath='/html/body/div[3]/div/div[2]/span') %>% html_text %>% sub('^[^ ]+ (.+).$','\\1',.) %>% ifelse(length(.),.,NA_character_)
          ,album=t %>% html_nodes(xpath='/html/body/div[3]/div/div[2]/div[12]/a') %>% html_text %>% sub('^\"(.+)+\" \\(([0-9]+)\\)$','\\1',.) %>% ifelse(length(.),.,NA_character_)
          ,year=t %>% html_nodes(xpath='/html/body/div[3]/div/div[2]/div[12]/a') %>% html_text %>% sub('^\"(.+)+\" \\(([0-9]+)\\)$','\\2',.) %>% as.integer %>% ifelse(length(.),.,NA_integer_)
          ,lyrics=t %>% html_nodes(xpath='/html/body/div[3]/div/div[2]/div[6]') %>% as.character %>% gsub('<br><br>','\n',.) %>% gsub('\n*<[^>]*>','',.) %>% sub('^[\r\n]+','',.) %>% ifelse(length(.),.,NA_character_)
        )
      }
    }
    ret<-rbindlist(ret)
    setattr(ret,'url',h)
    setattr(ret,'date.queried',Sys.Date())
    azl<-ret
    rm(ret)
  }
  pfs<-.Platform$file.sep
  if(!file.exists(dq)) dq<-getwd()
  dqf<-paste0(dq,pfs,slug,'-n',nrow(azl),'.RData') %>% gsub(paste0(pfs,'+'),pfs,.)
  if(file.exists(dd)) {
    dd<-paste0(dd,pfs) %>% gsub(paste0(pfs,'+'),pfs,.)
    cat('\nDumping text to',dd)
    azl[,names:=mapply(function(t,a,y) paste0(na.omit(c(y,t,a)) %>% paste(collapse='-') %>% gsub('[^A-Za-z0-9-]','',.)),t=title,a=artist,y=year)]
    azl[,mapply(function(l,n) writeLines(l,paste0(dd,n,'.txt')),l=lyrics,n=names)]
  } else {cat('\nNo data dump. Set cfso=T to load saved query and then specify a dd directory.')}
  if(!cfso) {
    cat('\nSaving',dqf)
    try(save(azl,file=dqf))
  }
  azl
}
