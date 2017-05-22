# R code for finding aberrantly-compartmentalized regions
require(GenomicRanges)
require(minfi)
require(biomaRt)

# construct a table of gene annotations
ensembl <- useMart(host="feb2014.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
geneInfo <-getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position","end_position","strand"),mart=ensembl)
geneInfo <- geneInfo[which(!geneInfo$hgnc_symbol==""),]
geneInfo <- geneInfo[!duplicated(geneInfo$hgnc_symbol),]
autosomes <- paste("chr",1:22,sep="")
geneInfo <- geneInfo[which(paste("chr",geneInfo$chromosome_name,sep="") %in% autosomes),]
geneInfo.gr <- GRanges(seqnames=paste("chr",geneInfo$chromosome_name,sep=""),ranges=IRanges(start=geneInfo$start_position,end=geneInfo$end_position),symbol=geneInfo$hgnc_symbol)

# load methylation data and call compartments: assumes all idat files are in working directory

basenames.t <- as.character(read.table("PRAD_Tumour_MethFilenames.txt")[[1]])
RGset.t <- read.450k(basenames.t)
MSet.t <- preprocessFunnorm(RGset.t)
MSet.t <- dropLociWithSnps(MSet.t,maf=0.01)
compartments.t <- lapply(autosomes,function(x)compartments(MSet.t, keep=FALSE, chr=x))
allcompartments.t <- compartments.t[[1]]
for(i in 2:length(compartments.t)){
        allcompartments.t <- c(allcompartments.t,compartments.t[[i]])
}

basenames.n <- as.character(read.table("PRAD_Normal_MethFilenames.txt")[[1]])
RGset.n <- read.450k(basenames.n)
MSet.n <- preprocessFunnorm(RGset.n)
MSet.n <- dropLociWithSnps(MSet.n,maf=0.01)
compartments.n <- lapply(autosomes,function(x)compartments(MSet.n, keep=FALSE, chr=x))
allcompartments.n <- compartments.n[[1]]
for(i in 2:length(compartments.n)){
        allcompartments.n <- c(allcompartments.n,compartments.n[[i]])
}

# alternatively, load pre-computed compartment objects
allcompartments.t <- load("PRADcompartmentsT.rda")
allcompartments.n <- load("PRADcompartmentsN.rda")

# load DNAse-, H3K27ac-IP-, H3K27ac-IP-seq data and compare open/closed compartments

require(rtracklayer)
lncap.dnase.gr <- import("LNCaP_ENCFF752YDY.bigWig")
allwindows.DHSscore <- sapply(allcompartments.t,function(x)sum(subsetByOverlaps(query=lncap.dnase.gr,subject=x)$score))
par(mfrow=c(1,2))
boxplot(list(open=allwindows.DHSscore[which(allcompartments.n$compartment=="open")],closed=allwindows.DHSscore[which(allcompartments.n$compartment=="closed")]),main="Normal Prostate Genome Compartments",ylab="LNCaP DNase-seq \n average normalized coverage")
boxplot(list(open=allwindows.DHSscore[which(allcompartments.t$compartment=="open")],closed=allwindows.DHSscore[which(allcompartments.t$compartment=="closed")]),main="PRAD Tumours Genome Compartments",ylab="LNCaP DNase-seq \n average normalized coverage")

h3k27ac.gff <- read.table("GSM1249448_LNCaP_+DHT_H3K27Ac_signifpeaks.gff",sep="\t",head=F)
h3k27ac.gr <- GRanges(seqnames=as.character(h3k27ac.gff[[1]]),ranges=IRanges(start=as.numeric(h3k27ac.gff[[4]]),end=as.numeric(h3k27ac.gff[[5]])),count=as.numeric(h3k27ac.gff[[6]]))
h3k27ac.cov <- coverage(h3k27ac.gr)
h3k27ac.cov <- h3k27ac.cov[seqlevels(allcompartments.t)]
h3k27ac.tAve <- binnedAverage(bins=allcompartments.t,numvar=h3k27ac.cov,varname="cov")
h3k27ac.nAve <- binnedAverage(bins=allcompartments.n,numvar=h3k27ac.cov,varname="cov")
par(mfrow=c(1,2))
boxplot(list(open=h3k27ac.nAve$cov[which(h3k27ac.nAve$compartment=="open")],closed=h3k27ac.nAve$cov[which(h3k27ac.nAve$compartment=="closed")]),main="Normal Prostate Genome Compartments",ylab="LNCaP H3K27ac ChIP-seq \n proportion of peak coverage")
boxplot(list(open=h3k27ac.tAve$cov[which(h3k27ac.tAve$compartment=="open")],closed=h3k27ac.tAve$cov[which(h3k27ac.tAve$compartment=="closed")]),main="Tumour Prostate Genome Compartments",ylab="LNCaP H3K27ac ChIP-seq \n proportion of peak coverage")

h3k27me3.bed <- read.table("GSE86532_Ctrl_MACS2_summits.bed",sep="\t",head=F)
h3k27me3.gr <- GRanges(seqnames=as.character(h3k27me3.bed[[1]]),ranges=IRanges(start=as.numeric(h3k27me3.bed[[2]]),end=as.numeric(h3k27me3.bed[[3]])),count=as.numeric(h3k27me3.bed[[5]]))
sum(countOverlaps(query=h3k27me3.gr,subject=allcompartments.t[allcompartments.t$compartment=="open"]))/sum(allcompartments.t$compartment=="open")
sum(countOverlaps(query=h3k27me3.gr,subject=allcompartments.t[allcompartments.t$compartment=="closed"]))/sum(allcompartments.t$compartment=="closed")
sum(countOverlaps(query=h3k27me3.gr,subject=allcompartments.n[allcompartments.n$compartment=="open"]))/sum(allcompartments.n$compartment=="open")
sum(countOverlaps(query=h3k27me3.gr,subject=allcompartments.n[allcompartments.n$compartment=="closed"]))/sum(allcompartments.n$compartment=="closed")

# compute smoothed window scores for plotting

smoothed.matrix.t <- cbind(as.numeric(allcompartments.t$pc)[1:(length(allcompartments.t$pc)-2)],as.numeric(allcompartments.t$pc)[2:(length(allcompartments.t$pc)-1)],as.numeric(allcompartments.t$pc)[3:length(allcompartments.t$pc)])
smoothed.matrix.t <- apply(smoothed.matrix.t,2,function(x)(2*(x-min(x))/(max(x)-min(x)))-1)
for(i in 1:ncol(smoothed.matrix.t)){
        smoothed.matrix.t[which(smoothed.matrix.t[,i]>0),i] <- (smoothed.matrix.t[which(smoothed.matrix.t[,i]>0),i]-min(smoothed.matrix.t[which(smoothed.matrix.t[,i]>0),i]))/(max(smoothed.matrix.t[which(smoothed.matrix.t[,i]>0),i])-min(smoothed.matrix.t[which(smoothed.matrix.t[,i]>0),i]))
        smoothed.matrix.t[which(smoothed.matrix.t[,i]<0),i] <- -(abs(smoothed.matrix.t[which(smoothed.matrix.t[,i]<0),i])-min(abs(smoothed.matrix.t[which(smoothed.matrix.t[,i]<0),i])))/(max(abs(smoothed.matrix.t[which(smoothed.matrix.t[,i]<0),i]))-min(abs(smoothed.matrix.t[which(smoothed.matrix.t[,i]<0),i])))
}
smoothed.means.t <- rowMeans(smoothed.matrix.t)
smoothed.matrix.n <- cbind(as.numeric(allcompartments.n$pc)[1:(length(allcompartments.n$pc)-2)],as.numeric(allcompartments.n$pc)[2:(length(allcompartments.n$pc)-1)],as.numeric(allcompartments.n$pc)[3:length(allcompartments.n$pc)])
smoothed.matrix.n <- apply(smoothed.matrix.n,2,function(x)(2*(x-min(x))/(max(x)-min(x)))-1)
for(i in 1:ncol(smoothed.matrix.n)){
        smoothed.matrix.n[which(smoothed.matrix.n[,i]>0),i] <- (smoothed.matrix.n[which(smoothed.matrix.n[,i]>0),i]-min(smoothed.matrix.n[which(smoothed.matrix.n[,i]>0),i]))/(max(smoothed.matrix.n[which(smoothed.matrix.n[,i]>0),i])-min(smoothed.matrix.n[which(smoothed.matrix.n[,i]>0),i]))
        smoothed.matrix.n[which(smoothed.matrix.n[,i]<0),i] <- -(abs(smoothed.matrix.n[which(smoothed.matrix.n[,i]<0),i])-min(abs(smoothed.matrix.n[which(smoothed.matrix.n[,i]<0),i])))/(max(abs(smoothed.matrix.n[which(smoothed.matrix.n[,i]<0),i]))-min(abs(smoothed.matrix.n[which(smoothed.matrix.n[,i]<0),i])))
}
smoothed.means.n <- rowMeans(smoothed.matrix.n)

par(mfrow=c(2,1))
plot(0,xlim=c(1,length(smoothed.means.t)+1),ylim=range(smoothed.means.t),col="white",xlab="genomic window rank",ylab="closed-ness score",main="PRAD tumour compartments")
rect(xleft=c(1:length(smoothed.means.t)),xright=c(1:length(smoothed.means.t))+1,ybottom=smoothed.means.t*(1-sign(smoothed.means.t))/2,ytop=smoothed.means.t*(1+sign(smoothed.means.t))/2,density=-1,col="black",border=NA)
abline(v=which(!duplicated(as.character(seqnames(allcompartments.t)))),col="black",lty=2)
plot(0,xlim=c(1,length(smoothed.means.t)+1),ylim=range(smoothed.means.n),col="white",xlab="genomic window rank",ylab="closed-ness score",main="normal prostate compartments")
rect(xleft=c(1:length(smoothed.means.n)),xright=c(1:length(smoothed.means.n))+1,ybottom=smoothed.means.n*(1-sign(smoothed.means.n))/2,ytop=smoothed.means.n*(1+sign(smoothed.means.n))/2,density=-1,col="black",border=NA)
abline(v=which(!duplicated(as.character(seqnames(allcompartments.t)))),col="black",lty=2)

# find regions with different compartment calls across the two cohorts (tumour vs normal)

getWindowSelectionVectors <- function(Tstatus="open",Nstatus="closed",pc.diff=0.1){
        tswitchopen.which = (allcompartments.t$compartment==Tstatus & allcompartments.n$compartment==Nstatus & abs(allcompartments.t$pc)>pc.diff & abs(allcompartments.n$pc)>pc.diff)
        tswitchclosed.which = (allcompartments.t$compartment==Nstatus & allcompartments.n$compartment==Tstatus & abs(allcompartments.t$pc)>pc.diff & abs(allcompartments.n$pc)>pc.diff)
        list(open=tswitchopen.which,closed=tswitchclosed.which)
}

tswitchwindows <- getWindowSelectionVectors()

# load Copy-Number data from TCGA and compare CNV rates for aberrantly-compartmentalized genes to all other genes
tcga.cn <- as.matrix(read.table("TCGA_PRAD_Gistic2ThreshCN.txt",sep="\t",head=T,row.names=1))
tswitchopen.genes <- unique(as.character(subsetByOverlaps(query=geneInfo.gr,subject=allcompartments.t[which(tswitchwindows$open)])$symbol))
tswitchclosed.genes <- unique(as.character(subsetByOverlaps(query=geneInfo.gr,subject=allcompartments.t[which(tswitchwindows$closed)])$symbol))
all.other.genes <- unique(as.character(subsetByOverlaps(query=geneInfo.gr,subject=allcompartments.t[which(!tswitchwindows$open)])$symbol))
cncomp.mat <- rbind(table(tcga.cn[intersect(tswitchopen.genes,rownames(tcga.cn)),]),table(tcga.cn[intersect(all.other.genes,rownames(tcga.cn)),]))

# testing pathway enrichments in aberrantly-compartmentalized genes
cpdb.pathways = read.table("CPDB_pathways_genes_2015.txt",sep="\t",head=T)
cpdb.genes = sapply(as.character(cpdb.pathways$hgnc_symbol_ids),function(x)strsplit(x,split=",",fixed=T)[[1]])
names(cpdb.genes) = as.character(cpdb.pathways$pathway)
cpdb.universe = setdiff(as.character(geneInfo.gr$symbol),c("",NA))
allpathways.open.pvals = sapply(cpdb.genes,function(x)1-phyper(length(intersect(x,tswitchopen.genes))-1,length(tswitchopen.genes),length(setdiff(cpdb.universe,tswitchopen.genes)),length(intersect(x,cpdb.universe))))
allpathways.closed.pvals = sapply(cpdb.genes,function(x)1-phyper(length(intersect(x,tswitchclosed.genes))-1,length(tswitchclosed.genes),length(setdiff(cpdb.universe,tswitchclosed.genes)),length(intersect(x,cpdb.universe))))
allpathways.ngenes <- sapply(cpdb.genes,length)
allpathways.opengenes <- sapply(cpdb.genes,function(x)length(intersect(x,tswitchopen.genes)))
allpathways.closedgenes <- sapply(cpdb.genes,function(x)length(intersect(x,tswitchopen.genes)))
cpdbtest.df = data.frame(pathway=names(cpdb.genes),nGenes=allpathways.ngenes,nTumourOpenGenes=allpathways.opengenes,open.p=allpathways.open.pvals,open.adj.p=p.adjust(allpathways.open.pvals,method="fdr"),nTumourClosedGenes=allpathways.closedgenes,closed.p=allpathways.closed.pvals,closed.adj.p=p.adjust(allpathways.closed.pvals,method="fdr"),TumourOpenGenes=sapply(cpdb.genes,function(x)paste(intersect(x,tswitchopen.genes),collapse=",")),TumourClosedGenes=sapply(cpdb.genes,function(x)paste(intersect(x,tswitchclosed.genes),collapse=",")))
cpdbtest.df = cpdbtest.df[order(cpdbtest.df$open.p,decreasing=F),]

# computing Observed/Expected peak overlaps for aberrantly-compartmentalized regions
# requires ENCODE ChIP-seq narrowpeak files to be in a folder called 'ENCODE', and AR ChIP-seq files to be downloaded and put into a folder called 'GSE65478_AR_ChIPseq'

peak.counts <- read.table("ENCODE_TFpeakCounts.txt", sep="\t",head=TRUE)
setwd("GSE65478_AR_ChIPseq")

ar.bedfiles <- c(list.files(pattern="wz25"),list.files(pattern="wz22"))
ar.bed <- lapply(ar.bedfiles,read.table,sep="\t",head=F)
ar.gr <- lapply(ar.bed,function(x)GRanges(seqnames=paste("chr",as.character(x[[1]]),sep=""),ranges=IRanges(start=x[[2]],end=x[[3]])))
allar.gr <- ar.gr[[1]]
for(i in 2:length(ar.gr)){
        allar.gr <- c(allar.gr,ar.gr[[i]])
}
setwd("..")

makeTFmatrix <- function(upstream=1000,downstream=1000){
        tfchip.matrix <- matrix(0,nrow=nrow(geneInfo),ncol=nrow(peak.counts)+1)
        rownames(tfchip.matrix) <- geneInfo$hgnc_symbol
        colnames(tfchip.matrix) <- c(as.character(peak.counts$TF),"AR.GSE65478")

        promoters <- geneInfo
        promoters[promoters$strand=="-1","start_position"] <- promoters[promoters$strand=="-1","end_position"]-upstream
        promoters[promoters$strand=="-1","end_position"] <- promoters[promoters$strand=="-1","end_position"]+upstream
        promoters[promoters$strand=="1","end_position"] <- promoters[promoters$strand=="1","start_position"]+downstream
        promoters[promoters$strand=="1","start_position"] <- promoters[promoters$strand=="1","start_position"]-downstream
        promoters.gr <- GRanges(seqnames=paste("chr",promoters$chromosome_name,sep=""),ranges=IRanges(start=promoters$start_position,end=promoters$end_position),strand=promoters$strand,symbol=promoters$hgnc_symbol)
        for(i in 1:nrow(peak.counts)){
                cat(paste("calculating overlaps for file",i,"of",ncol(tfchip.matrix),"\n"))
                this.chip <- read.table(paste("ENCODE/",peak.counts$file[i],sep=""),head=FALSE,sep="\t")
                chip.gr <- GRanges(seqnames=as.character(this.chip[[1]]),ranges=IRanges(start=this.chip[[2]],end=this.chip[[3]]))
                tfchip.matrix[subsetByOverlaps(query=promoters.gr,subject=chip.gr)$symbol,i] <- 1
        }
        tfchip.matrix[subsetByOverlaps(query=promoters.gr,subject=allar.gr)$symbol,ncol(tfchip.matrix)] <- 1
        tfchip.matrix
}
tfchip.matrix <- makeTFmatrix()

getTFObsvExp <- function(selectedOpen,selectedClosed){
        allchip.tswitchopen.FC <- rep(NA,ncol(tfchip.matrix))
        allchip.tswitchclosed.FC <- rep(NA,ncol(tfchip.matrix))
        allchip.tswitchopen.chisqp <- rep(NA,ncol(tfchip.matrix))
        allchip.tswitchclosed.chisqp <- rep(NA,ncol(tfchip.matrix))

        for(i in 1:nrow(peak.counts)){
                cat(paste("calculating enrichments for TF",i,"of",ncol(tfchip.matrix),"\n"))
                this.chip <- read.table(paste("ENCODE/",peak.counts$file[i],sep=""),head=FALSE,sep="\t")
                chip.gr <- GRanges(seqnames=as.character(this.chip[[1]]),ranges=IRanges(start=this.chip[[2]],end=this.chip[[3]]))
                allchip.tswitchopen.FC[i] = (sum(as.numeric(width(subsetByOverlaps(subject=allcompartments.t[which(selectedOpen)],query=chip.gr))))/sum(as.numeric(width(subsetByOverlaps(subject=allcompartments.t[which(!selectedOpen)],query=chip.gr)))))/(sum(as.numeric(width(allcompartments.t[which(selectedOpen)])))/sum(as.numeric(width(allcompartments.t[which(!selectedOpen)]))))
                allchip.tswitchclosed.FC[i] = (sum(as.numeric(width(subsetByOverlaps(subject=allcompartments.t[which(selectedClosed)],query=chip.gr))))/sum(as.numeric(width(subsetByOverlaps(subject=allcompartments.t[which(!selectedClosed)],query=chip.gr)))))/(sum(as.numeric(width(allcompartments.t[which(selectedClosed)])))/sum(as.numeric(width(allcompartments.t[which(!selectedClosed)]))))
                cat(paste("calculating Chi-squared p-values for TF",i,"of",ncol(tfchip.matrix),"\n"))
                allchip.tswitchclosed.chisqp[i] = chisq.test(matrix(c(length(subsetByOverlaps(subject=allcompartments.t[which(selectedClosed)],query=chip.gr)),length(subsetByOverlaps(subject=allcompartments.t[which(!selectedClosed)],query=chip.gr)),length(which(selectedClosed)),length(which(!selectedClosed))),nrow=2,ncol=2))$p.value
                allchip.tswitchopen.chisqp[i] = chisq.test(matrix(c(length(subsetByOverlaps(subject=allcompartments.t[which(selectedOpen)],query=chip.gr)),length(subsetByOverlaps(subject=allcompartments.t[which(!selectedOpen)],query=chip.gr)),length(which(selectedOpen)),length(which(!selectedOpen))),nrow=2,ncol=2))$p.value
        }
        allchip.tswitchopen.FC[ncol(tfchip.matrix)] = (sum(as.numeric(width(subsetByOverlaps(subject=allcompartments.t[which(selectedOpen)],query=allar.gr))))/sum(as.numeric(width(subsetByOverlaps(subject=allcompartments.t[which(!selectedOpen)],query=allar.gr)))))/(sum(as.numeric(width(allcompartments.t[which(selectedOpen)])))/sum(as.numeric(width(allcompartments.t[which(!selectedOpen)]))))
        allchip.tswitchclosed.FC[ncol(tfchip.matrix)] = (sum(as.numeric(width(subsetByOverlaps(subject=allcompartments.t[which(selectedClosed)],query=allar.gr))))/sum(as.numeric(width(subsetByOverlaps(subject=allcompartments.t[which(!selectedClosed)],query=allar.gr)))))/(sum(as.numeric(width(allcompartments.t[which(selectedClosed)])))/sum(as.numeric(width(allcompartments.t[which(!selectedClosed)]))))
        allchip.tswitchopen.chisqp[ncol(tfchip.matrix)] = chisq.test(matrix(c(length(subsetByOverlaps(subject=allcompartments.t[which(selectedOpen)],query=allar.gr)),length(subsetByOverlaps(subject=allcompartments.t[which(!selectedOpen)],query=allar.gr)),length(which(selectedOpen)),length(which(!selectedOpen))),nrow=2,ncol=2))$p.value
        allchip.tswitchclosed.chisqp[ncol(tfchip.matrix)] = chisq.test(matrix(c(length(subsetByOverlaps(subject=allcompartments.t[which(selectedClosed)],query=allar.gr)),length(subsetByOverlaps(subject=allcompartments.t[which(!selectedClosed)],query=allar.gr)),length(which(selectedClosed)),length(which(!selectedClosed))),nrow=2,ncol=2))$p.value

        tswitch.FC.df <- data.frame(ChIP=colnames(tfchip.matrix),T.open.log.FC=log(allchip.tswitchopen.FC,base=2),T.closed.log.FC=log(allchip.tswitchclosed.FC,base=2),T.open.p=allchip.tswitchopen.chisqp,T.closed.p=allchip.tswitchclosed.chisqp,T.open.adj.p=p.adjust(allchip.tswitchopen.chisqp,method="BH"),T.closed.adj.p=p.adjust(allchip.tswitchclosed.chisqp,method="BH"))
        tswitch.FC.df
}

tswitchFCs <- getTFObsvExp(selectedOpen=tswitchwindows$open,selectedClosed=tswitchwindows$closed)
filtereddf <- tswitchFCs[-c(grep("Ctcf",tswitchFCs$ChIP),grep("Pol2",tswitchFCs$ChIP)),]

# compute pair-wise Jaccard distances for TF ChIP bed files

tfbs.chip = lapply(peak.counts$file,function(x)read.table(paste("ENCODE/",x,sep=""),head=FALSE,sep="\t"))
tfbs.peaks = lapply(tfbs.chip,function(x)GRanges(seqnames=as.character(x[[1]]),ranges=IRanges(start=x[[2]],end=x[[3]])))
tfbs.overlaps = lapply(tfbs.peaks,function(x)subsetByOverlaps(query=allcompartments.t[which(tswitchwindows$open)],subject=x))
tfbs.overlaps = c(tfbs.overlaps,list(subsetByOverlaps(query=allcompartments.t[which(tswitchwindows$open)],subject=allar.gr)))
jdist = function(gr1,gr2){
        length(subsetByOverlaps(query=gr1,subject=gr2))/length(unique(c(gr1,gr2)))
}
pairdists = array(0,dim=rep(length(tfbs.overlaps),2))
for(i in 1:length(tfbs.overlaps)){
        cat(paste("screening dists for ChIP",i,"of",length(tfbs.overlaps),"\n"))
        for(j in 1:length(tfbs.overlaps)){
                pairdists[i,j] = jdist(tfbs.overlaps[[i]],tfbs.overlaps[[j]])
        }
}
rownames(pairdists) = colnames(tfchip.matrix)
colnames(pairdists) = colnames(tfchip.matrix)

# map TF binding sites to genes, load TCGA gene expression data, compute systematic overexpression in PRAD tumours relative to normal tissues

tcga.gx <- as.matrix(read.table("TCGA_PRAD_exp_HiSeqV2-2015-02-24/genomicMatrix",sep="\t",head=T,row.names=1))
tcga.clin <- read.table("TCGA_PRAD_exp_HiSeqV2-2015-02-24/clinical_data",sep="\t",head=T)
rownames(tcga.clin) = gsub(tcga.clin[[1]],pattern="-",replace=".")
sampleIDs = rownames(tcga.clin)[which(rownames(tcga.clin) %in% colnames(tcga.gx) & tcga.clin$sample_type %in% c("Primary Tumor","Solid Tissue Normal"))]

require(limma)
gxdiff.design <- cbind(intercept=1,tumour=as.numeric(tcga.clin[sampleIDs,"sample_type"]=="Primary Tumor"))
gxdiff.fit = lmFit(tcga.gx[,sampleIDs],design=gxdiff.design)
gxdiff.fit = eBayes(gxdiff.fit)
gxdiff.table = topTable(gxdiff.fit,coef=2,number=nrow(tcga.gx))
gxdiff.t = gxdiff.table$t
gxdiff.genes = rownames(gxdiff.table)

getTvNgxstats <- function(tfchip.matrix){
        tfgx.ids = intersect(rownames(tfchip.matrix),rownames(gxdiff.table))
        tf.Toverexpress = apply(tfchip.matrix[tfgx.ids,],2,function(x)geneSetTest(statistics=gxdiff.table[tfgx.ids,"t"],index=which(x==1),alternative="up"))
        tf.Tunderexpress = apply(tfchip.matrix[tfgx.ids,],2,function(x)geneSetTest(statistics=gxdiff.table[tfgx.ids,"t"],index=which(x==1),alternative="down"))
        tf.TvNgx.df = data.frame(ChIP=colnames(tfchip.matrix),targets.TvN.p=tf.Toverexpress,targets.NvT.p=tf.Tunderexpress,Thi.adj.p=p.adjust(tf.Toverexpress,method="fdr"),Tlo.adj.p=p.adjust(tf.Tunderexpress,method="fdr"))
        tf.TvNgx.df
}

TvNgx <- getTvNgxstats(tfchip.matrix)

TFtarget.df = cbind(filtereddf,TvNgx[as.character(filtereddf$ChIP),-1])

plot(x=gxdiff.table$logFC,y=-log(gxdiff.table$P.Value,base=10),xlab="Tumour v Normal median log2 Fold-Change",ylab="Tumour v Normal t-test -log10 p-value",main="PRAD vs Normal GX",col="grey",ylim=c(0,40),xlim=c(-4,4))
text(x=gxdiff.table[open.genes[open.genes %in% rownames(gxdiff.table)],"logFC"],y=-log(gxdiff.table[open.genes[open.genes %in% rownames(gxdiff.table)],"P.Value"],base=10),labels=open.genes[open.genes %in% rownames(gxdiff.table)],col="red")
text(x=gxdiff.table[closed.genes[closed.genes %in% rownames(gxdiff.table)],"logFC"],y=-log(gxdiff.table[closed.genes[closed.genes %in% rownames(gxdiff.table)],"P.Value"],base=10),labels=closed.genes[closed.genes %in% rownames(gxdiff.table)],col="blue")

