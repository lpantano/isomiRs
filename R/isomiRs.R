#' Do differential expression analysis with DESeq2
#'
#' This function does differential expression analysis with DESeq2.
#' It will return a DESeq2 object.
#' You can merge all isomiRs into miRNA by calling the function only
#' with the frist two parameters. You can get a table with isomiRs and
#' the reference miRBase sequence by calling the function with ref=T.
#' You can get a table with 5' trimming isomiRS, miRBase reference and
#' the rest by calling with ref=T,iso5=T.
#' If you set up all parameters to TRUE, you will get a table for
#' each differnt sequence mapping to a miRNA
#'
#' @param x object isomiDataSeq
#' @param formula used for DE analysis
#' @param ref differenciate reference miRNA from rest
#' @param iso5 differenciate trimming at 5 miRNA from rest
#' @param iso3 differenciate trimming at 3 miRNA from rest
#' @param add differenciate additions miRNA from rest
#' @param mism differenciate nt substitution miRNA from rest
#' @param seed differenciate changes in 2-7 nt from rest
#' @return DESeq object
#' @examples
#' library(DESeq2)
#' data(isomiRexp)
#' dds<-deIso(isomiRexp, formula=~condition)
#' @export
#' @import DESeq2
deIso<-function(x,formula,ref=FALSE,iso5=FALSE,iso3=FALSE,
                add=FALSE,mism=FALSE,seed=FALSE)
{
    if (ref | iso5 | iso3 | add | mism | seed){
        x<-countsIso(x,ref,iso5,iso3,add,mism,seed)
    }
    countData<-counts(x)
    dds<-DESeqDataSetFromMatrix(countData = countData,
                                colData = colData(x),
                                design = formula)
    dds <- DESeq(dds,quiet=TRUE)
    dds
}

#' Plot the amount of isomiRs in different samples divided by group
#'
#' @param x object isomiDataSeq
#' @param top number of isomiRs used
#' @export
#' @import ggplot2
#' @import gplots
#' @import RColorBrewer
plotTop<-function(x,top=20)
{
    # dds<-counts(x)
    # rld <- rlogTransformation(dds)
    select <- order(rowMeans(counts(x)),
                    decreasing=TRUE)[1:top]
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    heatmap.2(counts(dds,normalized=TRUE), col = hmcol,
            scale="none", 
            dendrogram="none", trace="none")
}

#' Plot the amount of isomiRs in different samples divided by group
#'
#' @param x object isomirDataSeq
#' @param type string (t5,t3,add,sub) to indicate what isomiR
#' change to use for the plot
#' @return ggplot2 figure showing the selected isomiR changes among samples
#' @export
#' @import ggplot2
#' @examples
#' data(isomiRexp)
#' plotIso(isomiRexp)
plotIso<-function(x, type="t5")
{
    freq=size=group=abundance=NULL
    codevn<-c(2,3,4,5)
    names(codevn)<-c("t5","t3","sub","add")
    ratiov<-c(1/6,1/6,1/23,1/3)
    names(ratiov)<-names(codevn)
    coden<-codevn[type]
    ratio<-ratiov[type]
    des<-colData(x)
    table<-data.frame()
    isoList <- isoinfo(x)
    for (sample in row.names(des)){
        
        uniq.dat <- as.data.frame( table(isoList[[sample]][[coden]]$size) )
        temp <- as.data.frame( isoList[[sample]][[coden]] %>%
                                group_by(size) %>%
                                summarise(freq=sum(freq)) )
        total <- sum(temp$freq)
        temp <- merge(temp,uniq.dat,by=1)
        Total <- sum(temp$Freq)
        temp$abundance <- temp$freq/total
        temp$unique <- temp$Freq/Total
        table <- rbind( table,data.frame( size=temp$size,abundance=temp$abundance,
                                        unique=temp$unique,
                                        sample=rep(sample,nrow(temp)),
                                        group=rep(des[sample,"condition"],
                                                nrow(temp)) ) )
    }
    isostats(x)[[type]]<-table
    p <- ggplot(table)+
        geom_jitter(aes(x=factor(size),y=unique,colour=factor(group),
                        size=abundance))+
        scale_colour_brewer("Groups",palette="Set1")+
        theme_bw(base_size = 14, base_family = "") +
        theme(strip.background=element_rect(fill="slategray3"))+
        labs(list(title=paste(type,"distribution"),y="# of unique sequences",
                x="position respect to the reference"))
    print(p)
    x
}

#' create count tables from isomirs
#' @param x IsomirDataSeq class
#' @param ref differenciate reference miRNA from rest
#' @param iso5 differenciate trimming at 5 miRNA from rest
#' @param iso3 differenciate trimming at 3 miRNA from rest
#' @param add differenciate additions miRNA from rest
#' @param mism differenciate nt substitution miRNA from rest
#' @param seed differenciate changes in 2-7 nt from rest
#' @return count table
#' @examples
#' library(DESeq2)
#' data(isomiRexp)
#' ma<-countsIso(isomiRexp)
#' @export
countsIso <- function(x, ref=FALSE,iso5=FALSE,iso3=FALSE,
                      add=FALSE,mism=FALSE,seed=FALSE)
    {
        counts <- IsoCountsFromMatrix(isoraw(x), colData(x), ref,
                                      iso5, iso3,
                                      add, mism, seed)
        se <- SummarizedExperiment(assays = SimpleList(counts=counts), 
                                   colData = colData(x))
        x <- IsomirDataSeq(se, isoraw(x), isoinfo(x), isostats(x))
        return(x)
    }


#' normalize count data
#'
#' @param x IsomirDataSeq object
#' @param formula formula that will be used for DE
#' library(DESeq2)
#' data(isomiRexp)
#' ma<-normIso(isomiRex, formula=~condition)
#' @export
#' @import DESeq2
normIso<-function(x,formula=~condition)
{
    dds<-DESeqDataSetFromMatrix(countData = counts(x),
                                colData = colData(x),
                                design = formula)
    rld<-rlogTransformation(dds,blind=FALSE)
    normcounts(x) <- assay(rld)
    x
}

