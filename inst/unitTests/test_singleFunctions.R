test_plotFunctions <-
    function()
    {
        data(isomiRexp)
        checkTrue(class(isoPlot(isomiRexp,type="t5"))[1]=="IsomirDataSeq")
        checkTrue(class(isoPlot(isomiRexp,type="t3"))[1]=="IsomirDataSeq")
        checkTrue(class(isoPlot(isomiRexp,type="add"))[1]=="IsomirDataSeq")
    }

test_dseFunctions <-
    function()
    {
        data(isomiRexp)
        checkTrue(class(isoDE(isomiRexp, formula=~condition, ref=TRUE,
                            iso5=TRUE))[1]=="DESeqDataSet")
    }

test_countsFunctions <-
    function()
    {
        data(isomiRexp)
        library(GenomicRanges)
        obj<-isoCounts(isomiRexp, ref=TRUE)
        checkTrue(class(counts(obj))[1]=="matrix")
        obj<-isoNorm(obj)
        checkTrue(class(counts(obj, norm=TRUE))[1]=="matrix")
    }