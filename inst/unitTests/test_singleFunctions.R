test_plotFunctions <-
    function()
    {
        data(mirData)
        checkTrue(class(isoPlot(mirData,type="iso5"))[1]=="IsomirDataSeq")
        checkTrue(class(isoPlot(mirData,type="iso3"))[1]=="IsomirDataSeq")
        checkTrue(class(isoPlot(mirData,type="add"))[1]=="IsomirDataSeq")
    }

test_dseFunctions <-
    function()
    {
        data(mirData)
        checkTrue(class(isoDE(mirData, formula=~condition, ref=TRUE,
                            iso5=TRUE))[1]=="DESeqDataSet")
    }

test_countsFunctions <-
    function()
    {
        data(mirData)
        library(GenomicRanges)
        obj<-isoCounts(mirData, ref=TRUE)
        checkTrue(class(counts(obj))[1]=="matrix")
        obj<-isoNorm(obj)
        checkTrue(class(counts(obj, norm=TRUE))[1]=="matrix")
    }