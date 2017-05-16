test_plotFunctions <-
    function()
    {
        data(mirData)
        checkTrue(class(isoPlot(mirData,type="iso5",column ="group"))[2]=="ggplot")
        checkTrue(class(isoPlot(mirData,type="iso3",column ="group"))[2]=="ggplot")
        checkTrue(class(isoPlot(mirData,type="add",column ="group"))[2]=="ggplot")
        checkTrue(class(isoPlotPosition(mirData, column="group"))[2]=="ggplot")
    }

test_dseFunctions <-
    function()
    {
        data(mirData)
        checkTrue(class(isoDE(mirData, formula=~group, ref=TRUE,
                            iso5=TRUE))[1]=="DESeqDataSet")
    }

test_countsFunctions <-
    function()
    {
        data(mirData)
        library(GenomicRanges)
        obj<-isoCounts(mirData, ref=TRUE)
        checkTrue(class(counts(obj))[1]=="matrix")
        obj<-isoNorm(obj,formula = ~group)
        checkTrue(class(counts(obj, norm=TRUE))[1]=="matrix")
    }