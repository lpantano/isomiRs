test_plotFunctions <-
    function()
    {
        data(isomiRex)
        checkTrue(class(plotIso(isomiRex,type="t5"))[1]=="IsomirDataSeq")
        checkTrue(class(plotIso(isomiRex,type="t3"))[1]=="IsomirDataSeq")
        checkTrue(class(plotIso(isomiRex,type="add"))[1]=="IsomirDataSeq")
    }

test_dseFunctions <-
    function()
    {
        data(isomiRex)
        checkTrue(class(deIso(isomiRex,formula=~condition,ref=TRUE,
                            iso5=TRUE))[1]=="DESeqDataSet")
    }

test_countsFunctions <-
    function()
    {
        data(isomiRex)
        obj<-countsIso(isomiRex,ref=TRUE)
        checkTrue(class(obj@counts)[1]=="matrix")
        obj<-normIso(obj)
        checkTrue(class(obj@normcounts)[1]=="matrix")
    }