#' Correct miRNA expression based on prior ligation bias information
#' 
#' this is the source file for fitting the linear quadratic normal family
#' 
#' @param train long data.frame to train model
#' @param data long data.frame to correct abundance
#' @param cycles number of cycles to reach convergency
#' @param long boolean if input is in long format instead of standard
#' wide format (rows:miRNAs, columns:samples)
#' @return data.frame with corrected expression
#' @examples
#' options(warn = -1) # this is only for tiny example
#' data(mirTritation)
#' ma <- isoCorrect(mirTritation[mirTritation$class=="train",],
#' mirTritation[mirTritation$class=="test",],cycles=5,long=TRUE)
#' library(ggplot2)
#' ggplot(ma,aes(y=log2(reads), x=Dilution)) + geom_jitter()
#' ggplot(ma,aes(y=m, x=Dilution)) + geom_jitter()
#' @author Christos Argyropoulos and Lorena Pantano
#' @details 
#' Methods adapted from:
#' 
#' Argyropoulos, Christos, et al. "Modeling bias and variation in 
#' the stochastic processes of small RNA sequencing." 
#' Nucleic Acids Research (2017).
#' @export
isoCorrect <- function(train, data, cycles=5000, long=FALSE){
    
    if (!class(train) %in% c("data.frame", "matrix"))
        stop("train data should be data.frame or matrix")
    if (!class(data) %in% c("data.frame", "matrix"))
        stop("target data should be data.frame or matrix")
    
    if (long == FALSE){
        train = reshape::melt(as.matrix(train))
        names(train) = c("miRNA", "SampleID", "reads")
        data = reshape::melt(as.matrix(data))
        names(data) = c("miRNA", "SampleID", "reads")
    }
    
    fLQNO.r <- gamlss(reads~ SampleID+random(miRNA),
                       sigma.formula=~ SampleID+random(miRNA),
                       data=train,family="LQNO",
                       control=gamlss.control(n.cyc=cycles,c.crit=0.1,
                                              mu.step=1,sigma.step=1),method=RS())

    
    corrfact <- function(fit,data) {
        pred <- predictAll(fit,type="terms",data=data)
        off.m <- pred$mu[,"random(miRNA)"]
        off.s <- pred$sigma[,"random(miRNA)"]
        off <- aggregate(cbind(off.m,off.s),by=list(data$miRNA),FUN=mean)
        names(off)[1] <- "miRNA"
        off
    }
    ## correction factors from LQNO
    cor.LQNO <- corrfact(fLQNO.r,train)
    ## correction factors from NBI

    ## in order to correct the development , merge the correction factors into it
    ## we will refit both models and see what we get
    data.LQNO <- merge(data,cor.LQNO, by="miRNA") ## correction factors from the LQNO model

    fLQNO.data.r <- gamlss(reads~ SampleID+random(miRNA)+offset(off.m),
                           sigma.formula=~ SampleID+random(miRNA)+offset(off.s),
                           data=data.LQNO,family="LQNO",
                           control=gamlss.control(n.cyc=cycles,c.crit=0.1,
                                                  mu.step=1,sigma.step=1,gd.tol=Inf), 
                           i.control=glim.control(cc = 0.001, cyc = cycles,  
                                                  glm.trace = FALSE, 
                                                  bf.cyc = 300, bf.tol = 0.1, 
                                                  bf.trace = FALSE),
                           method=RS())

    ## function to predict the effects of bias correction
    pred.fun <- function(x) {
        dummy <- predict(x[[1]],type="terms",se.fit=TRUE)
        d <- data.frame(m=dummy$fit[,"random(miRNA)"],
                        s=dummy$se.fit[,"random(miRNA)"])
        x[[2]]$m <- d$m
        x[[2]]$s <- d$s
        x[[2]]
    }

    
    fitSE.LQNO <- pred.fun(list(fLQNO.data.r, data.LQNO))

    if (long == TRUE)
        return(fitSE.LQNO)
    ma <- fitSE.LQNO[, c("miRNA", "SampleID", "m")] %>% spread(key="SampleID", value="m")
    row.names(ma) <- ma[,1]
    ma <- ma[,2:ncol(ma)]
}
