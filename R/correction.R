#' Calculate the correction factors from the equimolar 286 experiments
#' 
#' this is the source file for fitting the linear quadratic normal famile
#' change the string to the location you saved this file in your drive
#' 
#' @param train long data.frame to train model
#' @param data long data.frame to correct abundance
#' @param cycles number of cycles to reach convergency
#' @return data.frame with corrected expression
#' @examples
#' options(warn = -1)
#' data(mirTritation)
#' ma <- isoCorrect(mirTritation[mirTritation$class=="train",],
#' mirTritation[mirTritation$class=="test",],cycles=5)
#' library(ggplot2)
#' ggplot(ma,aes(y=log2(reads), x=Dilution)) + geom_jitter()
#' ggplot(ma,aes(y=m, x=Dilution)) + geom_jitter()
#' @export
isoCorrect <- function(train, data, cycles=5000){
    fLQNO.r <- gamlss(reads~ SampleID+random(miRNA),
                       sigma.formula=~ SampleID+random(miRNA),
                       data=train,family="LQNO",
                       control=gamlss.control(n.cyc=cycles,c.crit=0.1,
                                              mu.step=1,sigma.step=1),method=RS())
    ## fit the NBI to the 01 data.
    # fNBI.01.r<-gamlss(reads~ SampleID+random(miRNA),
    #                   sigma.fo=~ SampleID+random(miRNA),
    #                   data=train,family="NBI",
    #                   control=gamlss.control(n.cyc=5000,c.crit=0.001,
    #                                          mu.step=1,sigma.step=1,gd.tol=Inf),
    #                   method=RS())
    # 
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
    # cor.NBI<-corrfact(fNBI.01.r,train)
    
    ## in order to correct the development , merge the correction factors into it
    ## we will refit both models and see what we get
    data.LQNO <- merge(data,cor.LQNO, by="miRNA") ## correction factors from the LQNO model
    # datA.NBI<-merge(datA,cor.NBI, by="miRNA")
    
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
    # 
    # fNBI.1.r<-gamlss(reads~ SampleID+random(miRNA)+offset(off.m),
    #                  sigma.fo=~ SampleID+random(miRNA)+offset(off.s),
    #                  data=datA.NBI,family="NBI",
    #                  control=gamlss.control(n.cyc=5000,c.crit=0.001,
    #                                         mu.step=1,sigma.step=1,gd.tol=Inf), 
    #                  i.control=glim.control(cc = 0.001, cyc = 500,  glm.trace = FALSE, 
    #                                         bf.cyc = 300, bf.tol = 0.001, bf.trace = FALSE),
    #                  method=RS())
    # 
    
    ## function to predict the effects of bias correction
    pred.fun <- function(x) {
        dummy <- predict(x[[1]],type="terms",se.fit=TRUE)
        d <- data.frame(m=dummy$fit[,"random(miRNA)"],
                        s=dummy$se.fit[,"random(miRNA)"])
        x[[2]]$m <- d$m
        x[[2]]$s <- d$s
        # d<-ddply(x[[2]][,c("miRNA","m","s")],.(miRNA),summarize,m=mean(m)/log(10),s=mean(s)/log(10))
        x[[2]]
    }
    # 
    # 
    fitSE.LQNO <- pred.fun(list(fLQNO.data.r, data.LQNO))
    # fitSE.NBI<-pred.fun(list(fNBI.1.r,datA.NBI))
    # 
    fitSE.LQNO
}
