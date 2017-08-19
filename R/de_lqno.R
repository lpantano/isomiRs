#' Differential expression between two groups using LQNO model
#' 
#' 
#' @author Christos Argyropoulos and Lorena Pantano
#' @details 
#' Methods adapted from *Argyropoulos et al (2017)*.
#' 
#' @references 
#' Argyropoulos, Christos, et al. "Modeling bias and variation in 
#' the stochastic processes of small RNA sequencing." 
#' Nucleic Acids Research (2017).
#' @param counts Count matrix.
#' @param groups Character vector to indicate the group of each sample.
#' @param long Whether matrix is in long format. Default FALSE.
#' @return data.frame with estimates and p-values.
#' @examples
#' options(warn = -1) # this is only for tiny example
#' data(dat286)
#' datRat<-subset(dat286.long,(Series=="Equi" | Series =="RatioA") & Amount=="100 fmoles")
#' datRat$SampleID<-factor(datRat$SampleID)
#' datRat$Series<-factor(datRat$Series)
#' res <- isoLQNO(datRat, long=TRUE)
#' @export
isoLQNO <- function(counts, groups=NULL, long=FALSE){
    
    if (!class(counts) %in% c("data.frame", "matrix"))
        stop("counts data should be data.frame or matrix")
    
    if (long == FALSE){
        names(groups) = colnames(counts)
        counts = reshape::melt(as.matrix(counts))
        names(counts) = c("miRNA", "SampleID", "reads")
        counts[, "Series"] = groups[as.character(counts$SampleID)]
    }
    counts$SampleID<-factor(counts$SampleID)
    counts$Series<-factor(counts$Series)
    
    ## Sequential modeling strategy to ensure that the data are estimated robustly
    ## STEP 0: PREPARE THE DATA FOR THE FIT - THIS PART DOES NOT CHANGE
    u_X<-as.numeric(factor(counts$miRNA)) ## maps readings to the identity of the miRNA
    u_G<-as.numeric(factor(counts$Series)) ## maps counts to group
    y<-counts$reads  ## extract the actual counts
    X<-model.matrix(~Series,data=counts) ## design matrix (ANOVA) for group comparisons
    
    ## STEP 1: USE A POISSON MODEL TO OBTAIN ESTIMATES FOR THE MU SUBMODEL
    ##====================================================================
    ## fit the parameters for the mu submodel using the poisson GLM
    gl<-glmer(reads~Series+(0+Series|miRNA),data=counts,family="poisson")
    
    ## STEP 2: USE THE MU MODEL ESTIMATES TO FIT THE PHI SUBMODEL
    ##============================================================
    ## initializes standard deviation of RE for the mu submodel
    sigmu=sqrt(diag((summary(gl)[["varcor"]])[[1]]))  
    sigsig=rep(1,max(u_G)) ## initializes standard deviation of RE for the phi submodel
    b=fixef(gl) ## initial values for the overall group means (mean submodel) 
    ## initial values for the variation of miRNAs around their group mean (mean submodel)
    u_m=as.matrix(ranef(gl)$miRNA)
    ## Very rough initial values for the phi submodel parameters
    s_b=rep(0,ncol(X)) ## initial values for the overall group means (phi submodel) 
    ## initial values for the variation of miRNAs around their group mean (phi submodel)
    u_s= matrix(0,max(u_X),max(u_G))
    ## MAP list that allow us to fix some parameters to their values
    MAP<-NULL
    MAP[["b"]]<-factor(rep(NA,length(b)))
    MAP[["u_m"]]<-factor(rep(NA,length(c(u_m))))
    MAP[["sigmu"]]<-factor(rep(NA,length(sigmu)))
    ## construct the AD object - note that we fix the mu at their values while estimating the
    ## phi submodel
    obj.TMB<-MakeADFun(data=list(y=y,X=X,u_X=u_X,u_G=u_G),
                       parameters=list(b=b,s_b=s_b,u_m=u_m,u_s=u_s,
                                       sigmu=sigmu,sigsig=sigsig),
                       DLL="isomiRs",random=c("u_s"),hessian=FALSE,silent=TRUE,
                       method="BFGS",random.start=expression(last.par.best[random]),
                       ADReport=TRUE,map=MAP)
    ## parameter estimation - note errors may be generated during some iterations
    f.TMB<-suppressWarnings(nlminb(obj.TMB$par,obj.TMB$fn,obj.TMB$gr,
                  control=list(eval.max=10000,iter.max=10000),lower=-30,upper=30))
    ## obtain the report on the parameters to extract the fitted values of the gamlss model
    rep<-sdreport(obj.TMB) 
    u_s = matrix(summary(rep,"random",p.value=FALSE)[,1],ncol=max(u_G))
    dummy<-summary(rep,"fixed",p.value=FALSE)[,1] ## parameter estimates
    s_b=dummy[1:max(u_G)]
    sigsig=dummy[-(1:max(u_G))]
    
    ## STEP 3: REFIT THE MODEL WITHOUT FIXING ANY PARAMETERS
    ##=========================================================
    obj.TMB<-MakeADFun(data=list(y=y,X=X,u_X=u_X,u_G=u_G),
                       parameters=list(b=b,s_b=s_b,u_m=u_m,u_s=u_s,
                                       sigmu=sigmu,sigsig=sigsig),
                       DLL="isomiRs",random=c("u_m","u_s"),hessian=TRUE,silent=TRUE,
                       method="BFGS",random.start=expression(last.par.best[random]),
                       ADReport=TRUE)
    ## scale objective by the magnitude of the deviance of the fitted Poisson model
    f.TMB<-suppressWarnings(nlminb(obj.TMB$par,obj.TMB$fn,obj.TMB$gr,
                  control=list(eval.max=10000,iter.max=10000,scale=deviance(gl)),
                  lower=-30,upper=30))
    ## obtain the report on the parameters
    rep<-sdreport(obj.TMB) 
    ## differential expression ratios, standard errors z and p-values
    gamlssAD<-as.data.frame(summary(rep,"report",p.value=TRUE)[1:nlevels(counts$miRNA),])
    rownames(gamlssAD)<-levels(counts$miRNA)                     
    gamlssAD[,"FDR"] <- p.adjust(gamlssAD[,4], method = "fdr")
    gamlssAD
}

