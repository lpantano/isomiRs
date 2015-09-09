#' Partial Least Squares Discriminant Analysis for isomirSeqData
#' 
#' @aliases isoPLSDA
#' @usage isoPLSDA(ids, var, validation = NULL, learn = NULL, test = NULL,
#'  tol = 0.001, nperm = 400, refinment = FALSE, vip = 1.2)
#' @param ids IsomirDataSeq
#' @param var column name in design data.frame
#' @param validation type of validation, either NULL or "learntest". 
#' Default NULL
#' @param learn	optional vector of indices for a learn-set. 
#' Only used when validation="learntest". Default NULL
#' @param test	optional vector of indices for a test-set. 
#' Only used when validation="learntest". Default NULL
#' @param tol tolerance value based on maximum change of cumulative R-squared 
#' coefficient for each additional PLS component. Default tol=0.001
#' @param nperm	number of permutations to compute the PLD-DA p-value b
#' ased on R2 magnitude. Default nperm=400
#' @param refinment logical indicating whether a refined model, based on 
#' filtering out variables with low VIP values
#' @param vip Variance Importance in Projection threshole value when 
#' a refinement precess is considered. Default vip=1.2
#' @return PLS model
#' @details 
#' Partial Least Squares Discriminant Analysis (PLS-DA) is a technique specifically
#' appropriate for analysis of high dimensionality data sets and multicollineality
#' \cite{perezenciso}. PLS-DA is a supervised method (i.e. makes use of class
#' labels) with the aim to provide a dimension reduction strategy in a situation
#' where we want to relate a binary response variable (in our case young or old
#' status) to a set of predictor variables. Dimensionality reduction procedure is
#' based on orthogonal transformations of the original variables (isomiRs) into a
#' set of linearly uncorrelated latent variables (usually termed as components)
#' such that maximizes the separation between the different classes in the first
#' few components \cite{xia}. We used sum of squares captured by the model (R2) as
#' a goodness of fit measure. We implemented this method using the
#' \code{\link[DiscriMiner]{DiscriMiner-package}} into \code{\link{isoPLSDA}} function. The output
#' p-value of this function will tell about the statistical
#' significant of the group separation using miRNA/isomiR expression data.
#' @return 
#' \code{\link[base]{list}} with the followig elements: R2Matrix 
#' (R-squared coefficents of the PLS model),
#' components (of the PLS),
#' p.value obtained by the permutations and R2PermutationVector with all R2 from
#' the permutations.
#' 
#' If the option \code{refinment} is set to \code{TRUE}, then the following
#' elements will appear: 
#' R2RefinedMatrix and componentsRefinedModel (R-squared coefficients and components 
#'  (components of the PLS model only using the most important miRNAs). As well,
#' p.valRefined and R2RefinedPermutationVector with the p-value and the R2 of the
#' permutations shuffling individuals. And finally, 
#' p.valRefinedFixed and R2RefinedFixedPermutationVector with p-values and R2 of the
#' permutations shuffling the most important miRNAs.
#' 
#' @examples
#' library(DESeq2)
#' data(isomiRexp)
#' ids = isoCounts(isomiRexp, iso5=TRUE, iso3=TRUE, add=TRUE, ref=TRUE)
#' ids = isoNorm(ids)
#' pls.ids = isoPLSDA(ids, "condition", nperm = 10)
#' @export
isoPLSDA <- function(ids, var ,validation = NULL, learn = NULL, test = NULL,
                     tol = 0.001, nperm = 400, refinment = FALSE, vip = 1.2){
    tryCatch ({
        class(normcounts(ids))
    }, error = function(e){
        return("please, run first normIso.")
    })
    variables <- t(normcounts(ids))
    group <- colData(ids)[,var]
    if (length(group) < 6){
        return("this analysis only runs with group larger than 6.")
    }
    # Auxiliar data containing variable names and numeric ID per variable
    dataVariables <- data.frame( variable = colnames(variables),
                                 id = c(1:dim(variables)[2]))
    # PLS-DA model
    model.plsDA <- plsDA(variables, group, validation, learn, test)
    # R-squared coefficents (number of components are selected based on
    # cumulative R-squared coefficient changes)
    if (model.plsDA$R2[dim(model.plsDA$R2)[1],3] < tol){
        a <- which(model.plsDA$R2[,3] < tol)
        if (a[1] > 1)
            R2.mat <- model.plsDA$R2[1:(a[1] - 1),]
        if (a[1] == 1)
            R2.mat <- model.plsDA$R2[1,]
    }
    if (model.plsDA$R2[dim(model.plsDA$R2)[1],3] >= tol){
        R2.mat <- model.plsDA$R2
    }
    # model-based p-value (permutation analysis)
    R2.perm <- R2PermutationVector(variables, group, validation, learn,
                                   test, tol, nperm)
    p.val <- sum(R2.mat[dim(R2.mat)[1],4] <= R2.perm, na.rm=TRUE) / nperm
    if (refinment == TRUE){
        # refine model based on VIPs
        vip.max <- apply(model.plsDA$VIP[,1:dim(R2.mat)[1]], 1, max)
        dataVIP <- data.frame(variable = row.names(model.plsDA$VIP),
                              VIP = vip.max)
        dataVIPref <- dataVIP[dataVIP$VIP >= vip,]
        sel.vars <- merge(dataVariables , dataVIPref,
                          by.x="variable", by.y="variable")
        variables.ref <- variables[,sel.vars[,2]]
        # PLS-DA model with contributing variables
        model.plsDA.ref <- plsDA(variables.ref, group, validation, learn, test)
        if (model.plsDA.ref$R2[dim(model.plsDA.ref$R2)[1],3] < tol){
            a.ref <- which(model.plsDA.ref$R2[,3] < tol)
            if (a.ref[1] > 1)
                R2.mat.ref <- model.plsDA.ref$R2[1:(a.ref[1] - 1),]
            if (a.ref[1] == 1)
                R2.mat.ref <- model.plsDA.ref$R2[1,]
        }
        if (model.plsDA.ref$R2[dim(model.plsDA.ref$R2)[1],3] >= tol)
            R2.mat.ref <- model.plsDA.ref$R2
        # p-value using fixed vip variables while shuffling individuals
        R2.perm.ref1 <- R2PermutationVector(variables.ref, group,
                                            validation, learn, test, tol, nperm)
        p.val.ref1 <- sum(R2.mat.ref[dim(R2.mat.ref)[1],4] <= R2.perm.ref1,
                          na.rm=TRUE) / nperm
        # p-value based on performing the whole refinement process
        # in  each permutation iteration
        R2.perm.ref2 <- R2RefinedPermutationVector(variables, group,
                                                   validation, learn, test,
                                                   tol, nperm, vip)
        p.val.ref2 <- sum(R2.mat.ref[dim(R2.mat.ref)[1],4] <= R2.perm.ref2,
                          na.rm=TRUE) / nperm
        res <- list(R2.mat, model.plsDA$components[,1:dim(R2.mat)[1]],
                    p.val, R2.perm, R2.mat.ref,
                    model.plsDA.ref$components[,1:dim(R2.mat.ref)[2]],
                    p.val.ref2, R2.perm.ref2, p.val.ref1, R2.perm.ref1)
        names(res) <- c("R2Matrix", "components",
                        "p.val", "R2PermutationVector", "R2RefinedMatrix",
                        "componentsRefinedModel",
                        "p.valRefined", "R2RefinedPermutationVector",
                        "p.valRefinedFixed", "R2RefinedFixedPermutationVector")
        print(paste0("pval:",res$p.val))
        return(res)
    }
    if (refinment == FALSE){
        res <- list(R2.mat, model.plsDA$components[,1:dim(R2.mat)[1]],
                    p.val, R2.perm)
        names(res) <- c("R2Matrix", "components", "p.val",
                        "R2PermutationVector")
        # cat(paste0("pval:",res$p.val))
        return(res)
    }
}


# Function that computes p-values when only individuals 
# are considered to be permuted
R2PermutationVector <- function(variables, group, validation,
                                learn, test, tol, nperm){
    # intitialize R2 vector
    R2perm <- rep(NA,nperm)
    # number of individuals
    n <- length(group)
    # permutation loop
    for (p in 1:nperm){
        # shuffle individuals
        mysample <- group[sample(1:length(group), n, replace=FALSE)]
        # pls-da with groups shuffled
        model.perm <- plsDA(variables, mysample, validation, learn, test)
        # select components based on R2 contribution
        if (model.perm$R2[dim(model.perm$R2)[1],3] <= tol){
            a.perm <- which(model.perm$R2[,3] < tol)
            if (a.perm[1] > 1)
                R2perm[p] <- model.perm$R2[a.perm[1] - 1,4]
            if (a.perm[1] == 1)
                R2perm[p] <- model.perm$R2[1,4]
        }
        if (model.perm$R2[dim(model.perm$R2)[1],3] > tol)
            R2perm[p] <- model.perm$R2[dim(model.perm$R2)[1],4]
    }
    return(R2perm)
}


# Function that computes p-values when a refinment process is considered 
# per each permutation iteration
R2RefinedPermutationVector <- function(variables, group, validation, learn,
                                       test, tol, nperm, vip){
    dataVariables <- data.frame( variable = colnames(variables),
                                 id = c(1:dim(variables)[2]))
    n <- length(group)
    R2Refined.perm <- rep(NA,nperm)
    for (p in 1:nperm){
        mysample <- group[sample(1:length(group), n, replace=FALSE)]
        model.perm <- plsDA(variables, mysample, validation, learn, test)
        if (model.perm$R2[dim(model.perm$R2)[1],3] < tol){
            a.perm <- which(model.perm$R2[,3] < tol)
            if (a.perm[1] > 1)
                v <- (a.perm[1] - 1)
            if (a.perm[1] == 1)
                v <- 1
        }
        if (model.perm$R2[dim(model.perm$R2)[1],3] >= tol){
            v <- dim(model.perm$R2)[1]
        }
        vip.perm.max <- apply(model.perm$VIP[,1:v], 1, max)
        dataVIP.perm <- data.frame(variable = row.names(model.perm$VIP),
                                   VIP = vip.perm.max)
        dataVIP.perm.ref <- dataVIP.perm[dataVIP.perm$VIP >= vip,]
        sel.vars.perm <- merge(dataVariables , dataVIP.perm.ref,
                               by.x="variable", by.y="variable")
        variables.perm.ref <- variables[,sel.vars.perm[,2]]
        if(!is.null(dim(variables.perm.ref)[2])){
            model.perm.ref <- plsDA(variables.perm.ref, group,
                                    autosel = TRUE, comps = 2,
                                    validation, learn, test)
            if (model.perm.ref$R2[dim(model.perm.ref$R2)[1],3] <= tol){
                a.perm.ref <- which(model.perm.ref$R2[,3] < tol)
                if (a.perm.ref[1] > 1)
                    R2Refined.perm[p] <- model.perm.ref$R2[a.perm.ref[1] - 1,4]
                if (a.perm.ref[1] == 1)
                    R2Refined.perm[p] <- model.perm.ref$R2[1,4]
            }
            if (model.perm.ref$R2[dim(model.perm.ref$R2)[1],3] > tol)
                R2Refined.perm[p] <- model.perm.ref$R2[dim(model.perm.ref$R2)[1],4]
        }
    }
    return(R2Refined.perm)
}

#' Plot components from isoPLSDA(pairs plot)
#' 
#' @aliases isoPLSDAplot
#' @usage isoPLSDAplot(components, groups)
#' @param components PLS-DA components as it comes from isoPLSDA main function
#' @param groups	vector or factor with group memberships
#' @return plot
#' @details 
#' The function \code{isoPLSDAplot}
#' helps to visualize the results. It will plot the samples using the 
#' significant components (t1, t2, t3 ...) from the PLS-DA analysis and the 
#' samples distribution along the components.
#' @return \code{\link[ggplot2]{ggplot2-package}} object
#' @examples
#' data(isomiRexp)
#' ids = isoCounts(isomiRexp, iso5=TRUE, iso3=TRUE, add=TRUE, ref=TRUE)
#' ids = isoNorm(ids)
#' pls.ids = isoPLSDA(ids, "condition", nperm = 10)
#' isoPLSDAplot(pls.ids$component, colData(ids)[,"condition"])
#' @export
isoPLSDAplot <- function (components, groups){
    datacomponents <- data.frame(condition = groups, components)
    t <- dim(datacomponents)[2] - 1
    n <- length(levels(factor(groups)))
    ggplot <- function (...) ggplot2::ggplot(...) +
        scale_color_brewer(palette="Set1")
    unlockBinding("ggplot",parent.env(asNamespace("GGally")))
    assign("ggplot",ggplot,parent.env(asNamespace("GGally")))
    ggpairs(datacomponents, columns = 2:3, col="condition",
            upper="blank",legends=TRUE)
}
