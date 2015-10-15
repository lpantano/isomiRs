#' Partial Least Squares Discriminant Analysis for \code{\link{IsomirDataSeq}}
#' 
#' Use PLS-DA method with the normalized count data to detect the most
#' important features (miRNAs/isomiRs) that explain better
#' the group of samples given in the experimental design. It is a supervised
#' clustering method with permutations to calculate the significance
#' of analysis.
#' 
#' @aliases isoPLSDA
#' @usage isoPLSDA(ids, group, validation = NULL, learn = NULL, test = NULL,
#'  tol = 0.001, nperm = 400, refinment = FALSE, vip = 1.2)
#' @param ids object of class \code{\link{IsomirDataSeq}}
#' @param group column name in design data.frame
#' @param validation type of validation, either NULL or "learntest". 
#' Default NULL
#' @param learn	optional vector of indexes for a learn-set. 
#' Only used when validation="learntest". Default NULL
#' @param test	optional vector of indices for a test-set. 
#' Only used when validation="learntest". Default NULL
#' @param tol tolerance value based on maximum change of cumulative R-squared 
#' coefficient for each additional PLS component. Default tol=0.001
#' @param nperm	number of permutations to compute the PLD-DA p-value
#' based on R2 magnitude. Default nperm=400
#' @param refinment logical indicating whether a refined model, based on 
#' filtering out variables with low VIP values
#' @param vip Variance Importance in Projection threshold value when 
#' a refinement precess is considered. Default vip=1.2
#' @details 
#' Partial Least Squares Discriminant Analysis (PLS-DA) is a technique specifically
#' appropriate for analysis of high dimensionality data sets and multicollinearity
#' (\cite{Perez-Enciso, 2013}). PLS-DA is a supervised method (i.e. makes use of class
#' labels) with the aim to provide a dimension reduction strategy in a situation
#' where we want to relate a binary response variable (in our case young or old
#' status) to a set of predictor variables. Dimensionality reduction procedure is
#' based on orthogonal transformations of the original variables (miRNAs/isomiRs) into a
#' set of linearly uncorrelated latent variables (usually termed as components)
#' such that maximizes the separation between the different classes in the first
#' few components (\cite{Xia, 2011}). We used sum of squares captured by the model (R2) as
#' a goodness of fit measure. We implemented this method using the
#' \code{\link[DiscriMiner]{DiscriMiner-package}} into \code{\link{isoPLSDA}} function. The output
#' p-value of this function will tell about the statistical
#' significant of the group separation using miRNA/isomiR expression data.
#' @return 
#' \code{\link[base]{list}} with the following elements: \code{R2Matrix}
#' (R-squared coefficients of the PLS model),
#' \code{components} (of the PLS), \code{vip} (most important variables),
#' \code{group} (classification of the samples),
#' \code{p.value} and \code{R2PermutationVector} obtained by the permutations.
#' 
#' If the option \code{refinment} is set to \code{TRUE}, then the following
#' elements will appear: 
#' \code{R2RefinedMatrix} and \code{componentsRefinedModel} (R-squared coefficients
#' of the PLS model only using the most important miRNAs/isomiRs). As well,
#' \code{p.valRefined} and \code{R2RefinedPermutationVector} with p-value and R2 of the
#' permutations shuffling individuals. And finally, 
#' \code{p.valRefinedFixed} and \code{R2RefinedFixedPermutationVector} with p-value and R2 of the
#' permutations shuffling the most important miRNAs/isomiRs.
#' @references 
#' Perez-Enciso, Miguel and Tenenhaus, Michel. Prediction of clinical outcome with microarray data:
#' a partial least  squares discriminant analysis (PLS-DA) approach. Human
#' Genetics. 2003.
#' 
#' Xia, Jianguo and Wishart, David S. Web-based inference of biological patterns, functions and
#' pathways from metabolomic data using MetaboAnalyst. Nature Protocols. 2011.
#' @examples
#' library(DESeq2)
#' data(mirData)
#' ids = isoCounts(mirData, iso5=TRUE, iso3=TRUE, add=TRUE, ref=TRUE)
#' ids = isoNorm(ids)
#' pls.ids = isoPLSDA(ids, "condition", nperm = 10)
#' cat(paste0("pval:",pls.ids$p.val))
#' cat(paste0("components:",pls.ids$components))
#' @export
isoPLSDA <- function(ids, group ,validation = NULL, learn = NULL, test = NULL,
                     tol = 0.001, nperm = 400, refinment = FALSE, vip = 1.2){
    tryCatch ({
        class(normcounts(ids))
    }, error = function(e){
        return("please, run first normIso.")
    })
    variables <- t(normcounts(ids))
    group <- colData(ids)[,group]
    if (length(group) < 6){
        return("this analysis only runs with group larger than 6.")
    }
    # Auxiliary data containing variable names and numeric ID per variable
    dataVariables <- data.frame( variable = colnames(variables),
                                 id = c(1:dim(variables)[2]))
    # PLS-DA model
    model.plsDA <- plsDA(variables, group, validation, learn, test)
    # R-squared coefficients (number of components are selected based on
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
    vip.max <- apply(model.plsDA$VIP[,1:dim(R2.mat)[1]], 1, max)
    dataVIP <- data.frame(variable = row.names(model.plsDA$VIP),
                          VIP = vip.max)
    dataVIPref <- dataVIP[dataVIP$VIP >= vip,]
    sel.vars <- merge(dataVariables , dataVIPref,
                      by.x="variable", by.y="variable")
    variables.ref <- variables[,sel.vars[,2]]
    if (refinment == TRUE){
        # refine model based on VIPs
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
                    dataVIPref, p.val, R2.perm, R2.mat.ref,
                    model.plsDA.ref$components[,1:dim(R2.mat.ref)[2]],
                    p.val.ref2, R2.perm.ref2, p.val.ref1, R2.perm.ref1)
        names(res) <- c("R2Matrix", "components", "vip",
                        "p.val", "R2PermutationVector", "R2RefinedMatrix",
                        "componentsRefinedModel",
                        "p.valRefined", "R2RefinedPermutationVector",
                        "p.valRefinedFixed", "R2RefinedFixedPermutationVector")
        res[["group"]] = group
        return(res)
    }
    if (refinment == FALSE){
        res <- list(R2.mat, model.plsDA$components[,1:dim(R2.mat)[1]],
                    dataVIPref, p.val, R2.perm)
        names(res) <- c("R2Matrix", "components", "vip", "p.val",
                        "R2PermutationVector")
        res[["group"]] = group
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


# Function that computes p-values when a refinement process is considered 
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
#' Plot the most significant components that come from \code{\link{isoPLSDA}}
#' together with the density of the samples along those components.
#' 
#' @aliases isoPLSDAplot
#' @usage isoPLSDAplot(pls)
#' @param pls output from \code{\link{isoPLSDA}} function.
#' @return \code{\link[GGally]{ggpairs}} plot
#' @details 
#' The function \code{isoPLSDAplot}
#' helps to visualize the results from \code{\link{isoPLSDA}}.
#' It will plot the samples using the 
#' significant components (t1, t2, t3 ...) from the PLS-DA analysis and the 
#' samples distribution along the components. It uses \code{\link[GGally]{ggpairs}}
#' for the plot.
#' @return \code{\link[ggplot2]{ggplot2-package}} object
#' @examples
#' data(mirData)
#' ids = isoCounts(mirData, iso5=TRUE, iso3=TRUE, add=TRUE, ref=TRUE)
#' ids = isoNorm(ids)
#' pls.ids = isoPLSDA(ids, "condition", nperm = 10)
#' isoPLSDAplot(pls.ids)
#' @export
isoPLSDAplot <- function (pls){
    components = pls$component
    groups = pls$group
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
