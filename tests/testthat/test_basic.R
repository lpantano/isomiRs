data(mirData)
data(isoExample)

test_that("plots", {
    expect_true(class(isoPlot(mirData, type = "all", column = "condition"))[2] == "ggplot")
    expect_true(class(isoPlot(mirData, type = "iso5", column = "condition"))[2] == "ggplot")
    expect_true(class(isoPlot(mirData, type = "iso3", column = "condition"))[2] == "ggplot")
    expect_true(class(isoPlot(mirData, type = "add", column = "condition"))[2] == "ggplot")
    expect_true(class(isoPlotPosition(mirData, column = "condition"))[2] == "ggplot")
})

test_that("counts", {
    obj <- isoCounts(mirData, ref = TRUE)
    expect_true(class(counts(obj))[1] == "matrix")
    obj <- isoNorm(obj, formula = ~condition)
    expect_true(class(counts(obj, norm = TRUE))[1] == "matrix")
})  

test_that("psl-da", {
    expect_error(isoPLSDA(mirData, "condition", nperm = 2),
                 "please, run first isoNorm")
    pls.ids = isoPLSDA(isoNorm(mirData), "condition", nperm = 2)
    expect_output(str(pls.ids), "List of 6")
    expect_true(class(isoPLSDAplot(pls.ids)) == "data.frame")
})

test_that("target", {
    mirna_ma <- data.frame(gene = names(gene_ex_rse)[1:20],
                           mir = names(mirna_ex_rse))
    corMat <- findTargets(mirna_ex_rse, gene_ex_rse, mirna_ma)
    expect_output(str(corMat), "num [1:20, 1:20]", fixed = TRUE)
    expect_error(pairsMatrix(list()),
                 "Need a data.frame with 2 columns: gene, mir.")
    expect_equal(.detect_gene_symbol(names(gene_ex_rse)), "ENSEMBL")
    expect_false(.is_mapping_needed(names(gene_ex_rse), names(gene_ex_rse)))
})

test_that("accesor", {
    expect_output(str(design(mirData)), "formula")
    expect_true(class(counts(mirData)) == "matrix")
    expect_equal(nrow(isoSelect(mirData, "hsa-let-7a-5p")), 37)
})