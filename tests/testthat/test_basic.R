data(mirData)
data(isoExample)

test_that("plots", {
    expect_s3_class(isoPlot(mirData, type = "all", column = "condition"), "ggplot")
    expect_s3_class(isoPlot(mirData, type = "iso5", column = "condition"), "ggplot")
    expect_s3_class(isoPlot(mirData, type = "iso3", column = "condition"), "ggplot")
    expect_s3_class(isoPlot(mirData, type = "add", column = "condition"), "ggplot")
    expect_s3_class(isoPlotPosition(mirData, column = "condition"), "ggplot")
})

test_that("counts", {
    obj <- isoCounts(mirData, ref = TRUE)
    expect_is(counts(obj),"matrix")
    obj <- isoNorm(obj, formula = ~condition)
    expect_is(counts(obj, norm = TRUE), "matrix")
})  

# test_that("target", {
#     mirna_ma <- data.frame(gene = names(gene_ex_rse)[1:20],
#                            mir = names(mirna_ex_rse))
#     corMat <- findTargets(mirna_ex_rse, gene_ex_rse, mirna_ma)
#     expect_output(str(corMat), "num [1:20, 1:20]", fixed = TRUE)
#     expect_error(pairsMatrix(list()),
#                  "Need a data.frame with 2 columns: gene, mir.")
#     expect_equal(.detect_gene_symbol(names(gene_ex_rse)), "ENSEMBL")
#     expect_false(.is_mapping_needed(names(gene_ex_rse), names(gene_ex_rse)))
# })

test_that("accesor", {
    expect_output(str(design(mirData)), "formula")
    expect_is(counts(mirData), "matrix")
    expect_equal(nrow(isoSelect(mirData, "hsa-let-7a-5p")), 37)
})