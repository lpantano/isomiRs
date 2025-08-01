# context("calculus")
# data(isoExample)

# test_that("matrix", {
#     ma <- assay(gene_ex_rse)
#     g <- colData(gene_ex_rse)[["group"]]
#     maMean <- .apply_median(ma, g)
#     expect_equal(ncol(maMean), 6L)
#     geneGroup <- c(rep(1L, 15L), rep(2L, 10L))
#     names(geneGroup) <- rownames(ma)
#     expect_output(str(.plot_profiles(geneGroup, maMean)), "List of 2")
#     maRowMean <- mean(sapply(.scale(ma), mean))
#     expect_true(maRowMean > -1 & maRowMean < 1)
#     pairs <- as.matrix(ma_ex[1:2,1:3])
#     gene <- assay(gene_ex_rse)[rownames(pairs),]
#     mir <- assay(mirna_ex_rse)[colnames(pairs),]
#     mirTarMa <- .cor_matrix(mir, gene, pairs, -.54)
#     expect_equal(sum(mirTarMa), -0.55, tolerance = 0.01)
#     mirTarMa <- .cor_matrix(mir, gene, pairs, -.7)
#     expect_equal(sum(mirTarMa), 0, tolerance = 0.01)
#     expect_error(.cor_matrix(gene, mir, pairs))
#     expect_error(.cor_matrix(data.frame(gene), mir, pairs))
#     expect_equal(unique(.cluster_exp(gene)), c(1, 2))
# })

# test_that("mining", {
#     expect_match(.reduce_mirna("mmu-miR-181c-5p"), "miR-181")
#     df <- data.frame(gene = "ENSMUSG00000021822",
#                      mir = "mmu-miR-181c-5p,mmu-miR-30e-5p",
#                      f1 = "1",
#                      f2 = "2",
#                      term = "GO:fake")
#     expect_equal(nrow(.summary_mirna(df)), 2L)
#     expect_match(.change_seq("TTTTTT", "3TA"), "TTATTT")
# })