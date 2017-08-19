context("stats")

test_that("LQNO", {
    expect_equal(dLQNO(0.1, mu = 1, sigma = 1, log = FALSE), 0.230383)
    expect_equal(pLQNO(0.1, mu = 1, sigma = 1, log = FALSE), 0.2622591,
                 tolerance = 0.001)
    expect_equal(qLQNO(0.1, mu = 1, sigma = 1, log = FALSE), -0.8123876)
    expect_output(str(rLQNO(5, mu = 1, sigma = 1)), "num [1:5]", fixed = TRUE)
})