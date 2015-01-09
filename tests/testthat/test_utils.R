library(threshold_analysis)

context('Normalizers')
normalizers <- list(standardize,
                    max.normalizer,
                    max.min.normalizer,
                    q.min.normalizer)
test_that('Normalizers work for edge cases', {
    zero.vector <- c() # Empty input
    one.vector <- c(0) # Length one input
    two.vector <- c(0, 0) # Length n input with only one unique value
    zero.test <- lapply(normalizers,
                        function(f){
                            expect_equals(zero.vector, f(zero.vector))
                            expect_equals(one.vector, f(one.vector))
                            expect_equals(two.vector, f(two.vector))
                        })
})
test_that('Normalizers work for normal cases', {
    test.positive <- c(1, 1, 1, 2)
    test.negative <- c(-1, -2, -3)
    test.mixed <- c(1, -2)

    expect_equals(c(-1, -1, -1, 1), standardize(test.positive))
    expect_equals(c(0.5, 0.5, 0.5, 1), max.normalizer(test.positive))
    expect_equals(c(0, 0, 0, 1), max.min.normalizer(test.positive))
    expect_equals(c(0, 0, 0, 1), q.min.normalizer(test.positive))

    # Need more test cases
})

context('Miscellaneous Utils')

context('Distributions')
