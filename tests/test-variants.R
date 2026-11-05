# The supplied run_tests.R loads testthat and the functions from R/.
# Do not source original/messy_analysis.R here.

# A supplied ordinary example; keep it and add your own tests below.
test_that("MAF counts allele copies among observed diploid genotypes", {
    expect_equal(calculate_maf(c(0, 0, 1, 1)), 0.25)
})

#calculate_maf tests
test_that("MAF is correct with normal vectors", {
    #tests if the function handles a case where the counted allele is common
    expect_equal(calculate_maf(c(2, 2, 2, 1, 2, 0)), 0.25)
    
    expect_equal(calculate_maf(c(0, 0, 0)), 0)
    
    expect_equal(calculate_maf(c(2, 2, 2)), 0)
    
    expect_equal(calculate_maf(c(1, 1, 2, 2)), 0.75)
    
    #when there is one heterozygote the maf is correct
    expect_equal(calculate_maf(1), 0.5)
})

test_that("MAF calculation drops missing values correctly when na_rm is TRUE", {
    expect_equal(calculate_maf(c(0, 0, NA, 1, 1)), 0.25)
})

test_that("MAF returns NA when na_rm is FALSE"{
    expect_equal(calculate_maf(c(0, 0, NA, 1, 1), na_rm = FALSE), NA)
    
    #checks if this is done only when something is missing
    expect_equal(calculate_maf(c(0, 0, 1, 1), na_rm = FALSE), 0.25)
})

test_that("empty and completely missing input returns NA", {
    expect_true(is.na(calculate_maf(numeric(0))))
                 
    expect_true(is.na(calculate_maf(c(NA_real_, NA_real_))))
})

#summarize_variant tests
test_that("calculations are correct with full input", {
    results <- summarize_variant(c(0, 0, 1, 1))
    expect_equal(nrow(results), 1)
    expect_equal(results[1], 4)
    expect_equal(results[2], 0)
    expect_equal(results[3], 0.25)
    expect_true(results[4])
})

test_that("calculations are correct with some missing input", {
    results <- summarize_variant(c(0, 0, NA, 1, 0))
    expect_equal(nrow(results), 1)
    expect_equal(results[1], 4)
    expect_equal(results[2], 0.2)
    expect_equal(results[3], 0.125)
    expect_true(!results[4])
})

test_that("calculations are correct with empty and all missing input", {
    results <- summarize_variant(numeric(0))
    expect_equal(results[1], 0)
    expect_true(is.na(results[2]))
    expect_true(is.na(results[3]))
    expect_true(!results[4])
    
    results <- summarize_variant(c(NA_real_, NA_real_, NA_real_))
    expect_equal(results[1], 0)
    expect_equal(results[2], 1)
    expect_true(is.numeric(results[3])&&is.na(results[3]))
    expect_true(!results[4])
    
})

#threshold tests of summarize_variants
test_that("missing rate threshold functions correctly", {
    results <- summarize_variant(c(1, NA, NA, NA, NA))
    expect_equal(results[2], 0.8)
    expect_true(!results[4])
    
    results <- summarize_variant(c(2, 1, 2, 1, NA))
    expect_equal(results[2], 0.2)
    expect_true(results[4])
    
    results <- summarize_variant(c(2, 1, 2, 1, NA), max_missing = 0.1)
    expect_equal(results[2], 0.2)
    expect_true(!results[4])
    
    results <- summarize_variant(c(2, 1, 2, NA, NA), max_missing = 0.5)
    expect_equal(results[2], 0.4)
    expect_true(!results[4])
})

test_that("maf threshold functions correctly", {
    results <- summarize_variant(c(1, 1, 1, 1))
    expect_equal(results[3], 0.5)
    expect_true(results[4])
    
    results <- summarize_variant(c(0, 0, 0, 1))
    expect_equal(results[3], 0.125)
    expect_true(!results[4])
    
    results <- summarize_variant(c(1, 1, 1, 1), min_maf = 0.6)
    expect_equal(results[3], 0.5)
    expect_true(!results[4])
    
    results <- summarize_variant(c(0, 0, 0, 1), min_maf = 0.2)
    expect_equal(results[3], 0.125)
    expect_true(results[4])
})


#error checking
test_that("maf throws errors when given invalid genotypes", {
    expect_error(calculate_maf(c(0, 1, 2, 3)))
    
    expect_error(calculate_maf(c(0, -1)))
    
    expect_error(calculate_maf(c(0, 0.5)))
    
    expect_error(calculate_maf(c("a", "b")))
    
    #genotype is still validated even when na_rm = FALSE and an NA is present
    expect_error(calculate_maf(c(3, NA), na_rm = FALSE))
})


test_that("maf throws errors when given invalid na_rm", {
    expect_error(calculate_maf(c(0, 1, 2, 1), na_rm = 1))
    
    expect_error(calculate_maf(c(0, 1, 2, 1), na_rm = "a"))
    
    expect_error(calculate_maf(c(0, 1, 2, 1), na_rm = NA))
    
})


test_that("summarize_variants throws errors when invalid settings are given", {
    expect_error(summarize_variant(c(0, 1, 2, 1), min_maf = -.01))
    
    expect_error(summarize_variant(c(0, 1, 2, 1), min_maf = Inf))
    
    expect_error(summarize_variant(c(0, 1, 2, 1), min_maf = "a"))
    
    expect_error(summarize_variant(c(0, 1, 2, 1), max_missing = -.01))
    
    expect_error(summarize_variant(c(0, 1, 2, 1), max_missing = Inf))
    
    expect_error(summarize_variant(c(0, 1, 2, 1), max_missing = "a"))
    
})
