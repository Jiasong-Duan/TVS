context("TVS main functionality")

test_that("tvs() works with sample data", {
  data(data_tvs)  # Your example dataset

  # Basic smoke test
  result <-TVS(wrapper_beta, wrapper_beta0, wrapper_nu, wrapper_gamma, data_tvs)

  # Test output structure
  expect_true(is.list(result))
  expect_named(result, c("pvals"))  # Adjust as needed
})
