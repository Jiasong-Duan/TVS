test_that("tvs() works with sample data", {
  data(data_tvs)  # Your example dataset

  # Basic smoke test
  result <-TVS(data_tvs)

  # Test output structure
  expect_true(is.list(result))
  expect_named(result, c("selected_indices", "p_values"))
  expect_true(is.numeric(result$p_values))
  expect_true(is.integer(result$selected_indices) || is.numeric(result$selected_indices))
})
