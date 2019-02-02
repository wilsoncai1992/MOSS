library(testthat)
library(MOSS)

Sys.setenv(R_TESTS = "")
test_check("MOSS")
