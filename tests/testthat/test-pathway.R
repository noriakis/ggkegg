test_that("Generate test pathway without errors", {
	expect_error( create_test_pathway(), NA)
})
test_that("Pathway downloading without errors", {
	expect_error( pathway("hsa04110"), NA)
})
test_that("process_line without errors", {
	expect_error( create_test_pathway(line=TRUE) |> process_line(), NA)
})
test_that("process_reaction without errors", {
	expect_error( create_test_pathway(line=TRUE) |> process_reaction(), NA)
})
