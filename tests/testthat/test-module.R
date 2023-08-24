test_that("Module parsing without errors", {
	expect_error( create_test_module(), NA)
	mod <- create_test_module()
	expect_error( module_text(mod), NA)
	expect_error( module_text(mod) |> plot_module_text(), NA)
    mod <- module("M00004")
    expect_error( module_text(mod), NA)
    expect_error( module_text(mod) |> plot_module_text(), NA)
})
