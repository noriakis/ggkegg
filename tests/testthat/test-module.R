test_that("Module parsing to text and network without errors", {
	expect_error( create_test_module(), NA)
	mod <- create_test_module()
	
	## Text parsing
	expect_error( module_text(mod), NA)
	expect_error( module_text(mod) |> plot_module_text(), NA)
    mod <- module("M00004")
    expect_error( module_text(mod), NA)
    expect_error( module_text(mod) |> plot_module_text(), NA)
    
    ## Network parsing
    expect_error( obtain_sequential_module_definition(mod) |>
        plot_module_blocks(), NA)
})
