test_that("Module parsing without errors", {
	expect_error( create_test_network(), NA)
	net <- create_test_network()
	expect_error( net |> network_graph() |> plot_kegg_network(), NA)
})
