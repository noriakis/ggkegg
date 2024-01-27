test_that("Network parsing without errors", {
	expect_error( create_test_network(), NA)
	net <- create_test_network()
	expect_error( net %>% network_graph() %>% plot_kegg_network(), NA)
    net <- network("N00002")
    expect_error( net %>% network_graph() %>% plot_kegg_network(), NA)
})
