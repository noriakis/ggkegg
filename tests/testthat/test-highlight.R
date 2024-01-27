test_that("Highlight node functions without errors", {
	graph <- create_test_pathway()
	expect_error(graph %>% mutate(hl=highlight_set_nodes(c("hsa:51428"))), NA)
})
test_that("Highlight edge functions without errors", {
	graph <- create_test_pathway()
	expect_error(graph %>% activate("edges") %>%
        mutate(hl=highlight_set_edges(c("degradation"),
            name="subtype_name")), NA)
})
test_that("Highlight module functions without errors", {
	graph <- create_test_pathway()
    mo <- create_test_module()
	expect_error(graph <- graph %>% highlight_module(mo), NA)
})



