test_that("Highlight functions without errors", {
	graph <- create_test_pathway()
	expect_error( graph |> mutate(hl=highlight_set_nodes(c("hsa:51428"))), NA)
	expect_error( graph |> activate("edges") |>
        mutate(hl=highlight_set_edges(c("degradation"),
            name="subtype_name")), NA)
})


