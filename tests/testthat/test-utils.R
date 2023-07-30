test_that("Do utils without errors", {
    graph <- create_test_pathway()
    res <- data.frame(row.names="6737",log2FoldChange=1.2)
    expect_error( graph |>
        mutate(num=assign_deseq2(res, gene_type = "ENTREZID")),
        NA)
})
