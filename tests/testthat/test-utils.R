test_that("Do utils without errors", {
    graph <- create_test_pathway()
    res <- data.frame(row.names="6737",log2FoldChange=1.2)
    expect_error( graph %>%
        mutate(num=assign_deseq2(res, gene_type="ENTREZID")),
        NA)
    expect_error( graph %>% activate("edges") %>%
        mutate(num=edge_numeric_sum(c(1.2,-1.2) %>%
                                    setNames(c("TRIM21","DDX41")),
                                    name="graphics_name")),
        NA)
})

test_that("edge_matrix without errors", {
	graph <- create_test_pathway()
    num_df <- data.frame(row.names=c("6737","51428"),
                    "sample1"=c(1.1,1.2),
                    "sample2"=c(1.1,1.2),
                     check.names=FALSE)
    expect_error(graph <- graph %>% edge_matrix(num_df, gene_type="ENTREZID"),
    	NA)                     	
})

