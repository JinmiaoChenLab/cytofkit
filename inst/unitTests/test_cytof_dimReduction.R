test_cytof_dimReduction <- function() {
    data(iris)
    in_data <- iris[, 1:4]
    out_data <- cytof_dimReduction(in_data, method = "tsne")
    checkEquals(dim(out_data)[2], 2)
}