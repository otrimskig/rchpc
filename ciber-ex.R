library(CIBERSORT)
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
mixture_file <- system.file("extdata", "exampleForLUAD.txt", package = "CIBERSORT")
results <- cibersort(sig_matrix, mixture_file)
