#cat("What's your name?")
#x <- readLines(file("stdin"), 1)
#print(x)
x <- "test"
write.table(x, snakemake@output[[1]])
