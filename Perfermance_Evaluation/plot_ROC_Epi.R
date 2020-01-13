library(Epi)

data_A <- read.table("a1_MSI.txt", header=T, sep="\t")
x <- data_A$MSI_score
y <- data_A$MSI_status
z <- data_A$total_loci
rc <- ROC(form = y ~ x +z , plot="sp") 

## optimal combination
opt <- which.max(rowSums(rc$res[, c("sens", "spec")]))

## optimal cut-off point 
rc$res$lr.eta[opt]

ROC(form = y ~ x + z, plot = "ROC", MX = TRUE)


## ref https://stackoverflow.com/questions/23131897/how-can-i-get-the-optimal-cutoff-point-of-the-roc-in-logistic-regression-as-a-nu
## http://www.talkstats.com/threads/the-optimal-cutoff-score-in-the-classification-table.56212/

## https://smart-statistics.com/handling-roc-curves/
## http://ethen8181.github.io/machine-learning/unbalanced/unbalanced.html
