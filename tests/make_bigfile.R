dosages = matrix(sample((-1):2, 1000*100000, replace = TRUE, prob = c(0.1, 0.5, 0.3, 0.1)), ncol = 1000)
colnames(dosages) = paste("s", 1:ncol(dosages), sep = "")
output = data.frame(vid = paste("1", 1:100000, "A:T", sep = ":"), rsid = paste("rs", 1:100000, sep = ""))
vaf = (rowSums(dosages) + rowSums(dosages == -1)) / (2*rowSums(dosages != -1))
output$vaf = vaf
output = cbind(output, dosages)
write.table(output, file = "bigfile.dosages", row.names=FALSE, col.names=TRUE, quote=FALSE, sep = "\t")
