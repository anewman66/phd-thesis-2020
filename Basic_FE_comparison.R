FE.comp <- function (x) {
  comparisons <- x
  fisher.test(matrix(comparisons,nr=2))$p.value
}

#Put your data in for groups a and b
a.with <- 22
a.total <- 35
b.with <- 15
b.total <- 60

#Run the test
FE.comp(c(a.with,(a.total-a.with),b.with,(b.total-b.with)))