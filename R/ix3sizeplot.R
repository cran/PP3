ix3sizeplot <-
function(PP3manyobj, main="Projection Index Value Histogram", nbig=10){


hist(PP3manyobj$ix3, main=main, xlab="Projection Index", freq=FALSE)

lines(density(PP3manyobj$pseudp.vals), col=2)

pseudp.vals <- PP3manyobj$pseudp.vals

vvals <- c(quantile(pseudp.vals, prob=c(0.5, 0.75, 0.9)), max(pseudp.vals))

abline(v=vvals, lty=2, col=2)


bigix <- sort.list(PP3manyobj$ix3, decreasing=TRUE)[1:nbig]

cat("Big Projection Indices\n")
cat("Maximum Psuedo-index: ", max(pseudp.vals), "\n")
cat("Index Number and associated projection indices\n")
print(PP3manyobj$ix3[bigix])


}
