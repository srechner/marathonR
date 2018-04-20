
library("marathonR")

rowsums <- c(2, 1, 2)
colsums <- c(1, 1, 1, 2)

SampleFixed(rowsums, colsums, 2, method="exact")
SampleFixed(rowsums, colsums, 2, method="ktv-switch", steps=100)
SampleFixed(rowsums, colsums, 2, method="edge-switch", steps=100)
SampleFixed(rowsums, colsums, 2, method="curveball", steps=100)

rowsum.lower <- c(1, 0, 1)
rowsum.upper <- c(2, 1, 1)
colsum.lower <- c(0, 1, 1, 0)
colsum.upper <- c(1, 1, 2, 3)

SampleInterval(rowsum.lower, rowsum.upper, colsum.lower, colsum.upper, 2, method="exact")
SampleInterval(rowsum.lower, rowsum.upper, colsum.lower, colsum.upper, 2, method="simple", steps=1000)
SampleInterval(rowsum.lower, rowsum.upper, colsum.lower, colsum.upper, 2, method="informed", steps=100)