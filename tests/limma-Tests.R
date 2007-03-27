library(limma)

set.seed(0); u <- runif(100)

### strsplit2

x <- c("ab;cd;efg","abc;def","z","")
strsplit2(x,split=";")

### removeext

removeExt(c("slide1.spot","slide.2.spot"))
removeExt(c("slide1.spot","slide"))

### printorder
printorder(list(ngrid.r=4,ngrid.c=4,nspot.r=8,nspot.c=6),ndups=2,start="topright",npins=4)
printorder(list(ngrid.r=4,ngrid.c=4,nspot.r=8,nspot.c=6))

### merge.rglist

R <- G <- matrix(11:14,4,2)
rownames(R) <- rownames(G) <- c("a","a","b","c")
RG1 <- new("RGList",list(R=R,G=G))
R <- G <- matrix(21:24,4,2)
rownames(R) <- rownames(G) <- c("b","a","a","c")
RG2 <- new("RGList",list(R=R,G=G))
merge(RG1,RG2)
merge(RG2,RG1)

### background correction
RG <- new("RGList", list(R=c(1,2,3,4),G=c(1,2,3,4),Rb=c(2,2,2,2),Gb=c(2,2,2,2)))
backgroundCorrect(RG)
backgroundCorrect(RG, method="half")
backgroundCorrect(RG, method="minimum")
backgroundCorrect(RG, offset=5)

### normalizeWithinArrays

library(sma)
data(MouseArray)
MA <- normalizeWithinArrays(mouse.data, mouse.setup, method="robustspline")
MA$M[1:5,]
MA <- normalizeWithinArrays(mouse.data, mouse.setup)
MA$M[1:5,]

### normalizeBetweenArrays

MA <- normalizeBetweenArrays(MA,method="scale")
MA$M[1:5,]
MA$A[1:5,]
MA <- normalizeBetweenArrays(MA,method="quantile")
MA$M[1:5,]
MA$A[1:5,]

### unwrapdups

M <- matrix(1:12,6,2)
unwrapdups(M,ndups=1)
unwrapdups(M,ndups=2)
unwrapdups(M,ndups=3)
unwrapdups(M,ndups=2,spacing=3)

### trigammaInverse

trigammaInverse(c(1e-6,NA,5,1e6))

### lm.series, contrasts.fit, ebayes

M <- matrix(rnorm(10*6,sd=0.3),10,6)
M[1,1:3] <- M[1,1:3] + 2
design <- cbind(First3Arrays=c(1,1,1,0,0,0),Last3Arrays=c(0,0,0,1,1,1))
fit <- lm.series(M,design=design)
contrast.matrix <- cbind(First3=c(1,0),Last3=c(0,1),"Last3-First3"=c(-1,1))
fit2 <- contrasts.fit(fit,contrasts=contrast.matrix)
eb <- ebayes(fit2)

eb$t
eb$s2.prior
eb$s2.post
eb$df.prior
eb$lods
eb$p.value
eb$var.prior

### toptable

toptable(fit)

### topTable

fit <- lmFit(M,design)
fit2 <- eBayes(contrasts.fit(fit,contrasts=contrast.matrix))
topTable(fit2)
topTable(fit2,coef=3,resort.by="logFC")
topTable(fit2,coef=3,resort.by="p")
topTable(fit2,coef=3,sort="logFC",resort.by="t")
topTable(fit2,coef=3,resort.by="B")

### marray object

suppressMessages(suppressWarnings(gotmarray <- require(marray,quietly=TRUE)))
if(gotmarray) {
	data(swirl)
	snorm = maNorm(swirl)
	fit <- lmFit(snorm, design = c(1,-1,-1,1))
	fit <- eBayes(fit)
	topTable(fit,resort.by="AveExpr")
}

### duplicateCorrelation

cor.out <- duplicateCorrelation(M)
cor.out$consensus.correlation
cor.out$atanh.correlations

### gls.series

fit <- gls.series(M,design,correlation=cor.out$cor)
fit$coefficients
fit$stdev.unscaled
fit$sigma
fit$df.residual

### mrlm

fit <- mrlm(M,design)
fit$coef
fit$stdev.unscaled
fit$sigma
fit$df.residual

# Similar to Mette Langaas 19 May 2004
set.seed(123)
narrays <- 9
ngenes <- 5
mu <- 0
alpha <- 2
beta <- -2
epsilon <- matrix(rnorm(narrays*ngenes,0,1),ncol=narrays)
X <- cbind(rep(1,9),c(0,0,0,1,1,1,0,0,0),c(0,0,0,0,0,0,1,1,1))
dimnames(X) <- list(1:9,c("mu","alpha","beta"))
yvec <- mu*X[,1]+alpha*X[,2]+beta*X[,3]
ymat <- matrix(rep(yvec,ngenes),ncol=narrays,byrow=T)+epsilon
ymat[5,1:2] <- NA
fit <- lmFit(ymat,design=X)
test.contr <- cbind(c(0,1,-1),c(1,1,0),c(1,0,1))
dimnames(test.contr) <- list(1:3,c("alpha-beta","mu+alpha","mu+beta"))
fit2 <- contrasts.fit(fit,contrasts=test.contr)
eBayes(fit2)

### uniquegenelist

uniquegenelist(letters[1:8],ndups=2)
uniquegenelist(letters[1:8],ndups=2,spacing=2)

### classifyTests

tstat <- matrix(c(0,5,0, 0,2.5,0, -2,-2,2, 1,1,1), 4, 3, byrow=TRUE)
classifyTestsF(tstat)
FStat(tstat)
classifyTestsT(tstat)
classifyTestsP(tstat)
