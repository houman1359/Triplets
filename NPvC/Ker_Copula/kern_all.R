
library(kdecopula)
library(R.matlab)

args <- commandArgs(trailingOnly = TRUE)

temp.name <- args[[1]]

input.filename <- paste(temp.name, ".in.mat", sep="")
output.filename <- paste(temp.name, ".out.mat", sep="")

input <- readMat(input.filename)
u   <- input$e
ker <- input$ke
met <- input$method
bb  <- input$bw
dim <- input$dim

## bb = matrix(c(5, 1, 1, 5), nrow=2, ncol=2)
## bb <- bw_tll(ker,2)
## print(bb)


######## ino 19 july add kardam, bara test
## ker <- apply(ker, 2, function(x) rank(x)/(length(x)+1)) 

if (nrow(bb) == 1){
cop <- kdecop(ker, bw = NA, mult = 0.8, method = met, knots =100, info = FALSE, renorm.iter = 2L)
} else {
cop <- kdecop(ker, bw = NA, mult = 0.8, method = met, knots =100, info = FALSE, renorm.iter = 2L)
}

## logLK <- logLik(cop)
## aic <- AIC(cop)
## cdf <- pkdecop(u,cop)

pd <- dkdecop(u,cop)
ccdf <- hkdecop(u,cop,dim,inverse=FALSE)

pd_ker <- dkdecop(ker,cop)
ccdf_ker <- hkdecop(ker,cop,dim,inverse=FALSE)

## writeMat(output.filename,AIC=aic,loglk=logLK,pdf=pd,ccdf=ccdf,pdf_ker=pd_ker,ccdf_ker=ccdf_ker)

writeMat(output.filename,pdf=pd,ccdf=ccdf,pdf_ker=pd_ker,ccdf_ker=ccdf_ker)