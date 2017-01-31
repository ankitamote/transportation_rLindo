#------------------------------------------------------------------------#
# Transportation Problem
#------------------------------------------------------------------------#

#using gdata package to import ".xlsx" file in R
library(gdata)
library(rLindo)
library(pracma)
library(stringi)


#create LINDO enviroment object
rEnv <- rLScreateEnv()
#create LINDO model object
rModel <- rLScreateModel(rEnv)
#load LP data
data <- read.xls("TRANS_Lindo.xlsx", sheet=1, as.is=TRUE)
#str(data)
print(data)
nVars = 12                                            # define number of variables in your problem
nCons = length(data$RHS) - 1                         # find number of constraints
nDir = LS_MIN                                        # define objective of problem i.e MIN or MAX
dObjConst = 0.
a11 <- as.vector(data$x11[1])
b11 <- as.vector(data$x11[-1])
a12 <- as.vector(data$x12[1])
b12 <- as.vector(data$x12[-1])
a13 <- as.vector(data$x13[1])
b13 <- as.vector(data$x13[-1])
a14 <- as.vector(data$x14[1])
b14 <- as.vector(data$x14[-1])
a21 <- as.vector(data$x21[1])
b21 <- as.vector(data$x21[-1])
a22 <- as.vector(data$x22[1])
b22 <- as.vector(data$x22[-1])
a23 <- as.vector(data$x23[1])
b23 <- as.vector(data$x23[-1])
a24 <- as.vector(data$x24[1])
b24 <- as.vector(data$x24[-1])
a31 <- as.vector(data$x31[1])
b31 <- as.vector(data$x31[-1])
a32 <- as.vector(data$x32[1])
b32 <- as.vector(data$x32[-1])
a33 <- as.vector(data$x33[1])
b33 <- as.vector(data$x33[-1])
a34 <- as.vector(data$x34[1])
b34 <- as.vector(data$x34[-1])

adC <- as.vector(cbind(a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34))                       # coefficients of objective function
mat <- cbind(b11,b12,b13,b14,b21,b22,b23,b24,b31,b32,b33,b34)
adB <- as.vector(data$RHS[-1])              # define RHS coefficients of the constraints
acConTypes <- data$X[-1]                    # define constraints type Less than, greater than, or equal to
acConTypes <- paste(acConTypes[1],acConTypes[2],acConTypes[3],acConTypes[4],acConTypes[5],acConTypes[6],acConTypes[7],collapse = NULL)
acConTypes <- stri_replace_all_charclass(acConTypes, "\\p{WHITE_SPACE}","")

nNZ <- nnz(mat)                                      # number of non zeroes in the RHS of the constraints
adA <- as.vector(mat)
adA <- adA[-which(adA==0)]                       # non zero coefficients of constaint matrix by column
anRowX <- row(mat)[which(mat!=0)] - 1                # row indices of non zeroes in the constraints matrix by column
pdLower <- c(rep(0,12))                                    # lower limit of variables
pdUpper <- NULL # upper limit

mat1 <- NULL
for (i in c(1:nVars)){ 
  a <- as.vector(mat[,i])
  a <- a[-which(a==0)]
  mat1 <- cbind(mat1,a)
}
a <- matrix(which(mat1!=0)-1, 2, 12)
anBegCol <- c(a[1, ],nNZ)

#load data into model
rLSloadLPData(rModel , nCons, nVars, nDir, dObjConst, adC, adB, acConTypes,
              nNZ, anBegCol, NULL, adA, anRowX, pdLower, pdUpper)

#solve the model
rLSoptimize(rModel,LS_METHOD_FREE)
#rLSgetLPData(rModel)
message("Solution is:")
print(rLSgetSolution(rModel, 0))
message("Reduced costs are: ")
print(rLSgetReducedCosts(rModel))
#get primal solution
message("Primal Solution Is: ")
print(rLSgetPrimalSolution(rModel))
#get dual solution
message("Dual Solution Is: ")
print(rLSgetDualSolution(rModel))
#retrieve information
rLSgetDInfo(rModel,LS_DINFO_POBJ)
rLSgetIInfo(rModel,LS_IINFO_MODEL_STATUS)
#get basis
rLSgetBasis(rModel)

#delete enviroment and model objects
#free memory
rLSdeleteModel(rModel)
rLSdeleteEnv(rEnv)
