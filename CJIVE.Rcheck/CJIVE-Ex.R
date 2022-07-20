pkgname <- "CJIVE"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "CJIVE-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('CJIVE')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("GenToyDatBinRank")
### * GenToyDatBinRank

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: GenToyDatBinRank
### Title: Generate 'Toy' Data
### Aliases: GenToyDatBinRank

### ** Examples

ToyDat = GenToyDatBinRank(n = 200, p1 = 2000, p2 = 1000, JntVarEx1 = 0.05, JntVarEx2 = 0.05,
                           IndVarEx1 = 0.25, IndVarEx2 = 0.25, jnt_rank = 1, equal.eig = FALSE,
                           ind_rank1 = 2, ind_rank2 = 3, SVD.plots = TRUE, Error = TRUE,
                           print.cor = TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("GenToyDatBinRank", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MatVar")
### * MatVar

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MatVar
### Title: Matrix variation (i.e. Frobenius norm)
### Aliases: MatVar

### ** Examples

X = matrix(rnorm(10), 5,2)
MatVar(X)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MatVar", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MatVar2")
### * MatVar2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MatVar2
### Title: Alternative calculation - Matrix variation (i.e. Frobenius norm)
### Aliases: MatVar2

### ** Examples

X = matrix(rnorm(10), 5,2)
MatVar2(X)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MatVar2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cc.jive")
### * cc.jive

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cc.jive
### Title: Canonical (Correlation) JIVE
### Aliases: cc.jive

### ** Examples

#Assign sample size and the number of features in each dataset
n = 200 #sample size
p1 = 2000 #Number of features in data set X1
p2 = 1000 #Number of features in data set X2

# Assign values of joint and individual signal ranks
r.J = 1 #joint rank
r.I1 = 2 #individual rank for data set X1
r.I2 = 2 #individual rank for data set X2


# Simulate data sets
ToyDat = GenToyDatBinRank(n = 200, p1 = p1, p2 = p2, JntVarEx1 = 0.05, JntVarEx2 = 0.05,
                           IndVarEx1 = 0.25, IndVarEx2 = 0.25, jnt_rank = r.J, equal.eig = FALSE,
                           ind_rank1 = r.I1, ind_rank2 = r.I2, SVD.plots = TRUE, Error = TRUE,
                           print.cor = TRUE)
# Store simulated data sets in an object called 'blocks'
blocks <- ToyDat$'Data Blocks'

# Save Subject scores as R objects
JntScores = ToyDat[['Scores']][['Joint']]
IndivScore.X = ToyDat[['Scores']][["Indiv_1"]]
IndivScore.Y = ToyDat[['Scores']][["Indiv_2"]]

# Save joint variable loadings as R objects
JntLd.X = t(ToyDat$Loadings$Joint_1)
JntLd.Y = t(ToyDat$Loadings$Joint_2)

# Save individual variable loadings as R objects
IndivLd.X =t(ToyDat$Loadings$Indiv_1)
IndivLd.Y = t(ToyDat$Loadings$Indiv_2)

# Save joint, individual, and noise signal matrices as R objects
JX = ToyDat[[1]]$J1
JY = ToyDat[[1]]$J2
IX = ToyDat[[1]]$I1
IY = ToyDat[[1]]$I2
EX = ToyDat[[1]]$E1
EY = ToyDat[[1]]$E2


## Check that proportions of variation explained are (approximately) equal to intended values
JVE.X = MatVar(JX)/MatVar(blocks[[1]])
JVE.Y = MatVar(JY)/MatVar(blocks[[2]])

IVE.X = MatVar(IX)/MatVar(blocks[[1]])
IVE.Y = MatVar(IY)/MatVar(blocks[[2]])

TotVE.X = MatVar((JX + IX))/MatVar(blocks[[1]])
TotVE.Y = MatVar((JY + IY))/MatVar(blocks[[2]])


CJIVE.res = cc.jive(blocks, c(r.I1,r.I2)+r.J, r.J, perm.test = FALSE)
# CJIVE signal matrix estimates
J.hat = CJIVE.res$sJIVE$joint_matrices
I.hat = CJIVE.res$sJIVE$indiv_matrices

# CJIVE loading estimates
WJ = lapply(J.hat, function(x) x[['v']])
WI = lapply(I.hat, function(x) x[['v']])

# Plots of CJIVE estimates against true counterparts and include an estimate of their chordal norm
layout(matrix(1:6,2, byrow = TRUE))
plot(JntScores, CJIVE.res$CanCorRes$Jnt_Scores, xlab = "True Joint Scores",
    ylab = "CJIVE Joint Scores",
    sub = paste0("Chordal Norm = ",
                 round(chord.norm.diff(JntScores, CJIVE.res$CanCorRes$Jnt_Scores), 3)))
plot(JntLd.X, WJ[[1]][,1], xlab = "True Joint Loadings X", ylab = "CJIVE Joint Loadings X",
    sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, WJ[[1]][,1]), 3)))
plot(JntLd.Y, WJ[[2]][,1], xlab = "True Joint Loadings Y", ylab = "CJIVE Joint Loadings Y",
    sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, WJ[[2]][,1]), 3)))
plot.new(); legend("left", paste("Comp.", 1:2), pch = 1, col  = c("orange", "green"),bty = "n" )
plot(IndivLd.X, WI[[1]][,1:2], xlab = "True Individual Loadings X",
    ylab = "CJIVE Individual Loadings X",
    col = c(rep("orange",p1), rep("green",p2)),
    sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, WI[[1]][,1:2]), 3)))
plot(IndivLd.Y, WI[[2]][,1:2], xlab = "True Individual Loadings Y",
    ylab = "CJIVE Individual Loadings Y",
    col = c(rep("orange",p1), rep("green",p2)),
    sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, WI[[2]][,1:2]), 3)))
layout(1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cc.jive", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
