#############################################################################################################################
### Compare methods of implementing JIVE analysis on two different pairs of toy datasets ####################################
### Author: Raphiel J. Murden                                                            ####################################
### Supervised by Benjamin Risk                                                          ####################################
#############################################################################################################################
rep_number = 0
r.J = 3
r.I1 = 2
r.I2 = 2
#outdir = args[2]
n = 200
p1 = 200
p2 = 200 ####Note that we have p2 = 1000 here as opposed to p2 = 10,000 in simulations
JntVarEx1 = 0.5
JntVarEx2 = 0.5
#files = list.files(outdir)
IndVarEx1 = 0.25
IndVarEx2 = 0.25

###Load iJIVE package (Lock et al.)
###This package has the show.image function which is great for displaying the heatmaps that represent the dataset
#library(r.jive)
library(ggplot2)
library(xtable)
library(dplyr)
# library(ajive)
library(scales)
prog.dir = "C:/Users/rmurd/Documents/PhD Researach/P-JIVE/Programs/"
ajive.dir = "C:/Users/rmurd/Documents/PhD Researach/Applications/r_AJIVE/R"
files= list.files(ajive.dir)

source(file.path(prog.dir, "Functions_for_CJIVE.R"))
for (i in files) source(file.path(ajive.dir, i))

#######################################################################
set.seed(rep_number) 
###Construct Datasets

#######JIVE Implementations for Toy Data No ################################################################################
###Setup input parameters for JIVE implementations
true_signal_ranks.1 = r.J + c(r.I1,r.I2)                          ##true ranks of overall signals
equal.eig = F
ToyDat = GenToyDatBinRank(n, p1, p2, JntVarEx1, JntVarEx2, IndVarEx1, IndVarEx2, jnt_rank = r.J, equal.eig = equal.eig,
                           ind_rank1 = r.I1, ind_rank2 = r.I2, JntVarAdj = T, SVD.plots = T, Error = T, print.cor = F)
blocks <- lapply(ToyDat[[2]], scale)
blocks.t <- lapply(blocks,t)   ##iJIVE requires datasets organized as p-by-n matrices

JntScores = ToyDat[['Scores']][['Joint']]
D.J = (equal.eig + 1)*diag(r.J:1) + equal.eig*diag(rep(1,r.J))

IndivScore.X = ToyDat[['Scores']][["Indiv_1"]]
D.IX = (equal.eig + 1)*diag((r.I1:1) + 0.5) + equal.eig*diag(rep(1,r.I1) + 0.5)
IndivScore.Y = ToyDat[['Scores']][["Indiv_2"]]

JntLd.X = t(ToyDat[['Loadings']][["Joint_1"]])
JntLd.Y = t(ToyDat[['Loadings']][["Joint_2"]])
IndivLd.X =t(ToyDat[['Loadings']][["Indiv_1"]])
IndivLd.Y = t(ToyDat[['Loadings']][["Indiv_2"]])

JX = ToyDat[[1]]$J1
JY = ToyDat[[1]]$J2
IX = ToyDat[[1]]$I1
IY = ToyDat[[1]]$I2
EX = ToyDat[[1]]$E1
EY = ToyDat[[1]]$E2

AX = JX + IX
AY = JY + IY
X = blocks[[1]]
Y = blocks[[2]]

plot(svd(cbind(svd(AX, nu = true_signal_ranks.1[1])[['u']], svd(AY, nu = true_signal_ranks.1[2])[['u']]))$d^2)
plot(svd(cbind(svd(X, nu = true_signal_ranks.1[1])[['u']], svd(Y, nu = true_signal_ranks.1[2])[['u']]))$d^2)

JntBlock = cbind(JX,JY)

JVE.X = MatVar(JX)/MatVar(blocks[[1]])
JVE.Y = MatVar(JY)/MatVar(blocks[[2]])

IVE.X = MatVar(IX)/MatVar(blocks[[1]])
IVE.Y = MatVar(IY)/MatVar(blocks[[2]])

TotVE.X = MatVar((JX + IX))/MatVar(blocks[[1]])
TotVE.Y = MatVar((JY + IY))/MatVar(blocks[[2]])

####################################################################
################# Apply aJIVE package with correct ranks specified
####################################################################
ajive.oracle<-ajive(blocks,initial_signal_ranks = true_signal_ranks.1, joint_rank = r.J)
cctrue.jive.res = cc.jive(blocks, signal.ranks = true_signal_ranks.1, joint.rank = r.J, perm.test = FALSE)

a.rJ.c = ajive.oracle$joint_rank
cc.rJ = cctrue.jive.res[['CanCorRes']][['Jnt_Rank']]

a.rI.c1 = ajive.oracle$block_decomps[[1]]$individual$rank
a.rI.c2 = ajive.oracle$block_decomps[[2]]$individual$rank

## Joint Scores
a.c.JntScores.hat = ajive.oracle$joint_scores
cc.JntScores.hat = cctrue.jive.res$CanCorRes$Jnt_Scores
round(chord.norm.diff(a.c.JntScores.hat, scale(JntScores)), 4)
round(chord.norm.diff(cc.JntScores.hat, scale(JntScores)), 4)

a.c.W.JX = ajive.oracle$block_decomps[[1]]$joint$v
cc.W.JX = cctrue.jive.res$sJIVE$joint_matrices[[1]]$v
round(chord.norm.diff(scale(cc.W.JX), scale(JntLd.X)), 4)
round(chord.norm.diff(scale(a.c.W.JX), scale(JntLd.X)), 4)

a.c.W.JY = ajive.oracle$block_decomps[[2]]$joint$v
cc.W.JY = cctrue.jive.res$sJIVE$joint_matrices[[2]]$v
round(chord.norm.diff(scale(cc.W.JY), scale(JntLd.Y)), 4)
round(chord.norm.diff(scale(a.c.W.JY), scale(JntLd.Y)), 4)

a.c.bX.hat = ajive.oracle$block_decomps[[1]][["individual"]][["u"]]
cc.bX.hat = cctrue.jive.res$sJIVE$indiv_matrices[[1]]$u
round(chord.norm.diff(IndivScore.X, cc.bX.hat),4)
round(chord.norm.diff(IndivScore.X, a.c.bX.hat),4)

a.c.bY.hat = ajive.oracle$block_decomps[[2]][["individual"]][["u"]]
cc.bY.hat = cctrue.jive.res$sJIVE$indiv_matrices[[2]]$u
round(chord.norm.diff(IndivScore.Y, cc.bY.hat),4)
round(chord.norm.diff(IndivScore.Y, a.c.bY.hat),4)


plot(JntScores[,1], 7*a.c.JntScores.hat[,1], xlab = "True Joint Scores", ylab = "Estimated Joint Scores", main = "AJIVE Oracle")
# plot(JntScores[,2], 7*a.c.JntScores.hat[,2], xlab = "True Joint Scores", ylab = "Estimated Joint Scores", main = "AJIVE Oracle")
# plot(JntScores[,3], 7*a.c.JntScores.hat[,3], xlab = "True Joint Scores", ylab = "Estimated Joint Scores", main = "AJIVE Oracle")


a.c.W.IX.hat = ajive.oracle$block_decomps[[1]][["individual"]][["v"]]
cc.W.IX.hat = cctrue.jive.res$sJIVE$indiv_matrices[[1]]$v

a.c.W.IY.hat = ajive.oracle$block_decomps[[2]][["individual"]][["v"]]
cc.W.IY.hat = cctrue.jive.res$sJIVE$indiv_matrices[[2]]$v

############

########################
#Do aJIVE and ccJIVE joint subject scores match?
plot(cc.JntScores.hat, a.c.JntScores.hat, col = rep(2:4, each = n))
cor(cc.JntScores.hat, a.c.JntScores.hat)
chord.norm.diff(cc.JntScores.hat[,1], a.c.JntScores.hat[,1])
plot(cc.JntScores.hat2[,1], a.c.JntScores.hat[,1])
# plot(cc.JntScores.hat[,2], a.c.JntScores.hat[,2])

pdf(file = file.path(imgs.fldr, "CCJIVE_OracleLoadingExample.pdf"))
layout(matrix(1:4, ncol =2))
plot(JntScores, cc.JntScores.hat[,1], main = "CC JIVE Oracle Joint Scores", cex.sub = 0.7, xlim = c(-0.5,1.5),
     sub = paste("Norm = ", round(cc.jnt.subnorm, 3), sep = ""),
     xlab = "True Subject Scores", ylab = "Estimated Subject Scores", col = "blue", pch = 1)
# points(JntScores, cc.JntScores.hat[,2], col = "red", pch = 1)
legend("bottomright", c(paste("Comp. #", 1), paste("Comp. #", 2)), pch = c(1,1),
       col = c("blue", "red"), cex = 0.7)
plot(JntLd.X[,1], svd(cc.JX.hat)[['v']][,1], main = "CC JIVEOracle Joint Loadings", cex.sub = 0.7, ylim = c(-.55,.35), xlim = c(-3.5,3.5),
     sub = paste("Norm ", expression("X"[1]), " = ", round(cc.jnt.loadX.norm, 3), 
                 "; Norm ", expression("X"[2]), " = ", round(cc.jnt.loadY.norm, 3), sep = ""),
     xlab = "True Joint Loadings", ylab = "SVD of Estimated Joint Signal", col = "blue", pch = 1)
points(JntLd.Y[,1], svd(cc.JY.hat)[['v']][,1], col = "red", pch = 1)
legend("bottomright", c(paste(expression("X"[1]), " Comp. #", 1), paste(expression("X"[2]), " Comp. #", 1)), pch = c(1,1),
       col = c("blue", "red"), cex = 0.7)
# for(i in 2:3){
#   points(JntLd.X[,i], svd(cc.JX.hat)[['v']][,i], col = "blue", pch = i)
#   points(JntLd.Y[,i], svd(cc.JY.hat)[['v']][,i], col = "red", pch = i)
# }
# legend("bottomright", c(paste(expression("X"[1]), " Comp. #", 1:3), paste(expression("X"[2]), " Comp. #", 1:3)), pch = rep(1:3, 2), 
#        col = rep(c("blue", "red"), each = 3), cex = 0.56, bty = 'n')

plot(IndivScore.X[,1], svd(cc.IX.hat)[['u']][,1], main = "CC JIVEOracle Individual Subject Scores", 
     xlab = "True Subject Scores", ylab = "SVD of Estimated Indiv Signal", col = "blue", pch = 1, cex.sub = 0.7,
     sub = paste("Norm ", expression("X"[1]), " = ", round(cc.indiv.X.subnorm, 3), 
                 "; Norm ", expression("X"[2]), " = ", round(cc.indiv.Y.subnorm, 3), sep = ""))
points(IndivScore.X[,2], svd(cc.IX.hat)[['u']][,2], col = "blue", pch = 2)
points(IndivScore.Y[,1], svd(cc.IY.hat)[['u']][,1], col = "red", pch = 1)
points(IndivScore.Y[,2], svd(cc.IY.hat)[['u']][,2], col = "red", pch = 2)
legend("bottomright", c(paste(expression("X"[1]), " Comp. #", 1:2), paste(expression("X"[2]), " Comp. #", 1:2)), pch = c(1,2,1,2), 
       col = rep(c("blue", "red"), each = 2), cex = 0.7)

###Note in the r.J = 3 example, the first component is esimtated rather well. The second component here does not match the truth

plot(IndivLd.X[,1], svd(cc.IX.hat)[['v']][,1], main = "CC JIVEOracle Individual Loadings",
     xlab = "True Indiv. Loadings", ylab = "SVD of Estimated Indiv Signal", col = "blue", pch = 1, cex.sub = 0.7,
     sub = paste("Norm ", expression("X"[1]), " = ", round(cc.indiv.X.loadnorm, 3), 
                 "; Norm ", expression("X"[2]), " = ", round(cc.indiv.Y.loadnorm, 3), sep = ""))
points(IndivLd.X[,2], svd(cc.IX.hat)[['v']][,2], col = "blue", pch = 2)
points(IndivLd.Y[,1], svd(cc.IY.hat)[['v']][,1], col = "red", pch = 1)
points(IndivLd.Y[,2], svd(cc.IY.hat)[['v']][,2], col = "red", pch = 2)
legend("bottomright", c(paste(expression("X"[1]), " Comp. #", 1:2), paste(expression("X"[2]), " Comp. #", 1:2)), pch = c(1,2,1,2), 
       col = rep(c("blue", "red"), each = 2), cex = 0.7)
dev.off()

plot(IndivScore.X[,1], round(a.c.bX.hat[,1], 1), main = "AJIVE Oracle Individual Subject Scores", 
     xlab = "True Subject Scores", ylab = "SVD of Estimated Indiv Signal", col = alpha("blue", 0.75), pch = 0, cex.sub = 0.7,
     sub = paste("Norm ", expression("X"[1]), " = ", round(a.c.indiv.X.subnorm, 3), 
                 "; Norm ", expression("X"[2]), " = ", round(a.c.indiv.Y.subnorm, 3), sep = ""))
points(IndivScore.X[,2], round(a.c.bX.hat[,2],1), col = alpha("blue", 0.75), pch = 3)

plot(IndivLd.X[,1], a.c.W.IX.hat[,1], main = "AJIVE Oracle Individual Loadings",
     xlab = "True Indiv. Loadings", ylab = "SVD of Estimated Indiv Signal", col = alpha("blue", 0.75), pch = 0, cex.sub = 0.7,
     sub = paste("Norm ", expression("X"[1]), " = ", round(a.c.indiv.X.loadnorm, 3), 
                 "; Norm ", expression("X"[2]), " = ", round(a.c.indiv.Y.loadnorm, 3), sep = ""))
points(IndivLd.X[,2], a.c.W.IX.hat[,2], col = alpha("blue", 0.75), pch = 3)
plot(IndivLd.Y[,1], a.c.W.IY.hat[,1], col = alpha("red", 0.75), pch = 0)
plot(IndivLd.Y[,2], a.c.W.IY.hat[,2], col = alpha("red", 0.75), pch = 3)

###One thing to note: chordal norms seem to treat rank-1 subspaces differently than they do rank-2 estimates. 

#Show a plot that exhibits relationship between correlation and chordal norm
corrs = seq(0,1, len = 200)
norms  = sqrt((1 - corrs^2))
pdf(file = file.path(imgs.fldr, "PearsonCorr_v_ChordalNorm.pdf"))
plot(corrs, norms, xlab = expression(paste("Pearson Correlation (", rho, ")", sep = "")), 
     ylab = expression(paste("Chordal Norm (", delta, ")", sep = "")),
     main = "Relationship between measures")
chord.90= round(sqrt((1 - 0.9^2)),2)
chord.95= round(sqrt((1 - 0.95^2)),2)
chord.99= round(sqrt((1 - 0.99^2)),2)
lines(corrs, norms)
segments(x0 = 0, y0 = chord.90, x1 = 0.9, chord.90, col = 'red')
segments(x0 = 0.9, y0 = 0, x1 = 0.9, chord.90, col = 'red')
segments(x0 = 0, y0 = chord.95, x1 = 0.95, chord.95, col = 'blue')
segments(x0 = 0.95, y0 = 0, x1 = 0.95, chord.95, col = 'blue')
segments(x0 = 0, y0 = chord.99, x1 = 0.99, chord.99, col = 'green')
segments(x0 = 0.99, y0 = 0, x1 = 0.99, chord.99, col = 'green')
legend("bottomleft", legend = c(expression(paste(delta, '=', 0.44, ';', rho, '=', 0.90)),
                                expression(paste(delta, '=', 0.31, ';', rho, '=', 0.95)),
                                expression(paste(delta, '=', 0.14, ';', rho, '=', 0.99))),
       pch = '-', col = c('red', "blue", 'green'))
dev.off()

corrs = seq(0,1, len = 200)
cor2.v.chord.mat = matrix(NA, nrow = length(corrs), ncol = length(corrs))
for(x in 1:length(corrs)){
  for (y in 1:length(corrs)){
    sig = c(corrs[x],corrs[y])
    z = sqrt(sum((1 - sig^2))/2)
    cor2.v.chord.mat[x,y] = z
  }
}

pdf(file = file.path(imgs.fldr, "PearsonCorr2_v_ChordalNorm.pdf"))
image.plot(x = corrs, y = corrs, z = cor2.v.chord.mat, nlevel = 100,
           xlab = expression(paste("Comp. #1 (", rho[1], ")", sep = "")), 
           ylab = expression(paste("Comp. #2 (", rho[2], ")", sep = "")), 
           main = "Relationship between measures", 
           legend.lab = expression(paste("Chordal Norm (", delta, ")", sep = "")))
dev.off()

corrs = seq(0,1, len = 200)
norms2  = acos(corrs)
plot(corrs, norms2, xlab = expression(paste("Pearson Correlation (", rho, ")", sep = "")), 
     ylab = expression(paste("Alpha Norm (", delta, ")", sep = "")),
     main = "Relationship between measures")

corrs = seq(0,1, len = 200)
norms3  = log(corrs^-1)
plot(corrs, norms3, xlab = expression(paste("Pearson Correlation (", rho, ")", sep = "")), 
     ylab = expression(paste("Mu Norm (", delta, ")", sep = "")),
     main = "Relationship between measures")

# dat = data.frame(Est.JntLd.X = svd(a.c.JX.hat)[['v']][,1], Est.JndLd.Y = svd(a.c.JY.hat)[['v']][,1], 
#                   True.JntLd.X = svd(JX)[['v']][,1], True.JntLd.Y = svd(JY)[['v']][,1])
# 
# melt.dat = melt(dat)
# nm = as.character(melt.dat$variable)
# melt.dat[,"Var"] = factor(substr(nm, nchar(nm), nchar(nm)))
# melt.dat[,"Type"] = factor(substr(nm, 1, 3))
# 
# plot.1 = ggplot(data = melt.dat, aes())

# plot(svd(JX)[['v']][,1], svd(a.c.JX.hat)[['v']][,1], main = "AJIVE Oracle Joint Loadings", ylim = c(-.55, .3), xlim = c(-.55, .3),
#      xlab = "SVD of True Joint Signal", ylab = "SVD of Estimated Joint Signal", col = "blue", pch = 1)
# points(svd(JY)[['v']][,1], svd(a.c.JY.hat)[['v']][,1], col = "red", pch = 1)
# legend("bottomleft", c(paste(expression("X"[1]), " Comp. #", 1, paste("(", round(a.c.jnt.loadX.norm, 3), ")", sep = "")), 
#                        paste(expression("X"[2]), " Comp. #", 1, paste("(", round(a.c.jnt.loadX.norm, 3), ")", sep = ""))),
#                        pch = c(1,1), col = c("blue", "red"), cex = 0.7)
# 
# ###Does 'old definition' of chordal norm agree?
# old.chord.norm.diff(t(JX), t(a.c.JX.hat), 1)  # = 0.177; Not, it does not match
# old.chord.norm.diff(t(JX), t(a.c.JX.hat), 10)  # = 0.7692; Not, it does not match
# old.chord.norm.diff(t(JX), t(a.c.JX.hat), 20)  # = 0.7692; Not, it does not match
# old.chord.norm.diff(t(JX), t(a.c.JX.hat), 200)  # =~ 0 Now it matches

##Scatter Plots of oracle Loadings  with AJIVE-Oracle estimated loadings (estimated rank was correct)
# pdf(file = file.path(imgs.fldr, "Sim Scatter Plot - AJIVE-Oracle Loadings.pdf"))
# layout(matrix(1:3, 3))
# for (i in 1:3){
#   plot(V.X.Oracle[,i],V.X.Estimate[,i], xlab = paste("Oracle Loading #", i), ylab = paste("AJIVE-Oracle Loading #", i))
# }
# dev.off()
# 
