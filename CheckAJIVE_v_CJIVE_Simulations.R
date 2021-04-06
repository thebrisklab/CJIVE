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
round(chord.norm.diff(a.c.JntScores.hat, cc.JntScores.hat), 6)

# Joint loadings for X
a.c.W.JX = ajive.oracle$block_decomps[[1]]$joint$v
cc.W.JX = cctrue.jive.res$sJIVE$joint_matrices[[1]]$v
chord.norm.diff(a.c.W.JX, cc.W.JX)

# Joint loadings for Y
a.c.W.JY = ajive.oracle$block_decomps[[2]]$joint$v
cc.W.JY = cctrue.jive.res$sJIVE$joint_matrices[[2]]$v
chord.norm.diff(a.c.W.JY, cc.W.JY)

# Individual Scores for X
a.c.bX.hat = ajive.oracle$block_decomps[[1]][["individual"]][["u"]]
cc.bX.hat = cctrue.jive.res$sJIVE$indiv_matrices[[1]]$u
chord.norm.diff(a.c.bX.hat, cc.bX.hat)

# Individual Scores for Y
a.c.bY.hat = ajive.oracle$block_decomps[[2]][["individual"]][["u"]]
cc.bY.hat = cctrue.jive.res$sJIVE$indiv_matrices[[2]]$u
chord.norm.diff(a.c.bY.hat, cc.bY.hat)

# Individual Loadings for X
a.c.W.IX.hat = ajive.oracle$block_decomps[[1]][["individual"]][["v"]]
cc.W.IX.hat = cctrue.jive.res$sJIVE$indiv_matrices[[1]]$v
chord.norm.diff(cc.W.IX.hat, a.c.W.IX.hat)

# Individual Loadings for Y
a.c.W.IY.hat = ajive.oracle$block_decomps[[2]][["individual"]][["v"]]
cc.W.IY.hat = cctrue.jive.res$sJIVE$indiv_matrices[[2]]$v
chord.norm.diff(cc.W.IY.hat, a.c.W.IY.hat)
