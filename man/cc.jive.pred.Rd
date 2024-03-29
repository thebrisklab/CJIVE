% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CJIVE_Functions.R
\name{cc.jive.pred}
\alias{cc.jive.pred}
\title{CJIVE joint subject score prediction}
\usage{
cc.jive.pred(
  orig.dat.blocks,
  new.subjs,
  signal.ranks,
  cc.jive.loadings,
  can.cors
)
}
\arguments{
\item{orig.dat.blocks}{list of the two data matrices on which CJIVE was initially conducted}

\item{new.subjs}{list of two data matrices containing information on new subjects}

\item{signal.ranks}{a vector of length two which contains the rank for the signal within each data block. The rank corresponds to the number of principal
components (PCs) to be retained within each data block. If NULL, the ranks are determined by the parameter 'perc.var.' Default is NULL}

\item{cc.jive.loadings}{canonical loadings for the joint subspace}

\item{can.cors}{canonical correlations from the PCs of the data on which CJIVE was initially conducted - notated as rho_j in CJIVE manuscript}
}
\value{
matrix of joint subject score for new subjects
}
\description{
Predicts joint scores for new subjects based on CJIVE joint scores
}
\examples{
n = 200 #sample size
p1 = 100 #Number of features in data set X1
p2 = 100 #Number of features in data set X2
# Assign values of joint and individual signal ranks
r.J = 1 #joint rank
r.I1 = 2 #individual rank for data set X1
r.I2 = 2 #individual rank for data set X2
true_signal_ranks = r.J + c(r.I1,r.I2)
# Simulate data sets
ToyDat = GenerateToyData(n = n, p1 = p1, p2 = p2, JntVarEx1 = 0.05, JntVarEx2 = 0.05,
                          IndVarEx1 = 0.25, IndVarEx2 = 0.25, jnt_rank = r.J, equal.eig = FALSE,
                          ind_rank1 = r.I1, ind_rank2 = r.I2, SVD.plots = TRUE, Error = TRUE,
                          print.cor = TRUE)
# Store simulated data sets in an object called 'blocks'
blocks <- ToyDat$'Data Blocks'
# Split data randomly into two subsamples
rnd.smp = sample(n, n/2)
blocks.sub1 = lapply(blocks, function(x){x[rnd.smp,]})
blocks.sub2 = lapply(blocks, function(x){x[-rnd.smp,]})
# Joint scores for the two sub samples
JntScores.1 = ToyDat[['Scores']][['Joint']][rnd.smp]
JntScores.2 = ToyDat[['Scores']][['Joint']][-rnd.smp]
# Conduct CJIVE analysis on the first sub-sample and store the canonical loadings and canonical
# correlations
cc.jive.res_sub1 = cc.jive(blocks.sub1, signal.ranks = r.J+c(r.I1,r.I2), center = FALSE,
                           perm.test = FALSE, joint.rank = r.J)
cc.ldgs1 = cc.jive.res_sub1$CanCorRes$Loadings
can.cors = cc.jive.res_sub1$CanCorRes$Canonical_Correlations[1:r.J]
# Conduct CJIVE analysis on the second sub-sample. We will predict these joint scores using the
# results above
cc.jive.res_sub2 = cc.jive(blocks.sub2, signal.ranks = true_signal_ranks, center = FALSE,
                          perm.test = FALSE, joint.rank = r.J)
cc.jnt.scores.sub2 = cc.jive.res_sub2$CanCorRes$Jnt_Scores
cc.pred.jnt.scores.sub2 = cc.jive.pred(blocks.sub1, new.subjs = blocks.sub2,
                                      signal.ranks = true_signal_ranks,
                                      cc.jive.loadings = cc.ldgs1, can.cors = can.cors)
# Calculate the Pearson correlation coefficient between predicted and calculated joint scores
# for sub-sample 2
cc.pred.cor = diag(cor(cc.pred.jnt.scores.sub2, cc.jnt.scores.sub2))
print(cc.pred.cor)
# Plots of CJIVE estimates against true counterparts and include an estimate of their chordal
# norm
layout(matrix(1:2, ncol = 2))
plot(JntScores.2, cc.pred.jnt.scores.sub2, ylab = "Predicted Joint Scores",
    xlab = "True Joint Scores",
    col = rep(1:r.J, each = n/2),
    main = paste("Chordal Norm = ",
    round(chord.norm.diff(JntScores.2, cc.pred.jnt.scores.sub2),2)))
legend("topleft", legend = paste("Component", 1:r.J), col = 1:r.J, pch = 1)
plot(cc.jnt.scores.sub2, cc.pred.jnt.scores.sub2, ylab = "Predicted Joint Scores",
   xlab = "Estimated Joint Scores",
    col = rep(1:r.J, each = n/2),
    main = paste("Chordal Norm = ",
                  round(chord.norm.diff(cc.jnt.scores.sub2, cc.pred.jnt.scores.sub2),2)))
layout(1)

}
