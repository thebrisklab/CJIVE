###########          Use PC Scores to select joint rank based on cancor         ##################
#' Permutation Test for Joint Rank in CJIVE
#'
#' @description Conducts the permutation test for the number of joint components as described in CJIVE manuscript. Briefly, canonical correlations (CC) between principal component
#'              vectors of the data are obtained (PC). Then for 1:nperms, the rows of one data set are permuted and CCs between PC vectors are calculated, retaining
#'              the maximum CC. These maximum CCs form a null distribution against which the original CCs are tested. The number of original CCs exceeding the (1-alpha)^th
#'              percentile is the returned as the joint rank.
#'
#' @param dat.blocks a list of two matrices with samples along rows and features along columns, which contain data on the same n individuals/sampling units
#' @param signal.ranks a vector of length two which contains the rank for the signal within each data block. The rank corresponds to the number of principal
#'                     components (PCs) to be retained within each data block. If NULL, the ranks are determined by the parameter 'perc.var.' Default is NULL
#' @param nperms integer value indicating the number of permutations that should be performed
#' @param perc.var an alternative to signal.ranks that allows specification of ranks based on the desired proportion of total variation to be retained. F
#'                 For perc.var = p (where 0<p<1), rank is determined as the minimum number of eigenvalues whose cumulative sum is at least p*(total sum of eigenvalues).
#'                 Default is 0.95 (i.e. 95\% of total variation preserved for each data block)
#' @param alpha nominal type-I error rate
#' @param center logical (TRUE/FALSE) indicating whether data should be column-centered prior to testing. Default is TRUE
#' @return The Frobenius norm of the matrix X, calculated as the sum of square entries in X
#' @export
#'
perm.jntrank <- function(dat.blocks, signal.ranks = NULL, nperms = 500, perc.var = 0.95, alpha = 0.05, center = TRUE){
  n.r.1 = nrow(dat.blocks[[1]])
  n.r.2 = nrow(dat.blocks[[2]])
  if(n.r.1 != n.r.2){stop("The number of rows in each data matrix must match")}
  n = n.r.1

  K = length(dat.blocks)
  ##Column center data blocks
  if(center){
    cent.blocks = lapply(dat.blocks, function(D){scale(D, scale = FALSE)})
  }  else{
    cent.blocks = dat.blocks
  }

  if(is.null(signal.ranks)){
    all.singvals = lapply(cent.blocks, function(x) svd(x)$d)
    for(k in 1:K){
      d = all.singvals[[k]]
      signal.ranks = c(signal.ranks, which.min(cumsum(d^2) <= perc.var*sum(d^2)) )
    }
  }

  x.all.svd = list()
  for(k in 1:K){
    x.all.svd[[k]] = svd(cent.blocks[[k]], nu = signal.ranks[k], nv = signal.ranks[k])
  }

  U.all = lapply(x.all.svd, function(x) scale(x$u, scale = FALSE))

  U1tU2 = t(U.all[[1]])%*%U.all[[2]]
  orig.corrs = svd(U1tU2)$d

  perm.corrs = NULL

  for (i in 1:nperms){
    U1tU2.perm = t(U.all[[1]])%*%U.all[[2]][sample(n),]
    perm.svd = svd(U1tU2.perm)
    perm.corrs = rbind(perm.corrs, perm.svd$d)
  }

  test.val = stats::quantile(perm.corrs[,1],probs =  1-alpha)
  r.j = which.min(orig.corrs >= test.val) - 1
  ranks = c(r.j, signal.ranks)
  p.vals = NULL
  for(i in 1:r.j){
    p.vals = c(p.vals, sum(orig.corrs[i] <= perm.corrs[,1])/nperms)
    }
  names(ranks) = c("Joint", "Total Signal 1", "Total Signal 2")
  res = list(ranks, orig.corrs[1:r.j], p.vals)
  names(res) = c("Ranks", "Canonical Corrs", "P-values")
  return(res)
}

###########             CC.JIVE           ##################
###########   uses permutation test based on PC scores to find joint rank
###########   and estimates joint subject scores as scaled average of canonical variables
#' Canonical (Correlation) JIVE
#'
#' @description Performs Canonical JIVE as described in the CJVE manuscript. This method is equivalent to AJIVE for 2 data sets.
#'
#' @param dat.blocks a list of two matrices with samples along rows and features along columns, which contain data on the same n individuals/sampling units
#' @param signal.ranks a vector of length two which contains the rank for the signal within each data block. The rank corresponds to the number of principal
#'                     components (PCs) to be retained within each data block. If NULL, the ranks are determined by the parameter 'perc.var.' Default is NULL
#' @param joint.rank The rank of the joint subspace i.e., number of components in the joint subspace
#' @param perc.var an alternative to signal.ranks that allows specification of ranks based on the desired proportion of total variation to be retained. F
#'                 For perc.var = p (where 0<p<1), rank is determined as the minimum number of eigenvalues whose cumulative sum is at least p*(total sum of eigenvalues)
#'                 Default is 0.95 (i.e. 95\% of total variation preserved for each data block).
#' @param perm.test logical (TRUE/FALSE) of whether permutation test for joint rank should be performed. Overrides 'joint.rank' parameter if TRUE. Default is TRUE
#' @param center logical (TRUE/FALSE) indicating whether data should be column-centered prior to testing. Default is TRUE
#' @param nperms integer value indicating the number of permutations that should be performed. Default is 1000
#'
#' @return A list of two lists:
#'         1) 'CanCorRes' contains results from the canonical correlation of PC scores including, the joint rank, joint subject sores,
#'         canonical correlations (and their respective p-values if perm.test was used), canonical loadings for the joint subspace, and total signal ranks
#'         2) 'sJIVE', i.e. Simple JIVE results, correspond to the AJIVE when all ranks are known; includes the joint and individual signal matrices, concatenated PC scores,
#'         and the projection matrix used to project each data block onto the joint subspace
#' @examples
#' #Assign sample size and the number of features in each dataset
#' n = 200 #sample size
#' p1 = 2000 #Number of features in data set X1
#' p2 = 1000 #Number of features in data set X2
#'
#' # Assign values of joint and individual signal ranks
#' r.J = 1 #joint rank
#' r.I1 = 2 #individual rank for data set X1
#' r.I2 = 2 #individual rank for data set X2
#'
#'
#' # Simulate data sets
#' ToyDat = GenToyDatBinRank(n = 200, p1 = p1, p2 = p2, JntVarEx1 = 0.05, JntVarEx2 = 0.05,
#'                            IndVarEx1 = 0.25, IndVarEx2 = 0.25, jnt_rank = r.J, equal.eig = FALSE,
#'                            ind_rank1 = r.I1, ind_rank2 = r.I2, SVD.plots = TRUE, Error = TRUE,
#'                            print.cor = TRUE)
#' # Store simulated data sets in an object called 'blocks'
#' blocks <- ToyDat$'Data Blocks'
#'
#' # Save Subject scores as R objects
#' JntScores = ToyDat[['Scores']][['Joint']]
#' IndivScore.X = ToyDat[['Scores']][["Indiv_1"]]
#' IndivScore.Y = ToyDat[['Scores']][["Indiv_2"]]
#'
#' # Save joint variable loadings as R objects
#' JntLd.X = t(ToyDat$Loadings$Joint_1)
#' JntLd.Y = t(ToyDat$Loadings$Joint_2)
#'
#' # Save individual variable loadings as R objects
#' IndivLd.X =t(ToyDat$Loadings$Indiv_1)
#' IndivLd.Y = t(ToyDat$Loadings$Indiv_2)
#'
#' # Save joint, individual, and noise signal matrices as R objects
#' JX = ToyDat[[1]]$J1
#' JY = ToyDat[[1]]$J2
#' IX = ToyDat[[1]]$I1
#' IY = ToyDat[[1]]$I2
#' EX = ToyDat[[1]]$E1
#' EY = ToyDat[[1]]$E2
#'
#'
#' ## Check that proportions of variation explained are (approximately) equal to intended values
#' JVE.X = MatVar(JX)/MatVar(blocks[[1]])
#' JVE.Y = MatVar(JY)/MatVar(blocks[[2]])
#'
#' IVE.X = MatVar(IX)/MatVar(blocks[[1]])
#' IVE.Y = MatVar(IY)/MatVar(blocks[[2]])
#'
#' TotVE.X = MatVar((JX + IX))/MatVar(blocks[[1]])
#' TotVE.Y = MatVar((JY + IY))/MatVar(blocks[[2]])
#'
#'
#' CJIVE.res = cc.jive(blocks, c(r.I1,r.I2)+r.J, r.J, perm.test = FALSE)
#'# CJIVE signal matrix estimates
#'J.hat = CJIVE.res$sJIVE$joint_matrices
#'I.hat = CJIVE.res$sJIVE$indiv_matrices
#'
#'# CJIVE loading estimates
#'WJ = lapply(J.hat, function(x) x[['v']])
#'WI = lapply(I.hat, function(x) x[['v']])
#'
#'# Plots of CJIVE estimates against true counterparts and include an estimate of their chordal norm
#'layout(matrix(1:6,2, byrow = TRUE))
#'plot(JntScores, CJIVE.res$CanCorRes$Jnt_Scores, xlab = "True Joint Scores",
#'     ylab = "CJIVE Joint Scores",
#'     sub = paste0("Chordal Norm = ",
#'                  round(chord.norm.diff(JntScores, CJIVE.res$CanCorRes$Jnt_Scores), 3)))
#'plot(JntLd.X, WJ[[1]][,1], xlab = "True Joint Loadings X", ylab = "CJIVE Joint Loadings X",
#'     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.X, WJ[[1]][,1]), 3)))
#'plot(JntLd.Y, WJ[[2]][,1], xlab = "True Joint Loadings Y", ylab = "CJIVE Joint Loadings Y",
#'     sub = paste0("Chordal Norm = ", round(chord.norm.diff(JntLd.Y, WJ[[2]][,1]), 3)))
# plot(1,lwd=0,axes=F,xlab="",ylab="", "n")
#'plot.new(); legend("left", paste("Comp.", 1:2), pch = 1, col  = c("orange", "green"),bty = "n" )
#'plot(IndivLd.X, WI[[1]][,1:2], xlab = "True Individual Loadings X",
#'     ylab = "CJIVE Individual Loadings X",
#'     col = c(rep("orange",p1), rep("green",p2)),
#'     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.X, WI[[1]][,1:2]), 3)))
#'plot(IndivLd.Y, WI[[2]][,1:2], xlab = "True Individual Loadings Y",
#'     ylab = "CJIVE Individual Loadings Y",
#'     col = c(rep("orange",p1), rep("green",p2)),
#'     sub = paste0("Chordal Norm = ", round(chord.norm.diff(IndivLd.Y, WI[[2]][,1:2]), 3)))
#'layout(1)


#'
#' @export

cc.jive<-function(dat.blocks, signal.ranks = NULL, joint.rank = 1, perc.var = 0.95, perm.test = TRUE, center = FALSE, nperms = 1000){

  n.1 = nrow(dat.blocks[[1]])
  n.2 = nrow(dat.blocks[[2]])
  if(n.1 != n.2){stop("The number of rows in each data matrix must match")}
  n = n.1

  if(center){
    cent.blocks = lapply(dat.blocks, scale)
  } else {
    cent.blocks = dat.blocks
  }

  if(is.null(signal.ranks)){
    d1 = svd(cent.blocks[[1]])$d
    d2 = svd(cent.blocks[[2]])$d

    r.1 = which.min(cumsum(d1^2) <= perc.var*sum(d1^2))
    r.2 = which.min(cumsum(d2^2) <= perc.var*sum(d2^2))

    signal.ranks = c(r.1,r.2)
  }

  r.1 = signal.ranks[1]
  r.2 = signal.ranks[2]

  r.J = joint.rank

  if(perm.test){
    cca.perm.res = perm.jntrank(cent.blocks, signal.ranks = c(r.1,r.2), nperms = nperms)
    r.J = cca.perm.res$Ranks[1]
  }

  x1.svd = svd(cent.blocks[[1]], nu = signal.ranks[1], nv = signal.ranks[1])
  x2.svd = svd(cent.blocks[[2]], nu = signal.ranks[2], nv = signal.ranks[2])


  U.1 = x1.svd$u
  U.2 = x2.svd$u

  if(r.J > 0){
    U1tU2 = t(U.1)%*%U.2
    U1tU2.svd = svd(U1tU2, nu = r.J, nv = r.J)
    orig.corrs = U1tU2.svd$d

    if(r.J == 1){
      jnt.scores = (U.1%*%U1tU2.svd$u +U.2%*%U1tU2.svd$v)*(2*(1 + orig.corrs[r.J]))^-.5
    } else {
      jnt.scores = (U.1%*%U1tU2.svd$u +U.2%*%U1tU2.svd$v)%*%diag((2*(1 + orig.corrs[1:r.J]))^-.5)
    }
  } else {
    jnt.scores = rep(0, nrow(cent.blocks[[1]]))
    U1tU2 = t(U.1)%*%U.2
    U1tU2.svd = svd(U1tU2)
    orig.corrs = U1tU2.svd$d
  }
  sjive.out = sjive(cent.blocks, c(r.1, r.2), r.J, jnt.scores)

  if(perm.test){
    res = list(r.J, jnt.scores, list(orig.corrs,cca.perm.res$'P-values'), list(U1tU2.svd$u, U1tU2.svd$v), signal.ranks)
    names(res) = c("Jnt_Rank","Jnt_Scores", "Canonical_Correlations", "Loadings", "Signal Ranks")
  } else {
    res = list(r.J, jnt.scores, orig.corrs, list(U1tU2.svd$u, U1tU2.svd$v), signal.ranks)
    names(res) = c("Jnt_Rank","Jnt_Scores", "Canonical_Correlations", "Loadings", "Signal Ranks")
  }

  res.out = list(res, sjive.out)
  names(res.out) = c("CanCorRes","sJIVE")
  return(res.out)
}

#########   Compute predicted joint scores for new subjects using current CCJIVE Joint loadings #############
#' CJIVE joint subject score prediction
#'
#' @description Predicts joint scores for new subjects based on CJIVE joint scores
#'
#' @param orig.dat.blocks list of the two data matrices on which CJIVE was initially conducted
#' @param new.subjs list of two data matrices containing information on new subjects
#' @param signal.ranks a vector of length two which contains the rank for the signal within each data block. The rank corresponds to the number of principal
#'                     components (PCs) to be retained within each data block. If NULL, the ranks are determined by the parameter 'perc.var.' Default is NULL
#' @param cc.jive.loadings canonical loadings for the joint subspace
#' @param can.cors canonical correlations from the PCs of the data on which CJIVE was initially conducted - notated as rho_j in CJIVE manuscript
#'
#' @return matrix of joint subject score for new subjects
#' @export
cc.jive.pred<-function(orig.dat.blocks, new.subjs, signal.ranks, cc.jive.loadings, can.cors){

  r1 = signal.ranks[1]
  r2 = signal.ranks[2]
  X1.svd = svd(orig.dat.blocks[[1]], nu = r1, nv = r1)
  X2.svd = svd(orig.dat.blocks[[2]], nu = r2, nv = r2)

  U1.Ldngs = cc.jive.loadings[[1]]
  U2.Ldngs = cc.jive.loadings[[2]]

  Pred.CanVar.1 = new.subjs[[1]]%*%X1.svd$v%*%diag(X1.svd$d[1:r1]^-1)%*%U1.Ldngs
  Pred.CanVar.2 = new.subjs[[2]]%*%X2.svd$v%*%diag(X2.svd$d[1:r2]^-1)%*%U2.Ldngs
  pred.jnt.scores = sqrt(1/2*(1+can.cors))*(Pred.CanVar.1 + Pred.CanVar.2)

  return(pred.jnt.scores)
}

###########          sJIVE         ##################
#' Simple JIVE
#'
#' @description Conducts AJIVE estimation under the assumption that all ranks are known and no components are discarded
#'
#' @param blocks list of data blocks, i.e. matrices, all having the same number of rows, which correspond to the same sampling units (i.e. study participants, patients, etc.)
#' @param signal_ranks numerical vector of the same length as 'blocks' with each entry corresponding to the rank of the respective matrix in 'blocks'
#' @param joint.rank integer value corresponding to the rank of the joint signal subspace, i.e. number of components in the signal subspace
#' @param joint_scores numerical matrix containing joint subject scores if they were calculated by some other method, e.g. Canonical Correlation of PC scores.
#'                     Must have the same number of rows as each matrix in 'blocks' and number of columns equal to 'joint_rank'. If NULL, joint scores are calculated and
#'                     returned. Default is NULL.
#'
#' @return list of 4 items: 1) joint signal matrices, their SVDs, and the proportion of total variation in each matrix that is attributable to the joint signal
#'                          2) individual signal matrices, their SVDs, and the proportion of total variation in each matrix that is attributable to the individual signal
#'                          3) concatenated PC scores, used to determine joint subspace
#'                          4) projection matrix for joint subspace
#'                          5) joint subject scores (only returned if not provided initially)
#' @export
sjive<-function (blocks, signal_ranks, joint.rank, joint_scores = NULL)
{

  j=joint.rank
  K <- length(blocks)
  if (K < 2) {
    stop("ajive expects at least two data matrices.")
  }
  if (sum(sapply(blocks, function(X) any(is.na(X)))) > 0) {
    stop("Some of the blocks has missing data -- ajive expects full data matrices.")
  }

  block_svd <- list()
  scaled.centered <- list()
  total.var <- list()

  for (k in 1:K) {

    #Scale and center each block
    temp<-blocks[[k]]

    for (i in 1:dim(blocks[[k]])[2]) {
      temp[,i]<-blocks[[k]][,i]-mean(blocks[[k]][,i]) ##subtract the mean from each column
    }

    scaled.centered[[k]] <- temp
    block_svd[[k]] <- svd(scaled.centered[[k]],nu=signal_ranks[k]) ##Take SVD of each data block
    total.var[[k]] = sum(block_svd[[k]]$d^2)

    if (k==1) {
      stacked_mat<-block_svd[[k]]$u[,1:signal_ranks[k]] ##initialize matrix that results from stacking SVD matrices: contains first k vectors of U
    }
    else if (k>1){
      stacked_mat <- cbind(stacked_mat,block_svd[[k]]$u[,1:signal_ranks[k]]) ##stack the rest of the U matrices together with the first
    }
  }
  if(is.null(joint_scores)){
    if(j>0){
      joint_scores<-svd(stacked_mat,nu=j)$u ##take SVD of stacked matrix: nu=j indicates that the joint structure has rank=j
      joint_proj<-(joint_scores)%*%t(joint_scores) ##orthogonal projection matrix onto joint space based on stacked bases
    } else {
      joint_scores = rep(0, nrow(scaled.centered[[1]]))
      joint_proj<-(joint_scores)%*%t(joint_scores)
    }
  } else {
    joint_proj<-(joint_scores)%*%t(joint_scores)
  }
  indiv_structure <- list()
  joint_structure <- list()

  for (k in 1:K) {
    joint_structure[[k]] = list()
    joint_structure[[k]][['full']]<-joint_proj%*%as.matrix(scaled.centered[[k]])
    temp.svd = svd(joint_structure[[k]][['full']], nu = joint.rank, nv = joint.rank)

    if (j>0){
      joint_structure[[k]][['u']]<-temp.svd[['u']]
      joint_structure[[k]][['v']]<-temp.svd[['v']]
      joint_structure[[k]][['d']]<-temp.svd[['d']][1:joint.rank]
    } else {
      joint_structure[[k]][['u']]<-matrix(rep(0, nrow(scaled.centered[[1]])), ncol = 1)
      joint_structure[[k]][['v']]<-matrix(rep(0, ncol(scaled.centered[[k]])), ncol = 1)
      joint_structure[[k]][['d']]<-0
    }
    joint_structure[[k]][['full']] = joint_structure[[k]][['u']]%*%
      diag(joint_structure[[k]][['d']], nrow = max(joint.rank, 1), ncol = max(joint.rank, 1))%*%
      t(joint_structure[[k]][['v']])
    joint_structure[[k]][['VarEx']] = sum(joint_structure[[k]][['d']]^2)/total.var[[k]]

    dat_proj<-block_svd[[k]]$u%*%t(block_svd[[k]]$u)
    # indiv_structure[[k]]<-(dat_proj-joint_proj)%*%as.matrix(scaled.centered[[k]])
    temp = (diag(nrow(scaled.centered[[k]]))-joint_proj)%*%as.matrix(scaled.centered[[k]])
    indiv.rank = signal_ranks[k]-joint.rank
    temp.svd = svd(temp)

    indiv_structure[[k]] = list()

    indiv_structure[[k]][['u']] = temp.svd[['u']][,1:indiv.rank, drop = FALSE]
    indiv_structure[[k]][['v']] = temp.svd[['v']][,1:indiv.rank, drop = FALSE]
    indiv_structure[[k]][['d']] = temp.svd[['d']][1:indiv.rank]
    indiv_structure[[k]][['full']] = indiv_structure[[k]][['u']]%*%
      diag(indiv_structure[[k]][['d']], nrow = indiv.rank, ncol = indiv.rank)%*%
      t(indiv_structure[[k]][['v']])
    indiv_structure[[k]][['VarEx']] = sum(indiv_structure[[k]][['d']]^2)/total.var[[k]]
  }

  out = list(joint_structure,indiv_structure,stacked_mat,joint_proj)
  names(out) = c("joint_matrices", "indiv_matrices", "stacked_mat", "joint_projection")
  if(is.null(joint_scores)){out[[5]] = joint_scores; names(out)[5] = "joint_scores"}
  return(out)
}
