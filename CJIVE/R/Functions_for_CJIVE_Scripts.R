#' Generate 'Toy' Data
#'
#' @description Generates two Simulated Datasets that follow JIVE Model using binary subject scores
#'
#' @param n integer for sample size, i.e. number of subjects
#' @param p1 integer for number of features/variables in first data set
#' @param p2 integer for number of features/variables in second data set
#' @param JntVarEx1 numeric between (0,1) which describes proportion of variance in the first data set which is attributable to the joint signal
#' @param JntVarEx2 numeric between (0,1) which describes proportion of variance in the second data set which is attributable to the joint signal
#' @param IndVarEx1 numeric between (0,1) which describes proportion of variance in the first data set which is attributable to the individual signal
#' @param IndVarEx2 numeric between (0,1) which describes proportion of variance in the second data set which is attributable to the individual signal
#' @param jnt_rank integer for rank of the joint signal, i.e., number of joint components
#' @param equal.eig logical (TRUE/FALSE) for whether components should contribute equal variance to signal matrices - default is FALSE
#' @param ind_rank1 integer for rank of the individual signal in first data set, i.e., number of joint components
#' @param ind_rank2 integer for rank of the individual signal in second data set, i.e., number of joint components
#' @param SVD.plots logical (TRUE/FALSE)  for whether plots of singular values from signal should be produced - used to confirm number of components
#' @param Error logical (TRUE/FALSE)  final data sets should be noise contaminated - default is FALSE; use TRUE to obtain pure signal datasets
#' @param print.cor logical (TRUE/FALSE)  for whether to print matrix of correlations between subject scores)
#'
#' @return A 'list' object which contains 1) list of signal matrices which additively comprise the simulated data sets, i.e. joint, individual,
#'         and error matrices for each data set; 2) list of simulated data sets (each equal to the sum of the matrices in part 1);
#'         3) list of joint subject scores and individual subject scores for each data set, and 4) lsit of joint and individual loadings for each data set
#' @examples
#' ToyDat = GenToyDatBinRank(n = 200, p1 = 2000, p2 = 1000, JntVarEx1 = 0.05, JntVarEx2 = 0.05,
#'                            IndVarEx1 = 0.25, IndVarEx2 = 0.25, jnt_rank = 1, equal.eig = FALSE,
#'                            ind_rank1 = 2, ind_rank2 = 3, SVD.plots = TRUE, Error = TRUE,
#'                            print.cor = TRUE)
#' @export
#'
#'

GenToyDatBinRank <- function(n, p1, p2, JntVarEx1, JntVarEx2, IndVarEx1, IndVarEx2, jnt_rank = 1, equal.eig = FALSE, ind_rank1 = 2, ind_rank2 = 2,
                             SVD.plots = TRUE, Error = TRUE, print.cor = TRUE){

  #Write out both joint and indiv subject scores for both data sets first
  r.J = jnt_rank
  JntScores = matrix(stats::rbinom(n*r.J, size=1, prob=0.2), nrow = n, ncol = r.J)
  colnames(JntScores) = paste("Jnt Score", 1:r.J)

  r.I1 = ind_rank1
  r.I2 = ind_rank2

  b = stats::rbinom(n*(r.I1 + r.I2), size=1, prob=0.4)
  b = 1 - 2*b
  IndivScores = matrix(b, nrow = n, ncol = (r.I1 + r.I2))
  colnames(IndivScores) = c(paste("Ind X Score", 1:r.I1), paste("Ind Y Score", 1:r.I2))

  if(print.cor){print("The correlation between subject scores is given by")
    print(round(stats::cor(cbind(JntScores, IndivScores)),4))}

  ##############################Define X Dataset##############################
  ##Then write each of the 3 variable loading vectors for the first dataset
  ##Note that the joint signal has rank=1
  AdjJntLoad.X = matrix(stats::rnorm(r.J*p1), nrow = r.J, ncol = p1)

  #change relavent scaling of joint components
  D.J = ifelse(diag(rep(equal.eig, r.J)), diag(rep(1,r.J)), diag(r.J:1))
  if(equal.eig){D.J = diag(rep(1,r.J))}
  else{
    D.J = diag(r.J:1)
  }
  JX = JntScores%*%sqrt(D.J)%*%AdjJntLoad.X

  if(SVD.plots){
    plot(svd(JX)$d, ylab = "Singular Values")
    graphics::title("SVD of Joint Signal from X")
  }

  ##Note that the individual signal has rank = 2 as well
  IndScores.X = IndivScores[,1:r.I1]
  IndLoad.X = matrix(stats::rnorm(n = p1*r.I1), nrow = r.I1, ncol = p1)
  if(equal.eig){D.IX = diag(rep(1,r.I1))}
  else{
    D.IX = diag(r.I1:1)
  }

  IX = IndScores.X%*%sqrt(D.IX)%*%IndLoad.X

  if(SVD.plots){
    plot(svd(IX)$d, ylab = "Singular Values")
    graphics::title("SVD of Individual Signal from X")
  }

  AX = JX + IX

  ##############################Define Y Dataset##############################
  AdjJntLoad.Y = matrix(stats::rnorm(r.J*p2), nrow = r.J, ncol = p2)

  ##Note that the joint signal has rank = 3
  JY = JntScores%*%sqrt(D.J)%*%AdjJntLoad.Y

  if(SVD.plots){
    plot(svd(JY)$d, ylab = "Singular Values")
    graphics::title("SVD of Joint Signal from Y")
  }

  IndScores.Y = IndivScores[,(r.I1 + 1:r.I2)]
  IndLoad.Y = matrix(stats::rnorm(r.I2*p2), nrow = r.I2, ncol = p2)
  if(equal.eig){D.IY = diag(rep(1,r.I2))}
  else{
    D.IY = diag(r.I2:1)
  }

  ##Note that the individual signal has rank=2
  IY = IndScores.Y%*%sqrt(D.IY)%*%IndLoad.Y

  if(SVD.plots){
    plot(svd(IY)$d, ylab = "Singular Values")
    graphics::title("SVD of Individual Signal from Y")
  }

  ##Error matrix
  EX = matrix(stats::rnorm(n*p1), nrow=n, ncol=p1)*Error

  ##Error matrix
  EY = matrix(stats::rnorm(n*p2), nrow=n, ncol=p2)*Error

  Dat.X = AdjSigVarExp(JX, IX, EX, JntVarEx1, IndVarEx1)
  JX = Dat.X$J
  IX = Dat.X$I

  Dat.Y = AdjSigVarExp(JY, IY, EY, JntVarEx2, IndVarEx2)
  JY = Dat.Y$J
  IY = Dat.Y$I

  Blocks = list(Dat.X[["Data"]], Dat.Y[["Data"]])

  Dat.Comps = list(JX, JY, IX, IY, EX, EY)
  names(Dat.Comps) = c("J1", "J2", "I1", "I2", "E1", "E2" )

  Scores = list(JntScores, IndScores.X, IndScores.Y)
  names(Scores) = c("Joint", "Indiv_1", "Indiv_2")

  Loadings = list(AdjJntLoad.X, IndLoad.X, AdjJntLoad.Y, IndLoad.Y)
  names(Loadings) = c("Joint_1", "Indiv_1", "Joint_2", "Indiv_2")

  out = list(Dat.Comps, Blocks, Scores, Loadings)
  names(out) = c("Data Components", "Data Blocks", "Scores", "Loadings")

  out
}


###########   Adjust Dataset Components to get Desired R^2 Values  #################
#' Adjust Signal Variation Explained
#'
#' @description Adjusts the proportion of total variation attributable to each signal component to predetermined values
#'
#' @param J joint signal matrix of size n-by-p
#' @param I individual signal matrix of size n-by-p
#' @param N noise/error matrix of size n-by-p
#' @param JntVarEx desired proportion of total variation explained by the joint signal
#' @param IndVarEx desired proportion of total variation explained by the individual signal
#'
#' @return a list of 3 items: 1) adjusted joint signal matrix; 2) adjusted individual signal matrix; 3) data matrix additively comprised of the adjusted signal matrices
#'
#' @export
AdjSigVarExp <-function(J, I, N, JntVarEx, IndVarEx){
  simul.quads = function(x, parms){
    JJ = parms[1]
    II = parms[2]
    NN = parms[3]
    JN = parms[4]
    IN = parms[5]
    R_J = parms[6]
    R_I = parms[7]

    y1 = x[1]^2*II*(1 - R_I) - 2*x[1]*IN*R_I - R_I*(x[2]^2*JJ + 2*x[2]*JN + NN)
    y2 = x[2]^2*JJ*(1 - R_J) - 2*x[2]*JN*R_J - R_J*(x[1]^2*II + 2*x[1]*IN + NN)

    y = c(y1,y2)
    return(y)
  }

  JJ = MatVar2(J)
  II = MatVar2(I)
  NN = MatVar2(N)
  JN = sum(diag(J%*%t(N)))
  IN = sum(diag(I%*%t(N)))
  R_J = JntVarEx
  R_I = IndVarEx

  parms = c(JJ, II, NN, JN, IN, R_J, R_I)

  A = J + I
  AA = MatVar2(A)
  NN = MatVar2(N)
  AN = sum(diag(A%*%t(N)))

  ##Desired Total Variance Explained
  d0 = IndVarEx + JntVarEx
  a = AA*(1 - d0)
  b = -2*AN
  c = -d0*NN

  d.A = (-b+sqrt(b^2 - 4*a*c))/(2*a)

  start = c(0.5, 0.5)*d.A

  roots = rootSolve::multiroot(simul.quads, start, parms = parms)
  c = roots$root[1]
  d = roots$root[2]

  J.d = d*J; I.c = c*I;
  Dat = J.d + I.c + N

  res = list(J.d, I.c, Dat)
  names(res) = c("J", "I", "Data")
  return(res)
}


###########           Frobenius Norm of Data Matrix Values         #################
#' Matrix variation (i.e. Frobenius norm)
#'
#' @description Calculates the Frobenius norm of a matrix, which can be used as a measure of total variation
#'
#' @param X a matrix of any size
#'
#' @return The Frobenius norm of the matrix X, calculated as the sum of square entries in X
#'
#' @examples
#' X = matrix(rnorm(10), 5,2)
#' MatVar2(X)
#' @export

MatVar2 = function(X){
sum(X^2)
}

###########          Use PC Scores to calculate joint rank based on cancor         ##################
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
    all.singvals = sapply(cent.blocks, function(x) svd(x$d))
    for(k in 1:K){
      d = all.singvals
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

###########       Chordal norm for matrices with diff ranks         ##################
#' Chordal norm between column-subspaces of two matrices
#'
#' @description Calculates the chordal norm between the column subspaces of two matrices. Matrices must have the same number of rows. Let U_x and U_y represent the singular
#'              vectors of matrices X and Y, respectively. The chordal norm can be calculated as the square root of the sum of the singular values of t(U_x)%*%U_y
#'
#' @param X a matrix with the same number of rows as Y and any number of columns
#' @param Y a matrix with the same number of rows as X and any number of columns
#' @param tol threshold under which singular values of inner product are zeroed out
#'
#' @return (Numeric) Chordal norm between column-subspaces of X and Y, scaled to the interval [0,1]
#'
#' @export

chord.norm.diff = function(X,Y, tol = 1E-8){
  svd.X = svd(X)
  svd.Y = svd(Y)
  if(svd.Y$d[1]>0){
    Uy = svd.Y$u
  } else {
      Uy = matrix(0, dim(svd.Y$u))
  }
  if(svd.X$d[1]>0){
    Ux = svd.X$u
  } else {
    Ux = matrix(0, dim(svd.X$u))
  }
  k = max(min(sum(svd.X$d > tol), sum(svd.Y$d > tol)),1)
  inner.prod = t(Ux)%*%Uy
  sig = round(svd(inner.prod)$d[1:k], 8)
  sqrt(sum((1 - sig^2))/k)
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
#' @param joint.rank The rank of the joint subspace i.e., number of components in the joint subspace
#' @param signal.ranks a vector of length two which contains the rank for the signal within each data block. The rank corresponds to the number of principal
#'                     components (PCs) to be retained within each data block. If NULL, the ranks are determined by the parameter 'perc.var.' Default is NULL
#' @param cc.jive.loadings canonical loadings for the joint subspace
#'
#' @return matrix of joint subject score for new subjects
#' @export
cc.jive.pred<-function(orig.dat.blocks, new.subjs, joint.rank = 1, signal.ranks, cc.jive.loadings){

  r1 = signal.ranks[1]
  r2 = signal.ranks[2]
  X1.svd = svd(orig.dat.blocks[[1]], nu = r1, nv = r1)
  X2.svd = svd(orig.dat.blocks[[2]], nu = r2, nv = r2)

  U1.Ldngs = cc.jive.loadings[[1]]
  U2.Ldngs = cc.jive.loadings[[2]]

  Pred.CanVar.1 = new.subjs[[1]]%*%X1.svd$v%*%diag(X1.svd$d[1:r1]^-1)%*%U1.Ldngs
  Pred.CanVar.2 = new.subjs[[2]]%*%X2.svd$v%*%diag(X2.svd$d[1:r2]^-1)%*%U2.Ldngs
  pred.jnt.scores = sqrt(1/2)*(Pred.CanVar.1 + Pred.CanVar.2)

  return(pred.jnt.scores)
}

##############            Retrieve simulation results stored in a directory                ##################
#' Retrieve simulation results
#'
#' @description  Retrives and compiles results from simulation study which are stored in a directory. A directory should contain separate .csv files (one per replicate),
#'  each of which will include all evaluation metrics and most experimental settings for that particular replicate. For the CJIVE manuscript, a directory houses results of all
#'  100 replicates for each combination of experimental factors.
#'
#' @param sim.dir (character string) file path for the directory from which results will be retrieved
#' @param p1 number of features in data set 1
#' @param p2 number of features in data set 2
#' @param p1 (logical) do the replicate results contain correlations between predicted and true joint subject scores. Default is FALSE
#'
#' @return upper triangular p-by-p matrix
#' @export
GetSimResults_Dir = function(sim.dir, p1, p2, Preds=FALSE){
  files = list.files(sim.dir, pattern = ".csv")
  num_sims = length(files)

  JVEs = as.numeric(paste("0",strsplit(sim.dir, "[.]0")[[1]][-1], sep = "."))
  JntVarEx1 = JVEs[1]
  JntVarEx2 = JVEs[2]

  results = NULL
  for (nm in files){
    temp = utils::read.csv(file = file.path(sim.dir, nm), header = TRUE)
    varnames = temp[,1]; values = temp[,2]
    to.keep = !grepl("Pred",varnames)
    if(Preds==F & length(to.keep)>0){
      values = values[to.keep]
      varnames = varnames[to.keep]
      }
    results = rbind(results,c(values, p1, p2))
  }

  colnames(results) = c(varnames, "p1", "p2")

  results = data.frame(cbind(results, JntVarEx1, JntVarEx2))
}


############    Convert Vector to Network       ###############
#' Convert vector to network
#'
#' @description Converts a vector of size p choose 2 into a p-by-p upper triangular matrix
#'
#' @param invector numeric vector of size p choose 2
#'
#' @return upper triangular p-by-p matrix
#' @export

vec2net.u = function(invector) {
  #invector: p x 1, where p is the number of edges
  nNode = (1 + sqrt(1+8*length(invector)))/2
  outNet = matrix(0,nNode,nNode)
  outNet[upper.tri(outNet,diag = FALSE)] = invector
  dim(outNet) = c(nNode,nNode)
  outNet = outNet + t(outNet)
  #  diag(outNet) = 1
  outNet
}

#' Convert vector to network
#'
#' @description Converts a vector of size p choose 2 into a p-by-p lower triangular matrix
#'
#' @param invector numeric vector of size p choose 2
#'
#' @return lower triangular p-by-p matrix
#' @export
vec2net.l = function(invector) {
  #invector: p x 1, where p is the number of edges
  nNode = (1 + sqrt(1+8*length(invector)))/2
  outNet = matrix(0,nNode,nNode)
  outNet[lower.tri(outNet,diag = FALSE)] = invector
  outNet = t(outNet)
  dim(outNet) = c(nNode,nNode)
  outNet = outNet + t(outNet)
  #  diag(outNet) = 1
  outNet
}

############    Visually display heatmap of matrix: adapted from Erick Lock's show.image function employed in r.jive package   ###################
#' Display a heatmap of a matrix (adapted from Erick Lock's show.image funciton in the r.jive package)
#'
#' @description Visual display of a matrix as a heatmap with colors determined by entry values, and including a colorbar to aid interpretation of the heatmap
#'
#' @param Image matrix to display
#' @param ylab lab for y-axis of heatmap
#' @param xlab lab for x-axis of heatmap
#' @param net logical (TRUE/FALUSE) of whether entries correspond to edges between regions of interest in the Power-264 brain atlas. Default is FALSE
#' @param main main title for heatmap
#' @param sub subtitle for heatmap
#' @param colorbar logical (TRUE/FALUSE) of whether colorabar shouldl be included to aid interpretation. Default is TRUE
#'
#' @return graphical display of matrix as a heatmap
#' @export
show.image.2<-function (Image, ylab = "", xlab="",net=F,main="", sub="", colorbar = TRUE)
{
  #if(net){Image=Image-diag(Image)}

  #Image = Image*upper.tri(Image)
  lower = mean(Image) - 5 * stats::sd(Image)
  upper = mean(Image) + 5 * stats::sd(Image)
  Image[Image < lower] = lower
  Image[Image > upper] = upper
  if(colorbar){
    fields::image.plot(x = 1:dim(Image)[2], y = 1:dim(Image)[1], z = t(Image),
               zlim = c(lower, upper), axes = FALSE, col = gplots::bluered(100),
               xlab = xlab, ylab = ylab)
    graphics::title(main=main, sub = sub)
  }

  else if(!colorbar){
    graphics::image(x = 1:dim(Image)[2], y = 1:dim(Image)[1], z = t(Image),
               zlim = c(lower, upper), axes = FALSE, col = gplots::bluered(100),
               xlab = xlab, ylab = ylab)
    graphics::title(main=main, sub = sub)
  }

  if(net){
    mod = utils::read.table("C:/Users/Raphiel\'s PC/Dropbox/JIVE-PNC/PNC-Data/roi_id.txt")
    mod = mod[,1]

    Count.m = table(sort(mod))
    k = 0
    Grid.m = vector()
    for(i in 1:(length(Count.m)-1))
    {
      k = k + Count.m[i]
      Grid.m[i] = k
    }
    graphics::abline(v=c(0,Grid.m[-1]-32,232)+.5,
           h=c(0,232-Grid.m[-1]+32,232)+.5)
  }
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

    joint_structure[[k]][['u']]<-temp.svd[['u']]
    joint_structure[[k]][['v']]<-temp.svd[['v']]
    joint_structure[[k]][['d']]<-temp.svd[['d']][1:joint.rank]
    joint_structure[[k]][['full']] = joint_structure[[k]][['u']]%*%
      diag(joint_structure[[k]][['d']], nrow = joint.rank, ncol = joint.rank)%*%
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


#####################         Convert Simulation Results to a form that will allow ggplot       #######################
#' Convert simulation study results
#'
#' @description Convert results from simulation study into a form for graphing with ggplot
#'
#' @param AllSims matrix with each row representing results from a replicate in the simulation study described in CJIVE manuscript
#'
#' @return list of 2 items: 1) joint ranks determined by each method employed in the simulations study
#'                          2) chordal norms between true and estimated joint/individual loadings/scores for each method employed in the simulation study
#' @export
ConvSims_gg<-function(AllSims){
  sim.names = colnames(AllSims)
  #CompEvals.names = c(sim.names[grep("Scores", sim.names)], sim.names[grep("Loads", sim.names)])
  num.sims = dim(AllSims)[1]

  p1 = unique(AllSims$p1)
  p2 = unique(AllSims$p2)

  JVE_1 =  c(as.matrix(AllSims[,grep("JntVarEx1",sim.names)]))
  JVE1.labs = c(bquote("R"[J1]^2*"=0.05, p"[1]*"="*.(p1)), bquote("R"[J1]^2*"=0.5, p"[1]*"="*.(p1)))
  JVE_1 = factor(JVE_1, labels = JVE1.labs, levels = c(0.05, 0.5))

  JVE_2 = c(as.matrix(AllSims[,grep("JntVarEx2",sim.names)]))
  JVE2.labs = c(bquote("R"[J2]^2*"=0.05, p"[2]*"="*.(p2)), bquote("R"[J2]^2*"=0.5, p"[2]*"="*.(p2)) )
  JVE_2 = factor(JVE_2, labels = JVE2.labs, levels = c(0.05, 0.5))

  IVE_1 = c(as.matrix(AllSims[,grep("Indiv.Var.Exp.X", sim.names)]))
  IVE_1 = round(as.numeric(IVE_1), 2)

  IVE_2 = c(as.matrix(AllSims[,grep("Indiv.Var.Exp.Y", sim.names)]))
  IVE_2 = round(as.numeric(IVE_2), 2)

  if("True.Joint.Rank" %in% sim.names){
    sim.names.2 = sim.names[-grep("True.Joint.Rank", sim.names)]
    AllSims.2 = AllSims[,which(sim.names %in% sim.names.2)]
    Method = sim.names.2[grep("Joint.Rank", sim.names.2)]
    temp = colnames(AllSims.2)
    Rank = factor(c(as.matrix(AllSims.2[,temp[grep("Joint.Rank", temp)]])))
  } else {
    Method = sim.names[grep("Joint.Rank", sim.names)]
    Rank = factor(c(as.matrix(AllSims[,grep("Joint.Rank", sim.names)])))
  }

  Method = sub(".Joint.Rank", "", Method, fixed = TRUE);  Method = sub("aJI", "AJI", Method, fixed = TRUE)
  Method = sub("Correct", "", Method, fixed = TRUE); Method = sub(".Wrong", "-Over", Method, fixed = TRUE)
  Method = sub(".oracle", "-Oracle", Method, fixed = TRUE); Method = sub(".over", "-Over", Method, fixed = TRUE)
  Method = sub("Elbow", "Oracle", Method, fixed = TRUE); Method = sub(".95.", "-Over", Method, fixed = TRUE)
  Method = sub("iJI", "R.JI", Method, fixed = TRUE); Method = sub("AJIVE.", "AJIVE-", Method, fixed = TRUE)
  Method = sub("CC", "CJIVE", Method, fixed = TRUE); Method = sub("cJIVE", "CJIVE", Method, fixed = TRUE)
  Method[which(Method == "CJIVE.Oracle")] = "CJIVE-Oracle"
  Method[which(Method %in% c("CJIVE.Over", "cJIVE.Over"))] = "CJIVE-Over"
  Method = as.factor(Method)

  sim.ranks.gg = data.frame(Rank = Rank, Method = rep(Method, each = num.sims), JVE_1 = rep(JVE_1, each = length(levels(Method))),
                            JVE_2 = rep(JVE_2, each = length(levels(Method))),
                            IVE_1 = rep(IVE_1, each = length(levels(Method))),
                            IVE_2 = rep(IVE_2, each = length(levels(Method))))
  if("True.Joint.Rank" %in% sim.names){
    True_Rank = c(as.numeric(as.matrix(AllSims[,grep("True.Joint.Rank", sim.names)])))
    True_Rank = factor(True_Rank)
    sim.ranks.gg$True_Rank = rep(True_Rank, each = length(levels(Method)))
  }

  Norm = c(as.numeric(as.matrix(AllSims[,grep("Scores", sim.names)])))

  Method = rep(substr(sim.names[grep("Scores",sim.names)], 1,12), each = num.sims)
  Method = sub(".Joint.", "", Method, fixed = TRUE); Method = sub(".Indiv.", "", Method, fixed = TRUE);
  Method = sub("aJI", "AJI", Method, fixed = TRUE); Method = sub("iJI", "R.JI", Method, fixed = TRUE);
  Method = sub("r.J", "r", Method, fixed = TRUE); Method = sub("r.I", "r", Method, fixed = TRUE)
  Method[which(Method == "AJIVE")] = "AJIVE-Oracle"; Method = sub(".w", "-Over", Method, fixed = TRUE);
  Method[which(Method == "CC.Oracle")] = "CJIVE-Oracle"; Method[which(Method == "CC-Over")] = "CJIVE-Over"
  Method = sub("JIVE.", "JIVE-", Method, fixed = TRUE); Method = sub("[[:punct:]]$", "", Method, fixed = TRUE)
  Method = as.factor(Method)

  Type = rep(substring(sim.names[grep("Scores", sim.names)], 7), each = num.sims)
  Type = sub("w.", "", Type)
  Type = gsub("[.]", " ", Type)
  Type = gsub("X", "X[1]", Type)
  Type = gsub("Y", "X[2]", Type)
  Type = gsub("cle ", "", Type)
  Type = gsub("Over ", "", Type)
  Type = gsub("Ora", "", Type)
  Type = factor(Type)

  n.levs = nlevels(Type)*nlevels(Method)
  sim.score.norms.gg = data.frame(Norm = Norm, Method = Method, Type = Type, JVE_1 = rep(JVE_1, each = n.levs),
                                  JVE_2 = rep(JVE_2, each = n.levs),
                                  IVE_1 = rep(IVE_1, each = n.levs), IVE_2 = rep(IVE_2, each = n.levs))

  if("True.Joint.Rank" %in% sim.names){
    True_Rank = c(as.numeric(as.matrix(AllSims[,grep("True.Joint.Rank", sim.names)])))
    True_Rank = factor(True_Rank)
    sim.score.norms.gg$True_Rank = rep(True_Rank, each = length(levels(Method)))
  }

  Norm = c(as.numeric(as.matrix(AllSims[,grep("Loads", sim.names)])))

  Method = rep(substr(sim.names[grep("Loads",sim.names)], 1,12), each = num.sims)
  Method = sub(".Joint.", "", Method, fixed = TRUE); Method = sub(".Indiv.", "", Method, fixed = TRUE);
  Method = sub("aJI", "AJI", Method, fixed = TRUE); Method = sub("iJI", "R.JI", Method, fixed = TRUE);
  Method = sub("r.J", "r", Method, fixed = TRUE); Method = sub("r.I", "r", Method, fixed = TRUE)
  Method[which(Method == "AJIVE")] = "AJIVE-Oracle"; Method = sub(".w", "-Over", Method, fixed = TRUE);
  Method[which(Method == "CC.Oracle")] = "CJIVE-Oracle"; Method[which(Method == "CC-Over")] = "CJIVE-Over"
  Method = sub("JIVE.", "JIVE-", Method, fixed = TRUE); Method = sub("[[:punct:]]$", "", Method, fixed = TRUE)
  Method = as.factor(Method)

  Type = factor(rep(substring(sim.names[grep("Loads", sim.names)], 7), each = num.sims))
  Type = sub("w.", "", Type)
  Type = gsub("[.]", " ", Type)
  Type = gsub("X", "X[1]", Type)
  Type = gsub("Y", "X[2]", Type)
  Type = gsub("cle ", "", Type)
  Type = gsub("Over ", "", Type)
  Type = gsub("Ora", "", Type)
  Type = factor(Type)

  n.levs = nlevels(Type)*nlevels(Method)
  sim.load.norms.gg = data.frame(Norm = Norm, Method = Method, Type = Type, JVE_1 = rep(JVE_1, each = n.levs),
                                 JVE_2 = rep(JVE_2, each = n.levs), IVE_1 = rep(IVE_1, each = n.levs), IVE_2 = rep(IVE_2, each = n.levs))

  if("True.Joint.Rank" %in% sim.names){
    True_Rank = c(as.numeric(as.matrix(AllSims[,grep("True.Joint.Rank", sim.names)])))
    True_Rank = factor(True_Rank)
    sim.load.norms.gg$True_Rank = rep(True_Rank, each = length(levels(Method)))
  }

  sim.all.norms.gg = rbind(sim.score.norms.gg, sim.load.norms.gg)

  out = list(sim.ranks.gg, sim.all.norms.gg)
  names(out) = c("Ranks", "Subj and Ldg Norms")
  out
}


#####################     Convert square matrix/network for plotting      #######################
#' Function for plotting networks with ggplot
#'
#' @description Convert matrix representation of a network for graphical display via ggplot
#'
#' @param gmatrix square matrix of size p-by-p in which entries represent the strength of (un-directed) edges between the p nodes
#' @param sort_indices vector of length p by which nodes are sorted. If NULL, then nodes are not sorted. Default is NULL.
#'
#' @return a data frame of three variables: X1, which represents the row from which the edge comes;  X2, which represents the column from which the edge comes;
#'        3) value, matrix entry representing the strength of the edge between the nodes represented by X1 and X2
#' @export

create.graph.long = function(gmatrix,sort_indices=NULL) {
  nnode = nrow(gmatrix)
  X1 = c(1:nnode)%x%rep(1,nnode)
  X2 =  rep(1,nnode)%x%c(1:nnode)
  X1 = factor(X1, )
  if (!is.null(sort_indices)) {
    gmatrix = gmatrix[sort_indices,sort_indices]
  }
  value = as.vector(as.matrix(gmatrix))
  data.frame(X1,X2,value)
}

#####################         Plot CJIVE Norms from Simulation Study      #######################
#' Function for plotting chordal norms between estimated and true subspaces within the simulation study described in CJIVE manuscript
#'
#' @description Graphically displays the center and spread of chordal norms for joint/individual score/loading subspaces
#'
#' @param norm.dat data frame with at least the 5 following variables:
#'                   Norm - the value of the norm for a particular subspace;
#'                   Type - the subspace for which the norm is given (i.e., joint/individual score/loading for dataset X1 or X2 (except for joint scores))
#'                   Method - the method by which the subspace was estimated, e.g. CJIVE, AJIVE, R.JIVE
#'                   JVE_1 and JVE_2 - labels describing the proportion of joint variation explained in each dataset (and typically the number of variables in dataset X2)
#' @param cols a vector of colors, must have length equal to the number of methods used in the simulation
#' @param show.legend logical (TRUE/FALSE) for whether a legend should be included in the plot. Default is FALSE
#' @param text.size numeric value for the font size
#' @param lty linetype (see ggplot2). Default = 1
#' @param y.max maximum value for the horizontal axis of the plot
#'
#' @return graphical display (via ggplot2)
#' @export
gg.norm.plot<-function(norm.dat, cols, show.legend = FALSE, text.size, lty = 1, y.max = 1){
  labs = levels(norm.dat$Type)[c(3,6,7,1,2,4,5)]
  labs.ex = c("Joint Subj Scores", expression("Joint Loadings"*"X"[1]), expression("Joint Loadings"*"X"[2]),
              expression("Indiv Subj Scores"*"X"[1]), expression("Indiv Subj Scores"*"X"[2]),
              expression("Indiv Loadings"*"X"[1]), expression("Indiv Loadings"*"X"[2]))
  ggplot2::ggplot(data = norm.dat, ggplot2::aes(x = norm.dat$Type, y = norm.dat$Norm)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = norm.dat$Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.5) +
    ggplot2::labs(y = "Chordal Norm", x = "Type") +
    ggplot2::facet_grid(JVE_2 ~ JVE_1, labeller = ggplot2::label_parsed()) +
    ggplot2::scale_x_discrete(limits = levels(norm.dat$Type)[c(3,6,7,1,2,4,5)], labels = labs.ex) +
    ggplot2::scale_fill_manual(values=cols) +
    ggplot2::scale_colour_manual(values=cols) +
    ggplot2::theme_bw() + ggplot2::coord_cartesian(ylim = c(0, y.max)) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(face = "bold", hjust = 01, angle = 70, size = text.size-3),
          text = ggplot2::element_text(size = text.size))
}

#####################         Plot CJIVE Norms from Simulation Study      #######################
#' Function for plotting chordal norms between estimated and true subject score subspaces within the simulation study described in CJIVE manuscript
#'
#' @description Graphically displays the center and spread of chordal norms for joint/individual subject score subspaces
#'
#' @param norm.dat data frame with at least the 5 following variables:
#'                   Norm - the value of the norm for a particular subspace;
#'                   Type - the subspace for which the norm is given (i.e., joint and individual subject scores for dataset X1 or X2 (except joint scores, which are for both datasets))
#'                   Method - the method by which the subspace was estimated, e.g. CJIVE, AJIVE, R.JIVE
#'                   JVE_1 and JVE_2 - labels describing the proportion of joint variation explained in each dataset (and typically the number of variables in dataset X2)
#' @param cols a vector of colors, must have length equal to the number of methods used in the simulation
#' @param show.legend logical (TRUE/FALSE) for whether a legend should be included in the plot. Default is FALSE
#' @param text.size numeric value for the font size
#' @param lty linetype (see ggplot2). Default = 1
#' @param y.max maximum value for the horizontal axis of the plot
#'
#' @return graphical display (via ggplot2)
#' @export
gg.score.norm.plot<-function(norm.dat, cols, show.legend = FALSE, text.size, lty = 1, y.max = 1){
  labs = levels(norm.dat$Type)[grep("Score",levels(norm.dat$Type))][c(3,1,2)]
  labs.ex = c("Joint Subj Scores", expression("Indiv Subj Scores"*"X"[1]), expression("Indiv Subj Scores"*"X"[2]))
  ggplot2::ggplot(data = norm.dat, ggplot2::aes(x = norm.dat$Type, y = norm.dat$Norm)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = norm.dat$Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.5) +
    ggplot2::labs(y = "Chordal Norm", x = "Type") +
    ggplot2::facet_grid(norm.dat$JVE_2 ~ norm.dat$JVE_1, labeller = ggplot2::label_parsed()) +
    ggplot2::scale_x_discrete(limits = labs, labels = labs.ex) +
    ggplot2::scale_fill_manual(values=cols) +
    ggplot2::scale_colour_manual(values=cols) +
    ggplot2::theme_bw() + ggplot2::coord_cartesian(ylim = c(0, y.max)) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(face = "bold", hjust = 01, angle = 70, size = text.size-3),
          text = ggplot2::element_text(size = text.size))
}

#####################         Plot CJIVE Norms from Simulation Study      #######################
#' Function for plotting chordal norms between estimated and true variable loading subspaces within the simulation study described in CJIVE manuscript
#'
#' @description Graphically displays the center and spread of chordal norms for joint/individual variable loading subspaces
#'
#' @param norm.dat data frame with at least the 5 following variables:
#'                   Norm - the value of the norm for a particular subspace;
#'                   Type - the subspace for which the norm is given (i.e., joint/individual variable loadings for dataset X1 or X2)
#'                   Method - the method by which the subspace was estimated, e.g. CJIVE, AJIVE, R.JIVE
#'                   JVE_1 and JVE_2 - labels describing the proportion of joint variation explained in each dataset (and typically the number of variables in dataset X2)
#' @param cols a vector of colors, must have length equal to the number of methods used in the simulation
#' @param show.legend logical (TRUE/FALSE) for whether a legend should be included in the plot. Default is FALSE
#' @param text.size numeric value for the font size
#' @param lty linetype (see ggplot2). Default = 1
#' @param y.max maximum value for the horizontal axis of the plot
#'
#' @return graphical display (via ggplot2)
#' @export
gg.load.norm.plot<-function(norm.dat, cols, show.legend = FALSE, text.size, lty = 1, y.max = 1){
  labs = levels(norm.dat$Type)[grep("Load",levels(norm.dat$Type))][c(3,4,1,2)]
  labs.ex = c(expression("Joint Variable Loadings"*"X"[1]),expression("Joint Variable Loadings"*"X"[2]),
              expression("Indiv Variable Loadings"*"X"[1]), expression("Indiv Variable Loadings"*"X"[2]))
  ggplot2::ggplot(data = norm.dat, ggplot2::aes(x = norm.dat$Type, y = norm.dat$Norm)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = norm.dat$Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.5) +
    ggplot2::labs(y = "Chordal Norm", x = "Type") +
    ggplot2::facet_grid(norm.dat$JVE_2 ~ norm.dat$JVE_1, labeller = ggplot2::label_parsed()) +
    ggplot2::scale_x_discrete(limits = labs, labels = labs.ex) +
    ggplot2::scale_fill_manual(values=cols) +
    ggplot2::scale_colour_manual(values=cols) +
    ggplot2::theme_bw() + ggplot2::coord_cartesian(ylim = c(0, y.max)) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(face = "bold", hjust = 01, angle = 70, size = text.size-3),
          text = ggplot2::element_text(size = text.size))
}

