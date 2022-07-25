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
#' @param Preds (logical) do the replicate results contain correlations between predicted and true joint subject scores. Default is FALSE
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
  Method = sub(".oracle", "-r_k", Method, fixed = TRUE); Method = sub(".over", "-Over", Method, fixed = TRUE)
  Method = sub("Elbow", "r_k", Method, fixed = TRUE); Method = sub(".95.", "-Over", Method, fixed = TRUE)
  Method = sub("iJI", "R.JI", Method, fixed = TRUE); Method = sub("AJIVE.", "AJIVE-", Method, fixed = TRUE)
  Method = sub("CC", "CJIVE", Method, fixed = TRUE); Method = sub("cJIVE", "CJIVE", Method, fixed = TRUE)
  Method[which(Method == "CJIVE.Oracle")] = "CJIVE-r_k"
  Method[which(Method == "CJIVE.r_k")] = "CJIVE-r_k"
  Method[which(Method == "AJIVE-Oracle")] = "AJIVE-r_k"
  Method[which(Method == "R.JIVE.Free")] = "R.JIVE-Free"
  Method[which(Method %in% c("CJIVE.Over", "cJIVE.Over"))] = "CJIVE-Over"
  Method = factor(Method, levels = c("CJIVE-r_k","CJIVE-Over","AJIVE-r_k","AJIVE-Over","R.JIVE-Free"),
                  labels = c("CJIVE-r_k","CJIVE-Over","AJIVE-r_k","AJIVE-Over","R.JIVE-Free"))

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
  Method[which(Method == "AJIVE")] = "AJIVE-r_k"; Method = sub(".w", "-Over", Method, fixed = TRUE);
  Method[which(Method == "CC.Oracle")] = "CJIVE-r_k"; Method[which(Method == "CC-Over")] = "CJIVE-Over"
  Method[which(Method == "R.JIVE.Free")] = "R.JIVE-Free"
  Method = sub("JIVE.", "JIVE-", Method, fixed = TRUE); Method = sub("[[:punct:]]$", "", Method, fixed = TRUE)
  Method[which(Method == "AJIVE-Oracle")] = "AJIVE-r_k"; Method[which(Method == "CJIVE-Oracle")] = "CJIVE-r_k"
  Method = sub("R.JIVE-Oracl", "R.JIVE-Oracle", Method); Method = sub("R.JIVE-Free.", "R.JIVE-Free", Method)

  Method = factor(Method, levels = c("CJIVE-r_k","CJIVE-Over","AJIVE-r_k","AJIVE-Over","R.JIVE-Free","R.JIVE-Oracle"),
                  labels = c("CJIVE-r_k","CJIVE-Over","AJIVE-r_k","AJIVE-Over","R.JIVE-Free","R.JIVE-Oracle"))

  Type = rep(substring(sim.names[grep("Scores", sim.names)], 7), each = num.sims)
  Type = sub("w.", "", Type); Type = gsub("[.]", " ", Type)
  Type = gsub("X", "X[1]", Type); Type = gsub("Y", "X[2]", Type)
  Type = gsub("cle ", "", Type); Type = gsub("Over ", "", Type)
  Type = gsub("Ora", "", Type); Type = gsub("Free ", "", Type)
  Type = gsub(" Joint", "Joint", Type); Type = gsub(" Indiv", "Indiv", Type)
  Type = gsub(" Subj", "", Type); Type = gsub(" Variable", "", Type)
  Type = as.factor(Type)

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
  Method[which(Method == "AJIVE")] = "AJIVE-r_k"; Method = sub(".w", "-Over", Method, fixed = TRUE);
  Method[which(Method == "CC.Oracle")] = "CJIVE-r_k"; Method[which(Method == "CC-Over")] = "CJIVE-Over"
  Method = sub("JIVE.", "JIVE-", Method, fixed = TRUE); Method = sub("[[:punct:]]$", "", Method, fixed = TRUE)
  Method = sub("JIVE.", "JIVE-", Method, fixed = TRUE); Method = sub("[[:punct:]]$", "", Method, fixed = TRUE)
  Method[which(Method == "AJIVE-Oracle")] = "AJIVE-r_k"; Method[which(Method == "CJIVE-Oracle")] = "CJIVE-r_k"
  Method = sub("R.JIVE-Oracl", "R.JIVE-Oracle", Method); Method = sub("R.JIVE-Free.", "R.JIVE-Free", Method)
  Method = factor(Method, levels = c("CJIVE-r_k","CJIVE-Over","AJIVE-r_k","AJIVE-Over","R.JIVE-Free","R.JIVE-Oracle"),
                  labels = c("CJIVE-r_k","CJIVE-Over","AJIVE-r_k","AJIVE-Over","R.JIVE-Free","R.JIVE-Oracle"))

  Type = factor(rep(substring(sim.names[grep("Loads", sim.names)], 7), each = num.sims))
  Type = sub("w.", "", Type); Type = gsub("[.]", " ", Type)
  Type = gsub("X", "X[1]", Type); Type = gsub("Y", "X[2]", Type)
  Type = gsub("cle ", "", Type); Type = gsub("Over ", "", Type)
  Type = gsub("Ora", "", Type); Type = gsub("Free ", "", Type)
  Type = gsub(" Joint", "Joint", Type); Type = gsub(" Indiv", "Indiv", Type)
  Type = gsub(" Subj", "", Type); Type = gsub(" Variable", "", Type)
  Type = as.factor(Type)

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

#####################         Plot chordal norms from Simulation Study      #######################
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
#' @param x.lab.angle angle at which x-axis labels are tilted
#'
#' @return graphical display (via ggplot2)
#' @export
gg.norm.plot<-function(norm.dat, cols, show.legend = FALSE, text.size, lty = 1, y.max = 1, x.lab.angle = 70){
  labs = levels(norm.dat$Type)[c(3,6,7,1,2,4,5)]
  labs.ex = c(expression("Joint Subj Scores"), expression("Joint Loadings"*"X"[1]), expression("Joint Loadings"*"X"[2]),
              expression("Indiv Subj Scores"*"X"[1]), expression("Indiv Subj Scores"*"X"[2]),
              expression("Indiv Loadings"*"X"[1]), expression("Indiv Loadings"*"X"[2]))
  labs.rows = levels(norm.dat$JVE_2)
  labs.cols = levels(norm.dat$JVE_1)
  ggplot2::ggplot(data = norm.dat, ggplot2::aes(x = norm.dat$Type, y = norm.dat$Norm)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = norm.dat$Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.5) +
    ggplot2::labs(y = "Chordal Norm", x = "Type") +
    ggplot2::facet_grid(JVE_2 ~ JVE_1, labeller = ggplot2::label_parsed) +
    ggplot2::scale_x_discrete(limits = levels(norm.dat$Type)[c(3,6,7,1,2,4,5)], labels = labs.ex) +
    ggplot2::scale_fill_manual(values=cols, name = "Method") +
    ggplot2::scale_colour_manual(values=cols) +
    ggplot2::theme_bw() + ggplot2::coord_cartesian(ylim = c(0, y.max)) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(face = "bold", hjust = 01, angle = x.lab.angle, size = text.size-3),
          text = ggplot2::element_text(size = text.size))
}

#####################         Plot chordal norms from Simulation Study      #######################
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
#' @param x.lab.angle angle at which x-axis labels are tilted
#'
#' @return graphical display (via ggplot2)
#' @export
gg.score.norm.plot<-function(norm.dat, cols, show.legend = FALSE, text.size, lty = 1, y.max = 1, x.lab.angle = 70){
  norm.dat = norm.dat[grep("Score",norm.dat$Type),]
  labs = levels(norm.dat$Type)[grep("Score",levels(norm.dat$Type))][c(3,1,2)]
  labs.ex = c("Joint Scores", expression("Indiv Scores"*" X"[1]), expression("Indiv Scores"*" X"[2]))
  ggplot2::ggplot(data = norm.dat, ggplot2::aes(x = norm.dat$Type, y = norm.dat$Norm)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = norm.dat$Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.5) +
    ggplot2::labs(y = "Chordal Norm", x = "Type") +
    ggplot2::facet_grid(norm.dat$JVE_2 ~ norm.dat$JVE_1, labeller = ggplot2::label_parsed) +
    ggplot2::scale_x_discrete(limits = labs, labels = labs.ex) +
    ggplot2::scale_fill_manual(values=cols, name = "Method") +
    ggplot2::scale_colour_manual(values=cols) +
    ggplot2::theme_bw() + ggplot2::coord_cartesian(ylim = c(0, y.max)) +
  ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(face = "bold", hjust = 01, angle = x.lab.angle, size = text.size-3),
                 text = ggplot2::element_text(size = text.size))
}

#####################         Plot chordal norms from Simulation Study      #######################
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
#' @param x.lab.angle angle at which x-axis labels are tilted
#'
#' @return graphical display (via ggplot2)
#' @export
gg.load.norm.plot<-function(norm.dat, cols, show.legend = FALSE, text.size, lty = 1, y.max = 1, x.lab.angle = 70){
  norm.dat = norm.dat[grep("Load",norm.dat$Type),]
  labs = levels(norm.dat$Type)[grep("Load",levels(norm.dat$Type))][c(3,4,1,2)]
  labs.ex = c(expression("Joint Loadings"*" X"[1]),expression("Joint Loadings"*" X"[2]),
              expression("Indiv Loadings"*" X"[1]), expression("Indiv Loadings"*" X"[2]))
  ggplot2::ggplot(data = norm.dat, ggplot2::aes(x = norm.dat$Type, y = norm.dat$Norm)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = norm.dat$Method), position = "dodge", outlier.alpha = 0, show.legend = show.legend, linetype = lty,
                 fatten = 0.5) +
    ggplot2::labs(y = "Chordal Norm", x = "Type") +
    ggplot2::facet_grid(norm.dat$JVE_2 ~ norm.dat$JVE_1, labeller = ggplot2::label_parsed) +
    ggplot2::scale_x_discrete(limits = labs, labels = labs.ex) +
    ggplot2::scale_fill_manual(values=cols, name = "Method") +
    ggplot2::scale_colour_manual(values=cols) +
    ggplot2::theme_bw() + ggplot2::coord_cartesian(ylim = c(0, y.max)) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(face = "bold", hjust = 01, angle = x.lab.angle, size = text.size-3),
          text = ggplot2::element_text(size = text.size))
}

#####################         Plot Pearson correlations from Simulation Study      #######################
#' Function for plotting Pearson correlations between predicted and true subject scores within the simulation study described in CJIVE manuscript
#'
#' @description Graphically displays the center and spread of chordal norms for joint/individual subject score subspaces
#'
#' @param cor.dat data frame with at least the 5 following variables:
#'                   Norm - the value of the norm for a particular subspace;
#'                   Type - the subspace for which the norm is given (i.e., joint/individual score/loading for dataset X1 or X2 (except for joint scores))
#'                   Method - the method by which the subspace was estimated, e.g. CJIVE, AJIVE, R.JIVE
#'                   JVE_1 and JVE_2 - labels describing the proportion of joint variation explained in each dataset (and typically the number of variables in dataset X2)
#' @param cols a vector of colors, must have length equal to the number of methods used in the simulation
#' @param show.legend logical (TRUE/FALSE) for whether a legend should be included in the plot. Default is FALSE
#' @param text.size numeric value for the font size
#'
#' @return graphical display (via ggplot2)
#' @export
gg.corr.plot<-function(cor.dat, cols, show.legend = FALSE, text.size){
  ggplot2::ggplot(data = cor.dat, ggplot2::aes(x = cor.dat$Component, y = cor.dat$Prediction_Correlation)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = cor.dat$Type), position = "dodge", outlier.alpha = 0,
                          show.legend = show.legend, fatten = 0.5) +
    ggplot2::labs(y = "Absolute Pearson Correlation", x = "Component") +
    ggplot2::facet_grid(cor.dat$JVE_2 ~ cor.dat$JVE_1, labeller = ggplot2::label_parsed) +
    ggplot2::scale_fill_manual(values=cols, name = "Method") +
    ggplot2::scale_colour_manual(values=cols) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = text.size))
}

#####################         Plot Ranks CJIVE from Simulation Study      #######################
#' Function for plotting selected joint ranks
#'
#' @description Graphically displays the count of joint ranks selected by each method employed in the simulation study described in the CJIVE manuscript
#'
#' @param rank.dat data frame expected to be built with the functions dplyr::count and tidyr::complete, which should include the following variables
#'                   Rank - numeric values of the rank selected by each method in each replicate simulation
#'                   n - the number of times this value was selected as the rank
#'                   Type - the subspace for which the norm is given (i.e., joint/individual score/loading for dataset X1 or X2 (except for joint scores))
#'                   Method - the method by which the subspace was estimated, e.g. CJIVE, AJIVE, R.JIVE
#'                   JVE_1 and JVE_2 - labels describing the proportion of joint variation explained in each dataset (and typically the number of variables in dataset X2)
#' @param cols a vector of colors, must have length equal to the number of methods used in the simulation
#' @param show.legend logical (TRUE/FALSE) for whether a legend should be included in the plot. Default is FALSE
#' @param text.size numeric value for the font size
#' @param num.sims numeric value for the number of replicates evaluated in each full combination of experimental settings
#'
#' @return graphical display (via ggplot2)
#' @export
gg.rank.plot<-function(rank.dat, cols, show.legend = FALSE, text.size, num.sims){
  ggplot2::ggplot(data = rank.dat, ggplot2::aes(x = rank.dat$Rank, y = rank.dat$n/num.sims, fill=rank.dat$Method)) +
    ggplot2::geom_bar(position = ggplot2::position_dodge(), stat = "identity", width = 0.6, show.legend = show.legend) +
    ggplot2::labs(y = "Proportion", x = "Selected Rank") +
    ggplot2::facet_grid(rank.dat$JVE_2 ~ rank.dat$JVE_1, labeller = ggplot2::label_parsed) +
    ggplot2::scale_fill_manual(values=cols,  name = "Method") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = text.size))
}

#####################         'Melts' correlations from prediction portion of simulation study to prepare for plotting       #######################
#' Converts correlations of predicted to true joint subject scores to a format conducive to ggplot2
#'
#' @description Converts correlations of predicted to true joint subject scores into a format conducive to ggplot2
#'
#' @param sim.dat matrix with each row representing results from a replicate in the simulation study described in CJIVE manuscript
#' @param r.J (Numeric/integer) the joint rank, i.e. number of components in the joint subspace
#' @param p1 number of variables/features in data set X1
#' @param p2 number of variables/features in data set X2
#'
#' @return data frame with seven columns: one each for the joint variance explained in each data set,
#'         one column containing the method by which predictions were obtained,
#'         one column containing the component number (1,...,r.J),
#'
#' @export
##############   Prepare data from prediction portion of simulation study for graphical display      ################
Melt.Sim.Cors<-function(sim.dat,r.J,p1,p2){
  JVE1.labs = c(bquote("R"[J1]^2*"=0.05, p"[1]*"="*.(p1)), bquote("R"[J1]^2*"=0.5, p"[1]*"="*.(p1)))
  JVE2.labs = c(bquote("R"[J2]^2*"=0.05, p"[2]*"="*.(p2)), bquote("R"[J2]^2*"=0.5, p"[1]*"="*.(p2)))

  CC.PredCors.melt = reshape2::melt(sim.dat, id.vars = c("JntVarEx1", "JntVarEx2"),
                                    measure.vars = paste("CC_Pred_Corr", 1:r.J, sep = ""), variable_name = "Component",
                                    value.name = "Prediction_Correlation")
  colnames(CC.PredCors.melt)[which(colnames(CC.PredCors.melt)=="variable")]="Component"
  colnames(CC.PredCors.melt)[which(colnames(CC.PredCors.melt)=="value")]="Prediction_Correlation"

  levels(CC.PredCors.melt$Component) = paste("No.", 1:r.J)
  CC.PredCors.melt$Type = factor(1, levels = 1:3, labels = c("CJIVE-Prediction", "G-inverse Prediction", "R.JIVE-Prediction"))

  A.PredCors.melt = reshape2::melt(sim.dat, id.vars = c("JntVarEx1", "JntVarEx2"),
                                   measure.vars = paste("AJIVE_Pred_Corr", 1:r.J, sep = ""), variable_name = "Component")
  colnames(A.PredCors.melt)[which(colnames(A.PredCors.melt)=="variable")]="Component"
  colnames(A.PredCors.melt)[which(colnames(A.PredCors.melt)=="value")]="Prediction_Correlation"

  levels(A.PredCors.melt$Component) = paste("No.", 1:r.J)
  A.PredCors.melt$Type = factor(2, levels = 1:3, labels = c("CJIVE-Prediction", "G-inverse Prediction", "R.JIVE-Prediction"))

  R.JIVE.PredCors.melt = reshape2::melt(sim.dat, id.vars = c("JntVarEx1", "JntVarEx2"),
                                        measure.vars = paste("R.JIVE_Pred_Corr", 1:r.J, sep = ""), variable_name = "Component",
                                        value.name = "Prediction_Correlation")
  colnames(R.JIVE.PredCors.melt)[which(colnames(R.JIVE.PredCors.melt)=="variable")]="Component"
  colnames(R.JIVE.PredCors.melt)[which(colnames(R.JIVE.PredCors.melt)=="value")]="Prediction_Correlation"

  levels(R.JIVE.PredCors.melt$Component) = paste("No.", 1:r.J)
  R.JIVE.PredCors.melt$Type = factor(3, levels = 1:3, labels = c("CJIVE-Prediction", "G-inverse Prediction", "R.JIVE-Prediction"))

  Pred.Cors.melt = rbind(A.PredCors.melt, CC.PredCors.melt, R.JIVE.PredCors.melt)
  Pred.Cors.melt$Prediction_Correlation = abs(Pred.Cors.melt$Prediction_Correlation)
  Pred.Cors.melt$JVE_1 = factor(Pred.Cors.melt$JntVarEx1, labels = JVE1.labs)
  Pred.Cors.melt$JVE_2 = factor(Pred.Cors.melt$JntVarEx2, labels = JVE2.labs)
  Pred.Cors.melt
}

