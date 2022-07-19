#################################################################################################################################
#########################################        Conduct JIVE analyses of HCP data:        ######################################
##################################      Corr from fMRI and Streamline counts from DTI      ######################################
######################  Author: Raphiel J. Murden; Advised by Dr. Ben Risk                                               ########
######################
# Notes
# 01 JUNE 2022: Paper was submitted to JRSSC in April 2021. After 13.5 months we have not received reviews. 
#               Script is being reviewed and restructured for easier reproducibility. After revision the  
#               script will be posted to GitHub.
#################################################################################################################################


library(dplyr); library(mgcv); library(MASS); library(ICC); library(psych); library(lme4); library(arsenal)
library(RColorBrewer); library(grid); library(gridExtra); library(reshape); library(scales)

set.seed(0)
prog.dir = "H:/My Documents/FunctionalVsStructuralConnectivity/Programs/"
dat.dir = "H:/My Documents/FunctionalVsStructuralConnectivity/Data"
imgs.fldr = "H:/My Documents/FunctionalVsStructuralConnectivity/Results"

ajive.dir = "H:/My Documents/Applications2/r_AJIVE/R"
files= list.files(ajive.dir)
for (i in files) source(file.path(ajive.dir, i))

source(file.path(prog.dir, "Functions", "Functions_for_PJIVE.R"))

########Load SConn and FConn Data
sconn = read.csv(file.path(dat.dir, "SConn_Desikan.csv"), header = F) ##data are stored in p-by-n format
fconn = read.csv(file.path(dat.dir, "Z_TransCorrFC_Desikan.csv"), header = F) ##data are stored in p-by-n format

# load(file.path(dat.dir, "fconn.RData")) #data matrix is called fconn
# save(fconn, file = file.path(dat.dir, "fconn_cluster.RData"), version = 2) ##data would not load onto cluster bc of R.Version discrepancy
# fconn = fconn[FC_indices[,1],] ##Reduce FC matrix to only subjects with SC data

###Center Fconn, then transform and center Sconn data
fconn.cent = scale(fconn, scale = F)

sconn.sd = apply(sconn, 2, sd) #take SD of columns/edges
sconn.non0 = sconn[,sconn.sd!=0] #remove columns where SD = 0
# ncol(sconn) - ncol(sconn.non0) #reduces the number of edges by 1025
sconn.log = log(sconn.non0 + 1) #take log(sconn + 1)
sconn.cent = scale(sconn.log, scale = F)

###Center Fconn and Sconn data
fconn.mean = colMeans(fconn)
sconn.mean = colMeans(sconn.log)
ncol(fconn.cent)
ncol(sconn.cent)

FC.cent.svd = svd(fconn.cent)
SC.cent.svd = svd(sconn.cent)

FC.cent.eigvals = FC.cent.svd$d^2
SC.cent.eigvals = SC.cent.svd$d^2

cb.cols = c("#112e51", "#205493", "#0071bc", "#046b99", "#2e8540", "#4c2c92", "#212121",
            "#323a45", "#5b616b", "#494440", "#981b1e", "#cd2026", "#e31c3d", "#cc9900",
            "#ff9933", "#ff33cc", "#33ccff", "#cc0099", "#ff99cc", "#990099", "#33ffcc")

######Load Atlases
#mmp.yeo.atlas = read.csv(file.path(atlas.dir, "Yeo_MMP.csv"))
#dest.yeo.atlas = read.csv(file.path(atlas.dir, "Yeo_Dest.csv"))
load(file.path(atlas.dir, "Yeo_MMP.RData"))
load(file.path(atlas.dir, "Yeo_Dest.RData"))

r1.95 = which.max(cumsum(FC.cent.eigvals)>=0.95*sum(FC.cent.eigvals))
r1.elbow = 7
r1.elbow.2 = 25
r1.iJIVE = 53
FC.TotVarEx.iJIVE = sum(FC.cent.eigvals[1:r1.iJIVE])/sum(FC.cent.eigvals) # = 0.5740657
FC.TotVarEx = sum(FC.cent.eigvals[1:r1.elbow])/sum(FC.cent.eigvals) # = 0.5740657
FC.TotVarEx.2 = sum(FC.cent.eigvals[1:r1.elbow.2])/sum(FC.cent.eigvals) # = 0.5806622
pdf(file.path(imgs.fldr, "FC_HCP_ScreePlot.pdf"))
plot(log(FC.cent.eigvals[1:700]), xlab = "Index", ylab = "log(Eigenvalues)", main = "Scree Plot for FC") 
abline(v = r1.elbow + 0.5, col = cb.cols[3], lwd = 1.75)
abline(v = r1.elbow.2 + 0.5, col = cb.cols[14], lwd = 1.75)
abline(v = r1.95 + 0.5, col = cb.cols[11], lwd = 1.75)
abline(v = r1.iJIVE + 0.5, col = cb.cols[5], lwd = 1.75)
legend("topright", c(paste("Elbow (", 100*round(FC.TotVarEx, 2) ,"% Var.): r = ", r1.elbow, sep = ""), 
                     paste("Elbow (", 100*round(FC.TotVarEx.2, 2) ,"% Var.): r = ", r1.elbow.2, sep = ""), 
                     paste("Perm. Test (", 100*round(FC.TotVarEx.iJIVE, 2) ,"% Var.): r = ", r1.iJIVE, sep = ""),
                     paste("95 % Var.: r =", r1.95)), cex = 0.8,
       col = cb.cols[c(3,14,5,11)], text.col = cb.cols[c(3,14,5,11)])
dev.off()

#elbow of scree plots looks to occur around index 10 for FC data
r2.95 = which.max(cumsum(SC.cent.eigvals)>=0.95*sum(SC.cent.eigvals))
r2.elbow = 10
r2.elbow.2 = 25
r2.elbow.3 = 40
r2.iJIVE = 97
SC.TotVarEx.iJIVE = sum(SC.cent.eigvals[1:r2.iJIVE])/sum(SC.cent.eigvals) # = 0.5740657
SC.TotVarEx = sum(SC.cent.eigvals[1:r2.elbow])/sum(SC.cent.eigvals) # = 0.4133836
SC.TotVarEx.2 = sum(SC.cent.eigvals[1:r2.elbow.2])/sum(SC.cent.eigvals) # = 0.4133836
SC.TotVarEx.3 = sum(SC.cent.eigvals[1:r2.elbow.3])/sum(SC.cent.eigvals) # = 0.4133836
pdf(file.path(imgs.fldr, "SC_HCP_ScreePlot.pdf"))
plot(log(SC.cent.eigvals[1:700]), xlab = "Index", ylab = "log(Eigenvalues)", main = "Scree Plot for SC") 
abline(v = r2.elbow + 0.5, col = cb.cols[3], lwd = 1.75)
abline(v = r2.elbow.2 + 0.5, col = cb.cols[14], lwd = 1.75)
abline(v = r2.elbow.3 + 0.5, col = cb.cols[18], lwd = 1.75)
abline(v = r2.iJIVE + 0.5, col = cb.cols[5], lwd = 1.75)
abline(v = r2.95 + 0.5, col = cb.cols[11], lwd = 1.75)
legend("topright", c(paste("Elbow (", 100*round(SC.TotVarEx, 2) ,"% Var.): r = ", r2.elbow, sep = ""),
                     paste("Elbow (", 100*round(SC.TotVarEx.2, 2) ,"% Var.): r = ", r2.elbow.2, sep = ""),
                     paste("Elbow (", 100*round(SC.TotVarEx.3, 2) ,"% Var.): r = ", r2.elbow.3, sep = ""),
                     paste("Perm. Test (", 100*round(SC.TotVarEx.iJIVE, 2) ,"% Var.): r = ", r2.iJIVE, sep = ""),
                     paste("95 % Var.: r =", r2.95)), cex = 0.8,
       col = cb.cols[c(3,14,18,5,11)], text.col = cb.cols[c(3,14,18,5,11)])
dev.off()


