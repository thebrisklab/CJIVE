###########           Frobenius Norm of Data Matrix Values         #################
#' Matrix variation (i.e. Frobenius norm)
#'
#' @description Calculates the Frobenius norm of a matrix, which can be used as a measure of total variation
#'
#' @param X a matrix of any size
#'
#' @return The Frobenius norm of the matrix X, calculated as the square root of the sum of squared entries in X
#'
#' @examples
#' X = matrix(rnorm(10), 5,2)
#' MatVar(X)
#' @export

MatVar = function(X){
  sum(X^2)
}

###########           Frobenius Norm of Data Matrix Values         #################
#' Alternative calculation - Matrix variation (i.e. Frobenius norm)
#'
#' @description Calculates the Frobenius norm of a matrix, which can be used as a measure of total variation
#'
#' @param X a matrix of any size
#'
#' @return The Frobenius norm of the matrix X, calculated as the square root of the trace of t(X)%*%X
#'
#' @examples
#' X = matrix(rnorm(10), 5,2)
#' MatVar2(X)
#' @export

MatVar2 = function(X){
  sum(diag(t(X)%*%X))
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
#' Display a heatmap of a matrix (adapted from Erick Lock's show.image function in the r.jive package)
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

##############         Function to sign-correct and scale variable loadings       ################
#' Scale and sign-correct variable loadings to assist interpretation
#' @description Scale loadings for a joint or individual component by its largest absolute value resulting in loadings between -1 and 1.
#'              Loadings are also sign-corrected to result in positive skewness
#'
#' @param loading.comp numeric vector of variable loadings from a JIVE analysis
#'
#' @return numeric vector of loadings which have been scaled and sign-corrected
#' @export
scale_loadings = function(loading.comp){
  x = loading.comp
  x1 = x/max(abs(x), dna.rm = TRUE)
  pos = sign(psych::skew(x1, na.rm = TRUE))
  ((-1)^(pos))*x1
}
