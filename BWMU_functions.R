library(MASS)
library(EFA.dimensions)
library(CCA)

Principal_Axis <- function(A){
  temp <- eigen(A%*%t(A))
  rotA <- array(NA,dim=c(nrow(A),ncol(A)))
  for(d in 1:ncol(A)){
    rotA[,d] <- temp$vectors[,d]*sqrt(temp$values[d])
  }
  temp <- svd(t(rotA)%*%A)
  rotmat <- temp$v%*%t(temp$u)
  return(list(rotatedA=rotA,rotmat=rotmat))
}

plot_MDU<- function(delta,theta,phi,i=1,j=2,R=3){
  M <- length(phi)
  plot.circle <- function(x, y, r){
    theta <- seq(-pi, pi, length=100)
    points(x+r*cos(theta), y+r*sin(theta),type="l",col=rgb(0,0,1,alpha=0.5))
  }
  plot(c(-R,R), c(-R,R),type="n")
  text(delta[,i],delta[,j])
  for(m in 1:M){
    plot.circle(delta[m,i], delta[m,j], (phi[m]-min(phi)))
  }
  knl <- kde2d(theta[,i],theta[,j],n=50,lims=c(c(-R,R),c(-R,R)))
  contour(knl,add=T,drawlabels=T)
}

rotate_to_principal <- function(target_mat, other_mat) {
  svd_res <- svd(target_mat)
  rot_mat <- svd_res$v
  return(list(
    target_rot = target_mat %*% rot_mat,
    other_rot  = other_mat %*% rot_mat
  ))
}
