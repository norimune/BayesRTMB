library(MASS)
#library(EFA.dimensions)
#library(CCA)

map <- function(z) {
  if (abs(max(z) - min(z)) < 1e-10) return(z[1])
  d <- density(z)
  d$x[which.max(d$y)]
}


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

plot_MDU<- function(delta,theta,phi,i=1,j=2,R=3,minus=c(1,1), alpha = 0.2){
  D <- ncol(delta)
  minus2 <- numeric(D)
  for(d in 1:D) if(d==i) minus2[d] <- minus[1]
  for(d in 1:D) if(d==j) minus2[d] <- minus[2]
  M <- length(phi)
  plot.circle <- function(x, y, r){
    theta <- seq(-pi, pi, length=100)
    points(x+r*cos(theta), y+r*sin(theta),type="l",col=rgb(0,0,1,alpha=alpha))
  }
  plot(c(-R,R), c(-R,R),type="n")
  text(minus2[i]*delta[,i],minus2[j]*delta[,j])
  for(m in 1:M){
    plot.circle(minus2[i]*delta[m,i], minus2[j]*delta[m,j], (phi[m]-min(phi)))
  }
  knl <- kde2d(minus2[i]*theta[,i],minus2[j]*theta[,j],n=50,lims=c(c(-R,R),c(-R,R)))
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
