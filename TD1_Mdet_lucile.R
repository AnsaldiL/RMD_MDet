
library(scatterplot3d)
library(rgl)

#initialisation
x0 <- 1
K0 <- 1
lambda <- 1
sigmaa <- 1
r <-1

a <- function(x, y, sigmaa){
  return (exp((-0.5*(x-y)**2)/sigmaa**2))
}

K <- function(x, x0, lambda, sigmaa, K0){
  return (((K0-lambda*(x-x0)**2)/sigmaa**2))
}

s <- function(x,y){
  return(r*(1-a(x,y,sigmaa)*(K(y, x0, lambda, sigmaa, K0)/K(x, x0, lambda, sigmaa, K0))))
}


plot3d(s, xlim=c(0.1,1.9), ylim=c(0.1, 1.9), col="yellow")



planes3d(a = 0, b = 0, c = 1, d = 0, col = "red", alpha = 0.5, add = TRUE)

# Add a grid for better visualization
grid3d(c("x", "y", "z"))

# Optionally, you can set the aspect ratio to make the plot look better
aspect3d(1, 1, 1)

# View the plot
rglwidget()

contour(x, y, z, main="Contour Plot Example", xlab="X", ylab="Y") ## à tester !



### avec scatterplot3d, il faut avoir des vecteurs de la même taille... c'est compliqué...

x <- c(seq(0.01,1.99,0.01))
Lx=length(x)
x1 <- rep(0,Lx**2)
x2 <- rep(0,Lx**2)
y <- rep(0,Lx**2)

for (i in 1:Lx){
  for (j in 1:Lx){
    x1[i*j] = x[i]
    x2[i*j] = x[j]
    z[i*j] = s(x,y,sigmaa, x0, lambda, K0, r)
  }
}

scatterplot3d(x,y,z, color="blue", pch = 16, main = "3D Scatter Plot", type="p" )


y <- c(seq(0.01,1.99,0.01))
z <- c()
for (xi in x){
  for (xj in y){
    z<-c(z,s(x,y,sigmaa, x0, lambda, K0, r))
  }
}
length(x)
length(y)
length(z)

scatterplot3d(x,y,z, color="blue", pch = 16, main = "3D Scatter Plot", type="p" )




############################### 10/10/23 : il s'agit de simuler un modèle de LK pour N espèces.



#gnagngan bonnes pratiques : clear, commenter les lignes de code

library(deSolve)
library(tidyverse)



# valeurs initiales
K0 <- 1
lambda <- 1
x0 <- 1/2
xmin <- 0
xmax <- 1
N <- 8
sigma <- 0.5
r <- 0.5
epsilon <- 10^(-6)

temps <- seq(0,10000,1)
nbarre <- rep(K0/N, N)

#vecteur de paramètres
param <- c(K0, lambda, x0, xmin, xmax, N, sigma,r, epsilon)


matriceM <- function (param){
  #param
  K0 <- param[1]
  lambda <- param[2]
  x0 <- param[3]
  xmin <- param[4]
  xmax <- param[5]
  N <- param[6]
  sigma <- param[7]
  r <- param[8]
  epsilon <- param[9]
  K <- c()
  

  X = seq(xmin, xmax, length.out = N)

  
  for (x in X){
    K <- c(K, max((K0 -lambda*(x-x0)^2),0+epsilon)) #epsilon pour ne pas avoir de 0 ie, pas de pb quand on divise par 0
  }
  Kmoins <- K^(-1)
  B <- matrix(data = Kmoins, ncol=N, nrow=N, byrow=FALSE)
  B <- t(B) #ici on fait la transposée de B car le byrow n'a pas marché...
  
  D1 <- matrix(data=X, nrow= N, ncol = N, byrow= TRUE)
  D2 <- matrix(data=X, nrow= N, ncol = N, byrow= FALSE)
  C <- D2 - D1
  
  A <- exp(-(C^2)/(2*sigma^2))
  
  M <- A*B
  
  return(M)
}



EDO_LK <- function(t, y, param){
  #param
  r <- param[[1]][8]
  M <- param[[2]]

  dndt <- r*y*(1-y%*%M) #on écrit l'équadif
    
  return(list(dndt))
  }

MatM <- matriceM(param)

param <- list(param, MatM)

solution <- ode(y = nbarre, times = temps, func = EDO_LK, parms = param)


## évolution de la proportion du trait en fonction du temps
solution %>%
  as_tibble() %>%
  mutate_all(as.numeric) %>%
  pivot_longer(-time, names_to = "variable", values_to = "value") %>%
  ggplot() +
  aes(x = time, y = value, color = variable) +
  geom_line()



etatfinal = solution[nrow(solution),2:ncol(solution)]
X = seq(xmin, xmax, length.out = N)

plot(NULL, ylim = c(0, 100), xlim= c(xmin, xmax), xlab="valeur du trait", ylab="temps", type='l')
for (rgtrait in 1:length(etatfinal)){
  if (etatfinal[rgtrait] > 0.001){
    abline(v = X[rgtrait])
  }
}


########################################### 13/10/23 #############################################
## Question de recherche : comment est-ce que la valeur de sigma influe sur le nombre de phénotypes sélectionnés ?

valsigma <- seq(0.1,3,0.01)
NBtrait <- c() #on va stocker les nombre de trait dans ce vecteur

for (sig in valsigma){ #on simule pour différentes valeurs de sigma 
  sigma <- sig
  param <- c(K0, lambda, x0, xmin, xmax, N, sigma,r, epsilon)
  M <- matriceM(param)
  param <- list(param, M)
  
  
  solution <- ode(y = nbarre0, times = temps, func = EDO_LK, parms = param)
  
  etatfinal = solution[nrow(solution),2:ncol(solution)] # à la fin de la simulation...
  
  nb <- 0 #compte cb de traits sont fixés
  
  for (rgtrait in 1:length(etatfinal)){
    if (etatfinal[rgtrait] > 0.01){ # on considère le trait fixé si sa densité est suppérieure à 0.01
      nb <- nb + 1
    }
  }
  
  NBtrait <- c(NBtrait, nb)
}

# On représente le nombre de phénotype selectionné en fonction de la valeur de sigma
plot(x=valsigma, y=NBtrait)

##################### introduction de mutations
sigma = 0.5
N <- 20
param <- c(K0, lambda, x0, xmin, xmax, N, sigma, r, epsilon)
MatM <- matriceM(param)
param <- list(param, MatM)

M <- 100
p <- 0.1
X <- seq(xmin, xmax, length.out = N)

z <- rep(0,N)
z[1] <- 1

nbarre <- z #rep(K0/N, N)

nTmat = matrix(nrow=M, ncol=length(X))

for(m in 1:M){
  solution <- ode(y = nbarre, times = temps, func = EDO_LK, parms = param)
  nT = solution[nrow(solution),2:ncol(solution)]
  for (rgtrait in 1:length(nT)){
    if (nT[rgtrait] > 0.01){
      nTmute <- nT[1:length(nT)-1]*p
      a = runif(1,0,1) 
      if (a>0.5){ #mutation à droite
        nT[1:length(nT)-1] = nT[1:length(nT)-1] - nTmute
        nT[2:length(nT)] = nT[2:length(nT)] + nTmute
      }
      #else{ #mutation à gauche
        #nT[1:length(nT)-1] = nT[1:length(nT)-1] + nTmute
        #nT[2:length(nT)] = nT[2:length(nT)] - nTmute
      #}
    }
  }
  nbarre <- nT
  nTmat[m,] <- nT
}
 

 
#### et on plot !


colnames(nTmat) = c(1:N) # on renome pour acceder à l'indice ensuite

nTmat = as.data.frame(nTmat) 

temps_evol <- seq(1, M, 1)
nTmat$time <- temps_evol # ajout du temps


nTmat_long = nTmat%>%
  pivot_longer(-time, values_to = "dens", names_to = "phenotype")

nTmat_long$trait = X[ as.numeric(nTmat_long$phenotype)] # étiquette du trait



ggplot(nTmat_long)+
  geom_contour(aes(trait, time,  z=dens))+
  scale_fill_gradient2(low = "white" ,
                       high = "red")+
  main_theme

ggplot(nTmat_long)+
  geom_raster(aes(trait, time,  fill=dens))+
  scale_fill_gradient2(low = "white" ,
                       high = "red")+
  main_theme




main_theme = theme_bw()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=22),
        axis.text.y = element_text(colour = "black", size=22),
        legend.title = element_text(colour = "black", size=20),
        legend.title.align=0.5,
        legend.text = element_text(colour = "black", size=18),
        axis.title=element_text(size=28))
