## Nesse script vamos descrever alguns métodos comparativos, fazer simulações e
##    explorar as suas superfícies de likelihood. Vamos mais a frente voltar a
##    esses mesmos métodos para explorar a performance do MCMC.

set.seed(1234)

## Simulando o "andar do bebado."
bebado <- function(passos, sd, start){
  trace <- vector()
  trace[1] <- start
  for(i in 2:passos){
    trace[i] <- trace[i-1] + rnorm(n = 1, mean = 0, sd = sd)
  }
  return(trace)
}

## Usando a função para simular uma realização do processo.
res <- bebado(passos = 100, sd = 0.5, start = 10)
## Note que o trajeto é aleatório.
plot(x=1:100, y=res, type = "l", xlab="tempo", ylab="posição")
abline(h = 10, col="red")

## Repetindo várias vezes podemos notar que, em média, a expectativa
##    é que o valor não mude com o tempo.
## No processo de Brownian Motion a variância aumenta com o tempo
##    mas a média não muda.

many.trace <- sapply(X = 1:10, FUN = function(x) bebado(passos = 100
                                                        , sd = 0.5, start = 10) )
mx <- max( apply(many.trace, MARGIN = 2, FUN = max) ) # valor máximo
mn <- min( apply(many.trace, MARGIN = 2, FUN = min) ) # valor mínimo
plot(x=1:100, y=many.trace[,1], type = "l", ylim=c(mn, mx)
     , xlab="tempo", ylab="posição")
for(i in 2:ncol(many.trace)) lines(x=1:100, y=many.trace[,i])
abline(h = 10, col="red")

## Agora vamos demonstrar o efeito da taxa do modelo BM:
## Taxas altas aumentam o acúmulo de variância com o tempo.

par(mfrow= c(1,2))
low.rate <- sapply(X = 1:10, FUN = function(x) bebado(passos = 100
                                                        , sd = 0.02, start = 10) )
high.rate <- sapply(X = 1:10, FUN = function(x) bebado(passos = 100
                                                      , sd = 0.1, start = 10) )
plot(x=1:100, y=low.rate[,1], type = "l", ylim=c(8, 12)
     , xlab="tempo", ylab="posição")
for(i in 2:ncol(low.rate)) lines(x=1:100, y=low.rate[,i])
abline(h = 10, col="red")
plot(x=1:100, y=high.rate[,1], type = "l", ylim=c(8, 12)
     , xlab="tempo", ylab="posição")
for(i in 2:ncol(high.rate)) lines(x=1:100, y=high.rate[,i])
abline(h = 10, col="red")

## Nos vimos como funciona o movimento Browniano. Sabemos que temos um
##    valor inicial e depois este valor é modificado de acordo com realizações
##    de uma distribuição normal com média 0.
## Mas como esse processo pode ser visto em uma arvore filogenética?

if(!require(phytools)) {install.packages("phytools"); library(phytools)}
if(!require(geiger)) {install.packages("geiger"); library(geiger)}
if(!require(TreeSim)) {install.packages("TreeSim"); library(TreeSim)}

## Simulamos uma arvore:
phy <- sim.bd.taxa(n = 10, numbsim = 1, lambda = 1, mu = 0)[[1]]
plot.phylo(phy, direction="upwards", edge.color = "blue", edge.width = 6
           , show.tip.label = FALSE)
## Note que aqui o parâmetro é sigma^2 e não sigma como nos exemplos anteriores.
x <- bmPlot(tree = phy, type = "BM", sig2=0.02^2)

## Quais valores foram simulados para os tips?
plot(phy, direction = "upwards")
x$x[1:10]

## Likelihood function for the Brownian Motion model:

## Simulando arvore e dados.
phy <- rtree(100)
phy <- compute.brlen(phy)
dt <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM")[,,1]

singleBML<-function(phy, x, sigsq, mean) {
  ## phy = filogenia.
  ## x = dados.
  ## sigsq = sigma^2, a taxa do modelo.
  ## mean = valor da raiz.
  # count the number of tips in the tree
  n <- length(phy$tip.label)
  # calculate C, the phylogenetic variance-covariance matrix
  C <- vcv.phylo(phy)
  m <- match(rownames(C), names(x))
  x <- x[m]
  V <- sigsq*C
  lnlNum <- -0.5*(x-mean) %*% solve(V) %*% (x-mean)
  lnlDen <- log(sqrt((2*pi)^n*det(V)))
  L <- lnlNum-lnlDen
  L
}

## Making a likelihood surface for the BM model:
## These are x and y or x1 and x2 for the plot.
rep <- 20
rate <- seq(0.1, 10, length.out = rep)
root <- seq(-10, +10, length.out = rep)
z <- matrix(data = NA, nrow = rep, ncol = rep)
for(i in 1:rep){
  for(j in 1:rep){
   z[i,j] <- singleBML(phy = phy, x = dt, sigsq = rate[i], mean = root[j])
  }  
}

persp(x=rate, y=root, z=z, theta = -30, phi = 20, xlab = "sigma^2", ylab = "phylo.mean", zlab = "log(likelihood)", border = "black", col = "grey", r = 4)

###################################################################
###################################################################
## Modelo Ornstein-Uhlenbeck (OU).

## O modelo OU adiciona um parâmetro para o valo do ótimo e uma força de seleção que
##    puxa o valor da característica para mais próximo deste ótimo.
## OU ainda possuí um componente de BM, por isso, a simulação de um processo.
##    em um modelo OU será parecida com o processo em um modelo BM.

## Simulando o modelo OU.
simOU <- function(passos, sd, start, alpha, tetha){
  trace <- vector()
  trace[1] <- start
  for(i in 2:passos){
    trace[i] <- trace[i-1] + rnorm(n=1,mean=0,sd=sd) + (alpha * (tetha - trace[i-1]) )
  }
  return(trace)
}

OUtrace <- simOU(passos=300, sd=0.08, start=10, alpha=0.02, tetha=14)
plot(x=1:300, y=OUtrace, type="l", ylab="posição", xlab="tempo")
abline(h=10, col="red", lwd=2)
abline(h=14, col="blue", lwd=2)

par(mfrow = c(1,2))
many.OU <- sapply(X = 1:10, FUN = function(x) simOU(passos=300, sd=0.08, start=10, alpha=0.02, tetha=14) )
plot(x=1:300, y=many.OU[,1], type = "l", xlab="tempo", ylab="posição", main = "alpha = 0.02" )
for(i in 2:ncol(many.OU)) lines(x=1:300, y=many.OU[,i])
abline(h = 10, col="red")
abline(h = 14, col="blue")

many.OU <- sapply(X = 1:10, FUN = function(x) simOU(passos=300, sd=0.08, start=10, alpha=0.06, tetha=14) )
plot(x=1:300, y=many.OU[,1], type = "l", xlab="tempo", ylab="posição", main="alpha = 0.06")
for(i in 2:ncol(many.OU)) lines(x=1:300, y=many.OU[,i])
abline(h = 10, col="red")
abline(h = 14, col="blue")

par(mfrow = c(1,2))
many.OU <- sapply(X = 1:10, FUN = function(x) simOU(passos=300, sd=0.08, start=10, alpha=0.02, tetha=14) )
plot(x=1:300, y=many.OU[,1], type = "l", xlab="tempo", ylab="posição", main = "sigma = 0.08" )
for(i in 2:ncol(many.OU)) lines(x=1:300, y=many.OU[,i])
abline(h = 10, col="red")
abline(h = 14, col="blue")

many.OU <- sapply(X = 1:10, FUN = function(x) simOU(passos=300, sd=0.16, start=10, alpha=0.02, tetha=14) )
plot(x=1:300, y=many.OU[,1], type = "l", xlab="tempo", ylab="posição", main="sigma = 0.16")
for(i in 2:ncol(many.OU)) lines(x=1:300, y=many.OU[,i])
abline(h = 10, col="red")
abline(h = 14, col="blue")
