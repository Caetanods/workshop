## Nesse anexo vamos rodar todas as cadeias para gerar os resultados usados no tutorial.

## Delineamento das simulações:
## Número de tips na filogenia:
## 100 > 50 > 25 > 10
## Priors:
## --Informativo: Exponencial e normal.
## --Não-informativo: Uniforme.
## Estado inicial: Amostrado do prior.
## Número de cadeias para cada simulação: 2
## Número total de cadeias: 16.
## Número de gerações para cada simulação: 10000

## Funções para os priors:
prior.exp <- function(x){
  dexp(x = x, rate = 0.5, log = TRUE)
}
prior.norm <- function(x, mu){
  dnorm(x = x, mean = mu, sd = 10, log = TRUE)
}
prior.unif <- function(x){
  dunif(x = x, min=-100, max=100)
}

## Função de proposal.
sliding.window <- function(x, w){
  ## Sliding window proposal for unbounded trait.
  ## x = the current value.
  ## w = the width parameter of the proposal.
  y <- runif(1, min = x - (w/2), max = x + (w/2) )
  return(y)
}

## Definir a função de likelihood para o modelo BM:
singleBML <- function(phy, x, sigsq, mean){
  ## phy = filogenia.
  ## x = dados.
  ## sigsq = sigma^2, a taxa do modelo.
  ## mean = valor da raiz.
  n <- length(phy$tip.label)
  C <- vcv.phylo(phy)
  m <- match(rownames(C), names(x))
  x <- x[m]
  V <- sigsq*C
  lnlNum <- -0.5*(x-mean) %*% solve(V) %*% (x-mean)
  lnlDen <- log(sqrt((2*pi)^n*det(V)))
  L <- lnlNum-lnlDen
  return(L)
}

mcmc.BM <- function(phy, data, initial.rate, initial.root, prior.rate, prior.root, 
                    w=0.5, gen){
  ## phy = filogenia no formato 'ape'.
  ## data = vetor de dados para as espécies.
  ## initial.rate = valor inicial para rate.
  ## inital.root = valor inicial para root.
  ## prior.root = função de prior para root. Com unico argumento igual ao valor.
  ## prior.rate = função de prior para rate. Com unico argumento igual ao valor.
  ## w = a width (largura) da função uniforme de sliding window.
  ## gen = número de gerações da análise.
  
  ## Gerar vetores para guardar os resultados da cadeia.
  root.chain <- numeric(length = gen) ## Cadeia para root.
  rate.chain <- numeric(length = gen) ## Cadeia para rate.
  log.lik <- numeric(length = gen) ## A log(verossimilhança).
  accept <- numeric(length = gen) ## Para calcular o acceptance ratio.
  
  ## Primeiro estado é o initial state:
  root.chain[1] <- initial.root
  rate.chain[1] <- initial.rate
  accept[1] <- 1 ## Primeira geração é o valor inicial.
  log.lik[1] <- singleBML(phy=phy, x=data, sigsq=initial.rate, 
                          mean=initial.root)
  
  ## Fazer um loop do tamanho de gen. Esse é o MCMC.
  ## NOTE que o loop começa do número 2. Note também que uso 'i-1' para o estado
  ##    atual da cadeia e 'i' para o estado futuro.
  for( i in 2:gen ){
    
    ## Selecionar qual parâmetro fazer o update:
    par.update <- sample(x = c("root","rate"), size = 1, prob = c(0.5,0.5))
    
    ## Fazer o update, accept and reject e salvar o estado:
    if( par.update == "root"){
      ## Proposta de novo valor:
      prop.root <- sliding.window(x = root.chain[i-1], w = w)
      ## Calculo da odds ratio:
      lik.curr.root <- singleBML(phy=phy, x=data, sigsq=rate.chain[i-1], 
                                 mean=root.chain[i-1])
      lik.prop.root <- singleBML(phy=phy, x=data, sigsq=rate.chain[i-1], prop.root)
      prior.curr <- prior.rate(rate.chain[i-1]) + prior.root(root.chain[i-1])
      prior.prop <- prior.rate(rate.chain[i-1]) + prior.root(prop.root)
      ll <- lik.prop.root - lik.curr.root
      pp <- prior.prop - prior.curr
      r.root <- ll + pp ## ratio of log(lik) + ratio of log(prior)
      ## Passo de aceitar ou rejeitar:
      if(exp(r.root) > runif(1)){ ## Aceite.
        accept[i] <- 1
        root.chain[i] <- prop.root
        rate.chain[i] <- rate.chain[i-1]
        log.lik[i] <- lik.prop.root
      } else{ ## Rejeite.
        accept[i] <- 0
        root.chain[i] <- root.chain[i-1]
        rate.chain[i] <- rate.chain[i-1]
        log.lik[i] <- log.lik[i-1]
      }
    } else {
      ## Proposta de novo valor:
      prop.rate <- sliding.window(x = rate.chain[i-1], w = w)
      ## Refletir a proposta para valores positivos:
      ## Veja o slide da estratégia de refletir valores fora do limite.
      if(prop.rate < 0){ prop.rate <- -1 * prop.rate }
      ## Calculo da odds ratio:
      lik.curr.rate <- singleBML(phy=phy, x=data, sigsq=rate.chain[i-1], 
                                 mean=root.chain[i-1])
      lik.prop.rate <- singleBML(phy=phy, x=data, sigsq=prop.rate, 
                                 mean=root.chain[i-1])
      prior.curr <- prior.rate(rate.chain[i-1]) + prior.root(root.chain[i-1])
      prior.prop <- prior.rate(prop.rate) + prior.root(root.chain[i-1])
      ll <- lik.prop.rate - lik.curr.rate
      pp <- prior.prop - prior.curr
      r.rate <- ll + pp ## ratio of log(lik) + ratio of log(prior)
      ## Passo de aceitar ou rejeitar:
      if(exp(r.rate) > runif(1)){ ## Aceite.
        accept[i] <- 1
        root.chain[i] <- root.chain[i-1]
        rate.chain[i] <- prop.rate
        log.lik[i] <- lik.prop.rate
      } else{ ## Rejeite.
        accept[i] <- 0
        root.chain[i] <- root.chain[i-1]
        rate.chain[i] <- rate.chain[i-1]
        log.lik[i] <- log.lik[i-1]
      }
      
    }
    
  } ## Fim do for loop aqui.
  
  return(data.frame(log_lik = log.lik, root = root.chain, rate = rate.chain, 
                    accept = accept))
}

## Fixar o seed:
set.seed(1234)

library(geiger)

## 100 tips; informative priors:
## Nesse caso já temos as simulações prontas:
## Lembre-se de editar o caminho para os seus próprios arquivos.
info.100.mcmc.1 <- readRDS("~/Documents/Academicos/Harmon_Lab/MCMC - Workshop - Brasil/workshop/Slides/Work_6/mcmc.BM.1.rds")
info.100.mcmc.2 <- readRDS("~/Documents/Academicos/Harmon_Lab/MCMC - Workshop - Brasil/workshop/Slides/Work_6/mcmc.BM.2.rds")

## 100 tips; non-informative priors:
phy <- rtree(100); phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.100.mcmc.1 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                    initial.root = init.root, prior.rate = prior.unif, 
                    prior.root = prior.unif, w = 0.5, gen = 10000)
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.100.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                           initial.root = init.root, prior.rate = prior.unif, 
                           prior.root = prior.unif, w = 0.5, gen = 10000)

## 50 tips; informative priors:
phy <- rtree(50); phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
mu <- mean(data)
## Precisamos ajustar a função de prior.norm.
prior <- function(x) dnorm(x, mean = mu, sd = 10, log = TRUE)
init.rate <- rexp(n = 1, rate = 0.5)
init.root <- rnorm(n = 1, mean = mu, sd = 10)
inf.50.mcmc.1 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                           initial.root = init.root, prior.rate = prior.exp, 
                           prior.root = prior, w = 0.5, gen = 50000)
init.rate <- rexp(n = 1, rate = 0.5)
init.root <- rnorm(n = 1, mean = mu, sd = 10)
inf.50.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                           initial.root = init.root, prior.rate = prior.exp, 
                           prior.root = prior, w = 0.5, gen = 50000)

## 50 tips; uninformative priors:
phy <- rtree(50); phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.50.mcmc.1 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                           initial.root = init.root, prior.rate = prior.unif, 
                           prior.root = prior.unif, w = 0.5, gen = 50000)
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.50.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                           initial.root = init.root, prior.rate = prior.unif, 
                           prior.root = prior.unif, w = 0.5, gen = 50000)

## 25 tips; informative priors:
phy <- rtree(25); phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
mu <- mean(data)
## Precisamos ajustar a função de prior.norm.
prior <- function(x) dnorm(x, mean = mu, sd = 10, log = TRUE)
init.rate <- rexp(n = 1, rate = 0.5)
init.root <- rnorm(n = 1, mean = mu, sd = 10)
inf.25.mcmc.1 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                          initial.root = init.root, prior.rate = prior.exp, 
                          prior.root = prior, w = 0.5, gen = 50000)
init.rate <- rexp(n = 1, rate = 0.5)
init.root <- rnorm(n = 1, mean = mu, sd = 10)
inf.25.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                          initial.root = init.root, prior.rate = prior.exp, 
                          prior.root = prior, w = 0.5, gen = 50000)

## 25 tips; uninformative priors:
phy <- rtree(25); phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.25.mcmc.1 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                          initial.root = init.root, prior.rate = prior.unif, 
                          prior.root = prior.unif, w = 0.5, gen = 50000)
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.25.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                          initial.root = init.root, prior.rate = prior.unif, 
                          prior.root = prior.unif, w = 0.5, gen = 50000)

## 10 tips; informative priors:
phy <- rtree(10); phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
mu <- mean(data)
## Precisamos ajustar a função de prior.norm.
prior <- function(x) dnorm(x, mean = mu, sd = 10, log = TRUE)
init.rate <- rexp(n = 1, rate = 0.5)
init.root <- rnorm(n = 1, mean = mu, sd = 10)
inf.10.mcmc.1 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                          initial.root = init.root, prior.rate = prior.exp, 
                          prior.root = prior, w = 0.5, gen = 50000)
init.rate <- rexp(n = 1, rate = 0.5)
init.root <- rnorm(n = 1, mean = mu, sd = 10)
inf.10.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                          initial.root = init.root, prior.rate = prior.exp, 
                          prior.root = prior, w = 0.5, gen = 50000)

## 10 tips; uninformative priors:
phy <- rtree(10); phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.10.mcmc.1 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                          initial.root = init.root, prior.rate = prior.unif, 
                          prior.root = prior.unif, w = 0.5, gen = 50000)
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.10.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                          initial.root = init.root, prior.rate = prior.unif, 
                          prior.root = prior.unif, w = 0.5, gen = 50000)

## Agora vamos salvar os resultados. Mas antes vou transformar em um objeto
##   pronto para usar o pacote 'coda'.
library(coda)

inf.10 <- mcmc.list( mcmc(inf.10.mcmc.1), mcmc(inf.10.mcmc.2) )
unif.10 <- mcmc.list( mcmc(unif.10.mcmc.1), mcmc(unif.10.mcmc.2) )
inf.25 <- mcmc.list( mcmc(inf.25.mcmc.1), mcmc(inf.25.mcmc.2) )
unif.25 <- mcmc.list( mcmc(unif.25.mcmc.1), mcmc(unif.25.mcmc.2) )
inf.50 <- mcmc.list( mcmc(inf.50.mcmc.1), mcmc(inf.50.mcmc.2) )
unif.50 <- mcmc.list( mcmc(unif.50.mcmc.1), mcmc(unif.50.mcmc.2) )
inf.100 <- mcmc.list( mcmc(info.100.mcmc.1), mcmc(info.100.mcmc.2) )
unif.100 <- mcmc.list( mcmc(unif.100.mcmc.1), mcmc(unif.100.mcmc.2) )

setwd("~/Documents/Academicos/Harmon_Lab/MCMC - Workshop - Brasil/workshop/Slides/Work_8")
save(inf.10, unif.10, inf.25, unif.25, inf.50, unif.50, inf.100, unif.100
     , file = "data_long_sim_prior.Rdata")

## Todas as simulações até agora tinham o prior centrado em um valor razoável para
##    os parâmetros. No caso do valor de root, a média do prior normal é bem próxima
##    do valor que gerou os dados.

## O que acontece se o prior não for uma distribuição razoável para os valores dos
##    parâmetros?
## Agora vamos determinar uma distribuição de priors que não é boa para nossos dados.
## A posterior deveria, portanto, deslocar-se do prior. Por que esse é o resultado
##    esperado?

## Repetindo as simulações com prior diferente.
prior.1 <- function(x){
  dnorm(x = x, mean = 10, sd = 2, log = TRUE)
}
prior.2 <- function(x, mu){
  dnorm(x = x, mean = 20, sd = 2, log = TRUE)
}

## Como são os novos priors?
hist(rnorm(n = 1000, mean = 10, sd = 2))
hist(rnorm(n = 1000, mean = 20, sd = 2))

## Os valores para a variância estão muito maiores do que o valor que gerou os dados.
## Da mesma forma para o valor da raiz. O prior está centrado no dobro do valor que
##    gerou os dados.

## 100 tips; informative priors:
phy <- rtree(100); phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
init.rate <- rnorm(n = 1, mean = 10, sd = 2)
init.root <- rnorm(n = 1, mean = 20, sd = 2)
unif.100.mcmc.1 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                           initial.root = init.root, prior.rate = prior.1, 
                           prior.root = prior.2, w = 0.5, gen = 10000)
init.rate <- rnorm(n = 1, mean = 10, sd = 2)
init.root <- rnorm(n = 1, mean = 20, sd = 2)
unif.100.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                           initial.root = init.root, prior.rate = prior.1, 
                           prior.root = prior.2, w = 0.5, gen = 10000)

## 100 tips; non-informative priors:
phy <- rtree(100); phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.100.mcmc.1 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                           initial.root = init.root, prior.rate = prior.unif, 
                           prior.root = prior.unif, w = 0.5, gen = 10000)
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.100.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                           initial.root = init.root, prior.rate = prior.unif, 
                           prior.root = prior.unif, w = 0.5, gen = 10000)

## 50 tips; informative priors:
phy <- rtree(50); phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
init.rate <- rnorm(n = 1, mean = 10, sd = 2)
init.root <- rnorm(n = 1, mean = 20, sd = 2)
inf.50.mcmc.1 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                         initial.root = init.root, prior.rate = prior.1, 
                         prior.root = prior.2, w = 0.5, gen = 50000)
init.rate <- rnorm(n = 1, mean = 10, sd = 2)
init.root <- rnorm(n = 1, mean = 20, sd = 2)
inf.50.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                         initial.root = init.root, prior.rate = prior.1, 
                         prior.root = prior.2, w = 0.5, gen = 50000)

## 50 tips; uninformative priors:
phy <- rtree(50); phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.50.mcmc.1 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                          initial.root = init.root, prior.rate = prior.unif, 
                          prior.root = prior.unif, w = 0.5, gen = 50000)
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.50.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                          initial.root = init.root, prior.rate = prior.unif, 
                          prior.root = prior.unif, w = 0.5, gen = 50000)

## 25 tips; informative priors:
phy <- rtree(25); phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
mu <- mean(data)
## Precisamos ajustar a função de prior.norm.
prior <- function(x) dnorm(x, mean = mu, sd = 10, log = TRUE)
init.rate <- rnorm(n = 1, mean = 10, sd = 2)
init.root <- rnorm(n = 1, mean = 20, sd = 2)
inf.25.mcmc.1 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                         initial.root = init.root, prior.rate = prior.1, 
                         prior = prior.2, w = 0.5, gen = 50000)
init.rate <- rnorm(n = 1, mean = 10, sd = 2)
init.root <- rnorm(n = 1, mean = 20, sd = 2)
inf.25.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                         initial.root = init.root, prior.rate = prior.1, 
                         prior = prior.2, w = 0.5, gen = 50000)

## 25 tips; uninformative priors:
phy <- rtree(25); phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.25.mcmc.1 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                          initial.root = init.root, prior.rate = prior.unif, 
                          prior.root = prior.unif, w = 0.5, gen = 50000)
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.25.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                          initial.root = init.root, prior.rate = prior.unif, 
                          prior.root = prior.unif, w = 0.5, gen = 50000)

## 10 tips; informative priors:
phy <- rtree(10); phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
init.rate <- rnorm(n = 1, mean = 10, sd = 2)
init.root <- rnorm(n = 1, mean = 20, sd = 2)
inf.10.mcmc.1 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                         initial.root = init.root, prior.rate = prior.1, 
                         prior = prior.2, w = 0.5, gen = 50000)
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
inf.10.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                         initial.root = init.root, prior.rate = prior.1, 
                         prior = prior.2, w = 0.5, gen = 50000)

## 10 tips; uninformative priors:
phy <- rtree(10); phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.10.mcmc.1 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                          initial.root = init.root, prior.rate = prior.unif, 
                          prior.root = prior.unif, w = 0.5, gen = 50000)
init.rate <- runif(n = 1, min = 0, max = 100)
init.root <- runif(n = 1, min = -100, max = 100)
unif.10.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                          initial.root = init.root, prior.rate = prior.unif, 
                          prior.root = prior.unif, w = 0.5, gen = 50000)

inf.10 <- mcmc.list( mcmc(inf.10.mcmc.1), mcmc(inf.10.mcmc.2) )
unif.10 <- mcmc.list( mcmc(unif.10.mcmc.1), mcmc(unif.10.mcmc.2) )
inf.25 <- mcmc.list( mcmc(inf.25.mcmc.1), mcmc(inf.25.mcmc.2) )
unif.25 <- mcmc.list( mcmc(unif.25.mcmc.1), mcmc(unif.25.mcmc.2) )
inf.50 <- mcmc.list( mcmc(inf.50.mcmc.1), mcmc(inf.50.mcmc.2) )
unif.50 <- mcmc.list( mcmc(unif.50.mcmc.1), mcmc(unif.50.mcmc.2) )
inf.100 <- mcmc.list( mcmc(info.100.mcmc.1), mcmc(info.100.mcmc.2) )
unif.100 <- mcmc.list( mcmc(unif.100.mcmc.1), mcmc(unif.100.mcmc.2) )

save(inf.10, unif.10, inf.25, unif.25, inf.50, unif.50, inf.100, unif.100
     , file = "data_long_sim_bad_prior.Rdata")