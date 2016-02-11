## Na aula falamos mais sobre as distribuições de proposta e adicionamos a razão
##    de Hastings no algoritmo de MCMC, derivando o algoritmo de Metroplis-Hastings.

## Esse tutorial será a junção de muitos dos conceitos que vimos até o momento.
## Vamos desenvolver e aplicar uma função completa de MCMC que não deixa nada a desejar
##    a qualquer programa de análise Bayesiana que iria fazer tal análise.
## Contrário do que foi mostrado no tutorial anterior, vamos montar primeiro uma função
##    que vai rodar a cadeia de MCMC e então verificar os resultados.
## Se precisar entender o que está acontecendo em algum ponto em específico dessa
##    função de análise, retorne ao 'tutorial_5.R'.

## O DESAFIO deste tutorial é fazer uma análise de MCMC para o modelo OU.

## Vamos definir o sliding window:
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

## Pronto. Agora podemos testar nossa função.

## Primeiro vamos gerar os dados e filogenia.

set.seed(1234)
library(geiger)

phy <- rtree(100)
phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]

## Agora vamos criar as funções para o prior:
## Aqui usando os mesmos priors do tutorial anterior. Fique a vontade para colocar
##    seus próprios priors.
prior.rate <- function(x){
  dexp(x = x, rate = 0.5, log = TRUE)
}
mean.data <- mean(data)
prior.root <- function(x){
  dnorm(x = x, mean = mean.data, sd = 10, log = TRUE)
}

## Determinar valores iniciais. Desta vez escolhidos arbitrariamente.
init.rate <- 0.1
init.root <- 0

## Rodar a análise:
res.mcmc <- mcmc.BM(phy = phy, data = data, initial.rate = init.rate,
                    initial.root = init.root, prior.rate = prior.rate, 
                    prior.root = prior.root, w = 0.5, gen = 10000)

## Rode as linhas abaixo para salvar ou ler resultados prontos.
## Importe os resultados prontos caso o seu MCMC esteja demorando muito.
# saveRDS(res.mcmc, file = "mcmc.BM.1.rds")
# res.mcmc <- readRDS(file = "mcmc.BM.1.rds")


## O trace plot para a likelihood.
plot(1:10000, res.mcmc$log_lik, type = "l")
## Veja como o valor de likelihood sobe rapidamente e estabiliza. Note também que o
##   valor inicial do parâmetro era terrível!!

## O trace plot para os parâmetros.
plot(1:10000, res.mcmc$root, type = "l")
plot(1:10000, res.mcmc$rate, type = "l")
## Por que as primeiras gerações mostram um padrão diferente do resto?
## Marque no plot o valor inicial da cadeia e o valor de simulação ("valor real")
##    de cada parãmetro.

## Agora vamos plotar a densidade das distribuições posteriores.
## Também vamos plotar os priors, para comparação.
hist(x = rnorm(n = 10000, mean = mean.data, sd = 10), border = "blue", 
     main = "Posterior for root")
hist(x = res.mcmc$root, add = TRUE, col = "red")

hist(x = rexp(n = 10000, rate = 0.5), border = "blue", main = "Posterior for rate")
hist(x = res.mcmc$rate, add = TRUE, col = "red")

## Note como a distribuição posterior é distinta do prior. Mesmo o prior para o
##    valor do root ser bastante informativo, podemos ver que a informação contida
##    nos dados transformou notavelmente a distribuição do prior para a da
##    posterior.

## O plot abaixo mostra o trajeto feito pelo MCMC.
## A região onde as linhas se concentram é a distribuição posterior.
## O trajeto para chegar lá, ou seja, esse caminho que a cadeia fez até formar
##   esse "bolo" é chamado de 'burn-in'. Burn-in nada mais é do que as primeiras
##   gerações do MCMC quando a likelihood está subindo rapidamente e os parâmetros
##   mudando.
plot(res.mcmc$rate, res.mcmc$root, type = "l", xlab="Rate chain", ylab="Root chain",
     col="red")

## Vamos rodar mais duas cadeias com diferentes valores iniciais para mostrar a 
##   convergência. Convergência é quando diferentes cadeias começando de locais
##   distintos chegam em um mesmo resultado.

## Cadeia 2:
res.mcmc.2 <- mcmc.BM(phy = phy, data = data, initial.rate = 10,
                    initial.root = 8, prior.rate = prior.rate, 
                    prior.root = prior.root, w = 0.5, gen = 10000)
# saveRDS(res.mcmc.2, file="mcmc.BM.2.rds")
# res.mcmc.2 <- readRDS(file="mcmc.BM.2.rds")

## Cadeia 3:
res.mcmc.3 <- mcmc.BM(phy = phy, data = data, initial.rate = 1.5,
                    initial.root = 20, prior.rate = prior.rate, 
                    prior.root = prior.root, w = 0.5, gen = 10000)
# saveRDS(res.mcmc.3, file="mcmc.BM.3.rds")
# res.mcmc.3 <- readRDS(file="mcmc.BM.3.rds")

## Rodamos a mesma análise de MCMC três vezes. Em cada uma delas utilizamos diferentes
##   valores iniciais escolhidos para estarem espalhados na superfície de likelihood.
range.rate <- range( c(res.mcmc$rate, res.mcmc.2$rate, res.mcmc.3$rate) )
range.root <- range( c(res.mcmc$root, res.mcmc.2$root, res.mcmc.3$root) )
plot(res.mcmc$rate, res.mcmc$root, type = "l", xlab="Rate chain", ylab="Root chain",
     col="red", xlim=range.rate, ylim=range.root)
points(res.mcmc.2$rate, res.mcmc.2$root, type = "l", col="green")
points(res.mcmc.3$rate, res.mcmc.3$root, type = "l", col="blue")

## Que conclusão você pode tirar desse plot?

## No próximo tutorial vamos usar essas cadeias para explorar análises da posterior.
##    Portanto, lembre-se de salvar seus resultados.

## DESAFIO: Agora use a função de MCMC acima para fazer a inferência de dados
##    gerados sob um modelo de Ornstein–Uhlenbeck (OU).