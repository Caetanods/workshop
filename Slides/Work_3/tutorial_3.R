## Nesse tutorial vamos explorar os modelos comparativos de Brownian Motion e
##   Ornstein-Uhlenbeck.
## Ambos são modelos de evolução de características contínuas em uma filogenia.
## Para trabalhar com os modelos de BM ou de OU precisamos de uma filogenia
##   com comprimento de ramos e valores contínuos para as espécies.

set.seed(1234)

## Para esse tutorial vamos usar os pacotes 'geiger' versão 2.0 e 'phytools'.
## Embora o pacote 'phytools' seja muito útil, o número de funções é muito
##   grande e muitas vezes os exemplos na pagina de ajuda estão ausentes.
## Para entender melhor as utilidades do pacote 'phytools', visite o blog do
##   Dr. Liam Revell. -- http://blog.phytools.org/

library(geiger)
library(phytools)

## Simulando a filogenia:
phy <- rtree(100)
plot(phy)

## Note que essa filogenia não tem todas as espécies alinhadas em uma mesma linha
##   no tempo. Ou seja, essa filogenia não é ultramétrica.
## A grande maioria dos modelos que vamos explorar assumem que a filogenia seja
##   ultramétrica. É comum ver artigos onde a filogenia foi transformada para ter
##   todos os comprimentos de ramos iguais a 1. Essa é uma estratégia quando não temos
##   acesso aos comprimentos de ramos da filogenia, como no caso de uma análise de
##   parcimônia.
## Abaixo vamos plotar essa mesma filogenia em sua forma ultramétrica e com todos os
##   comprimentos iguais a 1.

par(mfrow=c(1,2))
phy <- compute.brlen(phy)
parc.phy <- phy
parc.phy$edge.length <- rep(1, times=length(parc.phy$edge.length) )
plot(phy, main = "Ultrametric tree", show.tip.label=FALSE)
plot(parc.phy, main = "Parsimony - all 1", show.tip.label=FALSE)
par(mfrow=c(1,1))

## Quando utilizamos métodos comparativos estamos muitas vezes usando modelos que
##   assumem que o tempo de divergência é proporcional ao comprimento dos ramos.
## Portanto, utilizar uma filogenia não ultramétrica como a gerada pela parcimônia
##   pode trazer resultados incorretos, visto que o pressuposto do modelo foi
##   de certa forma violado.
## Nessa situação é boa prática utilizar funções como "compute.brlen" para gerar
##   diferentes comprimentos de ramos para a filogenia e testar se estes influenciam
##   nos resutados.
## Por hora vamos utilizar a filogenia ultramétrica para fazer simulações.

## Simulando dados usando um modelo de Brownian motion.

## O modelo BM possuí dois parâmetros, a taxa (sigma^2) e o valor da raíz (root value).
## A taxa descreve o quão rápido a variância acumula com o tempo enquanto que o
##   valor da raíz determina o valor de inicio do processo de BM.

## Veja o help da função 'sim.char' para entender como simular os dados e quais
##   parâmetros foram escolhidos para essa simulação.
help(sim.char)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]
hist(data)

## Vamos plotar o valor das espécies na filogenia.
contMap(tree = phy, x = data, type = "fan")
contMap(tree = phy, x = data)

## Podemos observar o sinal filogenético dos dados acompanhando a mudança de cores
##   na filogenia. Note que linhagens com ancestrais comuns mais recentes apresentam
##   valores da característica mais próxima dos que espécies mais distantes.

## Pergunta: O que acontece com o sinal filogenético se aumentarmos o valor da
##   taxa do modelo BM? O que acontece se reduzirmos o valor da taxa?
## Faça simulações para verificar os resultados.

## Veja como os valores evoluem seguindo a filogenia. Vamos fazer uma série de
##   simulações para deixar claro o componente estocástico do modelo BM.
par(mfrow=c(2,2))
bmPlot(tree = phy, anc = 10, sig2 = 0.5)
bmPlot(tree = phy, anc = 10, sig2 = 0.5)
bmPlot(tree = phy, anc = 10, sig2 = 0.5)
bmPlot(tree = phy, anc = 10, sig2 = 0.5)
par(mfrow = c(1,1))

## Função de verossimilhança para o modelo BM:
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

## Agora vamos usar a função de likelihood para gerar a superfície de likelihood
##  para o modelo BM. Para tal, como no tutorial anterior, vamos verificar a
##  likelihood de diversas combinação de valores para os parâmetros.
rep <- 20
## Note que o valor da taxa não pode ser negativo.
rate <- seq(0.1, 10, length.out = rep)
root <- seq(1, 20, length.out = rep)
z <- matrix(data = NA, nrow = rep, ncol = rep)
for(i in 1:rep){
  for(j in 1:rep){
    z[i,j] <- singleBML(phy = phy, x = data, sigsq = rate[i], mean = root[j])
  }  
}

persp(x=rate, y=root, z=z, theta = -30, phi = 20, xlab = "sigma^2", ylab = "phylo.mean", zlab = "log(likelihood)", border = "black", col = "grey", r = 4)

## Pergunta: Veja que a superfície de verossimilhança do modelo BM é muito semelhante à
##   do modelo normal. Porque você acha que isso ocorre?

## Veja que há uma região quase plana na superfície de likelihood. Isso acontece
##   pois existe uma interação entre o valor de raíz e a taxa. Diferentes valores de
##   taxa e raíz no modelo BM podem gerar valores de verossimilhança não tão
##   distintos como espera-se.

## Modelo Ornstein-Uhlenbeck (OU).

## O modelo OU adiciona um parâmetro para o valo do ótimo e uma força de seleção que
##    puxa o valor da característica para mais próximo deste ótimo.
## OU ainda possuí um componente de BM, por isso, a simulação de um processo.
##    em um modelo OU será parecida com o processo em um modelo BM.

## Vamos simular o modelo OU como fizemos com o modelo BM:
data.ou <- fastBM(tree = phy, a = 10, sig2 = 0.5, alpha = 0.4, theta = 20)
## Note que o valor de 'theta' é o ótimo do regime OU e que 'alpha' é a força
##    que puxa os valores para o 'theta'.
## O valor da raíz é 10.

## Vamos plotar o valor das espécies na filogenia.
contMap(tree = phy, x = data.ou)

## Essa simulação é mais facilmente visualizada utilizando um phenogram.
phenogram(tree = phy, x = data.ou)

## Pergunta: O que aconteceu com o valor do ancestral nesses plots? Porque o valor
##    do ancestral não é o mesmo usado na simulação?

## Use a função 'fastBM' e os plots para explorar o que acontece com os valores
##    dos tips e da raíz no modelo OU.

## Phylogenetic half-life:
## A meia vida filogenética. Essa quantidade é o tempo que o processo de OU leva
##   para chegar até a metade da distância do ótimo. Ou seja, se o half-life for metade
##   da "altura" (tree depth) da filogenia, então esperamos que a característica
##   chegue no ótimo.
## Em outras palavras, quanto maior o alpha, mais rápido os valores convergem para
##   o ótimo.
## O cálculo do half-life é 'log(2)/alpha'.
## Lembre-se de verificar o half-life para saber se o processo OU é rápido ou
##   devagar.
## No exemplo anterior temos log(2)/0.3 = 2.31. Visto que a arvore tem altura 1.
##   O processo deve estar chegando na metade do caminho nas pontas. Veja usando
##   os gráficos se isso faz sentido.

## Processo OU sem uma filogenia.
## Agora vamos simular um processo de OU sem uma filogenia. Desta forma vai ser
##   mais fácil observar o efeito de theta e alpha.

## A função abaixo simula um processo de OU.
simOU <- function(passos, sd, start, alpha, tetha){
  trace <- vector()
  trace[1] <- start
  for(i in 2:passos){
    trace[i] <- trace[i-1] + rnorm(n=1,mean=0,sd=sd) + (alpha * (tetha - trace[i-1]) )
  }
  return(trace)
}

## Primeiro simulamos o processo de OU.
OUtrace <- simOU(passos=300, sd=0.08, start=10, alpha=0.02, tetha=14)
## Agora vamos fazer um plot do valor do trait em função do tempo.
## A linha azul mostra o valor de theta e a linha vermelha mostra o valor
##   inicial.
plot(x=1:300, y=OUtrace, type="l", ylab="posição", xlab="tempo")
abline(h=10, col="red", lwd=2)
abline(h=14, col="blue", lwd=2)

## Exercício: Mude os valores de alpha, tetha e sd. Veja como o modelo se
##   comporta. Note que:
## alpha = força de seleção puxando para o valor de ótimo;
## tetha = valor do ótimo.
## sd = desvio padrão (variância) do processo de BM dentro do modelo OU.

## Planeje os valores que vai utilizar antes. Lembre-se que se mudar todos
##   os parâmetros de uma só vez não será capaz de entender a função de cada
##   parâmetro no modelo.

## Função de verossimilhança do modelo OU.

## Vamos primeiro usar o pacote 'bayou'. Esse é um pacote com várias funções Bayesianas
##   para estimar diversos modelos OU.
library(bayou)

## A função abaixo usa elementos do 'bayou' para calcular o log(likelihood) do modelo
##    OU.
ou.log.lik <- function(pars){
  ## pars[1] = alpha
  ## pars[2] = sigma^2
  ## pars[3] = theta
  cache <- bayou:::.prepare.ou.univariate(phy, data.ou, SE=0)
  pars.list <- list( alpha=pars[1], sig2=pars[2], k=1, ntheta=1, theta=pars[3], sb=numeric(0), loc=numeric(0), t2=numeric(0) )
  ## pars.list é uma lista de argumentos para rodar a função de likelihood.
  ## alpha = valor de alpha.
  ## sig2 = valor de sigma^2
  ## k = número de regimes OU. No nosso caso temos somente um regime.
  ## ntheta = número de thetas, temos somente 1 regime em nosso exemplo.
  ## Os demais parâmetros são detalhes da posição das transições entre os regimes
  ##     e portanto são definidos como 0.
  bayou.lik(pars.list, cache, data.ou)$loglik[1,1]
}

## Note que para usar a função precisamos somente de um vetor com exatamente
##   3 elementos. c(alpha, sigma^2, theta).
ou.log.lik(pars = c(0.3, 0.5, 20))

## Diferente dos modelos anteriores agora temos mais de 2 parâmetros. Como visualizar
##    a superfície de likelihood?
## Uma forma de resolver o problema é manter um dos parâmetros fixo e modificar os
##    demais.

## Primeiro, vamos fixar sigma^2 no valor que gerou os dados.
## Vamos calcular a log(likelihood) do modelo com uma série de combinações de
##    parâmetros.
rep <- 20
## Note que o valor de alpha não pode ser negativo.
sigma <- 0.5
alpha <- seq(0.01, 1, length.out = rep)
theta <- seq(1, 20, length.out = rep)
z <- matrix(data = NA, nrow = rep, ncol = rep)
for(i in 1:rep){
  for(j in 1:rep){
    ## z[i,j] <- singleBML(phy = phy, x = data, sigsq = rate[i], mean = root[j])
    z[i,j] <- ou.log.lik(pars = c(alpha[i], sigma, theta[j]))
  }  
}

persp(x=rate, y=root, z=z, theta = -30, phi = 20, xlab = "alpha", ylab = "theta", zlab = "log(likelihood)", border = "black", col = "grey", r = 4)

## Curioso para saber como é a superfície de likelihood quando usamos outra
##    combinação de parâmetros? Então faça você mesmo!
## Escolha um parâmetro para manter fixo e crie combinação dos dois restantes.
## Use o exemplo acima como guia.

## Podemos também utilizar o pacote 'ouch'. Que é um crássico.
library(ouch)

## Lembrando dos parãmetros que usamos para gerar os dados:
## data.ou <- fastBM(tree = phy, a = 10, sig2 = 0.5, alpha = 0.2, theta = 20)

run.hansen <- function(phy, dat, alpha, sigma, maxit){
  ## Função para rodar a MLE do modelo OU usando o pacote 'ouch'.
  ## A função faz uma série de transformações de formatos de objetos para rodar
  ##    a MLE.
  ## phy = filogenia
  ## dat = dados
  ## alpha = valor inicial para alpha
  ## sigma = valor inicial para sigma
  ## maxit = número máximo de iterações.
  ouch.tree <- ape2ouch(phy)
  dat <- as.data.frame(dat)
  dat$labels <- rownames(dat)
  ouch.tree.dat <- as(ouch.tree,"data.frame")
  ouch.tree.dat <- merge(ouch.tree.dat, dat , by="labels", all=TRUE)
  rownames(ouch.tree.dat) <- ouch.tree.dat$nodes
  ouch.tree.2 <- with(ouch.tree.dat, ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))
  ouch.tree.dat$regimes <- as.factor("global")
  h1 <- hansen(tree=ouch.tree.2, data=ouch.tree.dat["dat"], regimes=ouch.tree.dat["regimes"], sqrt.alpha=alpha, sigma=sigma, maxit=maxit)
  return(h1)
}

## MLE usando a função 'hansen' do pacote 'ouch'.
run.hansen(phy=phy, dat=data.ou, alpha=1, sigma=1, maxit=10000)

## A estimativa de máxima verossimilhança usando a função do ouch faz sentido?
## Por que o valor de theta não é igual o valor da raíz e nem igual ao valor do ótimo
##    usado na simulação? O que esses resultados tem a ver com a phylogenetic
##    half-life (meia vida filogenética)?
