## Nessa unidade vimos a fórmula de Bayes e a odds ratio.
## A odds ratio vai ser um conceito fundamental do MCMC. Vamos usar essa razão para
##   comparar diferentes valores de parâmetros para os modelos.
## Nesse tutorial vamos utilizar o modelo de Brownian Motion para explorar como
##   utilizar a odds ratio para tomar decisões sobre valores de parâmetros.

set.seed(1234)

## Função de verossimilhança do modelo BM:
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

## Simulando dados. Filogenia e traits.
library(geiger)
phy <- rtree(100)
phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]

## Para calcularmos a 'odds ratio' precisamos de uma distribuição a priori (ou
##   prior distribution). A prior distribution é a distribuição de probabilidades
##   dos parâmetros independente dos dados. Por isso é chamada de distribuição a priori
##   pois sua forma não muda em função dos dados.

## Uma distribuição a priori informativa tem densidade maior para certos valores de
##   parâmetros e menor para outros. Veja os casos a seguir:
prior_1 <- rexp(n = 1000, rate = 0.5)
prior_2 <- rlnorm(n = 1000, meanlog = 0.5, sdlog = 0.25)
par(mfrow=c(1,2))
hist(prior_1, main = "Exponential prior")
abline(v=mean(prior_1), col="red")
hist(prior_2, main = "Log normal prior")
abline(v=mean(prior_2), col="red")
par(mfrow=c(1,1))

## O prior exponencial e log normal são reportados com frequência. Ambos são uma boa
##   escolha quando trabalhamos com parâmetros [0,+inf) como a taxa de um modelo BM.

## Uma distribuição a priori não-informativa é geralmente uma distribuição uniforme.
## No entanto, não é possível, na prática, determinar um prior completamente
##   não-informativo. Veja que a distribuição uniforme tem um valor de mínimo e outro
##   de máximo. Contrário das distribuições exponencial e log normal (ou normal) a
##   uniforme tem um 'hard bound' no máximo e no mínimo. Ou seja, se o valor sair
##   fora do prior uniforme pode criar um problema para o cálculo da odds ratio.
## Quando utilizamos uma distribuição uniforme como prior deve-se tomar cuidado com
##   os valores de máximo e de mínimo.

prior_3 <- runif(n = 1000, min = 0, max = 100)
## 100 me parece grande o bastante... Concorda?
hist(x = prior_3, main = "Uniform prior")

## Vamos usar ness exemplo a distribuição de prior exponencial.

## Como calcular P(theta)? Usamos a função 'dexp'.
dexp(x = 0:10, rate = 0.5, log = TRUE)
## Veja que a log(probabilidade) é inversamente proporcional ao valor nessa distribuição.
## Ou seja, estamos favorecendo os menores valores para a taxa do BM.

## Cálculo do odds ratio tem dois termos (lembrando que trabalhamos com o log):
## Razão do prior.
## Razão da likelihood.

## EXEMPLO:
## Odds ratio entre:
## Prop_1: mean = 10, sigma^2 = 0.1
## Prop_2: mean = 10, sigma^2 = 0.3
## odds ratio = log(Prop_1) - log(Prop_2)

lik.1 <- singleBML(phy = phy, x = data, sigsq = 0.1, mean = 10)
lik.2 <- singleBML(phy = phy, x = data, sigsq = 0.3, mean = 10)
prior.1 <- dexp(x = 0.1, rate = 0.5, log = TRUE)
prior.2 <- dexp(x = 0.3, rate = 0.5, log = TRUE)
odds <- (prior.1 - prior.2) + (lik.1 - lik.2)
odds

## Qual é o melhor vetor de parâmetros? Prop_1 ou Prop_2?

## Lembre-se que estamos trabalhando com log. Portanto:
(prior.1 - prior.1) + (lik.1 - lik.1)
## Ou seja, valores iguais de parâmetros tem odds ratio igual a 0.

## Teste a diferença de maiores valores de parâmetros. Tente construir uma função para
##   automatizar a comparação entre conjunto de parâmetros.

## Lembre-se que em um MCMC temos que calcular a odds ratio muitas vezes.

## Tente utilizar diferentes distribuições de priors, como o prior uniforme.
## Você nota alguma influência do prior nos resultados da 'odds ratio'?