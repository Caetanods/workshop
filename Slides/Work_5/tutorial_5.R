## Nesse tutorial vamos contruir passo a passo uma busca Bayesiana usando o algoritmo de
##    Metropolis.
## O fluxograma mostrado em aula pode ser uma boa referência para acompanhar cada um
##   dos passos.

## Nesse tutorial vamos aplicar o MCMC em um modelo de Brownian Motion. Podemos
##   facilmente mudar a função de verossimilhança e ajustar os parâmetros, priors e
##   distribuições de propostas para implementar uma busca de MCMC de outro modelo.

## Vou aproveitar para fazer comentários sobre a implementação da análise. Espero que
##   ajude a implementar seu próprio MCMC no futuro.

## Dica 1: Implemente cada passo do MCMC em funções separadas, dessa forma é mais fácil
##   modificar e corrigir erros.

## Distribuição de proposta:
sliding.window <- function(x, w){
  ## Sliding window proposal for unbounded trait.
  ## x = the current value.
  ## w = the width parameter of the proposal.
  y <- runif(1, min = x - (w/2), max = x + (w/2) )
  return(y)
}

## Função de verossimilhança:
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

## Por hora vamos definir somente estas duas funções. A distribuição de proposta
##   e a função de likelihood para o modelo BM.

set.seed(1234)

## Vamos novamente simular os dados e a filogenia:
library(geiger)
phy <- rtree(100)
phy <- compute.brlen(phy)
data <- sim.char(phy = phy, par = 0.5, nsim = 1, model = "BM", root = 10)[,,1]

## Distribuições a priori:
## Vou definir um prior exponencial para a taxa e uma normal para o valor da raíz.
## Antes mesmo de observar os dados eu sei que um bom chute para o valor da raíz é
##    a média dos valores das espécies, mas vamos deixar uma grande variância. Eu
##    sei também que a taxa é geralmente um número pequeno, muitas vezes entre 0 e 1.
prior.rate <- function(x){
  dexp(x = x, rate = 0.5, log = TRUE)
}
mean.data <- mean(data)
prior.root <- function(x){
  dnorm(x = x, mean = mean.data, sd = 10, log = TRUE)
}

## Valores iniciais.
## Vamos usar um valor inicial para os parâmetros tirado da distribuição dos priors.
initial.rate <- rexp(n = 1, rate = 0.5)
initial.root <- rnorm(n=1, mean = mean.data, sd = 10)

## Análise de MCMC usando o algoritmo de Metropolis:

## Valores iniciais da cadeia.
initial.rate; initial.root

## Update do parâmetro:
## Qual parâmetro vamos propor um novo valor? Podemos atualizar os dois ou um por vez.
## Vamos atualizar um parâmetro por vez. Vamos escolher um dos parãmetros. Para tal,
##    vamos sortear usando uma probabilidade simétrica.
par.update <- sample(x = c("rate","root"), size = 1, prob = c(0.5,0.5))

## Agora um if / else para fazer a proposta de cada um dos parâmetros:
if(par.update == "root"){
  prop.root <- sliding.window(x = initial.root, w = 0.5)
} else {
  prop.rate <- sliding.window(x = initial.rate, w = 0.5)
}

## Agora temos um valor de proposta para um dos parâmetros.
## Vamos então aceitar ou rejeitar. Novamente vamos usar um if / else para trabalhar
##   com o parâmetro correto.
if(par.update == "root"){
  ## Primeiro a odds ratio.
  lik.initial.root <- singleBML(phy=phy, x=data, sigsq=initial.rate, mean=initial.root)
  lik.prop.root <- singleBML(phy=phy, x=data, sigsq=initial.rate, mean=prop.root)
  prior.initial <- prior.rate(initial.rate) + prior.root(initial.root)
  prior.prop <- prior.rate(initial.rate) + prior.root(prop.root)
  ll <- lik.prop.root - lik.initial.root
  pp <- prior.prop - prior.initial
  r.root <- ll + pp ## ratio of log(lik) + ratio of log(prior)
} else {
  lik.initial.rate <- singleBML(phy=phy, x=data, sigsq=initial.rate, mean=initial.root)
  lik.prop.rate <- singleBML(phy=phy, x=data, sigsq=prop.rate, mean=initial.root)
  prior.initial <- prior.rate(initial.rate) + prior.root(initial.root)
  prior.prop <- prior.rate(prop.rate) + prior.root(initial.root)
  ll <- lik.prop.rate - lik.initial.rate
  pp <- prior.prop - prior.initial
  r.rate <- ll + pp ## ratio of log(lik) + ratio of log(prior)
}

## Agora temo a odds ratio. Com ela podemos fazer o passo de aceitar ou rejeitar.
## Novamente vamos fazer um if / else para ajustar de acordo com o parâmetro.

if(par.update == "root"){
  
  if(exp(r.root) > runif(1)){
    ## Aceite.
    ## Note que 'exp(r.root)' pois r foi calculado na escala log.
    print("Novo valor de 'root' aceito! :)")
    next.root <- prop.root
    next.rate <- initial.rate
  } else{
    print("Novo valor de 'root' rejeitado! :(")
    next.root <- initial.root
    next.rate <- initial.rate
  }
  
} else {
  
  if(exp(r.rate) > runif(1)){
    ## Aceite.
    ## Note que 'exp(r.rate)' pois r foi calculado na escala log.
    print("Novo valor de 'rate' aceito! :)")
    next.rate <- prop.rate
    next.root <- initial.root
  } else{
    print("Novo valor de 'rate' rejeitado! :(")
    next.rate <- initial.rate
    next.root <- initial.root
  }
  
}

## Vamos revisar o resultado dessa primeira geração do MCMC:

## Valores iniciais:
initial.root; initial.rate

## Proposta e valores propostos:
par.update
ifelse(test = par.update == "root", yes = paste("root", prop.root), no = paste("rate", prop.rate) )

## Próximo estado da cadeia, após aceitar ou rejeitar a proposta:
next.root; next.rate

## Esse foi somente um passo da cadeia de MCMC usando o algoritmo de Metropolis.
## Certifique-se de que entendeu cada um dos passos acima.
## Relacione o que foi feito aqui com o fluxograma mostrado na aula.

## Note que o passo de aceitar ou rejeitar não é exatamente igual ao que foi mostrado
##   na aula. Por quê?
## A forma que foi implementado aqui (em pseudo-código):
# if(exp(r) > runif(1)){
# # Aceite.
# { else{
# # Rejeite.
# }

## Qual a relação desta implementação com o que foi mostrado em aula?

## No próximo tutorial vamos construir uma função completa que fará a análise de
##    MCMC. Essa função terá os mesmo passos que foram mostrados acima. As diferenças
##    serão somente uma questão de implementação, visto que o código acima não
##    está preparado para rodar milhares de vezes, como é necessário no MCMC.

## Como você faria para rodar esses passos milhares de vezes?

## Quais informações (objetos/valores) devemos salvar ao longo das gerações?

## Sabemos que variância não aceita valores negativos. Existe algum problema com esse
##   esquema de propostas de novos valores? Qual é o problema?

## DESAFIO 1: Tente inventar uma maneira de não gerar ou não aceitar valores
##   negativos de rate.
## DESAFIO 2: Use o exemplo acima para fazer mais duas gerações de MCMC. Guarde os
##    resultados de cada uma das gerações.