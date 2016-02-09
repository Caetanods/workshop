## Nesse tutorial vamos explorar a utilização da função de verossimilhança.
## Vamos usar o modelo de distribuição Normal.

## Seed para o script, isso permite que todos os resultados sejam equivalentes.
set.seed(1234)

## Gerar uma distribuição normal e o plot.
nr <- rnorm(n = 2000, mean = 10, sd = 2)
hist(nr, main="", xlab="Data", col = "gray", border = "white")

## Função de verossimilhança para a distribuição normal.
## A função abaixo calcula a verossimilhança de um modelo normal.
lik.norm <- function(dt, mu, var){
  ## Calcula a verossimilhança de um modelo normal.
  ## dt = dados (class::vector).
  ## mu = média.
  ## var = variância.
  n <- length(dt)
  first <- (2*pi*var)^(-n/2)
  second <- (-1/2*var) * ( sum(dt-mu) )^2
  lik <- first * exp( second )
  return(lik)
}

## Calcule a verossimilhança de um par de valores para os parãmetros do modelo.
lik.norm(dt=nr, mu=11, var=4)
lik.norm(dt=nr, mu=10, var=4)
lik.norm(dt=nr, mu=10, var=5)

## Note que o valor de verossimilhança é muito pequeno (as vezes igual a 0). 
## Isso acontece pois a exponenciação com valor negativo dependente do número
## de amostra faz com que a quantidade calculada seja muito pequena.

## Vamos diminuir o número da amostra e ver qual o efeito sobre o cálculo da
##   verossimilhança.
nr.small <- rnorm(n = 10, mean = 10, sd = 2)
lik.norm(dt=nr.small, mu=10, var=5)

## Veja que agora conseguimos um valor de verossimilhança muito maior.
## Utilizar a função de likelihood sem transformar para o log pode provocar erros
##   nos cálculos, principalmente pois teremos que avaliar a função muitas vezes
##   quando implementamos uma análise usando MCMC.

## A função abaixo usa o log para transformar o cálculo da verossimilhança.
log.lik.norm <- function(dt, mu, var){
  ## dt = dados (class::vector).
  ## mu = média.
  ## var = variancia.
  n <- length(dt)
  first <- log(2)*(-n/2) + log(pi)*(-n/2) + log(var)*(-n/2)
  second <- -1/(2*var) * sum( (dt - mu)^2 )
  lik <- first + second
  return(lik)
}

## Vamos comparar os resultados produzidos pela função de likelihood em log e
##    sem a transformação:
lik.norm(dt=nr, mu=11, var=4)
log.lik.norm(dt=nr, mu=11, var=4)
lik.norm(dt=nr, mu=10, var=4)
log.lik.norm(dt=nr, mu=10, var=4)
lik.norm(dt=nr, mu=10, var=5)
log.lik.norm(dt=nr, mu=10, var=5)

## Note que os valores produzidos pela função em espaço log são negativos, no entanto,
##   o valor absoluto é muito maior. Cálculos com muitas casas decimais causam
##   problemas.
## Esse é o motivo porque muitas vezes a likelihood é mostrada em valores negativos.
## Pois na verdade está sendo mostrado os valores de log likelihood.

## Podemos usar a função de loglik para traçar a superfície de verossimilhança.
## O gráfico abaixo mostra somente o perfil de verossimilhança para a média
##     assumindo desvio padrão (sigma) = 2.
res_mu <- vector()
for(i in 1:20){
  res_mu[i] <- log.lik.norm(dt=nr, mu=i, var=4)
}
plot(res_mu~seq(1,20), xlab = "mean (mu)", ylab = "log(Likelihood)", pch = 16, type = "b")

## Podemos fixar o valor de mu e traçar o perfil de verossimilhança para sigma.
## O gráfico abaixo mostra somente o perfil de verossimilhança para a variância
##     assumindo média (mu) = 10.
res_sig <- vector()
for(i in 1:20){
  res_sig[i] <- log.lik.norm(dt=nr, mu=10, var=i)
}
plot(res_sig~seq(1,20), xlab = "variance (sigma^2)", ylab = "log(Likelihood)", pch = 16, type = "b")

## Agora vamos calcular a superfície completa de likelihood.
## Nesse caso temos dois parâmetros e, portanto, precisamos adicionar
##    mais um eixo no gráfico. O que vai tornar o gráfico 3D.
## Para tal vamos utilizar a função 'persp'.

x <- 1:20 ## valores para mu
y <- 1:20 ## valores para sigma
## Vamos usar z para criar uma tabela com a combinação dos valores de
##    log likelihood.
## Os valores de x estão nas linhas e os de y nas colunas.
z <- matrix(data = NA, nrow = 20, ncol = 20)
for(i in 1:20){
  for(j in 1:20){
    z[i,j] <- log.lik.norm(dt=nr, mu=x[i], var=y[j])
  }
}

persp(x, y, z, theta = 35, phi = 35, xlab = "mu", ylab = "sigma^2", zlab = "log(likelihood)", border = "black", col = "grey", r = 4)

## Para fazer esse plot também podemos usar o pacote 'lattice'.

library(lattice)
X <- as.matrix( expand.grid(x, y))
colnames(X) <- c("mu","var")
Z <- vector()
for(i in 1:nrow(X)){
  Z[i] <- log.lik.norm(dt=nr, mu=X[i,1], var=X[i,2])
}
df <- data.frame(X, Z)

wireframe(Z ~ mu*var, data=df, main="", shade=TRUE
          , screen = list(z = -50, x = -70) )

## Agora que temos uma visão da superfície de verossimilhança do modelo,
##   vamos fazer uma estimativa da Maximum Likelihood Estimate usando uma
##   função heurística de busca e depois as estimativas diretas.

## Estimativa heurística:
## Para tal vamos usar a função 'optim'.
## Essa é uma função geral e importante para fazer estimativas no R.
## Veja com calma o help dessa função, tente entender qual seu objetivo e como
##   ela opera.
help(optim)

## Note que precisamos criar uma nova função, pois 'optim' pede que o argumento 'par'
##   seja somente um vetor com os valores iniciais dos parâmetros da função.
## Ou seja, os dados precisam ser incluídos diretamente na função que será passada
##   para a função 'optim'.
to.optim <- function(par){
  log.lik.norm(dt=nr, mu=par[1], var=par[2])
}

## Aqui vamos escolher valores iniciais bem longe dos valores que geraram a
##   a distribuição. Vamos falar mais sobre valores iniciais nos próximos
##   tópicos.
## Precisamos mudar o valor de 'control$fnscale' para fazer uma maximização
##   ao invés de minimização usando 'optim'.
optim(c(1,1), to.optim, control=list(fnscale=-1))

## Agora podemos calcular o valor de MLE para a média e variância usando
##   as fórmulas derivadas.
mean(nr)
var(nr)

## Veja que existe uma diferença entre as duas estimativas.

## Essa diferença é importante?
## Qual os possíveis motivos da diferença?
## Por que nenhuma das duas estimativas converge com o valor que gerou os dados?

###########################################################
## Desafio:

## Use a função de likelihood do modelo linear (regressão) para explorar a
##    superfície de verossimilhança e estimar os melhores valores para A e B em
##    A*x + B = y

## Nessa regressão vamos usar uma distribuição normal para o erro, o erro nada
##    mais é do que o quão espalhado os pontos são da reta. Quando ajustamos
##    uma reta em uma regressão os nossos dados não se alinham completamente na
##    linha, essa distãncia das observações para a linha de regressão é explicada
##    por uma distribuição normal.

## Abaixo vou construir a função de verossimilhança para esse modelo.
## Note que a função é uma soma da log(likelihood) de uma distribuição normal
##    com média 'A * x + B' e desvio padrão 'sd' para cada ponto.
lm.log.lik <- function(a, b, sd){
  ## a = a term.
  ## b = b term.
  ## sd = standard deviation.
  pred <- a*x + b
  single.log.liks <- dnorm(y, mean = pred, sd = sd, log = T)
  sum(single.log.liks)
}

## Gerando os dados:
true.a <- 5
true.b <- 0
true.sd <- 10
sample <- 31

# create independent x-values 
x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
# create dependent values according to ax + b + N(0,sd)
y <-  trueA * x + trueB + rnorm(n=sampleSize,mean=0,sd=trueSd)
plot(x,y, main="Test Data")

## Dica: Copie e modifique os passos do tutorial acima.