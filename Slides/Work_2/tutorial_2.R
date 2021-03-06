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
  ## 'first' e 'second' estao somente para ajudar a ler a funcao.
  lik <- first * exp( second )
  return(lik)
}

## Calcule a verossimilhança de um par de valores para os parãmetros do modelo.
lik.norm(dt=nr, mu=11, var=4)
lik.norm(dt=nr, mu=10, var=4)
lik.norm(dt=nr, mu=10, var=5)

## Note que o valor de verossimilhança é muito pequeno (as vezes igual a 0).
## Isso acontece pois a exponenciação com valor negativo é dependente do número
## de amostra e faz com que a quantidade calculada seja muito pequena.
## Esse efeito já acontece com n = 2000!

## Vamos diminuir o número da amostra e ver qual o efeito sobre o cálculo da
##   verossimilhança. Menor valor de n deve facilitar o calculo da verossimilhança.
nr.small <- rnorm(n = 10, mean = 10, sd = 2)

lik.norm(dt=nr.small, mu=11, var=4)
lik.norm(dt=nr.small, mu=10, var=4)
lik.norm(dt=nr.small, mu=10, var=5)

## Veja que agora conseguimos um valor de verossimilhança diferente de 0.
## Utilizar a função de likelihood sem transformar para o log pode provocar erros
##   nos cálculos.
## Vamos ver que durante o MCMC temos que avaliar a função de verossimilhança com valores
##   de parametros "bem ruins". Ou seja, a verossimilhança do modelo é bem pequena. Nesse caso
##   seria praticamente impossivel proceder com a função de verossimilhança sem efetuar a
##   transformação para o log.

## Vamos repetir usando a transformação log:

## A função abaixo usa o log para transformar o cálculo da verossimilhança.
## Não basta simplesmente calcular o log da função anterior, pois o problema de calculo já acontece
##     durante as operações da função de verossimilhança.
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
##   o valor absoluto é muito maior!
## Cálculos com muitas casas decimais causam problemas para o computador.
## Esse é o motivo porque muitas vezes a likelihood é mostrada em valores negativos, pois
##   na verdade está sendo mostrado os valores de log likelihood.

## Podemos usar a função de loglik para traçar a superfície de verossimilhança.
## O gráfico abaixo mostra somente o perfil de verossimilhança para a média
##   assumindo desvio padrão (sigma^2) = 2 fixo. Em outras palavras, o gráfico
##   mostra o perfil de verossimilhança para a média condicionado em sigma^2 = 2.

res_mu <- vector()
for(i in 1:20){
  res_mu[i] <- log.lik.norm(dt=nr, mu=i, var=4)
}
plot(res_mu~seq(1,20), xlab = "mean (mu)", ylab = "log(Likelihood)", pch = 16, type = "b")

## Podemos, da mesma forma, fixar o valor de mu (média) e traçar o perfil de
##     verossimilhança para sigma^2 (variância).
## O gráfico abaixo mostra somente o perfil de verossimilhança para a variância
##     assumindo média (mu) = 10.

res_sig <- vector()
for(i in 1:20){
  res_sig[i] <- log.lik.norm(dt=nr, mu=10, var=i)
}
plot(res_sig~seq(1,20), xlab = "variance (sigma^2)", ylab = "log(Likelihood)", pch = 16, type = "b")

## Podemos comparar esses dois perfis de verossimilhança.
## Note como o efeito da média na verossimilhança é maior do que o efeito da variância.
## Aqui já temos um conceito importante: alguns parametros dos modelos são mais fáceis
##      de otimizar (i.e., estimar) do que outros.
par( mfrow = c(1,2) )
plot(res_mu~seq(1,20), xlab = "mean (mu)", ylab = "log(Likelihood)", pch = 16, type = "b"
     , ylim = c(-30000, -4000))
plot(res_sig~seq(1,20), xlab = "variance (sigma^2)", ylab = "log(Likelihood)", pch = 16, type = "b"
     , ylim = c(-30000, -4000))

## Agora vamos calcular a superfície completa de likelihood.
## Nesse caso temos dois parâmetros e, portanto, precisamos adicionar mais um eixo no gráfico.
## O que vai tornar o gráfico 3D. Para tal vamos utilizar a função 'persp'.

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

## A matrix 'z' guarda valores de log-likelihood para cada uma das 20x20 combinações de média e variãncia.
head(z)

## Agora o gráfico:
persp(x, y, z, theta = 35, phi = 35, xlab = "mu", ylab = "sigma^2", zlab = "log(likelihood)", border = "black", col = "grey", r = 4)
## Mude os valores de 'theta' e 'phi' para explorar diferentes perspectivas.
persp(x, y, z, theta = 20, phi = 15, xlab = "mu", ylab = "sigma^2", zlab = "log(likelihood)", border = "black", col = "grey", r = 4)
persp(x, y, z, theta = 200, phi = 15, xlab = "mu", ylab = "sigma^2", zlab = "log(likelihood)", border = "black", col = "grey", r = 4)

## Podemos também expandir a area desse perfil de likelihood.
x <- -40:40 ## valores para mu
y <- seq(from = 0.0, to = 10, by = 0.1) ## valores para sigma
## Vamos usar z para criar uma tabela com a combinação dos valores de
##    log likelihood.
## Os valores de x estão nas linhas e os de y nas colunas.
z <- matrix(data = NA, nrow = length(x), ncol = length(y))
for(i in 1:length(x)){
  for(j in 1:length(y)){
    z[i,j] <- log.lik.norm(dt=nr, mu=x[i], var=y[j])
  }
}
head( z )
persp(x, y, z, theta = 35, phi = 35, xlab = "mu", ylab = "sigma^2", zlab = "log(likelihood)"
      , border = "black", col = "grey", r = 4)
persp(x, y, z, theta = 20, phi = 15, xlab = "mu", ylab = "sigma^2", zlab = "log(likelihood)", border = "black", col = "grey", r = 4)
persp(x, y, z, theta = 200, phi = 15, xlab = "mu", ylab = "sigma^2", zlab = "log(likelihood)", border = "black", col = "grey", r = 4)

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
## (O que você vai fazer agora é muito similar com o que uma série de pacotes em R fazem
##   para encontrar o Maximum Likelihood Estimate de modelos!)
to.optim <- function(par){
  log.lik.norm(dt=nr, mu=par[1], var=par[2])
}

## Aqui vamos escolher valores iniciais bem longe dos valores que geraram a
##   a distribuição. Vamos falar mais sobre valores iniciais nos próximos
##   tópicos.
## Precisamos mudar o valor de 'control$fnscale' para fazer uma maximização
##   ao invés de minimização usando 'optim'.
## Faz sentido para você que precisamos fazer uma maximização dessa função? Esse conceito
##   é importante. Não deixe a dúvida passar.
model.fit <- optim(c(1,1), to.optim, control=list(fnscale=-1))
model.fit
names( model.fit )
model.fit$par ## O MLE dos parametros. Primero média e depois variancia.
model.fit$value ## O log-likelihood para o MLE.
model.fit$convergence ## Código para a convergência. Veja a pagina de ajuda de 'optim'.

## Note the a otimização anterior teve uma mensagem de Warning. Porque isso?
## Veja opção "lower" e "upper"
bounded.fit <- optim(par = c(1,1), fn = to.optim, method = "L-BFGS-B"
                     , lower = c(-100, 0.0000001), upper = c(100, 100)
                     , control = list(fnscale=-1))
model.fit$par ## O MLE dos parametros. Primero média e depois variancia.
model.fit$value ## O log-likelihood para o MLE.
model.fit$convergence ## Código para a convergência. Veja a pagina de ajuda de 'optim'.

bounded.fit <- optim(par = c(90,90), fn = to.optim, method = "L-BFGS-B"
                     , lower = c(-100, 0.0000001), upper = c(100, 100)
                     , control = list(fnscale=-1))
bounded.fit$par

## Agora podemos calcular o valor de MLE para a média e variância usando
##   as fórmulas analíticas.
mean(nr)
var(nr)

## Valores usados para gerar os dados:
abs( model.fit$par[1] - 10 )
abs( model.fit$par[2] - 4 ) ## SD = 2
abs( mean(nr) - 10 )
abs( var(nr) - 4 ) ## SD = 2

## Veja que existe uma diferença pequena entre as duas estimativas.

## Responda as perguntas:

## Essa diferença parece ser importante?
## Qual os possíveis motivos da diferença?
## Por que nenhuma das duas estimativas resultou no mesmo valor de média e variância
##    que gerou os dados?

library( nloptr )
## Note que vamos usar uma função de verossimilhança invertida pois "nloptr", como muitos outros
##      otimizadores, foca na minimização das funções:
to.nlopt <- function(par){
  ## Função que retorna o inverso do log( likelihood )
  -1 * log.lik.norm(dt=nr, mu=par[1], var=par[2])
}

sink(file = "trace.txt", type = "output")
S1 = nloptr(x0 = c(1,1), eval_f = to.nlopt
            , opts = list(algorithm = "NLOPT_LN_SBPLX", print_level = 3))
sink()
## Solução extraida do StackOverflow:
## https://stackoverflow.com/questions/33635174/getting-trace-to-store-a-value-in-r
# get the solutions from the output file
solutionPath = readLines(con = file("trace.txt"))
# extract the solution path data out of the raw output
solutionPathParamRaw = solutionPath[grepl("^\tx", solutionPath)]
solutionPathParamMatch = gregexpr("(-)?[0-9]+(\\.[0-9]+)?", solutionPathParamRaw, perl = TRUE)
solutionPathParam = as.data.frame( t( sapply( regmatches( solutionPathParamRaw, solutionPathParamMatch), as.numeric, simplify = TRUE) ) )
colnames(solutionPathParam) <- c("mean", "var")

## Podemos observar as ultimas iterações da busca. Veja que chegamos na estimativa esperada.
tail( solutionPathParam )

## Agora podemos fazer um plot com os valores visitados a cada iteração de busca.
plot(x = 1:nrow(solutionPathParam), y = solutionPathParam$mean, type = "l", xlab = "Iteration", ylab = "Parameter value", lwd = 2)
lines(x = 1:nrow(solutionPathParam), y = solutionPathParam$var, col = "red", lwd = 2)
legend(x = 100, y = max(solutionPathParam), legend = c("mean","var"), lty = 1
       , col = c("black", "red"), xjust = 1, lwd =2)

## Trace no espaço 2D:
plot(x = 1, y = 1, xlim = range(solutionPathParam$mean), ylim = range(solutionPathParam$var), type = "n", xlab = "mean", ylab = "var")
for( i in 1:nrow(solutionPathParam)){
lines(x = solutionPathParam$mean[i:(i+1)], y = solutionPathParam$var[i:(i+1)], type = "b")
}

## Podemos fazer uma animação desse plot.
## NOTA: ESSA PARTE DA ANIMAÇÃO PODE NÃO RODAR NO SEU SISTEMA.

## dir.create("animate")
## for( i in 1:nrow(solutionPathParam)){
##   png(filename = paste0("animate/Animation_MLE_search_", i, ".png"))
##   plot(x = 1, y = 1, xlim = range(solutionPathParam$mean), ylim = range(solutionPathParam$var), type = "n", xlab = "mean", ylab = "var", main = "Busca da média e variância")
##   lines(x = solutionPathParam$mean[i:(i+1)], y = solutionPathParam$var[i:(i+1)], type = "b")
##   dev.off()
## }
## library( magick )
## ## Vamos ler as imagens que geramos e construir a animação.
## img_list <- list()
## for(i in 1:nrow(solutionPathParam)) img_list[[i]] <- image_read( paste0("animate/Animation_MLE_search_", i, ".png") )
## ## join the images together
## img_joined <- image_join(img_list)
## ## animate at 2 frames per second
## img_animated <- image_animate(img_joined, fps = 2)
## ## view animated image
## img_animated
## ## save to disk
## image_write(image = img_animated, path = "Animation_MLE_search.gif")

###########################################################
## Função de verossimilhança para um modelo de diversificação

library( diversitree )
library( TreeSim )

## Primeiro geramos uma história de diversificação:
phy <- sim.bd.age(age = 40, numbsim = 1, lambda = 0.2, mu = 0.1
                  , complete = FALSE, mrca = TRUE)[[1]]
## Plote a filogenia e o LTT plot.
plot( phy ); axisPhylo()
ltt.plot(phy)

## Qual a função de verossimilhança para o modelo birth-death?
help( "make.bd" )
bd_fn <- make.bd(tree = phy)
bd_fn

## Função já retorna em log(likelihood)
bd_fn(pars = c(0.2, 0.1))

## O perfil de likelihood para a taxa de especiação:
res_lambda <- vector()
lambda_vals <- seq(from = 0.001, to = 2, by = 0.01)
for(i in 1:length(lambda_vals)){
  res_lambda[i] <- bd_fn(pars = c(lambda_vals[i], 0.1))
}
plot(res_lambda~lambda_vals, main = "", ylab = "log( likelihood )", xlab = "Lambda", type = "b")

## O perfil de likelihood para a taxa de extinção:
res_mu <- vector()
mu_vals <- seq(from = 0.0, to = 2, by = 0.01)
for(i in 1:length(mu_vals)){
  res_mu[i] <- bd_fn(pars = c(0.2, mu_vals[i]))
}
plot(res_mu~mu_vals, main = "", ylab = "log( likelihood )", xlab = "Mu", type = "b")

x <- seq(from = 0.001, to = 2, by = 0.1) ## lambda
y <- seq(from = 0.0, to = 2, by = 0.1) ## mu
z <- matrix(data = NA, nrow = length(x), ncol = length(y))
for(i in 1:length(x)){
  for(j in 1:length(y)){
    z[i,j] <- bd_fn(pars = c(x[i], y[j]))
  }
}
persp(x, y, z, theta = 35, phi = 35, xlab = "lambda", ylab = "mu", zlab = "log(likelihood)"
      , border = "black", col = "grey", r = 4)
persp(x, y, z, theta = 0, phi = 15, xlab = "lambda", ylab = "mu", zlab = "log(likelihood)"
      , border = "black", col = "grey", r = 4)
persp(x, y, z, theta = 0, phi = 0, xlab = "lambda", ylab = "mu", zlab = "log(likelihood)"
      , border = "black", col = "grey", r = 4)

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
lm.log.lik <- function(x, y, a, b, sd){
  ## x = dados x
  ## y = dados y
  ## a = a term.
  ## b = b term.
  ## sd = standard deviation.
  pred <- a*x + b
  single.log.liks <- dnorm(y, mean = pred, sd = sd, log = T)
  sum(single.log.liks)
}

## Gerando os dados.

## Valores para os parametros:
trueA <- 5
trueB <- 0
trueSd <- 10
sampleSize <- 30

# create independent x-values
x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
# create dependent values according to ax + b + N(0,sd)
y <-  trueA * x + trueB + rnorm(n=sampleSize,mean=0,sd=trueSd)
plot(x,y, main="Test Data")
