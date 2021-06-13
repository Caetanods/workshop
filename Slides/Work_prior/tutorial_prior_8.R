## Nesse tutorial vamos explorar o papel da distribuição a priori nas nossas estimativas.
## A busca feita pelo MCMC é guiada pelos valores da 'odds ratio'. Essa tem três
##   partes: likelihood ratio * prior ratio * Hastings ratio.
## A likelihood ratio é a razão da verossimilhança dos parâmetros, a prior ratio
##   é a razão da probabilidade dos parâmetros tirada da distribuição a priori ou 
##   prior distribution e, por fim, a Hastings ratio é uma correção no caso da
##   distribuição de proposta não ser simétrica.

## Nesse tutorial vamos focar na distribuição a priori que faz parte da 'odds ratio'.

## Existem dois tipos de priors, priors informativos (como a distribuição normal,
##   lognormal ou exponencial) e priors não-informativos (geralmente distribuição
##   uniforme). No caso de priors não-informativos, a distribuição uniforme é utilizada
##   pois ela garante que a probabilidade para qualquer valor dentro do prior
##   é o mesmo.

## Infelizmente precisamos dar um valor de máximo e minimo para a distribuição
##   uniforme. Ou seja, embora a ideia é ser 'naive' quanto à probabilidade dos
##   parâmetros, todo parâmetro que estiver fora da distribuição uniforme tem
##   densidade igual a 0.

## Exemplos de valores.
par <- c(-1, 1, 50, 99, 101)
## Probabilidades.
dunif(x = par, min = 0, max = 100)
## Probabilidades em log.
dunif(x = par, min = 0, max = 100, log = TRUE)

## Primeiro volte para o cálculo da 'odds ratio' e o passo de 'accept or reject'
##   (aceite ou rejeite) do MCMC.
## Qual é o efeito de log(prior) = -Inf? Esse valor de parâmetro vai ser aceito com
##   facilidade?

## Fica claro que a distribuição de prior influência no aceite dos valores de
##   parâmetros. Mas e os dados?

## Quando temos dados o bastante, o sinal das abservações, as evidências, a informação
##   em nossos dados é grande o bastante para desviar significativamente do prior.
## É comum que a posterior, após a análise de MCMC, tenha maior densidade super longe
##   da distribuição do prior. Nesses casos a distribuição do prior e da posterior em
##   um plot praticamente não se soprepõem.

## Vamos fazer alguns experimentos para mostrar o efeito do prior em nossas análises.
## Para tal vamos usar os mesmos priors do 'tutorial 6' e um prior uniforme.
## Vamos também variar a quantidade de observações. Para tal vamos reduzir o número
##   de espécies na filogenia. Quanto mais espécies temos em uma filogenia mais
##   confiável são nossas estimativas dos parâmetros. Quanto menos dados, maior é
##   a influência do prior na posterior resultante.

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

## Como o número de simulações para esse tutorial é muito grande, nós não vamos
##   rodar todas as cadeias aqui. Veja o arquivo 'anexo_tutorial_8.R' para rodar todas
##   as cadeias novamente.

library(coda)

## Carregue os resultados das simulações e verifique o 'traceplot'.
load(file = "./data_sim_prior.Rdata")
traceplot(inf.10)
traceplot(inf.25)
traceplot(inf.50)
traceplot(inf.100)

## Acima estão os traceplots das simulações com 10, 25, 50 e 100 espécies. Note que
##   somente a simulação com 100 espécies convergiu nesse número de gerações (10000).
## Por isso repeti o experimento com um número 5x maior de gerações (50000).

load(file = "./data_long_sim_prior.Rdata")
traceplot(inf.10) ## Converge ~20000 gen.
traceplot(inf.25) ## Converge ~10000 gen.
traceplot(inf.50) ## Converge ~5000 gen.
traceplot(inf.100) ## Rodou somente 10000, converge ~ 2000 gen.

## Pergunta: Note que o tempo para convergência parece ser inversamente proporcional
##    ao número de tips da filogenia. Por que isso?

## Agora vamos retirar o burn-in, o número de gerações depende do número de tips.
inf.10 <- window(x = inf.10, start = 30000)
unif.10 <- window(x = unif.10, start = 30000)
inf.25 <- window(x = inf.25, start = 20000)
unif.25 <- window(x = unif.25, start = 20000)
inf.50 <- window(x = inf.50, start = 15000)
unif.50 <- window(x = unif.50, start = 15000)
inf.100 <- window(x = inf.100, start = 3000)
unif.100 <- window(x = unif.100, start = 3000)

## Agora vamos checar o intervalo de 95% da Highest Posterior Density (HPD).
## Relembre que: rate = 0.5 e root = 10

summary(inf.10)[[2]][2:3,c(1,5)]
summary(unif.10)[[2]][2:3,c(1,5)]
summary(inf.25)[[2]][2:3,c(1,5)]
summary(unif.25)[[2]][2:3,c(1,5)]
summary(inf.50)[[2]][2:3,c(1,5)]
summary(unif.50)[[2]][2:3,c(1,5)]
summary(inf.100)[[2]][2:3,c(1,5)]
summary(unif.100)[[2]][2:3,c(1,5)]

## Vamos fazer o plot das posteriores. Primeiro do parãmetro 'rate':
par(mfrow = c(2,4) )
hist(x = inf.10[[1]][,3], main = "info 10")
hist(x = unif.10[[1]][,3], main = "unif 10")
hist(x = inf.25[[1]][,3], main = "info 25")
hist(x = unif.25[[1]][,3], main = "unif 15")
hist(x = inf.50[[1]][,3], main = "info 50")
hist(x = unif.50[[1]][,3], main = "unif 50")
hist(x = inf.100[[1]][,3], main = "info 100")
hist(x = unif.100[[1]][,3], main = "unif 100")

## Agora do valor da raíz, o 'root'.
hist(x = inf.10[[1]][,2], main = "info 10")
hist(x = unif.10[[1]][,2], main = "unif 10")
hist(x = inf.25[[1]][,2], main = "info 25")
hist(x = unif.25[[1]][,2], main = "unif 15")
hist(x = inf.50[[1]][,2], main = "info 50")
hist(x = unif.50[[1]][,2], main = "unif 50")
hist(x = inf.100[[1]][,2], main = "info 100")
hist(x = unif.100[[1]][,2], main = "unif 100")

## Note que nas simulações acima o prior causou uma diferença nas estimativas. No
##   entanto a distribuição de prior informativa incluí os valores reais da
##   simulação com alta probabilidade. O que aconteceria caso o prior fosse tão
##   informativo quanto, mas estivesse centrado em um valor de parâmetro longe
##   do valor esperado?

## Abaixo estão resultados da mesma simulação mas com os seguintes priors:
hist( rnorm(n = 1000, mean = 10, sd = 2) , main = "Prior for rate")
hist( rnorm(n = 1000, mean = 20, sd = 2) , main = "Prior for root")
## O prior para rate está muitas vezes maior do que o valor que gerou os dados,
##   o prior para o valor da raiz está 2x maior do que o valor que gerou os dados.

## Vamos carregar os resultados dessas novas simulações:

load(file = "./data_long_sim_bad_prior.Rdata")
traceplot(inf.10) ## Converge ~20000 gen.
traceplot(inf.25) ## Converge ~10000 gen.
traceplot(inf.50) ## Converge ~5000 gen.
traceplot(inf.100) ## Rodou somente 10000, converge ~ 2000 gen.

inf.10 <- window(x = inf.10, start = 30000)
unif.10 <- window(x = unif.10, start = 30000)
inf.25 <- window(x = inf.25, start = 20000)
unif.25 <- window(x = unif.25, start = 20000)
inf.50 <- window(x = inf.50, start = 15000)
unif.50 <- window(x = unif.50, start = 15000)
inf.100 <- window(x = inf.100, start = 3000)
unif.100 <- window(x = unif.100, start = 3000)

## Intervalo de 95% HPD
summary(inf.10)[[2]][2:3,c(1,5)]
summary(unif.10)[[2]][2:3,c(1,5)]
summary(inf.25)[[2]][2:3,c(1,5)]
summary(unif.25)[[2]][2:3,c(1,5)]
summary(inf.50)[[2]][2:3,c(1,5)]
summary(unif.50)[[2]][2:3,c(1,5)]
summary(inf.100)[[2]][2:3,c(1,5)]
summary(unif.100)[[2]][2:3,c(1,5)]

## Vamos fazer o plot das posteriores. Primeiro do parãmetro 'rate':
par(mfrow = c(2,4) )
hist(x = inf.10[[1]][,3], main = "info 10")
hist(x = unif.10[[1]][,3], main = "unif 10")
hist(x = inf.25[[1]][,3], main = "info 25")
hist(x = unif.25[[1]][,3], main = "unif 15")
hist(x = inf.50[[1]][,3], main = "info 50")
hist(x = unif.50[[1]][,3], main = "unif 50")
hist(x = inf.100[[1]][,3], main = "info 100")
hist(x = unif.100[[1]][,3], main = "unif 100")

## Agora do valor da raíz, o 'root'.
hist(x = inf.10[[1]][,2], main = "info 10")
hist(x = unif.10[[1]][,2], main = "unif 10")
hist(x = inf.25[[1]][,2], main = "info 25")
hist(x = unif.25[[1]][,2], main = "unif 15")
hist(x = inf.50[[1]][,2], main = "info 50")
hist(x = unif.50[[1]][,2], main = "unif 50")
hist(x = inf.100[[1]][,2], main = "info 100")
hist(x = unif.100[[1]][,2], main = "unif 100")