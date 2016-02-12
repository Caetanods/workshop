## No tutorial anterior completamos uma análise usando o algoritmo de Metropolis.
## Nesse tutorial vamos fazer análises da posterior.
## Caso tenha resultados de uma análise de alguma outra análise de MCMC você pode
##   usar seus próprios dados para rodar esse tutorial. Note que a formação dos
##   objetos e interpretação dos parãmetros fica por sua conta! ;)

## Para nossas análises vamos utilizar o pacote 'coda'.
library(coda)

## Esse pacote é dedicado a diversas análises da posterior. Veja a lista de funções
##   no pdf do pacote na pasta 'Work_7'.

## Primeiro precisamos carregar nossa cadeia.
mcmc.1 <- readRDS(file = "../Work_6/mcmc.BM.1.rds") ## Error here!
mcmc.2 <- readRDS(file = "../Work_6/mcmc.BM.2.rds")
mcmc.3 <- readRDS(file = "../Work_6/mcmc.BM.3.rds")

## Agora criamos um objeto da classe 'mcmc'.
## Note que a função tem opções para retirada do burn-in (informando a geração de
##   início) e aplicação do thinning. Mas, por enquanto, vamos usar a cadeia toda.
mcmc.2 <- mcmc(mcmc.2)
mcmc.3 <- mcmc(mcmc.3)

## Agora combinamos as cadeias em um só objeto. Fazendo isso informamos ao 'coda' que
##   essas cadeias são réplicas com diferentes pontos de partida.
help(mcmc.list)
mcmc <- mcmc.list(mcmc.2, mcmc.3)

## Primeiramente podemos ver o 'summary' das nossas cadeias:
summary(mcmc)
## Nesse summary podemos ver uma quantidade importante. O 95% HPD, a região de 95%
##   de maior densidade (i.e., probabilidade) dos parâmetros. Veja, por exemplo,
##   que o valor de 'root' (valor da raíz da filogenia) varia entre 15.06 e 44.08 no
##   95% HPD. Mas ainda temos outras análises a fazer.

## Antes de mais nada vamos retirar o burn-in. Para tal vamos observar o comportamento
##   das cadeias usando plots.
## Essa função vai produzir mais de um plot. Então talvez você terá que mover entre
##   as janelas dependendo do programa que esteja usando.
traceplot(x = mcmc)

## Note que as cadeias são colocadas no mesmo gráfico. Isso facilita a observação da
##   convergência. Podemos visualmente verificar se as diferentes cadeias estão
##   amostrando o mesmo ótimo, mesmo sendo inciadas em pontos diferentes.
## Lembre-se que os pontos iniciais foram escolhidos com o objetivo de estarem
##   espalhados pela superfície de verossimilhança.

## Pelos gráficos parece que a partir da geração '2000' já temos uma estabilização.
## Vamos retirar as gerações de 1:2000 e plotar novamente:
post.mcmc <- window(mcmc, start=2000, end=10000)
traceplot(post.mcmc)

## Já está com uma cara melhor, não?

## Agora que retiramos o burn-in vamos fazer alguns diagnósticos de convergência.
## Vamos usar o Gelman and Rubin R que verifica se a variação entre cadeias é menor
##   que a variação dentro das cadeias.
## Vamos colocar 'multivariate' pois queremos avaliar cada um dos parâmetros
##  separadamente. Veja outras opções no help.
help(gelman.diag)
gelman.diag(x = post.mcmc, transform = FALSE, autoburnin = FALSE, multivariate = FALSE)

## Queremos que o resultado seja o mais próximo possível de '1.00'. E... BINGO!

## OK. Retiramos o burn-in e sabemos que as cadeias convergiram. O que mais podemos
##  estar interessados?
## Antes de olhar diretamente para a distribuição posterior dos parâmetros vamos 
##  verificar a autocorrelação e calcular o tamanho efetivo da amostragem (ESS).
## 'realtive = FALSE' pois não aplicamos o thinning na cadeia.
## Essa função mostra tanto a autocorrelação entre parâmetros como de um parâmetro só.
## O valor 'Lag x' é o número de gerações que foram puladas para calcular a
##   autocorrelação. 'Lag 5' significa que foram somente avaliadas gerações de 5 em 5.
autocorr(x = post.mcmc, relative = FALSE)
