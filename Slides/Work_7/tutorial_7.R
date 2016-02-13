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
## Lembre-se de mudar o diretório abaixo para onde os seus arquivos estão:
setwd("~/Documents/Academicos/Harmon_Lab/MCMC - Workshop - Brasil/workshop/Slides/Work_7")
mcmc.1 <- readRDS(file = "../Work_6/mcmc.BM.1.rds")
mcmc.2 <- readRDS(file = "../Work_6/mcmc.BM.2.rds")
mcmc.3 <- readRDS(file = "../Work_6/mcmc.BM.3.rds")

## Agora criamos um objeto da classe 'mcmc'.
## Note que a função tem opções para retirada do burn-in (informando a geração de
##   início) e aplicação do thinning. Mas, por enquanto, vamos usar a cadeia toda.
mcmc.1 <- mcmc(mcmc.1)
mcmc.2 <- mcmc(mcmc.2)
mcmc.3 <- mcmc(mcmc.3)

## Agora combinamos as cadeias em um só objeto. Fazendo isso informamos ao 'coda' que
##   essas cadeias são réplicas com diferentes pontos de partida.
help(mcmc.list)
mcmc <- mcmc.list(mcmc.1, mcmc.2, mcmc.3)

## Primeiramente podemos ver o 'summary' das nossas cadeias:
summary(mcmc)
## Nesse summary podemos ver uma quantidade importante. O 95% HPD, a região de 95%
##   de maior densidade (i.e., probabilidade) dos parâmetros. Veja, por exemplo,
##   que o valor de 'root' (valor da raíz da filogenia) varia entre 9.41 e 11.47 no
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

## Vamos ver o novo intervalo de HPD das variáveis na posterior sem o burn-in:
summary(post.mcmc)

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

## A autocorrelação é bastante alta até um Lag de 10. Indicando que para diminuir
##   esse valor de autocorrelação deveriamos descartar valores da cadeia e reter somente
##   uma em cada 10 gerações. Não vamos fazer esse passo agora.

## Qual o número de Effective Sample Size?
effectiveSize(x = post.mcmc)

## Mesmo com a alta autocorrelação temos um número mínimo de amostragem igual a 407
##   combinando as diferentes cadeias. Ou seja, não há necessidade de rodas mais
##   gerações de MCMC para aumentar o número de ESS. Muito provavelmente ~ 400
##   pontos independentes de amostragem é o bastante para se ter uma idéia da
##   distribuição posterior dos parâmetros.

densplot(x = post.mcmc)

## Veja a densidade para 'root' a média do processo de Brownian Motion (valor da raíz)
##   e 'rate' a variância do processo (ou taxa).
## Essa distribuição é chamada de posterior distribution e é nosso resultado final.
## Agora podemos ustilizar tal distribuição para fazer conclusões sobre nosso modelo.

## Um ponto de discussão importante que é difícil de ser abordado é "Quantas
##  gerações de MCMC eu preciso rodar?". Uma prática comum é verificar artigos com
##  problemas similares e ver quanto é geralmente aceito. Para inferência filogenética
##  , por exemplo, geralmente são utilizados valores na casa dos milhões. No entanto,
##  problemas com modelos de 2 a 5 parâmetros acabam ficando na cada dos milhares.
## Existe um método que tenta estimar quantas gerações seriam necessárias para
##  atingir um certo número de samples. Esse método é bom para verificar
##  cadeias 'teste'.

raftery.diag(data = mcmc.1)

## DESAFIO: Veja que o número total de gerações para estimar o parâmetro 'rate'
##  é maior do que rodamos. Tente aumentar o número de gerações. Qual a diferença
##  na estimativa de convergência, ESS, e do 95% HPD?

## Nesse tutorial mostrei como acessar convergência e verificar burn-in quando
##  temos mais de uma cadeia de MCMC começando de pontos iniciais diferentes.
##  Essa é a situação ideal e deveria ser replicada sempre que possível.

## No entanto, algumas análises demoram dias, ou mesmo semanas, para serem completas
##  e pode ser bastante complicado rodar uma série de cadeias de MCMC com a mesma
##  análise. O contraponto para esse argumento é que no final das análises, após
##  verificar convergência você poderá COMBINAR AS CADEIAS DA POSTERIOR, ou seja,
##  o número de amostras vai aumentar além de ter mais certeza de que o ótimo
##  global foi atingido (pois as cadeias chegaram em um mesmo ponto vindo de
##  pontos diferentes).
## Bom, se eu vou rodar uma análise que leva semanas, eu gostaria de saber no final
##  que o resultado está aceitável, não é mesmo?

## Mas existem outras maneiras de verificar convergência sem rodar mais de uma
##  cadeia. Veja as funções abaixo e tente verificar a convergência de uma das
##  cadeias isoladamente:

help(geweke.diag)
help(heidel.diag)

## DESAFIO: Tente rodar esses diagnósticos e verificar a convergência das cadeias
##  uma a uma. Descreva as vantagens e desvantagens de cada diagnóstico.