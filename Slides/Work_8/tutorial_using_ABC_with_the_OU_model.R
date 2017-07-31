## Nesse tutorial vamos explorar o uso de MCMC no caso de que não temos a função de verossimilhança para
##  o modelo. Como vimos anteriormente, o MCMC usa a verossimilhança para comparar se o valor do parametro
##  da proposta com o valor atual do mesmo parametro na cadeia. No caso em que não temos a função de
##  verossimilhança precisamos usar simulações para comparar dados simulados com os dados observados.
## Esse método é conhecido como Approximate Bayesian Computation (ABC).

## Aqui vamos usar o modelo OU no lugar do modelo BM.

## Vamos definir o sliding window:
sliding.window <- function(x, w){
    ## Sliding window proposal for unbounded trait.
    ## x = the current value.
    ## w = the width parameter of the proposal.
    y <- runif(1, min = x - (w/2), max = x + (w/2) )
    return(y)
}

## Vamos usar uma distribuição de proposta para valores positivos.
## No caso do ABC a simulação dos dados acontece antes do passo de
##   accept and reject. Isso significa que o MCMC vai tentar simular
##   dados com a variância negativa, por exemplo, e vai gerar um erro.
## Isso é uma limitação do programa que simula os dados.
## Veja pagina de ajuda para 'fastBM' no pacote 'phytools'.
reflect.sliding.window <- function(x, w){
    ## Sliding window proposal for unbounded trait.
    ## x = the current value.
    ## w = the width parameter of the proposal.
    y <- abs( runif(1, min = x - (w/2), max = x + (w/2) ) )
    return(y)
}

## Abaixo está a função de MCMC. Modificada para o ABC.
## Compare essa função com a do tutorial anterior.

mcmc.ABC.OU <- function(phy, data, initial.root, initial.rate, initial.tetha, initial.alpha
                      , prior.root, prior.rate, prior.tetha, prior.alpha, eps=0.01, w=0.5, gen){
    ## phy = filogenia no formato 'ape'.
    ## data = vetor de dados para as espécies.
    ## initial.rate = valor inicial para rate.
    ## initial.root = valor inicial para root.
    ## initial.tetha = valor inicial para tetha.
    ## initial.alpha = valor inical para alpha. 
    ## prior.root = função de prior para root. Com unico argumento igual ao valor.
    ## prior.rate = função de prior para rate. Com unico argumento igual ao valor.
    ## prior.tetha = função de prior para tetha. Com unico argumento igual ao valor.
    ## prior.alpha = função de prior para alpha. Com unico argumento igual ao valor.
    ## eps = valor para a tolerância utilizada quando calculando a distância entre os
    ##     dados simulados e os dados observados.
    ## w = a width (largura) da função uniforme de sliding window.
    ## gen = número de gerações da análise.
    
    ## Gerar vetores para guardar os resultados da cadeia.
    tetha.chain <- numeric(length = gen) ## Cadeia para tetha.
    rate.chain <- numeric(length = gen) ## Cadeia para rate.
    alpha.chain <- numeric(length = gen) ## Cadeia para alpha.
    root.chain <- numeric(length = gen) ## Cadeia para root.
    accept <- numeric(length = gen) ## Para calcular o acceptance ratio.
    
    ## Primeiro estado é o initial state:
    tetha.chain[1] <- initial.tetha
    rate.chain[1] <- initial.rate
    alpha.chain[1] <- initial.alpha
    root.chain[1] <- initial.root
    accept[1] <- 1 ## Primeira geração é o valor inicial.

    ## Calcular a densidade do prior para o estado inicial.
    current.prior <- prior.root(initial.root) + prior.rate(initial.rate) + prior.tetha(initial.tetha) + prior.alpha(initial.alpha)

    ## Calcular a summary statistics para o estado inicial da cadeia.
    data.ou <- fastBM(tree = phy, a = initial.root, sig2 = initial.rate
                    , alpha = initial.alpha, theta = initial.tetha)
    
    ## Fazer um loop do tamanho de gen. Esse é o MCMC.
    ## NOTE que o loop começa do número 2. Note também que uso 'i-1' para o estado
    ##    atual da cadeia e 'i' para o estado futuro.
    for( i in 2:gen ){
        
        ## Selecionar qual parâmetro fazer o update:
        par.update <- sample(x = c("root","tetha","rate","alpha"), size = 1)
        
        ## Fazer o update, accept and reject e salvar o estado:
        if( par.update == "tetha"){
            ## Proposta de novo valor:
            prop.tetha <- sliding.window(x = tetha.chain[i-1], w = w)

            ## Simulamos novos dados para as espécies usando a mesma filogenia..
            data.ou <- fastBM(tree = phy, a = root.chain[i-1], sig2 = rate.chain[i-1]
                            , alpha = alpha.chain[i-1], theta = prop.tetha)

            ## Calculamos os indices de distancia (summary statistics):
            prop.sq.mean <- ( mean(data) - mean(data.ou) )^2
            prop.sq.sd <- ( sd(data) - sd(data.ou) )^2

            ## E o valor do prior para a proposta.
            prop.prior <- prior.root(root.chain[i-1]) + prior.rate(rate.chain[i-1]) + prior.tetha(prop.tetha) + prior.alpha(alpha.chain[i-1])

            ## Esse é o passo mais diferente do algoritmo.
            ## Primeiro comparamos a distância entre os dados simulados e os dados observados condicionados
            ##  nos parametros da proposta. Caso esteja dentro de 'eps' (nossa tolerância), nos passamos
            ##  para o passo de aceite ou rejeite. Do contrário a proposta é rejeitada.

            if( prop.sq.mean >= eps & prop.sq.sd >= eps){
                ## Distância da simulação está dentro da proposta.
                r <- min(1, exp(prop.prior - current.prior) ) ## Aceite ou rejeite depende somente do prior.
                if( r > runif(1) ){
                    ## Accept.
                    tetha.chain[i] <- prop.tetha
                    current.prior <- prop.prior
                    accept[i] <- 1
                } else{
                    ## Reject.
                    tetha.chain[i] <- tetha.chain[i-1]
                    accept[i] <- 0
                }            
            } else{
                tetha.chain[i] <- tetha.chain[i-1]
                accept[i] <- 0
            }
            rate.chain[i] <- rate.chain[i-1]
            root.chain[i] <- root.chain[i-1]
            alpha.chain[i] <- alpha.chain[i-1]
            
        }
        if( par.update == "root"){
            ## Proposta de novo valor:
            prop.root <- sliding.window(x = root.chain[i-1], w = w)

            ## Simulamos os dados novamente.
            data.ou <- fastBM(tree = phy, a = prop.root, sig2 = rate.chain[i-1]
                            , alpha = alpha.chain[i-1], theta = tetha.chain[i-1])

            ## Calculamos os indices de distancia (summary statistics):
            prop.sq.mean <- ( mean(data) - mean(data.ou) )^2
            prop.sq.sd <- ( sd(data) - sd(data.ou) )^2

            ## E o valor do prior para a proposta.
            prop.prior <- prior.root(prop.root) + prior.rate(rate.chain[i-1]) + prior.tetha(tetha.chain[i-1]) + prior.alpha(alpha.chain[i-1])

            ## Esse é o passo mais diferente do algoritmo.
            ## Primeiro comparamos a distância entre os dados simulados e os dados observados condicionados
            ##  nos parametros da proposta. Caso esteja dentro de 'eps' (nossa tolerância), nos passamos
            ##  para o passo de aceite ou rejeite. Do contrário a proposta é rejeitada.

            if( prop.sq.mean >= eps & prop.sq.sd >= eps){
                ## Distância da simulação está dentro da proposta.
                r <- min(1, exp(prop.prior - current.prior) ) ## Aceite ou rejeite depende somente do prior.
                if( r > runif(1) ){
                    ## Accept.
                    root.chain[i] <- prop.root
                    current.prior <- prop.prior
                    accept[i] <- 1
                } else{
                    ## Reject.
                    root.chain[i] <- root.chain[i-1]
                    accept[i] <- 0
                }            
            } else{
                root.chain[i] <- root.chain[i-1]
                accept[i] <- 0                
            }
            rate.chain[i] <- rate.chain[i-1]
            tetha.chain[i] <- tetha.chain[i-1]
            alpha.chain[i] <- alpha.chain[i-1]
        }
        if( par.update == "rate"){
            ## Proposta de novo valor:
            prop.rate <- reflect.sliding.window(x = rate.chain[i-1], w = w)

            ## Simulamos os dados novamente.
            data.ou <- fastBM(tree = phy, a = root.chain[i-1], sig2 = prop.rate
                            , alpha = alpha.chain[i-1], theta = tetha.chain[i-1])

            ## Calculamos os indices de distancia (summary statistics):
            prop.sq.mean <- ( mean(data) - mean(data.ou) )^2
            prop.sq.sd <- ( sd(data) - sd(data.ou) )^2

            ## E o valor do prior para a proposta.
            prop.prior <- prior.root(root.chain[i-1]) + prior.rate(prop.rate) + prior.tetha(tetha.chain[i-1]) + prior.alpha(alpha.chain[i-1])

            ## Esse é o passo mais diferente do algoritmo.
            ## Primeiro comparamos a distância entre os dados simulados e os dados observados condicionados
            ##  nos parametros da proposta. Caso esteja dentro de 'eps' (nossa tolerância), nos passamos
            ##  para o passo de aceite ou rejeite. Do contrário a proposta é rejeitada.

            if( prop.sq.mean >= eps & prop.sq.sd >= eps){
                ## Distância da simulação está dentro da proposta.
                r <- min(1, exp( prop.prior - current.prior ) ) ## Aceite ou rejeite depende somente do prior.
                if( r > runif(1) ){
                    ## Accept.
                    rate.chain[i] <- prop.rate
                    current.prior <- prop.prior
                    accept[i] <- 1
                } else{
                    ## Reject.
                    rate.chain[i] <- rate.chain[i-1]
                    accept[i] <- 0
                }            
            } else{
                rate.chain[i] <- rate.chain[i-1]
                accept[i] <- 0
            }
            root.chain[i] <- root.chain[i-1]
            tetha.chain[i] <- tetha.chain[i-1]
            alpha.chain[i] <- alpha.chain[i-1]
            
        }
        if( par.update == "alpha"){
            ## Proposta de novo valor:
            prop.alpha <- reflect.sliding.window(x = alpha.chain[i-1], w = w)

            ## Simulamos os dados novamente.
            data.ou <- fastBM(tree = phy, a = root.chain[i-1], sig2 = rate.chain[i-1]
                            , alpha = prop.alpha, theta = tetha.chain[i-1])

            ## Calculamos os indices de distancia (summary statistics):
            prop.sq.mean <- ( mean(data) - mean(data.ou) )^2
            prop.sq.sd <- ( sd(data) - sd(data.ou) )^2

            ## E o valor do prior para a proposta.
            prop.prior <- prior.root(root.chain[i-1]) + prior.rate(rate.chain[i-1]) + prior.tetha(tetha.chain[i-1]) + prior.alpha(prop.alpha)
            
            ## Esse é o passo mais diferente do algoritmo.
            ## Primeiro comparamos a distância entre os dados simulados e os dados observados condicionados
            ##  nos parametros da proposta. Caso esteja dentro de 'eps' (nossa tolerância), nos passamos
            ##  para o passo de aceite ou rejeite. Do contrário a proposta é rejeitada.

            if( prop.sq.mean >= eps & prop.sq.sd >= eps){
                ## Distância da simulação está dentro da proposta.
                r <- min(1, exp(prop.prior - current.prior) ) ## Aceite ou rejeite depende somente do prior.
                if( r > runif(1) ){
                    ## Accept.
                    alpha.chain[i] <- prop.alpha
                    current.prior <- prop.prior
                    accept[i] <- 1
                } else{
                    ## Reject.
                    alpha.chain[i] <- alpha.chain[i-1]
                    accept[i] <- 0
                }            
            } else{
                alpha.chain[i] <- alpha.chain[i-1]
                accept[i] <- 0
            }
            root.chain[i] <- root.chain[i-1]
            tetha.chain[i] <- tetha.chain[i-1]
            rate.chain[i] <- rate.chain[i-1]

        }
        
    }    

    return(data.frame(root = root.chain, rate = rate.chain, tetha = tetha.chain
                    , alpha = alpha.chain, accept = accept))
}

## Agora vamos gerar dados para testar a performance do ABC.
if(!require(geiger)) {install.packages("geiger"); library(geiger)}
if(!require(phytools)) {install.packages("phytools"); library(phytools)}

phy <- rtree(100)
phy <- compute.brlen(phy)
data.ou.obs <- fastBM(tree = phy, a = 10, sig2 = 0.5, alpha = 1.5, theta = 10)

## Veja na função definida acima 'mcmc.ABC.OU' que vamos utilizar a média e o
##   desvio padrão dos dados observados para comparar os dados simulados com os
##   dados observados.
## Veja a média e o desvio padrão dos dados:
mean( data.ou.obs )
sd( data.ou.obs )

## Agora vamos fazer a estimativa usando ABC.
## Primeiro definimos os priors. Como esse modelo tem mais parametros do que
##   o modelo BM, é uma boa idéia facilitar o trabalho da estimativa.
## Isso é importante especialmente pois estamos usando o método de ABC.
prior.rate <- function(x){ ## Usually a small number.
    dexp(x, rate=3, log=TRUE)
}
obs.mean <- mean(data.ou.obs)
prior.root <- function(x){ ## Centered in the mean of the observed data.
    dnorm(x, mean=obs.mean, sd=5, log=TRUE)
}
obs.range <- range( data.ou.obs )
prior.tetha <- function(x){ ## Something within the range of the data.
    dunif(x, min=obs.range[1], max=obs.range[2], log=TRUE)
}
prior.alpha <- function(x){ ## Have no idea! It is a positive number and 20 seems too large.
    dunif(x, min=0, max=20, log=TRUE)
}

## Vamos pegar os estados iniciais das distribuições prior.
initial.root <- rnorm(1, mean=obs.mean, sd=5)
initial.rate <- rexp(1, rate=3)
initial.tetha <- runif(1, min=obs.range[1], max=obs.range[2])
initial.alpha <- runif(1, min=0, max=20)

fit.abc <- mcmc.ABC.OU(phy=phy, data=data.ou.obs, initial.root, initial.rate, initial.tetha, initial.alpha
          , prior.root, prior.rate, prior.tetha, prior.alpha, eps=0.01, w=0.5, gen=100000)

## Podemos fazer gráficos para explorar a distribuição posterior.
## data.ou.obs <- fastBM(tree = phy, a = 10, sig2 = 0.5, alpha = 1.5, theta = 20)
hist( fit.abc$root, main="root")
abline( v=10, col="red")

hist( fit.abc$rate, main="rate")
abline( v=0.5, col="red")

hist( fit.abc$tetha, main="tetha")
abline( v=10, col="red")

hist( fit.abc$alpha, main="alpha")
abline( v=1.5, col="red")

## A distribuição posterior é diferente da distribuição prior nessa análise?
## E sobre as estimativas dos parametros? Parece que a cadeia de MCMC usando Approximate Bayesian Computation
##   está convergendo para o valor de parametros que usamos para criar os dados?

## Com base nos tutoriais anteriores, como você faria para checar a convergência dessa análise?
## Cheque se a análise está convergendo.

## O que está acontecendo aqui? Tente listar as diferenças entre o algoritmo usado para o ABC e o de
##   MCMC que aprendemos anteriormente. Quais destes passos parece poder criar problemas e por quê?
