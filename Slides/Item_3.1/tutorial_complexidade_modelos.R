## Para esse tutorial vamos utilizar novamente os modelos de diversificação.

library( diversitree )
library( TreeSim )

## Primeiro geramos uma história de diversificação:
phy <- sim.bd.age(age = 40, numbsim = 1, lambda = 0.2, mu = 0.1
                  , complete = FALSE, mrca = TRUE)[[1]]
## Plote a filogenia e o LTT plot.
plot( phy ); axisPhylo()
ltt.plot(phy)

## Crie a função de verossimilhança para o modelo:
bd_fn <- make.bd(tree = phy)
bd_fn

## Função já retorna em log(likelihood)
bd_fn(pars = c(0.2, 0.1))

## #############################################
## Modelo livre - two free parameters:
init <- runif(n = 1, min = 0.00001, max = 10)
( random_start <- c(init, init/2) )
## Estamos usando o método "L-BFGS-B" aqui para poder controlar limites da busca.
free_fit <- optim(par = random_start, fn = bd_fn, control=list(fnscale=-1)
                  , method = "L-BFGS-B", lower = 0.000001, upper = 10)
free_fit$par ## Estimativa não é ruim!

## #############################################
## Modelo reduzido - one free parameter:

## Vamos dizer que sabemos que a taxa de extinção é de 0.3 pois um trabalho anterior
##      estimou a taxa usando o registro fóssil.
## Nesse caso precisamos controlar a função de verossimilhança:
init <- runif(n = 1, min = 0.00001, max = 10)
known_mu_fn <- function(lambda){
    bd_fn( pars = c(lambda, 0.3) )
}
known_mu_fn( init )
known_mu_fit <- optim(par = init, fn = known_mu_fn, control=list(fnscale=-1)
                      , method = "L-BFGS-B", lower = 0.000001, upper = 10)
known_mu_fit$par

## #############################################
## Outro modelo reduzido - one free parameter:

## Agora sabemos que a especiação é 0.1 .
init <- runif(n = 1, min = 0.00001, max = 10)
known_lambda_fn <- function(mu){
    bd_fn( pars = c(0.1, mu) )
}
known_lambda_fn( init )
known_lambda_fit <- optim(par = init, fn = known_lambda_fn, control=list(fnscale=-1)
                          , method = "L-BFGS-B", lower = 0.000001, upper = 10)
known_lambda_fit$par ## Veja que a estimativa chegou no limite minimo!
known_lambda_fit

## Vamos repetir com um limite minimo diferente:
known_lambda_fit <- optim(par = init, fn = known_lambda_fn, control=list(fnscale=-1)
                          , method = "L-BFGS-B", lower = 0.0, upper = 10)
known_lambda_fit$par ## Extinção estimada como 0. Esse é o melhor que o modelo pode fazer.

## #############################################
## Outro modelo reduzido - one free parameter:

## Agora somente sabemos que as taxas de lambda e mu são iguais.
init <- runif(n = 1, min = 0.00001, max = 10)
equal_rates_fn <- function(x){
    bd_fn( pars = c(x, x) )
}
equal_rates_fn(init)
equal_rates_fit <- optim(par = init, fn = equal_rates_fn, control=list(fnscale=-1)
                         , method = "L-BFGS-B", lower = 0.000001, upper = 10)
equal_rates_fit$par

## #############################################
## Outro modelo reduzido - one free parameter:

## Agora sabemos que as taxas de especiação são o dobro das de extinção.
init <- runif(n = 1, min = 0.00001, max = 10)
double_lambda_fn <- function(x){
    bd_fn( pars = c(x, x/2) )
}
double_lambda_fn(init)
double_lambda_fit <- optim(par = init, fn = double_lambda_fn, control=list(fnscale=-1)
                           , method = "L-BFGS-B", lower = 0.000001, upper = 10)
double_lambda_fit$par

## ############################################
## ############################################

## Agora vamos ver que independente do modelo, o modelo livre sempre tem maior valor de verossimilhança:

log_lik_models <- c(free_fit$value, known_mu_fit$value, known_lambda_fit$value, equal_rates_fit$value, double_lambda_fit$value)
free_fit$value - log_lik_models

## ############################################
## Finalmente, vamos salvar os resultados para usar mais tarde:

model_results <- list(free = free_fit, known_mu = known_mu_fit
                      , known_lambda = known_lambda_fit
                      , equal_rates = equal_rates_fit
                      , double_lambda = double_lambda_fit)
saveRDS( object = model_results, file = "div_fit_resuls.rds" )
