## Nesse tutorial vamos calcular o AIC de alguns modelos.
## Vamos usar um modelo de evolução de caracteres discretos (o modelo MK).

library( corHMM )

## Importar os dados:
male_care <- readRDS("male_care.rds")
female_care <- readRDS("female_care.rds")
female_prec <- readRDS("female_prec.rds")

## #################################################
## Precursor de cuidado maternal.
## Temos 4 estados de precursores com frequencias distintas:
female_prec_mat <- data.frame(names(female_prec[[2]]), unname(female_prec[[2]]))

## Modelos:
## ER: Taxas iguais de transição de e para cada um dos estados:
female_prec_ER <- corHMM(phy = female_prec[[1]], data = female_prec_mat, rate.cat = 1
                    , model = "ER", n.cores = 2)
## SYM: Taxas diferentes entre os estados, mas taxas de ganhos e perdas são iguais.
female_prec_SYM <- corHMM(phy = female_prec[[1]], data = female_prec_mat, rate.cat = 1
                          , model = "SYM", n.cores = 2)
## ARD: Taxas diferentes para cada uma das transições:
female_prec_ARD <- corHMM(phy = female_prec[[1]], data = female_prec_mat, rate.cat = 1
                     , model = "ARD", n.cores = 2)
## BIO: Aqui temos uma hipótese biológica sobre as transições. Elas deveriam ser ordenadas, pois cada estado na verdade é um grupo de estados distintos.
female_prec_BIO <- corHMM(phy = female_prec[[1]], data = female_prec_mat
                          , rate.mat = female_prec[[3]], rate.cat = 1
                          , n.cores = 2)
## BIO1: Temos ainda mais um modelo. No caso vamos seguir as mesmas restrições do modelo BIO, mas todas as taxas de transição são iguais.
mat_BIO1 <- female_prec[[3]]
mat_BIO1[mat_BIO1 > 1] <- 1
female_prec_BIO1 <- corHMM(phy = female_prec[[1]], data = female_prec_mat
                          , rate.mat = mat_BIO1, rate.cat = 1
                          , n.cores = 2)

## Calculamos o delta AICc:
aic_vec <- setNames(object = c(female_prec_BIO1$AICc, female_prec_ARD$AICc, female_prec_ER$AICc, female_prec_BIO$AICc, female_prec_SYM$AICc), nm = c("BIO1", "ARD","ER","BIO", "SYM"))
aic_vec
delta_aic_vec <- aic_vec - min(aic_vec)
delta_aic_vec

## ##################################################################################
## AIC weights

## Agora podemos calcular o AICw:
rel_lik <- exp( -0.5 * delta_aic_vec ) ## Relative likelihood.

## Vamos primeiro criar um objecto para o AICw:
AICw <- vector(mode = "numeric", length = length( rel_lik ))
names( AICw ) <- names( rel_lik )
## Calculamos cada um dos AICw:
for( i in 1:length(rel_lik) ){
    AICw[i] <- rel_lik[i] / sum( rel_lik )
}

round(AICw, digits = 4)

## Modelo BIO tem suporte alto, com > 98% do peso dos modelos.

## NOTA: Usamos o AICc para calcular o peso dos modelos acima. Você pode repetir a operação usando o AIC para comparar os resultados.

## ###############################################################################
## Model averaging.

## Agora que temos o vetor de AICw, podemos usa-lo para efetuar a média entre esses modelos estimados. Para tal vamos simplesmente utilizar os valores dos parametros para cada modelo em função do seu AICw.

female_prec_BIO1$solution
female_prec_BIO$solution
female_prec_ARD$solution
female_prec_ER$solution
female_prec_SYM$solution

## Cada matrix acima é a taxa de transição entre os estados do modelo.
## Essas taxas são os parametros do modelo. Podemos usar o model averaging para calcular uma taxa média usando uma média ponderada em função do AICw dos modelos.

female_prec_BIO1_w <- female_prec_BIO1$solution * AICw["BIO1"]
female_prec_BIO_w <- female_prec_BIO$solution * AICw["BIO"]
female_prec_ARD_w <- female_prec_ARD$solution * AICw["ARD"]
female_prec_ER_w <- female_prec_ER$solution * AICw["ER"]
female_prec_SYM_w <- female_prec_SYM$solution * AICw["SYM"]

## Produzimos valores de parametro ajustados em função do AIC weight de cada modelo.
## Como a soma dos AICw é igual a 1, podemos criar um novo modelo que é o resultado da soma dos valores dos parametros multiplicados pelo AICw de cada modelo.

female_prec_model_ave <- female_prec_BIO1_w + female_prec_BIO_w + female_prec_ARD_w + female_prec_ER_w + female_prec_SYM_w

## Matriz de transição MK estimada usando model averaging:
female_prec_model_ave

## TEMOS UM PROBLEMA! Veja que acima temos várias transições que estão como NA!
## O problema é que usamos NA para excluir parametros de alguns desses modelos.
## Para conduzir o model averaging precisamos usar um grupo de modelos congruentes.

## #######################################################################
## Produzindo modelos congruentes.

## Nesse exemplo sabemos que alguns modelos estão com taxas de transição iguais a NA ao invés de 0.0 .
## Qual é o problema? -> O problema é que NA não é um valor válido de parametro. O NA é um conceito que usamos para informar que o parametro está fixo com valor 0.

## Solução:
female_prec_BIO1_w[is.na(female_prec_BIO1_w)] <- 0.0
diag( female_prec_BIO1_w ) <- NA
female_prec_BIO_w[is.na(female_prec_BIO_w)] <- 0.0
diag( female_prec_BIO_w ) <- NA
female_prec_ARD_w[is.na(female_prec_ARD_w)] <- 0.0
diag( female_prec_ARD_w ) <- NA
female_prec_ER_w[is.na(female_prec_ER_w)] <- 0.0
diag( female_prec_ER_w ) <- NA
female_prec_SYM_w[is.na(female_prec_SYM_w)] <- 0.0
diag( female_prec_SYM_w ) <- NA

## Verificando:
female_prec_BIO1_w
female_prec_BIO_w
female_prec_ARD_w
female_prec_ER_w
female_prec_SYM_w

## Ótimo! Agora podemos fazer o model averaging.
female_prec_model_ave <- female_prec_BIO1_w + female_prec_BIO_w + female_prec_ARD_w + female_prec_ER_w + female_prec_SYM_w

## Parametros estimados usando model averaging:
female_prec_model_ave
round(female_prec_model_ave, digits = 4)

## Qual a diferença entre a matriz estimada usando model averaging e a matrix do melhor modelo?
best_model <- female_prec_BIO$solution
best_model[is.na(best_model)] <- 0.0
diag( best_model ) <- NA

diff_par <- abs( female_prec_model_ave - best_model )
round( diff_par, digits = 4 )
