## Nesse tutorial vamos calcular o AIC de alguns modelos.
## Vamos usar um modelo de evolução de caracteres discretos (o modelo MK).

library( corHMM )

## Importar os dados:
male_care <- readRDS("male_care.rds")
female_care <- readRDS("female_care.rds")
female_prec <- readRDS("female_prec.rds")

## Plot da filogenia:
## NOTA: Esses são dados empiricos com um caso tipico. Falta de data nos nós e presença de politomias.
plot( ladderize(female_care[[1]]), show.tip.label = FALSE); axisPhylo()

## Prepare os dados para a função:
female_mat <- data.frame(names(female_care[[2]]), unname(female_care[[2]]))
head( female_mat )

## ###################################################
## Cuidado maternal em opiliões:

## ER: Taxas iguais de transição de e para cada um dos estados:
female_ER <- corHMM(phy = female_care[[1]], data = female_mat, rate.cat = 1
                  , model = "ER", n.cores = 2)
female_ER
## Podemos fazer o plot do modelo MK:
plotMKmodel( female_ER )

## ARD: Taxas diferentes para cada uma das transições:
female_ARD <- corHMM(phy = female_care[[1]], data = female_mat, rate.cat = 1
                     , model = "ARD", n.cores = 2)
female_ARD
plotMKmodel( female_ARD )

## Agora podemos comparar esses modelos:
female_ARD$AIC
female_ER$AIC

## O AICc é interessante, pois ele corrige para o viés causado pelo tamanho da amostra. Na verdade, a recomendação é sempre utilizar o AICc com dados empiricos, visto que os resultados da performance estatística do AIC são convergentes no limite assintótico.
female_ARD$AICc
female_ER$AICc

## Calculamos o delta AICc:
aic_vec <- setNames(object = c(female_ARD$AICc, female_ER$AICc), nm = c("ARD","ER"))
aic_vec
aic_vec - min(aic_vec) ## Qual o melhor modelo?

## Vamos fazer um plot da evolução do cuidado maternal usando o melhor modelo:
plotRECON(phy = female_ARD$phy, likelihoods = female_ARD$states)

## Não da para ver nada, daí a gente pode fazer um plot direto em um pdf:
pdf(file = "Recon_female_care_ARD.pdf", width = 7, height = 20)
plotRECON(phy = female_ARD$phy, likelihoods = female_ARD$states)
dev.off()

## #################################################
## Precursor de cuidado maternal.
## Vamos usar um caso mais complexo, com mais estados:

## Temos 4 estados de precursores com frequencias distintas:
table( female_prec[[2]] )

female_prec_mat <- data.frame(names(female_prec[[2]]), unname(female_prec[[2]]))
head( female_prec_mat )

## ER: Taxas iguais de transição de e para cada um dos estados:
female_prec_ER <- corHMM(phy = female_prec[[1]], data = female_prec_mat, rate.cat = 1
                    , model = "ER", n.cores = 2)
female_prec_ER
## Podemos fazer o plot do modelo MK:
plotMKmodel( female_prec_ER )

## SYM: Taxas diferentes entre os estados, mas taxas de ganhos e perdas são iguais.
female_prec_SYM <- corHMM(phy = female_prec[[1]], data = female_prec_mat, rate.cat = 1
                          , model = "SYM", n.cores = 2)
female_prec_SYM
plotMKmodel( female_prec_SYM )

## ARD: Taxas diferentes para cada uma das transições:
female_prec_ARD <- corHMM(phy = female_prec[[1]], data = female_prec_mat, rate.cat = 1
                     , model = "ARD", n.cores = 2)
female_prec_ARD
plotMKmodel( female_prec_ARD )

## BIO: Aqui temos uma hipótese biológica sobre as transições. Elas deveriam ser ordenadas, pois cada estado na verdade é um grupo de estados distintos.

print( female_prec[[3]] ) ## Cada um desses numeros indica uma transição possivel.
female_prec_BIO <- corHMM(phy = female_prec[[1]], data = female_prec_mat
                          , rate.mat = female_prec[[3]], rate.cat = 1
                          , n.cores = 2)
female_prec_BIO
plotMKmodel( female_prec_BIO )

## BIO1: Temos ainda mais um modelo. No caso vamos seguir as mesmas restrições do modelo BIO, mas todas as taxas de transição são iguais.
mat_BIO1 <- female_prec[[3]]
mat_BIO1[mat_BIO1 > 1] <- 1
female_prec_BIO1 <- corHMM(phy = female_prec[[1]], data = female_prec_mat
                          , rate.mat = mat_BIO1, rate.cat = 1
                          , n.cores = 2)
female_prec_BIO1
plotMKmodel( female_prec_BIO1 )

## Podemos listar os AICc para todos os modelos do set:
female_prec_ARD$AICc
female_prec_ER$AICc
female_prec_BIO$AICc
female_prec_BIO1$AICc
female_prec_SYM$AICc

## Calculamos o delta AICc:
aic_vec <- setNames(object = c(female_prec_BIO1$AICc, female_prec_ARD$AICc, female_prec_ER$AICc, female_prec_BIO$AICc, female_prec_SYM$AICc), nm = c("BIO1", "ARD","ER","BIO", "SYM"))
aic_vec
aic_vec - min(aic_vec) ## Qual o melhor modelo?

## Vamos fazer um plot da evolução do cuidado maternal usando o melhor modelo:
pdf(file = "Recon_female_prec_BIO.pdf", width = 7, height = 20)
plotRECON(phy = female_prec_BIO$phy, likelihoods = female_prec_BIO$states)
dev.off()
