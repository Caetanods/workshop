## Podemos usar a razão de verossimilhança (likelihood ratio) para testar se temos suporte
##      para aceitar o modelo mais complexo ao invés do modelo reduzido.
## A hipótese nula de um teste de likelihood ratio é que o modelo mais simples é o melhor modelo.

mod_res <- readRDS( "div_fit_resuls.rds" )
names( mod_res )

## ######################################
## Teste se o valor de mu está correto:
LR <- -2 * (mod_res$known_mu$value - mod_res$free$value)

## Verifique o p valor usando uma distribuição chi-square com degrees of freedom iguais a
##     diferença de numero de parametros entre o modelo nested e o modelo free.
pchisq(q = LR, df = 1, lower.tail = FALSE)
## Temos suporte para o valor de extinção estimado previamente?

## ######################################
## Teste se o valor de lambda está correto:
LR <- -2 * (mod_res$known_lambda$value - mod_res$free$value)
pchisq(q = LR, df = 1, lower.tail = FALSE)
## Temos suporte para o valor de especiação estimado previamente?

## ######################################
## Teste se mu e lambda tem valores diferentes:
LR <- -2 * (mod_res$equal_rates$value - mod_res$free$value)
pchisq(q = LR, df = 1, lower.tail = FALSE)
## Temos suporte para valores diferentes de lambda e mu?

## Podemos fazer uma representação gráfica do teste:
hist(x = rchisq(n = 10000, df = 1), main = "Chi-squared distribution", freq = FALSE
     , xlab = "LR test values", breaks = 50, xlim = c(0, 30))
abline(v = LR, col = "red", lwd = 2)
text(x = 20, y = 0.9, labels = paste0("LR test: free vs. equal rates \n p value = ", format(pchisq(q = LR, df = 1, lower.tail = FALSE))))

## ######################################
## Teste se lambda é diferente de 2*mu:
LR <- -2 * (mod_res$double_lambda$value - mod_res$free$value)
pchisq(q = LR, df = 1, lower.tail = FALSE)
## Temos suporte para um cenário em a especiação é duas vezes mais rápida que a extinção?

## PERGUNTE AOS UNIVERSITÁRIOS:
## Vocês conseguem listar algumas limitações do likelihood ratio test depois de ver esses exemplos?
