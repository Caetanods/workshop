## Alguns plots para os slides da introdução.

set.seed(1234) # Para reproducibilidade

## Um modelo de correlação.
## Exemplo de aplicação de um modelo linear.
a <- rnorm(n = 50, mean = 3, sd = 1)
c <- a + runif(n = 50, min = -2, max = 2)
plot(a~c, ylim = c(0,6), xlim = c(0,6), xlab = "Var A", ylab = "Var B", pch = 20)
abline(lm(a~c), col = "red", lty = 3, lwd = 3)

## Estimando os parâmetros do modelo.
lm(a~c)

## Exemplos de outros ajustes do model.
## Aqui vamos plotar diferentes valores para os parâmetros deste modelo e ver
##     como a linha de regressão ajusta nos dados.

## Relação negativa:
a.1 = (-0.354 * c) + 1.615
plot(a~c, ylim = c(0,6), xlim = c(0,6), xlab = "Var A", ylab = "Var B", pch = 20)
abline(lm(a.1~c), col = "blue", lty = 3, lwd = 3)
lm(a.1~c)

## Igual a zero:
a.2 = (0 * c) + 1.615
plot(a~c, ylim = c(0,6), xlim = c(0,6), xlab = "Var A", ylab = "Var B", pch = 20)
abline(lm(a.2~c), col = "blue", lty = 3, lwd = 3)
lm(a.2~c)
