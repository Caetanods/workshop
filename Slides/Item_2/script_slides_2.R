## Figures for the Item_2 presentation.

## Seed para o script:
set.seed(1234)

## Gerar uma distribuição normal e o plot.
nr <- rnorm(n = 2000, mean = 10, sd = 2)
hist(nr, main="", xlab="Data", col = "gray", border = "white")

## Log da verossimilhança de diferentes valores da média:
## Note que os valores são muito pequenos para calcular no espaço real.
## Para trabalhar com verossimilhanças é importante transformar para o log.
lik.norm <- function(dt, mu, var){
  ## dt = dados (vector).
  ## mu = média.
  ## var = variancia.
  n <- length(dt)
  first <- log(2)*(-n/2) + log(pi)*(-n/2) + log(var)*(-n/2)
  second <- -1/(2*var) * sum( (dt - mu)^2 )
  lik <- first + second
  return(lik)
}

## Podemos usar a função de loglik para traçar a superfície de verossimilhança.
## Note que esse gráfico mostra somente o perfil de verossimilhança para a média
##     assumindo sigma = 2.
res <- vector()
for(i in 1:20){
  res[i] <- lik.norm(dt=nr, mu=i, var=4)
}
plot(res~seq(1,20), xlab = "mean (mu)", ylab = "log(Likelihood)", pch = 16)

## Agora vamos calcular a superfície completa de likelihood.
## Nesse caso temos dois parâmetros e, portanto, precisamos adicionar
##    mais um eixo no gráfico. O que vai tornar o gráfico 3D.
## Para tal vamos utilizar a função 'persp'.

x <- 1:20 ## valores para mu
y <- seq(0, 10, length.out = 20) ## valores para sigma
## Vamos usar z para criar uma tabela com a combinação dos valores de
##    log likelihood.
## Os valores de x estão nas linhas e os de y nas colunas.
z <- matrix(data = NA, nrow = 20, ncol = 20)
for(i in 1:20){
  for(j in 1:20){
    z[i,j] <- lik.norm(dt=nr, mu=x[i], var=y[j])
  }
}

persp(x, y, z, theta = 35, phi = 35, xlab = "mu", ylab = "sigma^2", zlab = "log(likelihood)", border = "black", col = "grey", r = 4)
