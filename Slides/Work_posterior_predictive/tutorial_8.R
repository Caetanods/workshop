## Nesse tutorial vamos efetuar uma analise de posterior predictive check.

## O objetivo da posterior predictive check e verificar se o modelo descreve bem
##   o bastante as caracteristicas dos dados que estamos interessados.

## Usamos modelos para descrever padroes que encontramos nas nossas observacoes e
##   que estamos interessados. Por exemplo, o ajuste de um modelo normal para os
##   dados assume que a media e o desvio padrao seriam duas quantidades interessantes
##   para descrever nossos dados, pois, sao os parametros deste modelo.

## No entanto, como saber se o modelo captura (descreve) os padroes nos nossos dados
##   que estamos interessados? A inferencia Bayesiana por meio do MCMC nos permite usar
##   uma tecnica muito interessante para responder essa questao.

## Podemos utilizar a distribuicao posterior de valores dos parametros para gerar
##   novos dados que parecem com as nossas observacoes. Se o modelo e os parametros
##   estimados sao um bom fit para as nossas observacoes, entao os dados simulados
##   com base na posterior serao semelhantes aos dados reais.

## Vamos utilizar novamente os resultados da analise de MCMC do modelo BM feito nos
##   tutoriais anteriores.

set.seed(1234)

## Carregando os resultados da analise. Dessa vez precisamos tambem da arvore
##   filogenetica e dos dados observados.
load(file = "./post.data.tree.RData")

## Note que a posterior tem um formate de tabela (ou data.frame):
head(res.mcmc)

## Cada uma das linhas desta tabela é uma combinação de parametros avaliada e aceita
##   pela cadeia do MCMC. Primeiro vamos tirar a regiao do burn-in.
## Da mesma forma que no tutorial anterior vamos retirar 2.000 geracoes.

post <- res.mcmc[-(1:2000),]

## O plot da posterior para os parametros:
par(mfrow=c(1,2))
hist(x = post$root, xlab = "root", freq = FALSE)
hist(x = post$rate, xlab = "rate", freq = FALSE)
par(mfrow=c(1,1))

## Primeiro passo e selecionar um numero de amostras da posterior para fazer as
##   simulacoes. Vamos fazer 1000 simulações.
## Note que nosso objetivo e utilizar a combinacao de valores de parametros das
##   geracoes do MCMC, isso e importante pois essas sao as combinacoes que foram
##   avaliadas e aceitas pelo MCMC.

## Numero de linhas do data.frame:
dim(post)[1]

## O ideal seria fazer uma simulação de novos dados para cada uma das geracoes do MCMC
##   dessa forma os novos dados teriam uma densidade exatamente proporcional a
##   distribuicao dos parametros que constituem nossa posterior.
## No entanto, nesse exemplo simples teriamos que fazer 8000 simulacoes de novos dados.
## Vamos reduzir essa quantidade aplicando um 'thinning' na cadeia e deixando somente
##   uma entre cada 8 geracoes.

thin.post <- post[seq(from = 1, to = 8000, by = 8),]

## Vamos plotar novamente a posterior e comparar as distribuicoes.
par(mfrow=c(2,2))
hist(x = post$root, xlab = "root", freq = FALSE, main = "Complete")
hist(x = thin.post$root, xlab = "root", freq = FALSE, main = "Thinned")
hist(x = post$rate, xlab = "rate", freq = FALSE, main = "Complete")
hist(x = thin.post$rate, xlab = "rate", freq = FALSE, main = "Thinned")
par(mfrow=c(1,1))

## Embora o numero de amostras seja 8 vezes menor, podemos ver que o formato da
##   distribuicao posterior nao se altera.

## Vamos agora fazer uma simulacao para cada uma das geracoes.
## Pata tal vamos usar a mesma filogenia e comprimento de ramos que usamos para a
##   analise de MCMC que gerou a posterior do parametros. Lembre-se que para modelos
##   filogeneticos comparativos a filogenia e um dos parametros.

## Precisamos de um for loop para fazer o uso dos parametros na tabela.
## Vamos usar uma lista para guardar os resultados das simulacoes:
sims <- list()
for(i in 1:1000){
  sims[[i]] <- sim.char(phy = phy, par = thin.post$rate[i]
                        , nsim = 1, model = "BM", root = thin.post$root[i])[,,1]
}

## Um truque no R para transformar essa lista de vetores em uma tabela:
sim.table <- do.call(what = rbind, args = sims)
## Isso gera uma grande tabela com 1000 linhas e 100 colunas, as colunas são os
##   valores simulados para cada uma das especies.
dim(sim.table)

## Temos nossa distribuicao posterior, temos os dados simulados de acordo com os
##   parametros da posterior e tambem temos os nosso dados observados.
## Agora como podemos verificar se o modelo de BM estimado e um bom modelo para
##   explicar o padrao das observacoes?

## Precisamos de alguma medida de distancia. A ideia e saber se podemos gerar dados
##   parecidos com as nossas observacoes se simularmos os dados usando o mesmo
##   modelo e os parametros da posterior.

## Esse e um ponto bastante importante da posterior predictive check: Quais
##   estatisticas (ou summary statistics) vou utilizar para comparar as observacoes
##   reais com os dados simulados?

## As summary statistics devem estar relacionadas com as caracteristicas dos dados
##   que sao importantes para o seu estudo. Neste caso queremos saber sobre a
##   variancia dos valores das especies.

## Calculamos variancia para os dados observados (var.obs) e para as simulacoes
##   (var.sim):
var.obs <- var(data)
var.sim <- apply(sim.table, MARGIN = 1, FUN = var)

## Vamos verificar se podemos gerar um valor de variancia semelhante ao observado
##   nas simulacoes:

hist(x = var.sim, xlab = "Variancia simulada da posterior", freq = FALSE
     , col = "gray", border = "white")
abline(v = var.obs, lty = 2, col = "red", lwd = 4)

## Visualmente parece bom, certo? O valor observado da variancia (linha vermelha)
##   esta em uma regiao de alta densidade de valores simulados na posterior
##   predictive.

## Como temos uma distribuicao de valores de uma summary statistics e um valor
##   observado, podemos calcular um p valor. Esse p valor e a probabilidade de
##   gerar dados tao ou mais extremos que a quantidade que observamos.

## Esse p valor e semelhante ao calculado em analises usando o metodo de Monte Carlo.
## O metodo de Monte Carlo usa simulacoes para gerar um modelo nulo, depois
##   calculamos o valor de p baseado na proporcao de simulacoes maiores (ou menores)
##   que o valor observado.

## Para fazer esse calculo contamos o numero de observacoes maiores do que 'var.obs'.
## Depois dividimos pelo numero total de observacoes da posterior predictive
sum( var.sim > var.obs ) / 1000 ## 0.256

## O p valor de 0.256 nos indica que e plausivel pensar que o modelo estimado tem um
##   bom fit para nosso dados, pois a chance de gerar dados semelhantes aos dados
##   observados e relativamente grande.