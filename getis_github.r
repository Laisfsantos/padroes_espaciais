rm(list=ls()) # Limpa todas as variaveis
cat("\f") # Limpa o console
# chama as bibliotecas necessarias
library(rgdal) # Biblioteca para analise geoespacial -> sudo apt-get install libgdal-dev no terminal para instalar
library(spdep) # Biblioteca para criar matrizes de pesos espaciais -> sudo apt-get install -y libudunits2-0 libudunits2-dev no terminal para instalar
library(RColorBrewer)
library(classInt)

library(factoextra) # Carrega biblioteca para visualizar resultados de analise multivariada
library(cluster) # Carrega biblioteca para clusterizacao
library(fclust) # Carrega biblioteca para clusterizacao
library(ppclust)

library(readxl)# Importa biblioteca para ler arquivos .xlsx

#Funcao do GetisOrd

GetisOrd = function(X,w_cont)
{
n <- nrow(X); # numero de municipios
w <- diag(n)+w_cont;


G <- matrix(0,nrow(X),ncol(X));
zG <- matrix(0,nrow(X),ncol(X));
I <- matrix(0,nrow(X),ncol(X));

for(i in 1:nrow(X))
{
  for(j in 1:ncol(X))
  {
    G[i,j] <- sum(w[i,]*X[,j])/sum(X[,j]);
    xbar.ast <- mean(X[,j]);
    s.ast <- sqrt((sum(X[,j]^2))/n - xbar.ast^2);
    s2 <- sd(X[,j]); # *sqrt((n-1)/n))^2;
    zG[i,j] <- (sum(w[i,]*X[,j])-sum(w[i,])*xbar.ast)/(s2*sqrt(1/(n-1)*(n*sum(w[i,]^2) - sum(w[i,])^2)));
  }
}
return(zG)
}



X2006 <- read_excel("diretorio");
Xall <- data.matrix(X2006, rownames.force = NA); # Converte a lista para uma matriz
X <- Xall[,2:ncol(Xall)]; # Matriz de dados 
X <- cbind(X[,1]+X[,3], X[,2]+X[,4]);
X <-X*1.89
Variaveis <- names(X2006);
Variaveis <- Variaveis[2:ncol(Xall)];

MATOPIBA <- readOGR(dsn = "diretorio", layer = "matopiba_municipios2", encoding = "pt_BR.UTF-8")

## 1. Contiguidade
# poly2nb: Construção da lista de vizinhos
MATOPIBA_nbq <- poly2nb(MATOPIBA, queen=TRUE)
# 1.1 Construindo matriz de vizinhança
w_cont <-nb2mat(MATOPIBA_nbq, style="B");  # converte para matriz (style="B" binary e style="W" weighted)
## Calculando espaço de variáveis espaciais
zG <- GetisOrd(X,w_cont);



fviz_nbclust(zG, kmeans, method='silhouette')

# c <- 2;
J1 <- matrix(0,6,1);
J2 <- matrix(0,6,1);
R <- matrix(0,6,1);
NClust <- matrix(0,6,1);
i <- 1;
for(c in c(3))
{
   res.fcm <- fcm(zG, centers = c)

  lm.palette <- colorRampPalette(c("white","orange", "red"), space = "rgb")
  
  res.kmeans <- kmeans(zG, centers = res.fcm$v)
  J1[i] <- dunn(clusters = res.kmeans$cluster, Data=zG);
  J2[i] <- res.kmeans$betweenss;
  NClust[i]<-c;
  R[i] <- 1 - res.kmeans$tot.withinss/res.kmeans$totss;
  i <- i + 1;
}
res.kmeans$centers<-res.kmeans$centers[order(res.kmeans$centers[,1]),];
res.kmeans <- kmeans(zG, centers = res.kmeans$centers)

MATOPIBA$cluster <- factor(res.kmeans$cluster);
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
plot(MATOPIBA, col = res.kmeans$cluster + 1, axes = TRUE, xaxs="i", yaxs="i", cex.axis = 1.2)
legend(-44.,-6.5,
       legend = c("1", "2", "3"),
       fill = c("red", "green", "blue"),
       bty = "n",
       cex = 1.2,
       title = "Grupos")


plot(zG[,1:2], col = res.kmeans$cluster + 1, axes=TRUE, xaxt='n',  yaxt='n', xlab = "",  ylab = "", pch = 19, bty='l')

plot(X[,1:2], col = 1, axes=TRUE, xaxt='n',  yaxt='n', xlab = "",  ylab = "", pch = 19, bty='l')

par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
plot(MATOPIBA, col = res.kmeans$cluster + 1, axes = FALSE, xaxs="i", yaxs="i", cex.axis = 1.2)


