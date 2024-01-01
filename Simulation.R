#### This script can be used to reproduce the simulation plots

## setup
rm(list =ls())
library(igraph)
library(Matrix)
library(RColorBrewer)
source("SBP.R")

## generate empty graph
p <- 100
set.seed(123)
g <- make_lattice(c(10, 10))

plot(g,layout=layout.grid(g))

A <- as.matrix(get.adjacency(g,"both"))

d <- rowSums(A)

L <- t(t(A/sqrt(d))/sqrt(d))

U <- eigen(L)$vectors[,1:5]

n <- p
rownames(A) <- colnames(A) <- as.character(1:n)
g <- graph.adjacency(A,"undirected")
lo <- layout.grid(g)
#plot(g,layout=lo)



## generate separation enforced matrix 
R <- matrix(0,p,p)

R[1:10,11:20] <- 1
R[seq(8,98,by=10),seq(9,99,by=10)] <- 1
sum(R != 0)

R <- R+t(R)
sum(R != 0)

set.seed(123)
## lambda = 0
cols <- c(brewer.pal(8,"Set1"),brewer.pal(8,"Set2"),brewer.pal(8,"Dark2"))
fit1 <- SBP(U,R,K=15,max.iter=1,lambda=0)
V(g)$color <- cols[fit1$cluster]
V(g)$label.cex = 2
png(file="sim-lamba0.png",
    width=1200, height=1200,res=200)
par(mar = c(0,0,0,0)) 
plot(g,layout=lo,vertex.label=as.character(fit1$cluster),legend = F)
abline(h=-0.88, col="red", lwd = 2)
abline(v=0.67, col="red",lwd =2)
dev.off()


## lambda = 5
fit2 <- SBP(U,R,K=15,max.iter=20,lambda=5,tol=1e-5,nstart=10,trace=TRUE)
V(g)$color <- cols[fit2$cluster]
V(g)$label.cex = 2
png(file="sim-lamba5.png",
    width=1200, height=1200,res=200)
par(mar = c(0,0,0,0)) 
plot(g,layout=lo,vertex.label=as.character(fit2$cluster))
dev.off()


## lambda = 10
fit3 <- SBP(U,R,K=15,max.iter=20,lambda=10,tol=1e-5,nstart=10,trace=TRUE)
V(g)$color <- cols[fit3$cluster]
V(g)$label.cex = 2
png(file="sim-lamba10.png",
    width=1200, height=1200,res=200)
par(mar = c(0,0,0,0)) 
plot(g,layout=lo,vertex.label=as.character(fit3$cluster))
dev.off()



