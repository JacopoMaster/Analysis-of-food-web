library(igraph)
library(intergraph)
library(ergm)
library(ggplot2)
library(sbm)

load("C:/Users/manet/OneDrive/Desktop/Network analysis/progetto/floridaFoodWebs.RData")
set.seed(2107)
p_w<-asIgraph(bayWet)
p_d<-asIgraph(bayDry)
vertex.attributes(p_w)
V(p_w)
p_w<-delete_vertices(p_w,123:128)
p_w<-delete_vertices(p_w,12)
p_d<-delete_vertices(p_d,123:128)
p_d<-delete_vertices(p_d,12)
plot(p_d,size=10,vertex.color="orange",vertex.size=8,vertex.label=vertex.attributes(p_d)$vertex.names,vertex.label.color="black",edge.arrow.size=0.8,edge.arrow.width=0.8, main="Graph of the bay biome in the wet season")
V(p)
vertex.attributes(p)
gsize(p)
summary(V(p)$bio.masses...baywet.dat)
vertex.attributes(p)$vertex.names[order(V(p)$bio.masses...baywet.dat,decreasing = T)]

rho_w<-edge_density(p_w)
rho_d<-edge_density(p_d)
reciprocity(p)
odd_dens<-rho/(1-rho)
tr<-transitivity(p)
odd_tr<-tr/(1-tr)
odds_ratio<-odd_tr/odd_dens
odds_ratio
assortativity(p, V(p)$bio.masses...baywet.dat,directed = T)



degree(p)
st_in_deg<-degree(p,normalized = T,mode = "in")
vertex.attributes(p)$vertex.names[order(st_in_deg,decreasing = T)]
ord = order(st_in_deg, decreasing = T)
V(p)$color = "orange"
V(p)$color[ord[1:7]] = "lightblue"
V(p)$color[ord[118:120]] = "lightgray"
plot(p, vertex.size = st_in_deg*40,  
     vertex.label.cex = 1, main = "Wet Season In-degree centrality",vertex.label.color="black")
#aironi, rapaci  coccodrilli

st_out_deg<-degree(p,normalized = T,mode = "out")
ord = order(st_out_deg, decreasing = T)
vertex.attributes(p)$vertex.names[order(st_out_deg,decreasing = T)]
V(p)$color = "orange"
V(p)$color[ord[1:3]] = "lightblue"
V(p)$color[ord[118:120]] = "lightgray"
plot(p, vertex.size = st_out_deg*30,  
     vertex.label.cex = 1, main = "Wet Season Out-degree centrality",vertex.label.color="black")
vertex.attributes(p)

clos<-closeness(p,normalized = T)
clos[is.na(clos)] <- 0
vertex.attributes(p)$vertex.names[order(clos,decreasing = T)]
ord = order(clos, decreasing = T)
V(p)$color = "orange"
V(p)$color[ord[1:24]] = "lightblue"
plot(p, vertex.size = clos*15,  
     vertex.label.cex = 1, main = "Wet Season Closeness centrality",vertex.label.color="black")
vertex.attributes(p)


betw<-betweenness(p, directed = T, normalized = T)
vertex.attributes(p)$vertex.names[order(betw,decreasing = T)]
ord = order(betw, decreasing = T)
V(p)$color = "orange"
V(p)$color[ord[1:5]] = "lightblue"
plot(p, vertex.size = betw*700,  
     vertex.label.cex = 1, main = "Wet Season Betweenness centrality",vertex.label.color="black")


eig<-eigen_centrality(p, scale = T,directed = T)$vector
vertex.attributes(p)$vertex.names[order(eig,decreasing = T)]
ord = order(eig, decreasing = T)
V(p)$color = "orange"
V(p)$color[ord[1:3]] = "lightblue"
V(p)$color[ord[118:120]] = "lightgray"

plot(p, vertex.size = eig*30,  edge.arrow.size = 0.5,
     vertex.label.cex = 0.8, main = "Wet Season Eigenvector centrality")
vertex.attributes(p)
p[,119]


centr_degree(p, loops = F,mode="in")$centralization
centr_degree(p, loops = F,mode="out")$centralization
centr_clo(p)$centralization
centr_betw(p, directed = T)$centralization




rete_w<-asNetwork(p_w)
rete_d<-asNetwork(p_d)
rete_w
rete_d
mod0_w = ergm(rete_w ~ edges)
mod0_d= ergm(rete_d ~ edges)
summary(mod0_w)
summary(mod0_d)
exp(coef(mod0_w))/(1+exp(coef(mod0_w)))
exp(coef(mod0_d))/(1+exp(coef(mod0_d)))


sim0_d = simulate(mod0_d, nsim = 100, verbose = TRUE, seed = 1)
sim0_w = simulate(mod0_w, nsim = 100, verbose = TRUE, seed = 1)

fnc = function(xx){
  ig = asIgraph(xx)
  tr = transitivity(ig)
  ideg = sd(degree(ig, mode = "in"))
  odeg = sd(degree(ig, mode = "out"))
  return(c(tr, ideg, odeg))
}

#Sim HSRGM dry season
null.distr = matrix(,100,3)
for(b in 1:100){
  null.distr[b,]  = fnc(sim0_d[[b]])
}
dev.new()
par(mfrow = c(3,1))
hist(unlist(null.distr[,1]), xlab = "transitivity",xlim = c(0.2,0.28)); abline(v = transitivity(p_d), col = "red")
hist(unlist(null.distr[,2]), xlab = "Standard deviation in-degree",xlim = c(2,13)); abline(v = sd(degree(p_d, mode = "in")), col = "red")
hist(unlist(null.distr[,3]), xlab = "Standard deviation out-degree", ,xlim = c(2,15)); abline(v = sd(degree(p_d, mode = "out")), col = "red")

null.distr = matrix(,100,3)
for(b in 1:100){
  null.distr[b,]  = fnc(sim0_d[[b]])
}
dev.new()
par(mfrow = c(3,1))
hist(unlist(null.distr[,1]), xlab = "transitivity",xlim = c(0.2,0.28)); abline(v = transitivity(p_d), col = "red")
hist(unlist(null.distr[,2]), xlab = "Standard deviation in-degree",xlim = c(2,13)); abline(v = sd(degree(p_d, mode = "in")), col = "red")
hist(unlist(null.distr[,3]), xlab = "Standard deviation out-degree", ,xlim = c(2,15)); abline(v = sd(degree(p_d, mode = "out")), col = "red")

par(mfrow = c(3,1))
hist(unlist(null.distr[,1]), xlab = "transitivity",xlim = c(0.2,0.28)); abline(v = transitivity(p_d), col = "red")
hist(unlist(null.distr[,2]), xlab = "in-degree",xlim = c(2,13)); abline(v = sd(degree(p_d, mode = "in")), col = "red")
hist(unlist(null.distr[,3]), xlab = "out-degree", ,xlim = c(2,15)); abline(v = sd(degree(p_d, mode = "out")), col = "red")

null.distr = matrix(,100,3)
for(b in 1:100){
  null.distr[b,]  = fnc(sim0_w[[b]])
}
par(mfrow = c(3,1))
hist(unlist(null.distr[,1]), xlab = "transitivity",xlim = c(0.2,0.28)); abline(v = transitivity(p_w), col = "red")
hist(unlist(null.distr[,2]), xlab = "in-degree",xlim = c(2,13)); abline(v = sd(degree(p_w, mode = "in")), col = "red")
hist(unlist(null.distr[,3]), xlab = "out-degree", ,xlim = c(2,15)); abline(v = sd(degree(p_w, mode = "out")), col = "red")


mod_01_d = ergm(rete_d ~ edges + nodecov("bio.masses...baydry.dat") + 
               absdiff("bio.masses...baydry.dat"), control = control.ergm(seed = 1))

mod_01_w = ergm(rete_w ~ edges + nodecov("bio.masses...baywet.dat") + 
                absdiff("bio.masses...baywet.dat"), control = control.ergm(seed = 1))

summary(mod_01_d)
summary(mod_01_w)


mod1_d = ergm(rete_d ~ edges + receiver+sender)
summary(mod1_d)



mod2_d = ergm(rete_d ~ edges  + mutual, 
            control = control.ergm(seed = 1))
mod2_w = ergm(rete_w ~ edges  + mutual, 
              control = control.ergm(seed = 1))

summary(mod2_d)
summary(mod2_w)

sim2_d = simulate(mod2_d, nsim = 100, verbose = TRUE, seed = 1)
sim2_w = simulate(mod2_w, nsim = 100, verbose = TRUE, seed = 1)
null.distr = matrix(,100,3)
for(b in 1:100){
  null.distr[b,]  = fnc(sim2_d[[b]])
}
par(mfrow = c(3,1))
hist(unlist(null.distr[,1]), xlab = "transitivity",xlim = c(0.22,0.28)); abline(v = transitivity(p_d), col = "red")
hist(unlist(null.distr[,2]), xlab = "in-degree",xlim = c(2,13)); abline(v = sd(degree(p_d, mode = "in")), col = "red")
hist(unlist(null.distr[,3]), xlab = "out-degree", ,xlim = c(2,15)); abline(v = sd(degree(p_d, mode = "out")), col = "red")

for(b in 1:100){
  null.distr[b,]  = fnc(sim2_w[[b]])
}
par(mfrow = c(3,1))
hist(unlist(null.distr[,1]), xlab = "transitivity",xlim = c(0.22,0.28)); abline(v = transitivity(p_w), col = "red")
hist(unlist(null.distr[,2]), xlab = "in-degree",xlim = c(2,13)); abline(v = sd(degree(p_w, mode = "in")), col = "red")
hist(unlist(null.distr[,3]), xlab = "out-degree", xlim = c(2,15)); abline(v = sd(degree(p_w, mode = "out")), col = "red")

BIC(mod0,mod_01,mod2)

mod3_d = ergm(rete_d ~ edges + mutual + istar(2) + ostar(2) + triangle, 
                   control = control.ergm(seed = 1))

mod3_w = ergm(rete_w ~ edges + mutual + istar(2) + ostar(2) + triangle, 
              control = control.ergm(seed = 1))

mod4_w = ergm(rete_w ~ edges + mutual + istar(2) + ostar(2), 
              control = control.ergm(seed = 1))

mod5_d = ergm(rete_d ~ edges + mutual+ gwidegree(decay = 1, fixed = TRUE)
             + gwodegree(decay = 1, fixed = TRUE), 
             control = control.ergm(seed = 1))

mod6_d = ergm(rete_d ~ edges + mutual  + 
               gwesp(decay = 1, fixed = T) + 
               gwdsp(decay = 1, fixed = T), control = control.ergm(seed=1))
summary(mod6_d)

sim6_d = simulate(mod6_d, nsim = 1000, verbose = TRUE, seed = 1)
null.distr = matrix(,1000,3)
for(b in 1:1000){
  null.distr[b,]  = fnc(sim6_d[[b]])
}
par(mfrow = c(3,1))
hist(unlist(null.distr[,1]), xlab = "transitivity"); abline(v = transitivity(p_d), col = "red")
hist(unlist(null.distr[,2]), xlab = "in-degree"); abline(v = sd(degree(p_d, mode = "in")), col = "red")
hist(unlist(null.distr[,3]), xlab = "out-degree", xlim = c(10,15)); abline(v = sd(degree(p_d, mode = "out")), col = "red")
mean(transitivity(p_d) <= unlist(null.distr[,1]))
mean(sd(degree(p_d, mode = "in")) >= unlist(null.distr[,2]))
mean(sd(degree(p_d, mode = "out")) <= unlist(null.distr[,3]))

mod6_w = ergm(rete_w ~ edges + mutual  + 
                gwesp(decay = 1, fixed = T) + 
                gwdsp(decay = 1, fixed = T), control = control.ergm(seed=1))
summary(mod6_w)

sim6_w = simulate(mod6_w, nsim = 1000, verbose = TRUE, seed = 1)
for(b in 1:1000){
  null.distr[b,]  = fnc(sim6_w[[b]])
}
par(mfrow = c(3,1))
hist(unlist(null.distr[,1]), xlab = "transitivity"); abline(v = transitivity(p_w), col = "red")
hist(unlist(null.distr[,2]), xlab = "in-degree"); abline(v = sd(degree(p_w, mode = "in")), col = "red")
hist(unlist(null.distr[,3]), xlab = "out-degree" ,xlim = c(10,15)); abline(v = sd(degree(p_w, mode = "out")), col = "red")
mean(transitivity(p_w) <= unlist(null.distr[,1]))
mean(sd(degree(p_w, mode = "in")) >= unlist(null.distr[,2]))
mean(sd(degree(p_w, mode = "out")) <= unlist(null.distr[,3]))

mod7 = ergm(rete ~ edges  + mutual
               gwdsp(decay = 1,fixed = T), control = control.ergm(seed=1))
summary(mod14)



mod7_d = ergm(rete_d ~ edges  + 
               gwesp(decay = 1, fixed = T) + 
               gwdsp(decay = 1, fixed = T), control = control.ergm(seed=1))

mod8_d = ergm(rete_d ~ edges  + 
                gwesp(decay = 1, fixed = T) + 
                gwdsp(decay = 1, fixed = T), control = control.ergm(seed=1))


mod9_w = ergm(rete_w ~ edges + mutual  + gwidegree(decay = 1, fixed = TRUE)+ 
                gwesp(decay = 1, fixed = T) + 
                gwdsp(decay = 1, fixed = T), control = control.ergm(seed=1))

summary(mod9_w)
BIC(mod6_w,mod9_w)
sim9_w = simulate(mod9_w, nsim = 1000, verbose = TRUE, seed = 1)
null.distr = matrix(,1000,3)
for(b in 1:1000){
  null.distr[b,]  = fnc(sim9_w[[b]])
}
par(mfrow = c(3,1))
hist(unlist(null.distr[,1]), xlab = "transitivity"); abline(v = transitivity(p_w), col = "red")
hist(unlist(null.distr[,2]), xlab = "in-degree"); abline(v = sd(degree(p_w, mode = "in")), col = "red")
hist(unlist(null.distr[,3]), xlab = "out-degree" ,xlim = c(10,15)); abline(v = sd(degree(p_w, mode = "out")), col = "red")
mean(transitivity(p_w) <= unlist(null.distr[,1]))
mean(sd(degree(p_w, mode = "in")) >= unlist(null.distr[,2]))
mean(sd(degree(p_w, mode = "out")) <= unlist(null.distr[,3]))


BIC(mod0_d,mod_01_d,mod2_d,mod6_d)

Y_w<-as_adjacency_matrix(p_w)
Y_w<-as.matrix(Y_w)
Y_d<-as_adjacency_matrix(p_d)
Y_d<-as.matrix(Y_d)
plotMyMatrix(Y_w, dimLabels = list(row = 'organism', col = 'organism'))
plotMyMatrix(Y_d, dimLabels = list(row = 'organism', col = 'organism'))
sbm1_w = estimateSimpleSBM(Y_w, "bernoulli", dimLabels = 'organism', 
                         estimOptions = list(verbosity = 1))
sbm1_d = estimateSimpleSBM(Y_d, "bernoulli", dimLabels = 'organism', 
                           estimOptions = list(verbosity = 1))
sbm1_w$nbBlocks
sbm1_d$nbBlocks
# prior block probabilities
sbm1_w$blockProp #alpha prob a priori
sbm1_d$blockProp
# connectivity parameters
round(sbm1_w$connectParam$mean,3)
round(sbm1_d$connectParam$mean,3)
# a clearer representation
# ------------------------
plot(sbm1_d, type = "data")
plot(sbm1_w, type = "data")
# nodes are ordered wrt to the block they belong to and blocks are highlighted
plot(sbm1_d, type = "expected")
plot(sbm1_w, type = "expected")
# fitted connection probabilities
plot(sbm1_d, type = "meso")
plot(sbm1_w, type = "meso")

vertex.attributes(p_d)$vertex.names[sbm1_d$memberships==4]
vertex.attributes(p_w)$vertex.names[sbm1_d$memberships==4]
vertex.attributes(p_d)$vertex.names[sbm1_d$memberships==1]
vertex.attributes(p_w)$vertex.names[sbm1_w$memberships==1]
vertex.attributes(p_d)$vertex.names[sbm1_d$memberships==2]
vertex.attributes(p_w)$vertex.names[sbm1_w$memberships==2]
vertex.attributes(p_d)$vertex.names[sbm1_d$memberships==8]
vertex.attributes(p_w)$vertex.names[sbm1_w$memberships==8]
sbm1_w$memberships
plot(p_d, vertex.color = sbm1_d$memberships, vertex.size=10)
legend(x="topright", legend=c("1", "2","3","4","5","6","7","8"),
       fill = c("orange","lightblue","darkgreen","yellow","blue","red","magenta","grey"),
       title="Block", text.font=4, bg='white')

plot(p_w, vertex.color = sbm1_w$memberships, vertex.size=10)
legend(x="topright", legend=c("1", "2","3","4","5","6","7","8"),
       fill = c("orange","lightblue","darkgreen","yellow","blue","red","magenta","grey"),
       title="Block", text.font=4, bg='white')


sbm1_d$probMembership
sbm1_w$probMemberships
