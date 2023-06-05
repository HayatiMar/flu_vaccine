library(devtools)
install_git("http://github.com/bdearlove/treeImbalance.git")
library(treeImbalance)
#function to compute DR
diversificationRate <- function(tree) {
  paths <- nodepath(tree)
  n_paths <- length(paths)
  DR <- numeric(n_paths)
  for (i in 1:n_paths) {
    path <- paths[[i]]
    n <- length(path) - 1
    # Subset the edges based on the row indices
    edges <- which(tree$edge[, 1] %in% path[-(n + 1)] & tree$edge[, 2] %in% path[-1])
    edge_lengths <- tree$edge.length[edges]
    weights <- 1 / (2 ^ (n:1 - 1))
    s <- sum(edge_lengths*weights)
    DR[i] <- 1 / s
  }
  return(DR)
}
#test on a random tree
set.seed(17)
tr1= rtree(500)
dd = get_pairwise_distances(tr1,c(1:length(tr1$tip.label)),c(length(tr1$tip.label):1))
tree.lbi<-lbi(tr1,tau = 0.0625*mean(dd))
LBI_tr1 = tree.lbi[1:length(tr1$tip.label)]
DR_tr1 = diversificationRate(tr1)
plot(LBI_tr1,DR_tr1)

#compute LBI and DR for HA tree.
load("~/2020/flutree2020-2.Rdata")
#I got error when running lbi on the tree due to the size of the tree so I applied it on a subtree.
tr = extract.clade(tree,76876)
#mean pairwise distance in HA is 23.7
LBI_HA = lbi(tr ,tau = 0.0625*23.7)
LBI_HA_tip = LBI_HA[1:length(tr$tip.label)]
DR_HA_tip = diversificationRate(tr)

plot(LBI_HA_tip,DR_HA_tip)
ind = which(DR_HA_tip>= 100*median(DR_HA_tip))
LBI_HA_tip = LBI_HA_tip[-ind]
DR_HA_tip = DR_HA_tip[-ind]
plot(LBI_HA_tip,DR_HA_tip)
#compute LBI and DR for NA tree.

load("~/2020/flutree_NA_2020-2.Rdata")
tr = extract.clade(tree,59400)
#mean pairwise distance in NA is 24.3
LBI_NA = lbi(tr ,tau = 0.0625*24.3)
LBI_NA_tip = LBI_NA[1:length(tr$tip.label)]
DR_NA_tip = diversificationRate(tr)

plot(LBI_NA_tip,DR_NA_tip)

