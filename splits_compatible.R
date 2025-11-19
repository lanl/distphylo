### Goal: Is the information in tree2 consistent with the information in tree1?
### Assume that tree1 is fully resolved/bifurcating, and that polytomies in tree2
### represent a lack of information.  And that both trees are unrooted.

### Approach: See whether all the splits in tree2 are also in tree1 (pruned to
### the tips in tree2).
args <- commandArgs(trailingOnly = TRUE)
tree1_newick <- args[1]
tree2_newick <- args[2]

library(TreeTools)
suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})

splits_compatible <- function(tree1, tree2)
{
    tree1 <- read.tree(text = tree1)
    tree1 <- unroot(tree1)
    tree2 <- read.tree(text = tree2)
    tree2 <- unroot(tree2)
    stopifnot(all(tree2$tip.label %in% tree1$tip.label))

    # prune tree1 if necessary
    dropme <- which(!(tree1$tip.label %in% tree2$tip.label))
    if (length(dropme) > 0)
        tree1 <- drop.tip(tree1, dropme)

    # compare splits
    x <- tree1$tip.label[1]
    spl1 <- PolarizeSplits(as.Splits(tree1), x)
    spl2 <- PolarizeSplits(as.Splits(tree2), x)
    result <- all(spl2 %in% spl1)
    num_splits <- length(WithoutTrivialSplits(as.Splits(unroot(tree2))))
    
    # Print the result and number of splits
    cat(as.character(result), as.character(num_splits), sep=" ")
}

splits_compatible(tree1_newick, tree2_newick)
library(ape)
library(phytools)
library(TreeDist)
# Function to resolve all polytomies in tree1 and check compatibility with tree2
check_all_resolutions <- function(tree1, tree2) {
  # Read the trees
  #tree1 <- read.tree(text = tree1)
  #tree2 <- read.tree(text = tree2)
  
  # Get all possible fully resolved versions of tree1 using resolveAllNodes
  resolutions <- resolveAllNodes(tree1)
  
  # Compare each resolution with tree2
  for (resolved_tree in resolutions) {
    result <- splits_compatible(write.tree(resolved_tree), tree2)
    if (result == TRUE) {
      return(TRUE)  # If any resolution is compatible, return TRUE
    }
  }
  
  return(FALSE)  # If no resolution is compatible, return FALSE
}

