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
### Returns TRUE if tree2 is consistent with tree1, FALSE if not.
splits_compatible <- function(tree1, tree2)
{
   # stopifnot(is.binary(tree1))
    tree1 <- read.tree(text = tree1)
    tree1 <- unroot(tree1)
    tree2 <- read.tree(text = tree2)
    tree2 <- unroot(tree2)
    stopifnot(all(tree2$tip.label %in% tree1$tip.label))

    # prune tree1 if necessary
    dropme <- which(!(tree1$tip.label %in% tree2$tip.label))
    if (length(dropme) > 0)
        tree1 <- drop.tip(tree1, dropme)
    # order tree2 labels to match tree1
    trees <- .compressTipLabel(c(tree1, tree2))
    tree2 <- trees[[2]]
    stopifnot(all.equal(tree1$tip.label, tree2$tip.label))
    # get all splits on each tree
    spl1 <- as.splits(unroot(tree1))
    spl2 <- as.splits(unroot(tree2))
    stopifnot(all.equal(attr(spl1, "labels"), attr(spl2, "labels")))
    # ensure the splits are in the same orientation
    polarize_splits <- function(spl)
    {
        all_tips <- seq(length(attributes(spl)$labels))
        for (i in seq_along(spl))
            if (!(1 %in% spl[[i]]))
                spl[[i]] <- setdiff(all_tips, spl[[i]])
        return(spl)
    }
    spl1 <- polarize_splits(spl1)
    spl2 <- polarize_splits(spl2)
    # see if all the tree2 splits are in the tree1 splits
    ans <- sapply(seq_along(spl2), function(i) spl2[i] %in% spl1)
    false_splits = NULL 
    if (all(ans))
    {
        result = TRUE

    } else {
        result = FALSE
        false_splits = spl2[which(!ans)]
    }
    num_splits <- length(WithoutTrivialSplits(as.Splits(unroot(tree2))))
    return(list(istrue = as.character(result), numsplits = as.character(num_splits), falsesplits = length(false_splits)))
}

splits_compatible(tree1_newick, tree2_newick)
