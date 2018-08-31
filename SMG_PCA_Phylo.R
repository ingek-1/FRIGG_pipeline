library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
library("devtools")
library("hierfstat")

library("stats")
library("ade4")
library("phangorn")
library("hierfstat")
library("seqinr")


################ FUNCTION MAKING THE PCA ANALYSIS - OUTPUT ONE PDF WITH FIGURES #############################
PCA_analysis <- function(alignment_file) {
  # getting the alignment from fasta file
  Alignment <- read.alignment(alignment_file, "fasta")
  
  # getting the hfam name
  namefile <- strsplit(alignment_file,"/")
  namefile_short <- tail(namefile[[1]], n=1)
  name_split <- strsplit(namefile_short,"_")
  hfam <- name_split[[1]][1]
  print(hfam)
  Alignment <- read.alignment(alignment_file, "fasta")

  # making a genind object from the Alignment
  dat <- alignment2genind(Alignment)
  
  # making population based on the section
  temp <- indNames(dat)
  temp1 <- lapply(strsplit(temp,"_"), rev)
  section <- sapply(temp1, function(vec) vec[4])
  pop(dat) <- section
  
  # making pop based on clust
  clust_info <- sapply(temp1, function(vec) vec[1])
  
  # saving both clust and section info in the S4 data frame of dat - using strata
  strata_df <- data.frame(section, clust_info) # Create a data frame from the list.
  strata(dat) <- strata_df
  
  # making a matrix with the data for the pca
  x <- tab(dat, freq=TRUE, NA.method="mean")
  
  # doing the pca
  pca.x <- dudi.pca(x, center=TRUE, scale=FALSE, scannf=FALSE, nf=3)
  end <- length(pca.x$eig)
  
  # Making a pdf file with the resulting plots
  filename <- paste("PCA_Hfam_", hfam, ".pdf")
  filename <- gsub(" ", "", filename, fixed = TRUE)
  pdf(file = filename, paper="a4")
  par(mfrow=c(3,2), "oma"=c(0,0,5,0))

  col_flag <-c("#87CEEB", "#EEB422","#43CD80", "#CD6889")
  
  # draw the plots
  # PCA 1+2 section
  setPop(dat) <- ~section
  s.class(pca.x$li, fac=pop(dat),  col=transp(funky(22),.6), axesel=FALSE, cstar=0, cpoint=3, sub ="hfam test")
  add.scatter.eig(pca.x$eig[1:end],posi="bottomright", 3,1,2, ratio=.3)

  # PCA 1+2 clust
  setPop(dat) <- ~clust_info
  s.class(pca.x$li, fac=pop(dat),  col=transp(col_flag,.6), axesel=FALSE, cstar=0, cpoint=3)
  add.scatter.eig(pca.x$eig[1:end],posi="bottomright", 3,1,2, ratio=.3)
  
  # PCA 1+3 section
  setPop(dat) <- ~section
  s.class(pca.x$li, fac=pop(dat), xax=1, yax=3, col=transp(funky(22),.6), axesel=FALSE, cstar=0, cpoint=3)
  add.scatter.eig(pca.x$eig[1:end],posi="topleft", 3,1,3, ratio=.3)

  # PCA 1+3 clust
  setPop(dat) <- ~clust_info
  s.class(pca.x$li, fac=pop(dat), xax=1, yax=3, col=transp(col_flag,.6), axesel=FALSE, cstar=0, cpoint=3)
  add.scatter.eig(pca.x$eig[1:end],posi="topleft", 3,1,3, ratio=.3)
  
  # PCA 2+3 section
  setPop(dat) <- ~section
  s.class(pca.x$li, fac=pop(dat), xax=2, yax=3, col=transp(funky(22),.6), axesel=FALSE, cstar=0, cpoint=3)
  add.scatter.eig(pca.x$eig[1:end],posi="topleft", 3,2,3, ratio=.3)
  
  # PCA 2+3 clust
  setPop(dat) <- ~clust_info
  s.class(pca.x$li, fac=pop(dat), xax=2, yax=3, col=transp(col_flag,.6), axesel=FALSE, cstar=0, cpoint=3)
  add.scatter.eig(pca.x$eig[1:end],posi="topleft", 3,2,3, ratio=.3)
  
  # Adding title to the page
  mtext(text = paste("Hfam ", hfam), side = 3, line=2, outer = TRUE)
  # close the device to do the drawing
  dev.off()
  
}


################ FUNCTION MAKING THE PHYLO ANALYSIS - OUTPUT ONE PDF WITH FIGURES #############################
ML_phylo <- function(alignment_file) {
    # getting the alignment from fasta file
  Alignment <- read.alignment(alignment_file, "fasta")
  
  # getting the hfam number 
  namefile <- strsplit(alignment_file,"/")
  namefile_short <- tail(namefile[[1]], n=1)
  name_split <- strsplit(namefile_short,"_")
  hfam <- name_split[[1]][1]
  print(hfam)
  #################### ML Tree #######################
  # based on example from:
  # https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.pdf
  
  #phyAA <- read.phyDat(alignment_file, format="fasta", type="AA")
  phyAA <- as.phyDat(Alignment, type="AA")
  
  n_phyAA <- names(phyAA)

  # testing which model to choose
  print("Testing which model to use of JTT, LG and WAG")
  mt <- modelTest(phyAA, model=c("JTT", "LG", "WAG"), multicore=TRUE)
  
  env <- attr(mt, "env")
  print("Making initial tree")
  fitStart = eval(get(mt$Model[which.min(mt$BIC)], env), env)
  print("getting the maximum likelihood tree")
  fit = optim.pml(fitStart, rearrangement = "stochastic", optInv=TRUE, optGamma=TRUE)
  
  tree_output_initial <- plot(fit$tree, "phylogram")
  filename2 <- paste("Phylo_", hfam, "_1.nwk")
  filename2 <- gsub(" ", "", filename2, fixed = TRUE)
  write.tree(fit$tree, file = filename2)
  
  # Only using the neighbout joining as start og model JTT - so not testing which model fits best
  # dm  <- dist.ml(phyAA, model = "JTT")
  # treeNJ  <- NJ(dm)
  # fitNJ = pml(treeNJ, phyAA, model="JTT", k=4, inv=.2)
  # fit = optim.pml(fitNJ, rearrangement = "stochastic", optInv=TRUE, optGamma=TRUE)
  
  # calculating bootstrap values  
  print("calculating bootstrap")
  bs = bootstrap.pml(fit, bs=500, optNni=TRUE, multicore=TRUE)
  
  ############## Plotting ###############
  
  # making groups for coloring for section and clust_flag
  temp <- fit$tree$tip.label 
  temp1 <- lapply(strsplit(temp,"_"), rev)
  section <- sapply(temp1, function(vec) vec[4])
  clust_flag <- sapply(temp1, function(vec) vec[1])
  
  print("making pdf file")
  # making pdf file
  filename <- paste("Phylo_Hfam_", hfam, ".pdf")
  filename <- gsub(" ", "", filename, fixed = TRUE)
  pdf(file = filename, width=15 , height=23, pointsize=12) ### Change paper size here for phylo tree pdf
  
  col_flag <-c("#87CEEB", "#EEB422","#43CD80", "#CD6889")
  names(col_flag) <- c("0", "Clust", "outsideSC", "StrictClust")
  col_flag_all <- col_flag[clust_flag]
  
  # Plotting with legend color of section
  plotBS(fit$tree, bs, p = 50, "phylogram") 
  add.scale.bar()
  text1 = paste("Maximum Likelihood Tree Hfam ", hfam, " section")
  title(text1)
  tiplabels(frame = "none", pch=20, col=transp(fac2col(section),0.7), cex=2, fg="transparent")
  legend("bottomleft", leg =unique(section), fill= fac2col(unique(section) ) )
  
  # Plotting with legend color of clust_flag   
  plotBS(fit$tree, bs, p = 50, type="p")
  add.scale.bar()
  text2 = paste("Maximum Likelihood Tree Hfam ", hfam, " clust flag")
  title(text2)
  tiplabels(frame = "none", pch=20, col=col_flag_all, cex=2, fg="transparent")
  legend("bottomleft", leg =unique(clust_flag), fill= unique(col_flag_all) )
  
  dev.off()
  
  print("making nwk file")
  # write the tree into nwk
  tree_output <- plotBS(fit$tree, bs, p = 50, "phylogram")
  filename2 <- paste("Phylo_Hfam_", hfam, ".nwk")
  filename2 <- gsub(" ", "", filename2, fixed = TRUE)
  write.tree(tree_output, file = filename2)
  
}



###### For fasta_hfam_final
# set working directory
setwd("/path/to/the/outputDirectory/fasta_best_Pipeline_aligned/")

# get the fasta filenames of the aligned protein families
filenames <- list.files(pattern="*.fasta$", full.names=TRUE)

# run the PCA on all the protein families
PCA <- lapply(filenames, PCA_analysis)
# run the ML_phylo to create phylogenetic trees with 500 bootstrap - output in pdf and nwk format
test2 <- lapply(filenames, ML_phylo)
