### Calculate summary statistics from 2D-mSFS
### Number of row is number of Pop2 + 1
### Number of column is number of Pop1 + 1
### At first, make FST matrix for saving the computational costs and then, calculate FST and the other statistics.
### Output S, PI, Tajima's D and pairwise FST.
### PI = average number of pairwise differences / seq length.

calc.fst.mat <- function(infile.header, n.pop) {
    ## Make 2 population combinations
    pop.comb <- combn(1:n.pop, 2) # Matrix of 2 x n.pop
    n.comb <- ncol(pop.comb)
    infile <- paste(infile.header, "_jointMAFpop",
                    pop.comb[2, ]-1, "_", pop.comb[1, ]-1, ".obs",
                    sep = "")

    ## Calculate FST matrix
    out.list <- list() # List for output

    for (i in 1:n.comb) {
        d <- read.table(infile[i], skip = 2, header = FALSE)[, -1]
        d <- as.matrix(d)

        ## Number of genecopies
        N1 <- dim(d)[2] - 1 # Pop1
        N2 <- dim(d)[1] - 1 # Pop2
        NT <- N1 + N2 # Total

        ## Number of minor alleles (matrix)
        A1 <- matrix(rep(0:N1, times = N2+1), nrow = N2+1, byrow = TRUE)
        A2 <- matrix(rep(0:N2, times = N1+1), ncol = N1+1, byrow = FALSE)
        AT <- A1 + A2

        ## Minor allele frequency (matirx)
        P1 <- A1/N1
        P2 <- A2/N2
        PT <- AT/NT

        ## Expected heterozygosity (matrix)
        H1 <- 2*P1*(1 - P1)
        H2 <- 2*P2*(1 - P2)
        HT <- 2*PT*(1 - PT)
    
        HS <- (H1+H2)/2
        
        FST <- (HT-HS)/HT
        FST[is.na(FST)] <- 0 # Non-variable sites

        out.list <- c(out.list, list(FST))
    }

    ## Prepare output list
    out.list <- c(out.list, list(pop.comb), list(n.pop))
    names(out.list) <- c(paste("FST_", pop.comb[2, ]-1, "_",
                               pop.comb[1, ]-1, sep = ""),
                         "PopComb", "NPop")

    
    return(out.list)
}

sfs2sumstat <- function(infile.header, fst.mat) {
    n.comb <- ncol(fst.mat$PopComb)
    infile <- paste(infile.header, "_jointMAFpop",
                    fst.mat$PopComb[2, ]-1, "_", fst.mat$PopComb[1, ]-1,
                    ".obs", sep = "")

    ## 1. Calculate single population statistics
    S <- PI <- D <- rep(NA, fst.mat$NPop)

    ## Pop1
    d <- read.table(infile[1], skip = 2, header = FALSE)[, -1]
    d <- as.matrix(d) # For speed up
    N <- dim(d)[2] - 1 # Number of genecopies
    
    msfs <- apply(d, 2, sum) # Minor allele SFS
    msfs <- msfs[1:(0.5*N + 1)] + c(rev(msfs[(0.5*N + 2):(N+1)]), 0) # Fold
    msfs <- as.vector(msfs)

    S[1] <- sum(msfs[-1]) # Number of segregating sites

    n.diff <- 0 # Number of pairwise differences
    for (i in 1:(length(msfs)-1)) {
        n.diff <- n.diff + i*(N-i)*msfs[i+1]
    }
    n.pair <- N*(N-1)/2
    seq.length <- sum(msfs)
    k <- n.diff/n.pair
    PI[1] <- k/seq.length

    a1 <- sum(1/(1:(N-1))) # Tajima's D
    a2 <- sum(1/((1:(N-1))^2))
    b1 <- (N+1)/3/(N-1)
    b2 <- 2*(N^2 + N + 3)/9/N/(N-1)
    c1 <- b1 - 1/a1
    c2 <- b2 - (N+2)/a1/N + a2/(a1^2)
    e1 <- c1/a1
    e2 <- c2/(a1^2 + a2)

    D[1] <- (k - S[1]/a1)/(sqrt(e1*S[1] + e2*S[1]*(S[1] - 1)))
    
    ## Pop2 or later
    for (i in 2:fst.mat$NPop) {
        d <- read.table(infile[i-1], skip = 2, header = FALSE)[, -1]
        d <- as.matrix(d) # For speed up
        N <- dim(d)[1] - 1 # Number of genecopies

        msfs <- apply(d, 1, sum) # Minor allele SFS
        msfs <- msfs[1:(0.5*N + 1)] + c(rev(msfs[(0.5*N + 2):(N+1)]), 0) # Fold
        msfs <- as.vector(msfs)
        
        S[i] <- sum(msfs[-1]) # Number of segregating sites

        n.diff <- 0 # Number of pairwise differences
        for (j in 1:(length(msfs)-1)) {
            n.diff <- n.diff + j*(N-j)*msfs[j+1]
        }
        n.pair <- N*(N-1)/2
        seq.length <- sum(msfs)
        k <- n.diff/n.pair
        PI[i] <- k/seq.length

        a1 <- sum(1/(1:(N-1))) # Tajima's D
        a2 <- sum(1/((1:(N-1))^2))
        b1 <- (N+1)/3/(N-1)
        b2 <- 2*(N^2 + N + 3)/9/N/(N-1)
        c1 <- b1 - 1/a1
        c2 <- b2 - (N+2)/a1/N + a2/(a1^2)
        e1 <- c1/a1
        e2 <- c2/(a1^2 + a2)
        
        D[i] <- (k - S[i]/a1)/(sqrt(e1*S[i] + e2*S[i]*(S[i] - 1)))
    }

    ## 2. Pairwise FST (for variable sites only)
    FST <- rep(NA, n.comb)
    
    for (i in 1:n.comb) {
        d <- read.table(infile[i], skip = 2, header = FALSE)[, -1]
        d <- as.matrix(d) # For speed up

        n.variable.site <- sum(d) - d[1, 1] - d[dim(d)[1], dim(d)[2]]
        
        FST[i] <- sum(fst.mat[[i]]*d)/n.variable.site
    }

    ## Output
    out <- c(S, PI, D, FST)
    names(out) <- c(paste("S", 0:(fst.mat$NPop - 1), sep = ""),
                    paste("PI", 0:(fst.mat$NPop - 1), sep = ""),
                    paste("D",  0:(fst.mat$NPop - 1), sep = ""),
                    names(fst.mat)[1:n.comb])

    return(out)
}

#fst.mat <- calc.fst.mat("./output/haha_40inds", 4)
#out <- sfs2sumstat("./output/haha_40inds", fst.mat) 


