### Make 2D-mSFS from VCF and popmap of stacks.
### Refer to vcf2sfs.r (https://github.com/shenglin-liu/vcf2sfs).
### In vcf, Reference allele = 0, Derived allele = 1  and missing = .
### In the output, ref homo = 0, hetero = 1, derived homo = 2 and missing = NA
### Populaitons must be numbered from 1, for example, 1, 2, 3...
### Missing numbers for population ID are not acceptable.

library(ggplot2)

### Read VCF and convert genotype (0, 1, 2)
vcf2gt <- function (infile.vcf, infile.popmap) {
    oldw <- getOption("warn")
    options(warn = -1)

    ## Read VCF file and popmap files.
    vcf.gt <- read.table(infile.vcf, sep = "\t", stringsAsFactors = FALSE)
    locus <- vcf.gt[, 1] # Extract RAD-locus information
    vcf.gt <- as.matrix(vcf.gt[, -c(1:9)]) # Remove headers

    popmap <- read.table(infile.popmap, sep = "\t", stringsAsFactors = FALSE)
    popmap <- popmap[, 2]

    nrow.vcf <- nrow(vcf.gt)
    ncol.vcf <- ncol(vcf.gt)
    
    allele1 <- substring(vcf.gt, 1, 1) # 1st allele
    allele2 <- substring(vcf.gt, 3, 3) # 2nd allele
    genotype <- matrix(as.integer(allele1) + as.integer(allele2),
                       nrow = nrow.vcf, ncol = ncol.vcf)
    
    options(warn = oldw)
    list(popmap = popmap, genotype = genotype, locus = locus)
}


### Fill missings (NA) with bootstraped allele.
### This prosidure is conducted within each populasion.
fill.missing <- function(gt, seed = 123) {
    set.seed(seed)
    d <- out <- gt$genotype
    pop <- gt$popmap
    n.pop <- unique(length((pop)))
    filter_out_locus <- c() #set locus ID for excluding
    for (i in 1:nrow(d)) {
        if (any(is.na(d[i, ]))) {
            for (j in 1:n.pop) {
                temp <- d[i, pop == j] # Extract genotypes of jth pop
                if(any(is.na(temp))) {
                    na_count_temp <- sum(is.na(temp))
                    n_ind_witin_pop <- length(temp)
                    if(na_count_temp==n_ind_witin_pop){
                        filter_out_locus <- c(filter_out_locus,i)
                    }else{
                        n.missing.ind <- sum(is.na(temp))
                        temp <- na.omit(temp) # Remove NA
                        n.allele <- 2*length(temp)

                        # Allele frequency of ref allele
                        freq.ref <- (sum(temp == 0)*2 + sum(temp == 1))/n.allele
                    
                        allele1 <- sample(c(0, 1), n.missing.ind, replace = TRUE,
                                      prob = c(freq.ref, 1-freq.ref))
                        allele2 <- sample(c(0, 1), n.missing.ind, replace = TRUE,
                                      prob = c(freq.ref, 1-freq.ref))
                
                        out[i, pop == j & is.na(out[i, ])] <- allele1 + allele2
                    }
                }
            }
        }
    }
    gt$genotype <- out
    gt$genotype <- gt$genotype[-filter_out_locus,]
    gt$locus <- gt$locus[-filter_out_locus]

    return(gt)
}


### Calculate 2D-mSFS
### All missings must be filled.
gt22dmsfs <- function(gt, n.all.sites) {
    ## Make 2 population combinations
    d <- gt$genotype
    popmap <- gt$popmap
    n.pop <- length(unique(popmap))
    pop.comb <- combn(1:n.pop, 2) # Matrix of 2 x n.pop
    n.comb <- ncol(pop.comb)
    
    ## Determin minor allele
    ## 0, ref allele; 1, derived allele; 0.5, allele freqs are equal 
    minor.allele <- sapply(1:nrow(d), function(i) {
        temp <- d[i, ]
        n.allele0 <- sum(temp == 0)*2 + sum(temp == 1)
        n.allele1 <- ncol(d)*2 - n.allele0
        if(n.allele0 < n.allele1) return(0)
        if(n.allele0 > n.allele1) return(1)
        if(n.allele0 == n.allele1) return(0.5)
    })
    
    ## Calculate 2D-mSFS for each 2 pop combination
    out.list <- list() # List for output
    
    for (i in 1:n.comb) {
        pop1 <- pop.comb[1, i]
        pop2 <- pop.comb[2, i]
        
        d2 <- d[, popmap == pop1 | popmap == pop2]
        popmap2 <- popmap[popmap == pop1 | popmap == pop2]
        
        n.all.ind <- ncol(d2)
        n.allele.pop1 <- 2*sum(popmap2 == pop1)
        n.allele.pop2 <- 2*sum(popmap2 == pop2)
        out <- matrix(0, nrow = n.allele.pop2+1, ncol = n.allele.pop1+1)
        
        for (j in 1:nrow(d2)) {
            temp <- d2[j, ]

            ## Count minor allele in each pop
            temp1 <- temp[popmap2 == pop1] # genotype of pop1
            temp2 <- temp[popmap2 == pop2] # genotype of pop2
        
            if(minor.allele[j] == 0) {
                n.minor.allele.pop1 <- sum(temp1 == 0)*2 + sum(temp1 == 1)
                n.minor.allele.pop2 <- sum(temp2 == 0)*2 + sum(temp2 == 1)
            }
        
            if(minor.allele[j] == 1) {
                n.minor.allele.pop1 <- sum(temp1 == 1) + sum(temp1 == 2)*2
                n.minor.allele.pop2 <- sum(temp2 == 1) + sum(temp2 == 2)*2
            }
            
            out[n.minor.allele.pop2+1, n.minor.allele.pop1+1] <-
                out[n.minor.allele.pop2+1, n.minor.allele.pop1+1] + 1

            if(minor.allele[j] == 0.5) {
                ## Ref allele freq + 0.5
                n.ref.allele.pop1 <- sum(temp1 == 0)*2 + sum(temp1 == 1)
                n.ref.allele.pop2 <- sum(temp2 == 0)*2 + sum(temp2 == 1)
                out[n.ref.allele.pop2+1, n.ref.allele.pop1+1] <-
                    out[n.ref.allele.pop2+1, n.ref.allele.pop1+1] + 0.5
                ## Derived allele freq + 0.5
                n.der.allele.pop1 <- sum(temp1 == 1) + sum(temp1 == 2)*2
                n.der.allele.pop2 <- sum(temp2 == 1) + sum(temp2 == 2)*2
                out[n.der.allele.pop2+1, n.der.allele.pop1+1] <-
                    out[n.der.allele.pop2+1, n.der.allele.pop1+1] + 0.5
            }
        }

        colnames(out) <- paste("d", pop1-1, "_", 0:n.allele.pop1, sep = "")
        rownames(out) <- paste("d", pop2-1, "_", 0:n.allele.pop2, sep = "")

        out[1, 1] <- out[1, 1] + n.all.sites - sum(out)
        #print(out)
        out.list <- c(out.list, list(out))
    }
    names(out.list) <- paste("jointMAFpop", pop.comb[2, ]-1, "_",
                             pop.comb[1, ]-1, sep = "")
    return(out.list)
}


### Output 2D-mSFS with fastsimcoal format
write.2dmsfs.fsc <- function(sfs, outfile.base) {
    oldw <- getOption("warn")
    options(warn = -1)

    n.sfs <- length(sfs)
    outfile <- paste(outfile.base, "_", names(sfs), ".obs", sep = "")

    for (i in 1:n.sfs) {
        sink(outfile[i], append = FALSE)
        cat("1 observations\n")
        cat("\t")
        sink()
        write.table(sfs[[i]], file = outfile[i], append = TRUE, quote = FALSE,
                    col.names = TRUE, row.names = TRUE, sep = "\t")
    }
    options(warn = oldw)
}


### Make 2D-mSFS fig with ggplot
plot.2dmsfs <- function(sfs, outfile.base, width = 6, height = 5) {
    library(ggplot2)
    n.sfs <- length(sfs)
    outfile <- paste(outfile.base, "_", names(sfs), ".pdf", sep = "")

    for (i in 1:n.sfs) {
        n.pop1 <- ncol(sfs[[i]])-1
        n.pop2 <- nrow(sfs[[i]])-1
        d <- data.frame(X = rep(0:n.pop1, each = n.pop2+1),
                        Y = rep(0:n.pop2, times = n.pop1+1),
                        Count = c(sfs[[i]]))
        d$Count[1] <- 0 # Fill 0 for (0, 0)
        d$Count[d$Count == 0] <- NA

        temp <- names(sfs)[i]
        pop1.name <- paste("Pop",
                           substring(temp, nchar(temp), nchar(temp)),
                           sep = "")
        pop2.name <- paste("Pop",
                           substring(temp, nchar(temp)-2, nchar(temp)-2),
                           sep = "")
        g <- ggplot(data = d, mapping = aes(x = X, y = Y, fill = Count)) +
            geom_raster() +
            labs(x = pop1.name, y = pop2.name) +
            theme_bw() +
            scale_fill_gradientn(colors = c("yellow", "red"))
        ggsave(outfile[i], width = width, height = height)
    }
}


### Bootstraping of SFS
boot.sfs <- function(gt, n.all.sites, outfile.base, n.sim) {
    for (i in 1:n.sim) {
        use.genotype <- sample(1:nrow(gt$genotype), replace = TRUE)
        boot.gt <- list(popmap = gt$popmap,
                        genotype = gt$genotype[use.genotype, ])
        boot.2dmsfs <- gt22dmsfs(boot.gt, n.all.sites)
        write.2dmsfs.fsc(boot.2dmsfs,
                         paste(outfile.base, "_", i, sep = ""))
    }
}

#infile.vcf <- "~/works/works21/setsuko/data/stacks/haha_40inds/populations.snps.vcf"
#infile.popmap <- "~/works/works21/setsuko/data/stacks/popmap/haha_40inds.txt"
#mygt <- vcf2gt(infile.vcf, infile.popmap)
#mygt.filled <- fill.missing(mygt, seed = 46)
#my2dmsfs <- gt22dmsfs(mygt.filled, 781321)
#write.2dmsfs.fsc(my2dmsfs, "haha_40inds") 
#plot.2dmsfs(my2dmsfs, "haha_40inds")
#boot.sfs(mygt.filled, 781321, "output/boot", 3)


