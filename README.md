Calculate 2D-msfs from VCF and popmap files of stacks, and summary statistics from the 2D-msfs using R
======
* When I made these scripts, I refered to vcf2sfs.r  <https://github.com/shenglin-liu/vcf2sfs>
* In vcf, Reference allele = 0, Derived allele = 1  and missing = .
* In the output, ref homo = 0, hetero = 1, derived homo = 2 and missing = NA
* Populaitons must be numbered from 1, for example, 1, 2, 3...
* Missing numbers for population ID are not acceptable.
* vcf must be created by stacks with __-p X option__. X must be the same as number of populations. If X < number of populations, some population may have 0 numbers of genotyped individuals. In such case, this script will stop
* 2 dimentional minor allel site frequency spectrum (2D-msfs) calculated using this R script can be used for fastsimcoal2  <http://cmpg.unibe.ch/software/fastsimcoal27/>
* Missing data are compensated by bootstrapping within the same population
* Utilize observed summary statistics (S, Pi, Tajima's D and pairwise FST) for confirmation of data set. These summary statistics can be used for approximate Bayesian computation.


Example
------
Example for calculation of 2D-msfs

    source("r211223vcf2sfs.r")
    dir.create("./figs", showWarnings = FALSE, recursive = TRUE)
    dir.create("./output", showWarnings = FALSE, recursive = TRUE)

    infile.vcf <- "~/works/works22/Lee/data/stacks/220722Abelia/Abelia2pops187inds/populations.snps.vcf" # Path for input VCF file
    infile.popmap <- "~/works/works22/Lee/data/stacks/220722Abelia/popmap/Abelia2pops187inds.txt" # Path for input popmap file
    n.all.sites <- 242107 # Length of the sequence including non-variable sites
    outfile.base.sfs <- "./output/abelia2pops187inds" # Base name of output sfs files
    outfile.base.plot <- "./figs/abelia2pops187inds" # Base name of putput sfs plot files

    gt <- vcf2gt(infile.vcf, infile.popmap)
    gt.filled <- fill.missing(gt, seed = 46)
    msfs <- gt22dmsfs(gt.filled, n.all.sites)
    write.2dmsfs.fsc(msfs, outfile.base.sfs)
    plot.2dmsfs(msfs, outfile.base.plot)

Example for calculation of sumstats

    source("r220105sfs2sumstat.r")
    
    n.pop <- 2 # Number of populations
    infile.header <- "./output/abelia2pops187inds" # Base name of input sfs files
    outfile <- ./output/observed_sumstat_abelia2pops187inds220726.csv # Name of output csv file
    
    fst.mat <- calc.fst.mat(infile.header, n.pop)
    out <- sfs2sumstat(infile.header, fst.mat)

    out <- as.data.frame(t(out))
    write.table(out, file = outfile, sep = ",", row.names = FALSE, quote = FALSE)

