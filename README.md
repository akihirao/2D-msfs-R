Calculate 2D-msfs from VCF and popmap files of stacks using R
=======
* I refered to vcf2sfs.r  <https://github.com/shenglin-liu/vcf2sfs>
* In vcf, Reference allele = 0, Derived allele = 1  and missing = .
* In the output, ref homo = 0, hetero = 1, derived homo = 2 and missing = NA
* Populaitons must be numbered from 1, for example, 1, 2, 3...
* Missing numbers for population ID are not acceptable.
* vcf should be created by stacks with -p X option. X must be the same as number of populations. If X < number of populations, some population may have 0 numbers of genotyped individuals. In such case, this script will stop with errors.








