# Refactor the calculations in original/messy_analysis.R into these interfaces.
# Preserve the required behavior described in README.md.
# This file should define functions only: no reads, writes, or analysis run.

calculate_maf <- function(genotypes, na_rm = TRUE) {
    
    if (!(is.logical(na_rm) && length(na_rm) == 1 && !is.na(na_rm))) { #Validate na_rm 
        stop("na_rm must be one non-missing logical value")
    }
    
    if(missing(genotypes) || length(genotypes) == 0 || all(is.na(genotypes))){ #In case of missing or empty input
        return(NA_real_)
    }
    
    if (!is.numeric(genotypes) || !is.null(dim(genotypes))) { #Validate genotype
        stop("genotypes must be a numeric vector")
    }
    
    g <- genotypes[!is.na(genotypes)]
    
    for(value in g){
        if(!value %in% c(0,1,2)){ #Check on input vector values
            stop("genotypes must contain only 0, 1, 2, or NA")
        }
    }
    
    if(anyNA(genotypes) && na_rm == FALSE){
        return(NA_real_)
    }
    
    #Once all checks have passed calculate MAF
    p <- sum(g)/(2*length(g))
    maf <- min(p, 1 - p)
    return(maf)   
}

summarize_variant <- function(genotypes, min_maf = 0.25,
                              max_missing = 0.25) {
    if(!is.numeric(min_maf) || length(min_maf) != 1 || !(is.finite(min_maf)) || 
        min_maf > 0.5 || min_maf < 0){
        stop("min_maf must be a finite numeric value between 0 and 0.5 inclusive")
    }

    if(!is.numeric(max_missing) ||  length(max_missing) != 1 || !(is.finite(max_missing)) || 
       max_missing > 1 || max_missing < 0){
        stop("max_missing must be a finite numeric value between 0 and 1 inclusive")
    }
    
    if(missing(genotypes) || length(genotypes) == 0){ #In case of empty input
        n_called <- 0; missing_rate = NA_real_; maf = NA_real_; keep = FALSE
        return(data.frame(n_called, missing_rate, maf, keep))
    }
    
    if(all(is.na(genotypes))){ #In case of non-empty missing input
        n_called <- 0; missing_rate = 1; maf = NA_real_; keep = FALSE
        return(data.frame(n_called, missing_rate, maf, keep))
    }
    
    
    #Number of non-missing genotypes
    n_called <- length(genotypes[!is.na(genotypes)])
    
    #Number missing divided by total length
    missing_rate <- (n_called)/length(genotypes) #changed
    
    #Minor allele frequency
    maf <- calculate_maf(genotypes, na_rm = TRUE)
    
    #Keep is TRUE when all conditions are met
    keep <- (n_called >= 1 && 
                 !(is.na(maf)) &&
                 maf >= min_maf && 
                 missing_rate <= max_missing)
    
    return(data.frame(n_called, missing_rate, maf, keep))
}
