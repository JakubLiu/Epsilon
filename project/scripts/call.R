library(data.table)
library(MASS)
library(stringr)
library(glmmTMB)


# =========================================================== setup command line arguments ========================================================================
args <- commandArgs(trailingOnly = TRUE)

# Required: fitted model RDS
if (length(args) >= 1) {
    fitted_model_rds <- args[1]
} else {
    stop("Missing required argument: fitted_model_rds")
}

# Required: input data
if (length(args) >= 2) {
    input_data <- args[2]
} else {
    stop("Missing required argument: input_data")
}

# Required: output VCF
if (length(args) >= 3) {
    output_vcf <- args[3]
} else {
    stop("Missing required argument: output_vcf")
}

# Optional: variant calling mode
if (length(args) >= 4) {
    variant_calling_mode <- args[4]
} else {
    variant_calling_mode <- "binomial_test"
}

if (variant_calling_mode == "binomial_test") {

    if (length(args) >= 5) {
        multiple_testing_correction_method <- args[5]
    } else {
        multiple_testing_correction_method <- "BH"
    }

    if (length(args) >= 6) {
        significance_level_alpha <- as.numeric(args[6])
    } else {
        significance_level_alpha <- 0.05
    }

} else if (variant_calling_mode == "bayes_posterior") {

    if (length(args) >= 5) {
        prior <- as.numeric(args[5])
    } else {
        prior <- 0.3
    }

    if (length(args) >= 6) {
        posterior_cutoff <- as.numeric(args[6])
    } else {
        posterior_cutoff <- 0.5
    }

} else {
    stop("variant_calling_mode must be 'binomial_test' or 'bayes_posterior'")
}

# Optional: alt mode
if (length(args) >= 7) {
    alt_mode <- args[7]
} else {
    alt_mode <- "generic"
}

if (!(alt_mode %in% c("specific", "generic"))) {
    stop("alt_mode must be 'specific' or 'generic'")
}

# =======================================================================================================================================================


# ========================================================== variant calling function ===================================================================

variant_calling_generic_binomial <- function(alt_base, alt_base_count){

    EPS_PRED <- predict(model, newdata=data, allow.new.levels=TRUE, type='response')  # the predicted error rate
    observed_coverage <- data$coverage
    observed_vaf <- alt_base_count/observed_coverage
    PVALUES <- pbinom(alt_base_count - 1, size = observed_coverage, prob = EPS_PRED, lower.tail = FALSE)
    PVALUES_ADJUSTED <- p.adjust(PVALUES, method = multiple_testing_correction_method, n = length(PVALUES))
    significant_idx <- which(PVALUES_ADJUSTED < significance_level_alpha)
    sig_data <- data[significant_idx,]


    vcf <- data.table(
    CHROM = sig_data$chrom,
    POS = sig_data$pos,
    ID = sig_data$sampleID,
    REF = sig_data$ref,
    ALT = rep(alt_base, length(significant_idx)),
    PVAL_ADJ = PVALUES_ADJUSTED[significant_idx],
    PVAL_RAW = PVALUES[significant_idx],
    SIG_LEVEL = rep(significance_level_alpha, length(significant_idx)),
    INFO = sprintf(
        "DP=%d;AD=%d;AF=%.6f",
        sig_data$coverage,
        alt_base_count[significant_idx],
        alt_base_count[significant_idx] / sig_data$coverage
     )
    )

    return(vcf)
}




variant_calling_generic_bayesian <- function(alt_base, alt_base_count){

    observed_coverage <- data$coverage
    observed_vaf <- alt_base_count/observed_coverage
    EPS_PRED <- predict(model, newdata=data, allow.new.levels=TRUE, type='response')
    P_H0 <- dbinom(x=alt_base_count, size=observed_coverage, prob=EPS_PRED, log=TRUE)
    P_H1 <- dbinom(x=alt_base_count, size=observed_coverage, prob=observed_vaf, log=TRUE)
    BF <- exp(P_H1 - P_H0)
    POSTERIOR <- BF * prior / (BF * prior + (1 - prior))
    significant_idx <- which(POSTERIOR > posterior_cutoff)
    sig_data <- data[significant_idx,]

    vcf <- data.table(
    CHROM = sig_data$chrom,
    POS = sig_data$pos,
    ID = sig_data$sampleID,
    REF = sig_data$ref,
    ALT = rep(alt_base, length(significant_idx)),
    PRIOR = rep(prior, length(significant_idx)),
    POSTERIOR = POSTERIOR[significant_idx],
    POSTERIOR_CUTOFF = rep(posterior_cutoff, length(significant_idx)),
    INFO = sprintf(
        "DP=%d;AD=%d;AF=%.6f",
        sig_data$coverage,
        alt_base_count[significant_idx],
        alt_base_count[significant_idx] / sig_data$coverage
     )
    )

    return(vcf)
    
}


variant_calling_specific_binomial <- function(){

    EPS_PRED <- predict(model, newdata=data, allow.new.levels=TRUE, type='response')  # the predicted error rate
    observed_coverage <- data$coverage
    observed_alt_count <- data$specific_alt_counts
    observed_vaf <- observed_alt_count/observed_coverage
    alt <- data$alt
    PVALUES <- pbinom(observed_alt_count - 1, size = observed_coverage, prob = EPS_PRED, lower.tail = FALSE)
    PVALUES_ADJUSTED <- p.adjust(PVALUES, method = multiple_testing_correction_method, n = length(PVALUES))
    significant_idx <- which(PVALUES_ADJUSTED < significance_level_alpha)
    sig_data <- data[significant_idx,]
    

    vcf <- data.table(
    CHROM = sig_data$chrom,
    POS = sig_data$pos,
    ID = sig_data$sampleID,
    REF = sig_data$ref,
    ALT = alt[significant_idx],
    PVAL_ADJ = PVALUES_ADJUSTED[significant_idx],
    PVAL_RAW = PVALUES[significant_idx],
    SIG_LEVEL = rep(significance_level_alpha, length(significant_idx)),
    INFO = sprintf(
        "DP=%d;AD=%d;AF=%.6f",
        sig_data$coverage,
        observed_alt_count[significant_idx],
        observed_alt_count[significant_idx] / sig_data$coverage
     )
    )

    
    return(vcf)

}



variant_calling_specific_bayesian <- function(){

    observed_coverage <- data$coverage
    observed_alt_count <- data$specific_alt_counts
    observed_vaf <- observed_alt_count/observed_coverage
    alt <- data$alt
    EPS_PRED <- predict(model, newdata=data, allow.new.levels=TRUE, type='response')
    P_H0 <- dbinom(x=observed_alt_count, size=observed_coverage, prob=EPS_PRED, log=TRUE)
    P_H1 <- dbinom(x=observed_alt_count, size=observed_coverage, prob=observed_vaf, log=TRUE)
    BF <- exp(P_H1 - P_H0)
    POSTERIOR <- BF * prior / (BF * prior + (1 - prior))
    significant_idx <- which(POSTERIOR > posterior_cutoff)
    sig_data <- data[significant_idx,]

    vcf <- data.table(
    CHROM = sig_data$chrom,
    POS = sig_data$pos,
    ID = sig_data$sampleID,
    REF = sig_data$ref,
    ALT = alt[significant_idx],
    PRIOR = rep(prior, length(significant_idx)),
    POSTERIOR = POSTERIOR[significant_idx],
    POSTERIOR_CUTOFF = rep(posterior_cutoff, length(significant_idx)),
    INFO = sprintf(
        "DP=%d;AD=%d;AF=%.6f",
        sig_data$coverage,
        observed_alt_count[significant_idx],
        observed_alt_count[significant_idx] / sig_data$coverage
     )
    )

    return(vcf)
    
}



# =======================================================================================================================================================
model <- readRDS(fitted_model_rds)
data <- fread(input_data, head = TRUE, sep = ',')
data <- data[data$general_alt_counts > 0,]   # this stops potential 0 alt calling with a strong prior

if (nrow(data) == 0) {
    cat("No alternative reads found. Writing empty VCF.\n")

    header <- c(
        "##fileformat=VCFv4.3",
        "##source=Epsilon",
        "#CHROM\tPOS\tID\tREF\tALT"
    )

    writeLines(header, output_vcf)

    quit(save = "no", status = 0)
}

# optionally scale the predictors
script_path <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())])
script_dir <- dirname(normalizePath(script_path))
project_dir <- normalizePath(file.path(script_dir, ".."))

scaling_file <- file.path(project_dir, "fitted_error_model", "scaling_params.rds")

if (file.exists(scaling_file)) {
    
  scaling_params <- readRDS(scaling_file)
  data$GC_content <- (data$GC_content - scaling_params$GC_content$center) / scaling_params$GC_content$scale
  data$mismatches <- (data$mismatches - scaling_params$mismatches$center) / scaling_params$mismatches$scale
  data$median_mapQ <- (data$median_mapQ - scaling_params$median_mapQ$center) / scaling_params$median_mapQ$scale
  data$median_baseQ <- (data$median_baseQ - scaling_params$median_baseQ$center) / scaling_params$median_baseQ$scale
}



if(alt_mode == 'specific'){
    
    if(variant_calling_mode == 'binomial_test'){
        
        res_vcf <- variant_calling_specific_binomial()
        header <- c(
            "##fileformat=VCFv4.3",
            "##source=Epsilon Binomial test mode (specific alternative reads)",
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
            '##INFO=<ID=AD,Number=1,Type=Integer,Description="Alternate Allele Count">',
            '##INFO=<ID=AF,Number=1,Type=Float,Description="Observed Variant Allele Fraction">',
            "#CHROM\tPOS\tID\tREF\tALT\tPVAL_ADJ\tPVAL_RAW\tSIG_LEVEL\tINFO")

        writeLines(header, output_vcf)

        alt_count <- as.numeric(str_extract(res_vcf$INFO, "(?<=AD=)\\d+"))
        mask1 <- alt_count > 0
        mask2 <- res_vcf$REF != res_vcf$ALT
        mask3 <- nchar(res_vcf$REF) == 1
        mask4 <- nchar(res_vcf$ALT) == 1
        mask <- mask1 & mask2 & mask3 & mask4
        res_vcf <- res_vcf[mask,]

        fwrite(res_vcf,file = output_vcf,append = TRUE,sep = "\t",quote = FALSE,col.names = FALSE)

    }else if(variant_calling_mode == 'bayes_posterior'){
        
        res_vcf <- variant_calling_specific_bayesian()
        header <- c(
            "##fileformat=VCFv4.3",
            "##source=Epsilon Bayesian posterior mode (specific alternative reads)",
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
            '##INFO=<ID=AD,Number=1,Type=Integer,Description="Alternate Allele Count">',
            '##INFO=<ID=AF,Number=1,Type=Float,Description="Observed Variant Allele Fraction">',
            "#CHROM\tPOS\tID\tREF\tALT\tPRIOR\tPOSTERIOR\tPOSTERIOR_CUTOFF\tINFO")
            
        writeLines(header, output_vcf)

        alt_count <- as.numeric(str_extract(res_vcf$INFO, "(?<=AD=)\\d+"))
        mask1 <- alt_count > 0
        mask2 <- res_vcf$REF != res_vcf$ALT
        mask3 <- nchar(res_vcf$REF) == 1
        mask4 <- nchar(res_vcf$ALT) == 1
        mask <- mask1 & mask2 & mask3 & mask4
        res_vcf <- res_vcf[mask,]

        fwrite(res_vcf,file = output_vcf,append = TRUE,sep = "\t",quote = FALSE,col.names = FALSE)
    }

   

}else if (alt_mode == 'generic'){
    
    if(variant_calling_mode == 'binomial_test'){

        res_vcf_A <- variant_calling_generic_binomial(alt_base = 'A', alt_base_count = data$alt_A)
        res_vcf_T <- variant_calling_generic_binomial(alt_base = 'T', alt_base_count = data$alt_T)
        res_vcf_G <- variant_calling_generic_binomial(alt_base = 'G', alt_base_count = data$alt_G)
        res_vcf_C <- variant_calling_generic_binomial(alt_base = 'C', alt_base_count = data$alt_C)
        res_vcf <- data.frame(rbind(res_vcf_A, res_vcf_C, res_vcf_G, res_vcf_T))

         header <- c(
            "##fileformat=VCFv4.3",
            "##source=Epsilon Binomial test mode (generic alternative reads)",
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
            '##INFO=<ID=AD,Number=1,Type=Integer,Description="Alternate Allele Count">',
            '##INFO=<ID=AF,Number=1,Type=Float,Description="Observed Variant Allele Fraction">',
            "#CHROM\tPOS\tID\tREF\tALT\tPVAL_ADJ\tPVAL_RAW\tSIG_LEVEL\tINFO"
            )

        writeLines(header, output_vcf)

        alt_count <- as.numeric(str_extract(res_vcf$INFO, "(?<=AD=)\\d+"))
        mask1 <- alt_count > 0
        mask2 <- res_vcf$REF != res_vcf$ALT
        mask3 <- nchar(res_vcf$REF) == 1
        mask4 <- nchar(res_vcf$ALT) == 1
        mask <- mask1 & mask2 & mask3 & mask4
        res_vcf <- res_vcf[mask,]

        fwrite(res_vcf,file = output_vcf,append = TRUE,sep = "\t",quote = FALSE,col.names = FALSE)



    }else if (variant_calling_mode == 'bayes_posterior'){

        res_vcf_A <- variant_calling_generic_bayesian(alt_base = 'A', alt_base_count = data$alt_A)
        res_vcf_T <- variant_calling_generic_bayesian(alt_base = 'T', alt_base_count = data$alt_T)
        res_vcf_G <- variant_calling_generic_bayesian(alt_base = 'G', alt_base_count = data$alt_G)
        res_vcf_C <- variant_calling_generic_bayesian(alt_base = 'C', alt_base_count = data$alt_C)
        res_vcf <- data.frame(rbind(res_vcf_A, res_vcf_C, res_vcf_G, res_vcf_T))

        header <- c(
            "##fileformat=VCFv4.3",
            "##source=Epsilon Bayesian posterior mode (generic alternative reads)",
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
            '##INFO=<ID=AD,Number=1,Type=Integer,Description="Alternate Allele Count">',
            '##INFO=<ID=AF,Number=1,Type=Float,Description="Observed Variant Allele Fraction">',
            "#CHROM\tPOS\tID\tREF\tALT\tPRIOR\tPOSTERIOR\tPOSTERIOR_CUTOFF\tINFO"
            )

        writeLines(header, output_vcf)

        alt_count <- as.numeric(str_extract(res_vcf$INFO, "(?<=AD=)\\d+"))
        mask1 <- alt_count > 0
        mask2 <- res_vcf$REF != res_vcf$ALT
        mask3 <- nchar(res_vcf$REF) == 1
        mask4 <- nchar(res_vcf$ALT) == 1
        mask <- mask1 & mask2 & mask3 & mask4
        res_vcf <- res_vcf[mask,]

        fwrite(res_vcf,file = output_vcf,append = TRUE,sep = "\t",quote = FALSE,col.names = FALSE)

    }
}





