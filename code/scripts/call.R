#!/usr/bin/Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(glmmTMB)
})


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(paste(
    "Usage: Rscript epsilon_call.R model.rds input.csv output.vcf[.gz]",
    "[binomial_test|bayes_posterior] [BH|prior] [alpha|posterior_cutoff]",
    "[generic|specific]"
  ))
}

model_path <- args[1]
input_path <- args[2]
output_path <- args[3]
calling_mode <- if (length(args) >= 4L) args[4] else "binomial_test"

if (!calling_mode %in% c("binomial_test", "bayes_posterior")) {
  stop("calling_mode must be 'binomial_test' or 'bayes_posterior'")
}

if (calling_mode == "binomial_test") {
  correction_method <- if (length(args) >= 5L) args[5] else "BH"
  alpha <- if (length(args) >= 6L) as.numeric(args[6]) else 0.05
  if (!correction_method %in% p.adjust.methods) {
    stop("Unsupported p-value adjustment method: ", correction_method)
  }
  if (is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1")
  }
} else {
  prior <- if (length(args) >= 5L) as.numeric(args[5]) else 0.3
  posterior_cutoff <- if (length(args) >= 6L) as.numeric(args[6]) else 0.5
  if (is.na(prior) || prior <= 0 || prior >= 1) stop("prior must be between 0 and 1")
  if (is.na(posterior_cutoff) || posterior_cutoff <= 0 || posterior_cutoff >= 1) {
    stop("posterior_cutoff must be between 0 and 1")
  }
}

alt_mode <- if (length(args) >= 7L) args[7] else "generic"
if (!alt_mode %in% c("specific", "generic")) {
  stop("alt_mode must be 'specific' or 'generic'")
}

if (!file.exists(model_path)) stop("Model file not found: ", model_path)
if (!file.exists(input_path)) stop("Input file not found: ", input_path)

model <- readRDS(model_path)
data <- fread(input_path, header = TRUE, sep = ",")

required <- c("chrom", "pos", "sampleID", "ref", "coverage")
missing <- setdiff(required, names(data))
if (length(missing)) stop("Missing required input columns: ", paste(missing, collapse = ", "))


data[, pos := suppressWarnings(as.integer(pos))]
data[, coverage := suppressWarnings(as.numeric(coverage))]
data <- data[
  !is.na(chrom) & nzchar(chrom) &
  !is.na(pos) & pos > 0 &
  !is.na(ref) & nzchar(ref) &
  !is.na(coverage) & coverage > 0
]


script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg)) {
  script_path <- sub("^--file=", "", script_arg[1])
  script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
  project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  scaling_file <- file.path(project_dir, "fitted_error_model", "scaling_params.rds")
  if (file.exists(scaling_file)) {
    scaling_params <- readRDS(scaling_file)
    scale_cols <- c("GC_content", "mismatches", "median_mapQ", "median_baseQ")
    for (nm in scale_cols) {
      if (nm %in% names(data) &&
          !is.null(scaling_params[[nm]]$center) &&
          !is.null(scaling_params[[nm]]$scale) &&
          is.finite(scaling_params[[nm]]$scale) &&
          scaling_params[[nm]]$scale != 0) {
        data[[nm]] <- (data[[nm]] - scaling_params[[nm]]$center) /
          scaling_params[[nm]]$scale
      }
    }
  }
}

call_one_alt <- function(d, alt_base, alt_counts, mode, model,
                         correction_method = NULL, alpha = NULL,
                         prior = NULL, posterior_cutoff = NULL) {
  alt_counts <- suppressWarnings(as.numeric(alt_counts))
  valid <- !is.na(alt_counts) & alt_counts > 0 &
    alt_counts <= d$coverage & alt_base != d$ref
  if (!any(valid)) return(data.table())

  dd <- copy(d[valid])
  ac <- alt_counts[valid]
  eps <- as.numeric(predict(
    model, newdata = as.data.frame(dd),
    allow.new.levels = TRUE, type = "response"
  ))
  eps <- pmin(pmax(eps, 0), 1)

  ok <- is.finite(eps) & is.finite(dd$coverage) & dd$coverage > 0
  if (!any(ok)) return(data.table())
  dd <- dd[ok]
  ac <- ac[ok]
  eps <- eps[ok]

  if (mode == "binomial_test") {
    p_raw <- pbinom(ac - 1, size = dd$coverage, prob = eps,
                    lower.tail = FALSE)
    p_adj <- p.adjust(p_raw, method = correction_method)
    keep <- which(is.finite(p_adj) & p_adj < alpha)
    if (!length(keep)) return(data.table())

    out <- data.table(
      CHROM = as.character(dd$chrom[keep]),
      POS = as.integer(dd$pos[keep]),
      ID = as.character(dd$sampleID[keep]),
      REF = as.character(dd$ref[keep]),
      ALT = rep(alt_base, length(keep)),
      QUAL = ".",
      FILTER = "PASS",
      INFO = sprintf(
        "DP=%d;AD=%d;AF=%.8g;PVAL_RAW=%.8g;PVAL_ADJ=%.8g;SIG_LEVEL=%.8g",
        as.integer(round(dd$coverage[keep])),
        as.integer(round(ac[keep])),
        ac[keep] / dd$coverage[keep],
        p_raw[keep], p_adj[keep], alpha
      )
    )
  } else {
    vaf <- ac / dd$coverage
    
    log_bf <- dbinom(ac, size = dd$coverage, prob = vaf, log = TRUE) -
      dbinom(ac, size = dd$coverage, prob = eps, log = TRUE)
    log_odds <- log_bf + log(prior) - log1p(-prior)
    posterior <- plogis(log_odds)
    keep <- which(is.finite(posterior) & posterior > posterior_cutoff)
    if (!length(keep)) return(data.table())

    out <- data.table(
      CHROM = as.character(dd$chrom[keep]),
      POS = as.integer(dd$pos[keep]),
      ID = as.character(dd$sampleID[keep]),
      REF = as.character(dd$ref[keep]),
      ALT = rep(alt_base, length(keep)),
      QUAL = ".",
      FILTER = "PASS",
      INFO = sprintf(
        "DP=%d;AD=%d;AF=%.8g;PRIOR=%.8g;POSTERIOR=%.8g;POSTERIOR_CUTOFF=%.8g",
        as.integer(round(dd$coverage[keep])),
        as.integer(round(ac[keep])),
        vaf[keep], prior, posterior[keep], posterior_cutoff
      )
    )
  }

  out
}


if ("general_alt_counts" %in% names(data)) {
  data <- data[!is.na(general_alt_counts) & general_alt_counts > 0]
}

results <- list()
if (nrow(data) > 0L) {
  if (alt_mode == "specific") {
    needed <- c("alt", "specific_alt_counts")
    missing <- setdiff(needed, names(data))
    if (length(missing)) stop("Specific mode requires columns: ", paste(missing, collapse = ", "))

    
    alts <- unique(as.character(data$alt))
    alts <- alts[!is.na(alts) & alts %in% c("A", "C", "G", "T")]
    for (a in alts) {
      d <- data[as.character(alt) == a]
      results[[length(results) + 1L]] <- call_one_alt(
        d, a, d$specific_alt_counts, calling_mode, model,
        correction_method, alpha,
        if (calling_mode == "bayes_posterior") prior else NULL,
        if (calling_mode == "bayes_posterior") posterior_cutoff else NULL
      )
    }
  } else {
    count_cols <- c(A = "alt_A", C = "alt_C", G = "alt_G", T = "alt_T")
    missing <- setdiff(unname(count_cols), names(data))
    if (length(missing)) stop("Generic mode requires columns: ", paste(missing, collapse = ", "))

    for (a in names(count_cols)) {
      results[[length(results) + 1L]] <- call_one_alt(
        data, a, data[[count_cols[[a]]]], calling_mode, model,
        correction_method, alpha,
        if (calling_mode == "bayes_posterior") prior else NULL,
        if (calling_mode == "bayes_posterior") posterior_cutoff else NULL
      )
    }
  }
}

if (length(results)) {
  calls <- rbindlist(results, use.names = TRUE, fill = TRUE)
} else {
  calls <- data.table(
    CHROM = character(), POS = integer(), ID = character(), REF = character(),
    ALT = character(), QUAL = character(), FILTER = character(), INFO = character()
  )
}


if (nrow(calls)) {
  calls <- calls[
    !is.na(CHROM) & nzchar(CHROM) &
    !is.na(POS) & POS > 0 &
    !is.na(REF) & nzchar(REF) &
    !is.na(ALT) & nzchar(ALT) &
    REF != ALT
  ]
  setorder(calls, CHROM, POS, REF, ALT)
  
  calls[, ID := gsub("[[:space:];,]", "_", fifelse(is.na(ID) | ID == "", ".", ID))]
}

header <- c(
  "##fileformat=VCFv4.3",
  "##source=Epsilon",
  '##FILTER=<ID=PASS,Description="All filters passed">',
  '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">',
  '##INFO=<ID=AD,Number=1,Type=Integer,Description="Alternate allele read count">',
  '##INFO=<ID=AF,Number=1,Type=Float,Description="Observed alternate allele fraction">'
)
if (calling_mode == "binomial_test") {
  header <- c(
    header,
    '##INFO=<ID=PVAL_RAW,Number=1,Type=Float,Description="Raw one-sided binomial p-value">',
    '##INFO=<ID=PVAL_ADJ,Number=1,Type=Float,Description="Multiple-testing adjusted p-value">',
    '##INFO=<ID=SIG_LEVEL,Number=1,Type=Float,Description="Significance threshold">'
  )
} else {
  header <- c(
    header,
    '##INFO=<ID=PRIOR,Number=1,Type=Float,Description="Prior probability of a variant">',
    '##INFO=<ID=POSTERIOR,Number=1,Type=Float,Description="Posterior probability of a variant">',
    '##INFO=<ID=POSTERIOR_CUTOFF,Number=1,Type=Float,Description="Posterior calling threshold">'
  )
}
header <- c(header, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")


tmp_vcf <- tempfile(fileext = ".vcf")
writeLines(header, tmp_vcf)
if (nrow(calls)) {
  fwrite(calls, file = tmp_vcf, append = TRUE, sep = "\t", quote = FALSE,
         col.names = FALSE, na = ".")
}

if (grepl("\\.vcf\\.gz$", output_path, ignore.case = TRUE)) {
  
  if (!nzchar(Sys.which("bgzip"))) stop("Output requested as .vcf.gz, but bgzip is not on PATH")
  if (!nzchar(Sys.which("tabix"))) stop("Output requested as .vcf.gz, but tabix is not on PATH")
  status <- system2("bgzip", args = c("-c", shQuote(tmp_vcf)), stdout = output_path)
  if (status != 0L) stop("bgzip failed with exit status ", status)
  if (nrow(calls) > 0L) {
    status <- system2("tabix", args = c("-f", "-p", "vcf", shQuote(output_path)))
    if (status != 0L) stop("tabix indexing failed with exit status ", status)
  } else {
    message("No variants passed the calling threshold; wrote a valid header-only VCF.gz (no index created).")
  }
} else {
  if (!grepl("\\.vcf$", output_path, ignore.case = TRUE)) {
    warning("Output filename does not end in .vcf or .vcf.gz; writing plain-text VCF.")
  }
  if (!file.copy(tmp_vcf, output_path, overwrite = TRUE)) stop("Could not write output VCF: ", output_path)
}

unlink(tmp_vcf)
message("Wrote ", nrow(calls), " variant records to ", output_path)
if (grepl("\\.vcf\\.gz$", output_path, ignore.case = TRUE)) {
  message("Created tabix index: ", output_path, ".tbi")
}
