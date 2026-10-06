library(data.table)

args <- commandArgs(trailingOnly = TRUE)


if (length(args) >= 1) {
    tumor_vcf <- args[1]
} else {
    stop("Missing required argument: tumor_vcf")
}

args <- commandArgs(trailingOnly = TRUE)


if (length(args) >= 2) {
    normal_vcf <- args[2]
} else {
    stop("Missing required argument: normal_vcf")
}


if (length(args) >= 3) {
    output_vcf <- args[3]
} else {
    stop("Missing required argument: output_vcf")
}


vcf_header <- readLines(tumor_vcf)
vcf_header <- vcf_header[grepl("^##|^#CHROM", vcf_header)]


tumor_vcf <- fread(tumor_vcf, skip = "#CHROM", sep = "\t")
normal_vcf <- fread(normal_vcf, skip = "#CHROM", sep = "\t")

colnames(tumor_vcf)[c(1,2,4,5)] <- c('CHROM', 'POS', 'REF', 'ALT')
colnames(normal_vcf)[c(1,2,4,5)] <- c('CHROM', 'POS', 'REF', 'ALT')

tumor_vcf$CHROM <- as.character(tumor_vcf$CHROM)
normal_vcf$CHROM <- as.character(normal_vcf$CHROM)


tumor_filtered <- tumor_vcf[
    !normal_vcf,
    on = .(CHROM, POS, REF, ALT)
]




writeLines(vcf_header, output_vcf)

fwrite(
    tumor_filtered,
    output_vcf,
    sep = "\t",
    append = TRUE,
    col.names = FALSE,
    quote = FALSE
)