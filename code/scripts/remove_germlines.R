#!/usr/bin/Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: Rscript remove_germlines.R tumor.vcf.gz normal.vcf.gz output.vcf.gz")
}

tumor_path  <- args[1]
normal_path <- args[2]
output_path <- args[3]

if (!file.exists(tumor_path)) stop("Tumor VCF not found: ", tumor_path)
if (!file.exists(normal_path)) stop("Normal VCF not found: ", normal_path)

read_vcf_header <- function(path) {
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        gzfile(path, open = "rt")
    } else {
        file(path, open = "rt")
    }
    on.exit(close(con))

    header <- character()
    repeat {
        line <- readLines(con, n = 1L, warn = FALSE)
        if (length(line) == 0L) break
        if (startsWith(line, "#")) {
            header <- c(header, line)
        } else {
            break
        }
    }

    if (!any(startsWith(header, "#CHROM"))) {
        stop("Could not find #CHROM header in ", path)
    }

    header
}


read_vcf <- function(path) {
  
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)

  header_lines <- character()
  header_line_number <- 0L
  col_header <- NULL

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) break

    header_line_number <- header_line_number + 1L

    if (startsWith(line, "#")) {
      header_lines <- c(header_lines, line)
    }

    if (startsWith(line, "#CHROM\t")) {
      col_header <- line
      break
    }
  }

  if (is.null(col_header)) {
    stop("Could not find the #CHROM header in: ", path)
  }

  col_names <- strsplit(col_header, "\t", fixed = TRUE)[[1]]

  data <- data.table::fread(
    file = path,
    sep = "\t",
    header = FALSE,
    skip = header_line_number,
    na.strings = NULL,
    quote = "",
    check.names = FALSE,
    showProgress = FALSE
  )

  if (ncol(data) != length(col_names)) {
    stop(
      "Column count mismatch in ", path,
      ". Header has ", length(col_names),
      " columns, but data has ", ncol(data), " columns."
    )
  }

  data.table::setnames(data, col_names)


  data
}



vcf_header <- read_vcf_header(tumor_path)

tumor <- read_vcf(tumor_path)
normal <- read_vcf(normal_path)

# Confirm required VCF columns exist.
key_cols <- c("#CHROM", "POS", "REF", "ALT")

if (!all(key_cols %in% names(tumor))) {
    stop("Tumor VCF is missing one or more required columns: ",
         paste(key_cols, collapse = ", "))
}

if (!all(key_cols %in% names(normal))) {
    stop("Normal VCF is missing one or more required columns: ",
         paste(key_cols, collapse = ", "))
}

# Normalize key-column types.
for (dt in list(tumor, normal)) {
    dt[, `#CHROM` := as.character(`#CHROM`)]
    dt[, POS := as.integer(POS)]
    dt[, REF := as.character(REF)]
    dt[, ALT := as.character(ALT)]
}

# Remove tumor records whose exact CHROM/POS/REF/ALT combination
# occurs in the normal VCF.
setkeyv(normal, key_cols)

tumor_filtered <- tumor[
    !normal,
    on = key_cols
]

message("Tumor records: ", nrow(tumor))
message("Normal records: ", nrow(normal))
message("Remaining records: ", nrow(tumor_filtered))

# Write an uncompressed VCF first.
tmp_vcf <- tempfile(fileext = ".vcf")
writeLines(vcf_header, tmp_vcf)

fwrite(
    tumor_filtered,
    file = tmp_vcf,
    sep = "\t",
    append = TRUE,
    col.names = FALSE,
    quote = FALSE,
    na = "."
)

# Compress and index if output is .vcf.gz.
if (grepl("\\.vcf\\.gz$", output_path, ignore.case = TRUE)) {
    if (!nzchar(Sys.which("bgzip"))) {
        stop("bgzip is required to create a BGZF-compressed VCF")
    }
    if (!nzchar(Sys.which("tabix"))) {
        stop("tabix is required to index the output VCF")
    }

    status <- system2(
        "bgzip",
        args = c("-c", shQuote(tmp_vcf)),
        stdout = output_path
    )
    if (status != 0L) stop("bgzip failed")

    status <- system2(
        "tabix",
        args = c("-f", "-p", "vcf", shQuote(output_path))
    )
    if (status != 0L) stop("tabix indexing failed")
} else {
    if (!file.copy(tmp_vcf, output_path, overwrite = TRUE)) {
        stop("Failed to write output VCF")
    }
}

unlink(tmp_vcf)

message("Wrote filtered VCF: ", output_path)
