library(data.table)
library(MASS)
library(glmmTMB)


# ========================================= function definitions ==================================================

fit_standard <- function(data,alt,ref,output_rds){
    
    print('Fitting mixed effects model with all unscaled predictors...')

    model <- glmmTMB(
    cbind(alt, ref) ~ mismatches + median_baseQ  + GGC_upstream + CGG_downstream + GGT_upstream + TGG_downstream + GC_content + (1 | sampleID) + median_mapQ,
    data = data,
    ziformula = ~ 1,
    family = binomial(),
    verbose = TRUE)

    print(summary(model))
    saveRDS(model, file = output_rds)
    print(paste0('Model saved to: ', output_rds))

    return(list('convergence' = model$fit$convergence))
}






fit_scaled <- function(data,alt,ref,output_rds,scaling_params_rds){

    print('Due to convergence problems reverting to scaled mixed effects model...')
    
    data$GC_content <- scale(data$GC_content)
    gc_center <- attr(data$GC_content, "scaled:center")
    gc_scale  <- attr(data$GC_content, "scaled:scale")

    data$mismatches <- scale(data$mismatches)
    mismatch_center <- attr(data$mismatches, "scaled:center")
    mismatch_scale  <- attr(data$mismatches, "scaled:scale")

    data$median_mapQ <- scale(data$median_mapQ)
    mapq_center <- attr(data$median_mapQ, "scaled:center")
    mapq_scale  <- attr(data$median_mapQ, "scaled:scale")

    data$median_baseQ <- scale(data$median_baseQ)
    baseq_center <- attr(data$median_baseQ, "scaled:center")
    baseq_scale  <- attr(data$median_baseQ, "scaled:scale")



    scaling_parameters <- list(
        GC_content = list(
            center = gc_center,
            scale = gc_scale
        ),
        mismatches = list(
            center = mismatch_center,
            scale = mismatch_scale
        ),
        median_mapQ = list(
            center = mapq_center,
            scale = mapq_scale
        ),
        median_baseQ = list(
            center = baseq_center,
            scale = baseq_scale
        )
    )

    saveRDS(
        scaling_parameters,
        file = scaling_params_rds
    )


    model <- glmmTMB(
    cbind(alt, ref) ~ mismatches + median_baseQ  + GGC_upstream + CGG_downstream + GGT_upstream + TGG_downstream + GC_content + (1 | sampleID) + median_mapQ,
    data = data,
    ziformula = ~ 1,
    family = binomial(),
    verbose = TRUE)

    print(summary(model))
    saveRDS(model, file = output_rds)
    print(paste0('Model saved to: ', output_rds))
    print(paste0('Scaling parameters saved to: ', scaling_params_rds))

    return(list('convergence' = model$fit$convergence))
}







fit_small_sample <- function(data,alt,ref,output_rds,scaling_params_rds){

    print('Due to small sample size and/or convergence problems reverting to limited scaled fixed effects model...')
    
    data$GC_content <- scale(data$GC_content)
    gc_center <- attr(data$GC_content, "scaled:center")
    gc_scale  <- attr(data$GC_content, "scaled:scale")

    data$mismatches <- scale(data$mismatches)
    mismatch_center <- attr(data$mismatches, "scaled:center")
    mismatch_scale  <- attr(data$mismatches, "scaled:scale")

    scaling_parameters <- list(
        GC_content = list(
            center = gc_center,
            scale = gc_scale
        ),
        mismatches = list(
            center = mismatch_center,
            scale = mismatch_scale
        )
    )

    saveRDS(
        scaling_parameters,
        file = scaling_params_rds
    )

    model <- glmmTMB(
    cbind(alt, ref) ~ mismatches + GGC_upstream + CGG_downstream + GGT_upstream + TGG_downstream + GC_content,
    data = data,
    ziformula = ~ 1,
    family = binomial(),
    verbose = TRUE)

    print(summary(model))
    saveRDS(model, file = output_rds)
    print(paste0('Model saved to: ', output_rds))
    print(paste0('Scaling parameters saved to: ', scaling_params_rds))

    return(list('convergence' = model$fit$convergence))
}
# =================================================================================================================


output_dir <- "fitted_error_model"

if (dir.exists(output_dir)) {
    # Remove previous fitted models and scaling parameters
    old_files <- list.files(output_dir, full.names = TRUE)
    if (length(old_files) > 0) {
        file.remove(old_files)
    }
} else {
    dir.create(output_dir, showWarnings = FALSE)
}


args <- commandArgs(trailingOnly = TRUE)
negative_control_data <- args[1]
noise_level <- as.numeric(args[2])
alt_mode <- args[3]   # 'specific' or 'generic'
fitted_model_rds <- args[4]

# read in the data
data <- fread(negative_control_data, header = TRUE, sep = ',')
# clean the data
data <- na.omit(data)
data <- data[data$ref_counts > 0,]


# choose the alt mode
if(alt_mode == 'specific') {
    
    background_error_rate <- data$specific_alt_counts/(data$specific_alt_counts + data$ref_counts)
    data_cleaned <- data[background_error_rate <= noise_level,]
    alt <- as.numeric(data_cleaned$specific_alt_counts)
    data_cleaned$alt <- alt
    ref <- as.numeric(data_cleaned$ref_counts)
    data_cleaned$ref <- ref

    
    if(length(unique(data$sampleID)) <= 5){  #if the sample size is <= 5 run fit_small_sample()

        res <- fit_small_sample(data = data_cleaned,
                    alt = alt,
                    ref = ref,
                    output_rds = paste0('fitted_error_model/', fitted_model_rds, '.rds'),
                    scaling_params_rds = 'fitted_error_model/scaling_params.rds')
    }else{  # if the sample size is > 5, run fit_standard()

        res <- fit_standard(data = data_cleaned,
                    alt = alt,
                    ref = ref,
                    output_rds = paste0('fitted_error_model/', fitted_model_rds, '.rds'))

        if(res$convergence != 0){  # if the standard model did not converge, run fit_scaled()

            res <- fit_scaled(data = data_cleaned,
                            alt = alt,
                            ref = ref,
                            output_rds = paste0('fitted_error_model/', fitted_model_rds, '.rds'),
                            scaling_params_rds = 'fitted_error_model/scaling_params.rds')
        }
    }
    



}else if (alt_mode == 'generic') {
   background_error_rate <- data$general_alt_counts/(data$general_alt_counts + data$ref_counts)
   data_cleaned <- data[background_error_rate <= noise_level,]
   alt <- as.numeric(data_cleaned$general_alt_counts)
   data_cleaned$alt <- alt
   ref <- as.numeric(data_cleaned$ref_counts)
   data_cleaned$ref <- ref

   # now here based on the samplesize check invoke the corresponding functions
   if(length(unique(data$sampleID)) <= 5){  #if the sample size is <= 5 run fit_small_sample()

        res <- fit_small_sample(data = data_cleaned,
                    alt = alt,
                    ref = ref,
                    output_rds = paste0('fitted_error_model/', fitted_model_rds, '.rds'),
                    scaling_params_rds = 'fitted_error_model/scaling_params.rds')
    }else{  # if the sample size is > 5, run fit_standard()

        res <- fit_standard(data = data_cleaned,
                    alt = alt,
                    ref = ref,
                    output_rds = paste0('fitted_error_model/', fitted_model_rds, '.rds'))

        if(res$convergence != 0){  # if the standard model did not converge, run fit_scaled()

            res <- fit_scaled(data = data_cleaned,
                            alt = alt,
                            ref = ref,
                            output_rds = paste0('fitted_error_model/', fitted_model_rds, '.rds'),
                            scaling_params_rds = 'fitted_error_model/scaling_params.rds')
        }
    }
}


cat('








')


if(res$convergence == 0){
    print('Model fitted successfully!')
}else{
    print('Model failed to converge. Check, negative control data.')
}

