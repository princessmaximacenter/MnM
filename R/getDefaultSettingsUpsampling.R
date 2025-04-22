getDefaultSettingsUpsampling <- function() {
  upsimpler_init_kwargs <- list(dataDF = NULL, # TO BE DEFINED AT RUNTIME
                                metadataDF = NULL, # TO BE DEFINED AT RUNTIME
                                class_col = NULL, # TO BE DEFINED AT RUNTIME
                                feature_list = NULL, # Use all features
                                metric = "manhattan",
                                precomputed_dists = FALSE,
                                verbosity = 0L)

  upsimpler_upsimple_kwargs <- list(target_sample_id = NULL, # TO BE DEFINED ON EACH OF THE N RUNS
                                    n_neighbors = 10L,
                                    n_samples = 10L,
                                    exclude_same_class = FALSE,
                                    allow_interclass_dists = FALSE,
                                    subsample_frac = 0.7,
                                    with_replacement = FALSE,
                                    rseed = 1L,
                                    scaling_strategy = "nnmaxrepel",
                                    nnmax_factor = 0.45,
                                    nnrepel_factor = 1.1,
                                    neg_strategy = "clip")

  upsamplerArgs <- list(init = upsimpler_init_kwargs,
                        upsimple = upsimpler_upsimple_kwargs,
                        maj_strategy = "het+")

  return(upsamplerArgs)
}
