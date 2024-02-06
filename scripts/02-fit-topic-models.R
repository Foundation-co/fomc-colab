### This script creates a chain of topic models of FOMC minutes ----
# for the FOMC data received from Wonseong Kim and Jan Spoerer

set.seed(8675309)

# load libraries
library(tidyverse)
library(Matrix)
library(tidylda)

library(patchwork)

# prep parallel framework
library(furrr)

plan(multisession, workers = parallel::detectCores() - 1)

### load data ----
dtm_list <- 
  read_rds(
    "data-derived/dtm-list.rds"
  )

### Get an estimate of the number of topics per year ----
estimate_num_topics <- function(dtm, k_range) {
  
  est_k <- 
    map(
      k_range, 
      function(k) {
        
        # use an asymmetric alpha?
        alph <- rgeom(n = k, prob = 0.05) + .Machine$double.eps
        
        alph <- alph / sum(alph) * 0.1 * k # keep the sum same as default
        
        alph <- alph |> sort(decreasing = TRUE) # most prevalent prior first
        
        # use an asymmetric beta?
        bet <- colSums(dtm) / sum(dtm) * 0.05 * ncol(dtm) # sum same as default
        
        m <- 
          tidylda(
            data = dtm,
            k = k,
            iterations = 200,
            burnin = 150,
            alpha = alph,
            beta = bet,
            calc_r2 = TRUE
          )
        
        tibble(
          k = k,
          mean_coherence = try(mean(m$summary$coherence)),
          var_coherence = try(var(m$summary$coherence)),
          pct_high_coherence = try(sum(m$summary$coherence > 0.01) / nrow(m$summary)),
          r2 = try(m$r2),
          likelihood = try(mean(m$log_likelihood$log_likelihood[51:200])),
          pct_high_skewness = sum(
            m$beta |>
              apply(1, function(x){
                out <- moments::skewness(x)
                
                out[is.nan(out)] <- 0
                
                out
              }) > 0 
          )/ nrow(m$summary)
          
        )
      }
    ) |> 
    bind_rows()
  
  # use loess to get optimals by variable
  
  optimal_k <- try({
    
    opt_c <- loess(mean_coherence ~ k, data = est_k)
    
    opt_r2 <- loess(r2 ~ k, data = est_k)
    
    opt_ll <- loess(likelihood ~ k, data = est_k)
    
    
    tibble(
      coherence = opt_c$x[which.max(opt_c$fitted)],
      r2 = opt_r2$x[which.max(opt_r2$fitted)],
      likelihood = opt_ll$x[which.max(opt_ll$fitted)]
    )
  })
  
  
  
  list(
    est_k = est_k,
    optimal_k = optimal_k
  )
  
}


topic_estimates <- dtm_list |> 
  future_map(function(x){
    
    estimate_num_topics(dtm = x, k_range = seq(5, 200, by = 5))
    
  }, .options = furrr_options(
    globals = "estimate_num_topics", 
    packages = c("tidylda", "tidyverse", "Matrix"), 
    seed = TRUE)
  )

for (j in seq_along(topic_estimates)) {
  
  topic_estimates[[j]]$optimal_k$year <- names(topic_estimates)[j] |> as.numeric()
  
  topic_estimates[[j]]$est_k$year <- names(topic_estimates)[j] |> as.numeric()
  
}

optimal_k <- 
  do.call(
    rbind, 
    map(topic_estimates, 
        function(x){
          if (inherits(x$optimal_k, "tbl_df")) {
            x$optimal_k
          } else {
            NULL
          }
        }
    )
  )


optimal_k <- 
  optimal_k |> 
  mutate(
    nrow = sapply(dtm_list, nrow), 
    ncol = sapply(dtm_list, ncol)
  ) 

topic_estimates <-
  map(topic_estimates, function(x) x$est_k) |>
  bind_rows()

# save objects in case we need them later
write_rds(
  topic_estimates,
  file = "data-derived/topic-number-estimates.rds"
)

write_rds(
  optimal_k,
  file = "data-derived/topic-num-optimum.rds"
)

# plot optimal number of topics 
p1 <- 
  topic_estimates |>
  ggplot(aes(x = k, y = mean_coherence)) + 
  geom_point(aes(color = year |> factor()), alpha = 0.5) + 
  geom_smooth() + theme_minimal() + theme(legend.position = "none") + 
  geom_vline(aes(xintercept = mean(optimal_k$coherence)), color = "red") +
  xlab("")

p2 <- 
  topic_estimates |>
  ggplot(aes(x = k, y = r2)) + 
  geom_point(aes(color = year |> factor()), alpha = 0.5) + 
  geom_smooth() + theme_minimal() + theme(legend.position = "none") + 
  geom_vline(aes(xintercept = mean(optimal_k$r2)), color = "blue") +
  xlab("")

p3 <- 
  topic_estimates |>
  ggplot(aes(x = k, y = pct_high_skewness)) + 
  geom_point(aes(color = year |> factor()), alpha = 0.5) + 
  geom_smooth() + theme_minimal() + theme(legend.position = "none") + 
  geom_vline(aes(xintercept = mean(optimal_k$r2)), color = "blue") + 
  geom_vline(aes(xintercept = mean(optimal_k$coherence)), color = "red")


plot(p1 / p2 / p3)


### Fit topic model list ----

# I am going to choose a fixed number of topics, because it makes for an easier
# analysis. The number, as shown in p1, is the mean number of topics across
# all years based on coherence. 

# Note for the future: we could do dynamic topics and then use skewness to determine
# whether to keep newly added topics. IDK about retiring topics though... That's
# an analytical mess
choose_k <- 
  mean(optimal_k$coherence) |>
  round()

# build models unique date
# Note that we have info bleeding from the future b/c vocab was established
# annually. In practice, we'd want to have a per-document vocab procedure

tidy_docs <- read_rds("data-derived/tidy-docs.rds")

model_list <- 
  tidy_docs |>
  by(INDICES = tidy_docs$date, function(x){
    
    ids <- x$id |> 
      unique() |> 
      sort() |>
      as.character()
    
    # get dtm of just the ids at that date
    dtm <- 
      dtm_list[[x$year[1] |> as.character()]][ids, ]
    
    # prune vocab not in that document to speed computation time
    dtm <- dtm[, colSums(dtm) > 0]
    
    # remove columns with blank colnames
    dtm <- dtm[, ! colnames(dtm) %in% c(" ", "")]
    
    dtm
    
  })

# use an asymmetric alpha
alph <- rgeom(n = choose_k, prob = 0.05) + .Machine$double.eps

alph <- alph / sum(alph) * 0.1 * choose_k # keep the sum same as default

alph <- alph |> sort(decreasing = TRUE) # most prevalent prior first

# use an asymmetric beta
bet <- colSums(model_list[[1]]) / sum(model_list[[1]]) * 0.05 * ncol(model_list[[1]]) # sum same as default


model_list[[1]] <-
  tidylda(
    data = model_list[[1]],
    k = choose_k,
    iterations = 200,
    burnin = 150,
    alpha = alph,
    beta = bet,
    calc_r2 = TRUE
  )

model_list[[1]]$summary <- 
  model_list[[1]]$summary |>
  mutate(
    top_terms_lambda = model_list[[1]]$lambda |>
      apply(1, function(x){
        names(x)[order(x, decreasing = TRUE)][1:5] |>
          paste(collapse = ", ") |>
          paste0(", ...")
      }),
    skewness = model_list[[1]]$beta |>
      apply(1, function(x){
        out <- moments::skewness(x)
        
        out[is.nan(out)] <- 0
        
        out
      })
  )

for (j in 2:length(model_list)) {
  
  model_list[[j]] <-
    refit(
      object = model_list[[j - 1]],
      new_data = model_list[[j]],
      iterations = 200,
      burnin = 150,
      prior_weight = 0.7,
      calc_likelihood = TRUE,
      calc_r2 = TRUE
    )
  
  model_list[[j]]$summary <- 
    model_list[[j]]$summary |>
    mutate(
      top_terms_lambda = model_list[[j]]$lambda |>
        apply(1, function(x){
          names(x)[order(x, decreasing = TRUE)][1:5] |>
            paste(collapse = ", ") |>
            paste0(", ...")
        }),
      skewness = model_list[[j]]$beta |>
        apply(1, function(x){
          out <- moments::skewness(x)
          
          out[is.nan(out)] <- 0
          
          out
        })
    )
  
}

write_rds(
  model_list,
  file = "data-derived/model-list.rds"
)


### Construct time series matrices of prevalence and word popularity ----

model_summaries <- 
  model_list |>
  map(function(x) x$summary)

for (j in seq_along(model_summaries)) {
  model_summaries[[j]]$date <- names(model_summaries)[j] |> as.Date()
}

model_summaries <-
  model_summaries |>
  bind_rows()

model_summaries_monthly <-
  model_summaries |>
  mutate(
    date = paste(year(date), month(date) + 1, "01", sep = "-") |> 
      as.Date() - 1
  ) |>
  group_by(date, topic) |>
  summarize(
    prevalence = mean(prevalence, na.rm = TRUE),
    coherence = mean(coherence, na.rm = TRUE),
    top_terms = top_terms[1],
    top_terms_lambda = top_terms_lambda[1],
    skewness = mean(skewness, na.rm = TRUE)
  ) |>
  ungroup() |>
  group_by(topic) |>
  mutate(
    prev_roll3m = map(0:2, function(n) lag(prevalence, n)) |>
      bind_cols() |>
      rowMeans(),
    prev_roll6m = map(0:5, function(n) lag(prevalence, n)) |>
      bind_cols() |>
      rowMeans(),
    prev_roll12m = map(0:11, function(n) lag(prevalence, n)) |>
      bind_cols() |>
      rowMeans(),
    prev_pct3m = (prevalence / lag(prevalence, 3) - 1) * 100,
    prev_pct6m = (prevalence / lag(prevalence, 6) - 1) * 100,
    prev_pct12m = (prevalence / lag(prevalence, 12) - 1) * 100,
  ) |>
  ungroup()


write_rds(
  model_summaries,
  file = "data-derived/model-summaries.rds"
)

write_csv(
  model_summaries,
  file = "data-derived/model-summaries.csv",
  append = FALSE,
  col_names = TRUE
)

write_rds(
  model_summaries_monthly,
  file = "data-derived/model-summaries-monthly.rds"
)

write_csv(
  model_summaries_monthly,
  file = "data-derived/model-summaries-monthly.csv",
  append = FALSE,
  col_names = TRUE
)


# get a tidy beta and lambda for each model
# to keep things a reasonable size, only do it for the top 10 words per topic
tidy_beta <-
  model_list |>
  map(
    function(m) {
      tb <-
        tidy(m, "beta") |>
        group_by(topic) |>
        slice_max(beta, n = 10) |>
        ungroup()
      
      tb
    }
  )

for(j in seq_along(tidy_beta)) {
  tidy_beta[[j]]$date = names(tidy_beta)[j] |> as.Date()
}

tidy_beta <- 
  tidy_beta |> 
  bind_rows()

write_rds(
  tidy_beta,
  "data-derived/tidy-beta.rds"
)

write_csv(
  tidy_beta,
  file = "data-derived/tidy-beta.csv",
  append = FALSE,
  col_names = TRUE
)


tidy_lambda <-
  model_list |>
  map(
    function(m) {
      tb <-
        tidy(m, "lambda") |>
        group_by(topic) |>
        slice_max(lambda, n = 10) |>
        ungroup()
      
      tb
    }
  )

for(j in seq_along(tidy_lambda)) {
  tidy_lambda[[j]]$date = names(tidy_lambda)[j] |> as.Date()
}

tidy_lambda <- 
  tidy_lambda |> 
  bind_rows()

write_rds(
  tidy_lambda,
  "data-derived/tidy-lambda.rds"
)

write_csv(
  tidy_lambda,
  file = "data-derived/tidy-lambda.csv",
  append = FALSE,
  col_names = TRUE
)

### Calculate aggregate model summaries

agg_model_summary <-
  model_list |>
  map(function(m){
    m$summary |>
      mutate(
        coherence = case_when(
          is.na(coherence) ~ 0,
          TRUE ~ coherence
        )
      ) |>
      summarize(
        mean_coherence = mean(coherence),
        var_coherence = var(coherence),
        skew_coherence = moments::skewness(coherence),
        var_prevalence = var(prevalence),
        skew_prevalence = moments::skewness(prevalence),
        r2 = m$r2
      )
  })

for (j in seq_along(agg_model_summary)){
  agg_model_summary[[j]]$date <- names(agg_model_summary)[j] |> as.Date()
}

agg_model_summary <- bind_rows(agg_model_summary)

write_rds(
  agg_model_summary,
  file = "data-derived/agg-model-summary.rds"
)

write_csv(
  agg_model_summary,
  file = "data-derived/agg-model-summary.csv",
  append = FALSE,
  col_names = TRUE
)

