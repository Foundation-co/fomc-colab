### This script creates a list of document term matrices, one per year, ----
# for the FOMC data received from Wonseong Kim and Jan Spoerer

# load libraries
library(tidyverse)
library(tidytext)
library(Matrix)

# prep parallel framework
library(furrr)

plan(multisession, workers = parallel::detectCores() - 1)

### Import data ----
fomc_sent_tag <- readxl::read_xlsx(
  "data-raw/FOMC-sentences.xlsx", 
  sheet = 1,
  col_types = c("numeric", "date", "text", "text", "text", "text")
)

fomc_sent_untagged <- readxl::read_xlsx(
  "data-raw/FOMC-sentences.xlsx", 
  sheet = 2,
  col_types = c("numeric", "date", "text", "text", "text", "text")
)

# verify that I only need untagged data to continue
if (sum(fomc_sent_tag$id %in% fomc_sent_untagged$id) != nrow(fomc_sent_tag)) {
  stop("There are sentences in `fomc_sent_tag` not in `fomc_sent_untagged`")
}

### Filter vocabulary ----

# To construct our FOMC vocabulary we are going to:
#   1. *not* remove stop words at first
#   2. stem each word
#   3. go by year and...
#     a) remove words appearing in half or more sentences
#     b) remove words appearing in [2] or fewer sentences
#   4. construct bigrams and filter based on...
#     a) remove bigrams that don't start with a retained stem
#     b) remove bigrams that appear in half or more sentences
#     c) remove bigrams that appear in [5] or fewer sentences

fomc_tidy_docs <- 
  fomc_sent_untagged |>
  mutate(
    year = year(date)
  ) |>
  select(
    id,
    year,
    date,
    sentence
  ) |>
  unnest_tokens(
    word,
    sentence
  ) |>
  mutate( # remove any punctuation 
    word = word |> str_replace_all("[^\\w\\s]+", "")
  ) |> 
  mutate( # Stem using Porter's word stemmer
    stem = tokenizers::tokenize_word_stems(word, language = "porter") |> unlist()
  ) |>
  filter( # Remove any words/stems that don't contain alphabetical character
    ! str_detect(stem, "^[^a-z|A-Z]+$") 
  ) |>
  filter( # remove any "words" that are blanks or spaces
    stem != "" | str_detect(stem, "^\\s*$")
  ) |> 
  mutate( # remove any numbers (e.g. "8per cent" should just be "per cent")
    stem = str_replace_all(stem, "[0-9]+", "")
  )

fomc_tidy_docs <-
  fomc_tidy_docs |>
  group_by(id) |>
  mutate(
    bigram = paste(stem, lead(stem), sep = " ")
  ) |> 
  ungroup()

# filter #3, unigrams, above
# filter #4, bigrams, above
dtm_list <- 
  by(fomc_tidy_docs, INDICES = fomc_tidy_docs$year, function(x) x) |>
  future_map(
    function(x) {
      
      # construct unigram dtm
      unigram_dtm <- x |>
        count(
          id,
          stem
        ) |> 
        cast_sparse(
          id,
          stem,
          n
        )
      
      tf <-
        tibble(
          word = colnames(unigram_dtm),
          term_freq = colSums(unigram_dtm),
          doc_freq = colSums(unigram_dtm > 0)
        )
      
      vocab <- 
        tf |> 
        filter(
          doc_freq > 2 & doc_freq < (nrow(unigram_dtm) / 2)
        ) |>
        select(
          word
        ) 
      
      unigram_dtm <- unigram_dtm[, vocab[[1]]]
      
      # construct bigram DTM
      bigram_dtm <- 
        x |>
        filter(
          stem %in% vocab[[1]] & ! str_detect(bigram, "NA")
        ) |>
        count(
          id,
          bigram
        ) |> 
        cast_sparse(
          id,
          bigram,
          n
        )
      
      tf <-
        tibble(
          word = colnames(bigram_dtm),
          term_freq = colSums(bigram_dtm),
          doc_freq = colSums(bigram_dtm > 0)
        )
      
      vocab <- 
        tf |> 
        filter(
          doc_freq > 5 & doc_freq < (nrow(bigram_dtm) / 2)
        ) |>
        select(
          word
        ) 
      
      bigram_dtm <- bigram_dtm[, vocab[[1]]]
      
      
      # combine dtms into one
      ids <- x$id |> 
        unique() |> 
        sort() |>
        as.character() # id is numeric and below will break without this
      
      if ((setdiff(ids, rownames(unigram_dtm)) |> length()) > 0) {
        
        missing_ids <- setdiff(ids, rownames(unigram_dtm))
        
        rows_to_add <- 
          Matrix(
            0,
            nrow = length(missing_ids),
            ncol = ncol(unigram_dtm),
            sparse = TRUE
          )
        
        rownames(rows_to_add) <- missing_ids
        
        unigram_dtm <- 
          rbind(unigram_dtm, rows_to_add)
        
      }
      
      if ((setdiff(ids, rownames(bigram_dtm)) |> length()) > 0) {
        
        missing_ids <- setdiff(ids, rownames(bigram_dtm)) 
        
        rows_to_add <- 
          Matrix(
            0,
            nrow = length(missing_ids),
            ncol = ncol(bigram_dtm),
            sparse = TRUE
          )
        
        rownames(rows_to_add) <- missing_ids
        
        bigram_dtm <- 
          rbind(bigram_dtm, rows_to_add)
        
      }
      
      dtm <- cbind(unigram_dtm[ids, ], bigram_dtm[ids, ])
      
      dtm
      
    }
  ) 

### Plot diagnostics... just in case :) ----
dtm_diagnostics <- 
  tibble(
    year = names(dtm_list) |> as.numeric(),
    vocab = sapply(dtm_list, ncol),
    sentences = sapply(dtm_list, nrow),
    volume = sapply(dtm_list, sum),
    empty_docs = sapply(dtm_list, function(x){
      (rowSums(x) == 0) |> sum()
    }),
    empty_vocab = sapply(dtm_list, function(x){
      (colSums(x) == 0) |> sum()
    })
  )
  

dtm_diagnostics |> 
  ggplot(aes(x = year, y = vocab)) + 
  geom_line() +
  ggtitle("Number of unique tokens per year")

dtm_diagnostics |> 
  ggplot(aes(x = year, y = sentences)) + 
  geom_line() +
  ggtitle("Number of sentences per year")

if (sum(dtm_diagnostics$empty_docs) > 0) {
  warning("found empty sentences. Check your dtms")
}

if (sum(dtm_diagnostics$empty_vocab) > 0) {
  warning("found empty vocab. Check your dtms")
}

### Save output for later ----
write_rds(
  dtm_list,
  "data-derived/dtm-list.rds"
)

write_rds(
  fomc_tidy_docs,
  "data-derived/tidy-docs.rds"
)

write_rds(
  dtm_diagnostics,
  "data-derived/dtm-diagnostics.rds"
)
