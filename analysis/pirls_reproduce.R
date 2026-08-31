## =====================================================================
## analysis/pirls_reproduce.R
##
## Reproduces the PIRLS 2011 Indonesia application reported in the
## copyfair manuscript.
##
## DATA NOT INCLUDED: the PIRLS 2011 International Database must be
## obtained separately (non-commercial research use) from
## https://timssandpirls.bc.edu/pirls2011/international-database.html
## Download "P11_SPSSData_pt2.zip (IDN-ZAF)" and extract asaidnr3.sav.
##
## Usage:
##   Rscript analysis/pirls_reproduce.R path/to/asaidnr3.sav
## =====================================================================

args <- commandArgs(trailingOnly = TRUE)
sav_path <- if (length(args) >= 1) args[1] else "asaidnr3.sav"

if (!file.exists(sav_path)) {
  stop("File not found: ", sav_path, "\n",
       "Download the PIRLS 2011 International Database (Indonesia) from\n",
       "https://timssandpirls.bc.edu/pirls2011/international-database.html\n",
       "and pass the path to asaidnr3.sav as the first argument.")
}

library(copyfair)
library(haven)

cat("=== copyfair applied to PIRLS 2011 Indonesia ===\n\n")

df <- read_sav(sav_path)
cat("Loaded", nrow(df), "students,", ncol(df), "columns\n\n")

## --- 1. Identify raw multiple-choice variables and their answer keys ---
mc_vars <- grep("^R[0-9]{2}[A-Z][0-9]{2}M$", names(df), value = TRUE)

get_correct_code <- function(var) {
  labs <- attr(df[[var]], "labels")
  opt_labs <- labs[names(labs) %in% c("A", "B", "C", "D", "A*", "B*", "C*", "D*")]
  correct_name <- names(opt_labs)[grepl("\\*$", names(opt_labs))]
  if (length(correct_name) == 0) return(NA_real_)
  unname(opt_labs[correct_name])
}
correct_key_list <- sapply(mc_vars, get_correct_code)
mc_vars <- mc_vars[!is.na(correct_key_list)]
correct_key_list <- correct_key_list[!is.na(correct_key_list)]

cat("Multiple-choice items with identified answer key:", length(mc_vars), "\n\n")

## --- 2. Build a wide response matrix per booklet (rotated-booklet design) ---
booklets <- sort(unique(as.character(df$IDBOOK)))
class_map <- setNames(as.character(df$IDCLASS), as.character(df$IDSTUD))

all_pairs <- list()
for (bk in booklets) {
  sub <- df[as.character(df$IDBOOK) == bk, ]
  if (nrow(sub) < 10) next

  ids <- paste0("s", as.character(sub$IDSTUD))
  R <- matrix(NA_real_, nrow = nrow(sub), ncol = length(mc_vars),
              dimnames = list(ids, mc_vars))
  for (v in mc_vars) {
    vals <- as.numeric(sub[[v]])
    vals[!vals %in% 1:4] <- NA
    R[, v] <- vals
  }

  key <- setNames(correct_key_list, mc_vars)
  pd <- copyfair_pairs(R, key, min_common = 5, n_ability_bins = 5)
  if (nrow(pd) > 0) {
    pd$booklet <- bk
    all_pairs[[bk]] <- pd
  }
}

pair_df <- do.call(rbind, all_pairs)
cat("Total directed pairs tested:", nrow(pair_df), "\n\n")

## --- 3. Flag and test ---
flagged <- copyfair_flag(pair_df, q_threshold = 0.05, min_denom = 5)
cat("Pairs flagged (q < .05):", sum(flagged$flagged), "\n\n")

## --- 4. External validity check: same-classroom enrichment ---
raw_id <- function(x) sub("^s", "", x)
flagged$source_class <- class_map[raw_id(flagged$source)]
flagged$copier_class <- class_map[raw_id(flagged$copier)]
flagged$same_class <- flagged$source_class == flagged$copier_class

cat("Baseline same-classroom rate (all tested pairs):\n")
print(prop.table(table(flagged$same_class)))

top_flagged <- flagged[flagged$flagged, ]
cat("\nSame-classroom rate among flagged pairs (q<.05):\n")
print(prop.table(table(top_flagged$same_class)))
