clanc_fit <- function(expression, class_data, classes) {
  class_means <- collapse::fmean(expression, classes, na.rm = TRUE)
  overall_means <- colMeans(expression, na.rm = TRUE)

  pooled_sds <- calculate_pooled_sd(expression, class_data, classes)

  class_stats <- calculate_class_stats(
    classes, class_means, overall_means, pooled_sds
  )

  abs_stats <- abs(class_stats)
  sign_stats <- sign(class_stats)
  # up_only   <- ifelse(class_stats > 0,  abs(class_stats), -Inf)
  # down_only <- ifelse(class_stats < 0,  abs(class_stats), -Inf)

  # class_data_half <- class_data
  # class_data_half$active <- as.integer(floor(class_data$active / 2))

  # ranks_up <- t(apply(up_only,   1, \(x) rank(-x, ties.method = "min")))
  # ties_up  <- t(apply(ranks_up,  1, duplicated))
  
  # ranks_dw <- t(apply(down_only, 1, \(x) rank(-x, ties.method = "min")))
  # ties_dw  <- t(apply(ranks_dw,  1, duplicated))
  
  # selected_up <- select_genes(up_only,   ranks_up, ties_up, class_data_half)
  # selected_dw <- select_genes(down_only, ranks_dw, ties_dw, class_data_half)
  
  # selected <- rbind(selected_up, selected_dw)
  
  ranks <- t(apply(abs_stats, 1, \(x) rank(-x, ties.method = "min")))
  ties <- t(apply(ranks, 1, duplicated))
  selected <- select_genes(abs_stats, ranks, ties, class_data,sign_stats)

  centroids <- create_centroids(selected, class_means, overall_means)

  pooled_sds <- pooled_sds[match(rownames(centroids), names(pooled_sds))]

  tall_centroids <- matrix(unlist(centroids), ncol = 1) |> as.data.frame()
  tall_centroids$gene <- rownames(centroids)
  tall_centroids$class <- rep(colnames(centroids), each = nrow(centroids))
  tall_centroids <- tall_centroids[, c(3, 2, 1)]
  w_sd <- merge(tall_centroids, pooled_sds, by.x = "gene", by.y = "row.names")
  out <- merge(w_sd, class_data, by = "class")
  out$class <- factor(out$class, levels = levels(class_data$class))
  colnames(out) <- c(
    "class", "gene", "expression", "pooled_sd", "active", "prior"
  )

  list(centroids = out)
}

calculate_pooled_sd <- function(expression, class_data, classes) {
  df <- nrow(expression) - nrow(class_data)
  sum_sq_diffs <- tapply(data.frame(expression), classes, sum_sq_diff)
  sum_sq_diffs <- do.call(rbind, sum_sq_diffs)
  out <- sqrt(colSums(sum_sq_diffs / df))
  names(out) <- colnames(expression)
  out
}

sum_sq_diff <- function(exp) {
  colSums(scale(exp, center = TRUE, scale = FALSE)^2)
}

#' @importFrom collapse %c/% %r-% %r/%
calculate_class_stats <- function(classes,
                                  class_means,
                                  overall_means,
                                  class_pooled_sds) {
  mks <- sqrt(1 / table(classes) - 1 / length(classes))
  class_means <- as.matrix(class_means)
  class_means %r-%
    overall_means %r/%
    class_pooled_sds %c/%
    as.data.frame(mks)$Freq
}

select_genes <- function(abs_stats, ranks, ties, class_data, sign_stats) {
  df <- data.frame(
    gene = colnames(abs_stats),
    abs = matrix(t(abs_stats), ncol = 1),
    sign = matrix(t(sign_stats), ncol = 1),
    rank = matrix(t(ranks), ncol = 1),
    tie = matrix(t(ties), ncol = 1),
    class = rep(class_data$class, each = ncol(abs_stats)),
    reserved = FALSE,
    n_win = 0
  )

  # class_data is where `active` lives
  class_data_half <- class_data 
  class_data_half$active <- as.integer(floor(class_data$active/2))
  df_pos <- merge(df, class_data_half, by = "class") |>
    selection_recurse(sign_int=1)

  df_neg <- merge(df, class_data_half, by = "class") |>
    selection_recurse(sign_int=-1)

  df <- rbind(df_pos,df_neg)
  df <- df[df$gets, ]

  data.frame(
    class = factor(df$class, levels = levels(class_data$class)),
    gene = df$gene
  )
}



#' @importFrom rlang .data
selection_recurse <- function(df,sign_int = 1) {
  if (all(df$reserved) || all(df$n_win >= df$active)) {
    return(NULL)
  }

  df <- df |>
    dplyr::filter(!.data$reserved, .data$n_win < .data$active) |>
    dplyr::arrange(.data$gene, .data$rank, dplyr::desc(.data$abs)) |>
    dplyr::mutate(
      n = seq_len(dplyr::n()),
      win = (.data$n == 1)&(.data$sign == sign_int), # win only when the sign matches
      .by = "gene"
    ) |>
    dplyr::arrange(.data$class, .data$rank) |>
    dplyr::mutate(
      n_win = .data$n_win + cumsum(.data$win),
      gets = .data$win & (.data$n_win <= .data$active),
      .by = "class"
    ) |>
    dplyr::mutate(reserved = any(.data$gets), .by = "gene") |>
    dplyr::mutate(n_win = max(.data$n_win), .by = "class")
  
  rbind(df, selection_recurse(df))
}

create_centroids <- function(selected, class_means, overall_means) {
  lvls <- levels(selected$class)
  n_levels <- length(lvls)
  class <- class_means[, match(selected$gene, colnames(class_means))]
  overall <- overall_means[match(selected$gene, names(overall_means))]
  mm <- matrix(rep(overall, n_levels))

  update_col <- function(class_means, winner, overall_means) {
    idx <- which(winner == lvls)
    new <- rep(overall_means, length(class_means))
    new[idx] <- class_means[idx]
    new
  }

  out <- mapply(update_col, class, selected$class, overall)
  out <- t(out)
  colnames(out) <- lvls
  out
}

# Predictors to use are only 'discovered' *after* fitting
# Therefore the default way of dealing with predictors isn't sufficient
make_new_ptypes <- function(fit, processed) {
  genes <- unique(fit$centroids$gene)
  # make a dummy tibble w proper dims
  preds <- matrix(1.0, nrow = 1, ncol = length(genes))
  colnames(preds) <- genes
  preds <- tibble::as_tibble(preds)
  preds <- preds[-1, ]
  list(predictors = preds, outcomes = processed$blueprint$ptypes$outcomes)
}
