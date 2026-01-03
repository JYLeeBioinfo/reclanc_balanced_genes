make_discrim_linear_clanc <- function() {
  parsnip::set_model_engine("discrim_linear", "classification", eng = "clanc")
  parsnip::set_dependency("discrim_linear", eng = "clanc", pkg = "reclancbalgene") ###
  parsnip::set_fit(
    model = "discrim_linear",
    eng = "clanc",
    mode = "classification",
    value = list(
      interface = "formula",
      protect = c("formula", "data"),
      func = c(pkg = "reclancbalgene", fun = "clanc"), ###
      defaults = list()
    )
  )
  parsnip::set_encoding(
    model = "discrim_linear",
    eng = "clanc",
    mode = "classification",
    options = list(
      predictor_indicators = "none",
      compute_intercept = TRUE,
      remove_intercept = TRUE,
      allow_sparse_x = FALSE
    )
  )
  parsnip::set_pred(
    model = "discrim_linear",
    eng = "clanc",
    mode = "classification",
    type = "class",
    value = list(
      pre = NULL,
      post = NULL,
      func = c(fun = "predict"),
      args = list(
        object = quote(object$fit),
        new_data = quote(new_data),
        type = "class"
      )
    )
  )
}
