# #### ERROR CONSTRUCTORS ####
#
# ## The purpose of these functions is to create error messages within our main functions.

#' @noRd
check_contains_cols <- function(col_names, dat, arg_name) {
  missing <- col_names[!col_names %in% colnames(dat)]
  if (length(missing) > 0) {
    cli_abort(
      "{.arg {arg_name}} column(s) not found in data: {.val {missing}}"
    )
  }
}

#' @noRd
check_positive <- function(x, arg = rlang::caller_arg(x),
                           call = rlang::caller_env()) {
  if (x <= 0) {
    cli_abort("{.arg {arg}} must be a positive value above 0",
              call = call)
  }
}

#' @noRd
check_positive0 <- function(x, arg = rlang::caller_arg(x),
                            call = rlang::caller_env()) {
  if (x < 0)
    cli_abort("{.arg {arg}} must equal 0 or a positive value", call = call)
}

#' @noRd
check_positive <- function(x, arg = rlang::caller_arg(x),
                           call = rlang::caller_env()) {
  if (x <= 0){
    cli_abort("{.arg {arg}} must be a positive value above 0", call = call)
  }
}

#' @noRd
check_negative0 <- function(x, arg = rlang::caller_arg(x),
                            call = rlang::caller_env()) {
  if (x > 0){
    cli_abort("{.arg {arg}} must equal 0 or a negative value", call = call)
  }
}

#' @noRd
check_negative <- function(x, arg = rlang::caller_arg(x),
                           call = rlang::caller_env()) {
  if (x >= 0){
    cli_abort("{.arg {arg}} must be a negative value below 0", call = call)
  }
}

#' @noRd
check_NAs <- function(x, arg = rlang::caller_arg(x),
                      call = rlang::caller_env(), threshold) {

  if(threshold == "any"){

    if (any(is.na(x))){
      cli_abort("{.arg {arg}} cannot contain missing (NA) values", call = call)
    }

  }

  if(threshold == "all"){

    if (all(is.na(x))){
      cli_abort("{.arg {arg}} cannot contain only missing (NA) values", call = call)
    }

  }
}

#' @noRd
check_class <- function(x, class, arg = rlang::caller_arg(x),
                        call = rlang::caller_env()) {
  ok <- switch(class,
               numeric   = is.numeric(x),
               integer   = is.integer(x),
               character = is.character(x),
               factor    = is.factor(x),
               logical   = is.logical(x),
               Date      = lubridate::is.Date(x),
               data.frame = is.data.frame(x),
               list      = is.list(x),
               cli_abort("Unknown class '{class}' passed to check_class()",
                         call = call))
  if (!ok){
    cli_abort("{.arg {arg}} must be of class {.cls {class}}", call = call)
  }
}

#' @noRd
check_contains <- function(a, b, arg_a = rlang::caller_arg(a),
                           arg_b = rlang::caller_arg(b),
                           call = rlang::caller_env()) {
  if (any(!a %in% b)){
    cli_abort("Values in {.arg {arg_a}} cannot be found in {.arg {arg_b}}",
              call = call)
  }
}

#' @noRd
check_empty <- function(obj, type, arg = rlang::caller_arg(obj),
                        call = rlang::caller_env()) {

  empty <- switch(type,
                  data.frame = any(dim(obj) == c(0, 0)) || is.null(obj),
                  vector     = ,
                  list       = length(obj) == 0 || is.null(obj),
                  FALSE)

  if (empty){
    cli_abort("{.arg {arg}} is empty", call = call)
  }
}

#' @noRd
check_length <- function(obj, obj_len, arg = rlang::caller_arg(obj),
                         call = rlang::caller_env()) {

  if(length(obj) != obj_len){

    cli_abort("{.arg {arg}} should be of length {obj_len}", call = call)

  }

}

#' @noRd
# check_single_col <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
#   if (length(x) != 1)
#     cli_abort("{.arg {arg}} must be a single column name, not {length(x)}", call = call)
# }

check_single_col <- function(x, arg = rlang::caller_arg(x),
                             call = rlang::caller_env()) {
  if (!rlang::is_symbol(x)) {
    cli_abort(
      "{.arg {arg}} must be a single column name, not {.code {rlang::expr_deparse(x)}}",
      call = call
    )
  }
}
