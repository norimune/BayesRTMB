#' Wrapper function to create an RTMB_Model instance
#'
#' @param data A list of data used for modeling.
#' @param par_list A list of parameters defined by Dim().
#' @param log_prob A function to calculate the log-posterior probability.
#' @param transform An optional function to calculate transformed parameters.
#' @param generate An optional function to calculate generated quantities.
#'
#' @return An instance of the RTMB_Model class.
#' @export
rtmb_model <- function(data, par_list, log_prob,
                       transform = NULL,
                       generate = NULL) {
  model_instance <- RTMB_Model$new(
    data         = data,
    par_list     = par_list,
    log_prob     = log_prob,
    transform    = transform,
    generate     = generate
  )
  return(model_instance)
}
