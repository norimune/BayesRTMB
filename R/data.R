#' Beverage Preference Data
#'
#' @description
#' A dataset containing preference ratings for various beverages.
#' This data was used to demonstrate a random effect model in metric
#' multidimensional unfolding.
#'
#' @format A data frame (or matrix) with [INSERT NUMBER OF ROWS] rows and [INSERT NUMBER OF COLUMNS] columns.
#' \describe{
#'   \item{AJ}{Apple Juice preference rating}
#'   \item{BT}{Black Tea preference rating}
#'   \item{CF}{Coffee preference rating}
#'   \item{CL}{Cola preference rating}
#'   \item{GT}{Green Tea preference rating}
#'   \item{OJ}{Orange Juice preference rating}
#'   \item{OT}{Oolong Tea preference rating}
#' }
#'
#' @source
#' Adachi, K. (2000). A Random Effect Model in Metric Multidimensional Unfolding.
#' \emph{The Japanese Journal of Behaviormetrics}, 27(1), 12-23.
#' \url{https://doi.org/10.2333/jbhmk.27.12}
#'
#' @keywords datasets
"beverage"


#' Big Five Personality Traits Data
#'
#' @description
#' A dataset containing responses to 20 items measuring the Big Five personality traits.
#' This dataset was originally collected by the package author.
#' The items are arranged in a cyclical order representing the five factors:
#' Extraversion, Neuroticism, Openness, Agreeableness, and Conscientiousness.
#'
#' @format A data frame with 10 rows and 20 columns:
#' \describe{
#'   \item{BF1}{Extraversion item 1}
#'   \item{BF2}{Neuroticism item 1}
#'   \item{BF3}{Openness item 1}
#'   \item{BF4}{Agreeableness item 1}
#'   \item{BF5}{Conscientiousness item 1}
#'   \item{BF6}{Extraversion item 2}
#'   \item{BF7}{Neuroticism item 2}
#'   \item{BF8}{Openness item 2}
#'   \item{BF9}{Agreeableness item 2}
#'   \item{BF10}{Conscientiousness item 2}
#'   \item{BF11}{Extraversion item 3}
#'   \item{BF12}{Neuroticism item 3}
#'   \item{BF13}{Openness item 3}
#'   \item{BF14}{Agreeableness item 3}
#'   \item{BF15}{Conscientiousness item 3}
#'   \item{BF16}{Extraversion item 4}
#'   \item{BF17}{Neuroticism item 4}
#'   \item{BF18}{Openness item 4}
#'   \item{BF19}{Agreeableness item 4}
#'   \item{BF20}{Conscientiousness item 4}
#' }
#'
#' @source
#' Originally collected by the package author.
#'
#' @keywords datasets
"BigFive"

#' Group Discussion Simulation Data
#'
#' @description
#' A simulated dataset containing individual satisfaction scores and related
#' variables from a 3-person group discussion experiment. This dummy data was
#' created by the package author to demonstrate hierarchical (multilevel) modeling.
#'
#' @format A data frame with 300 rows and 6 columns:
#' \describe{
#'   \item{group}{Group identifier. Each group consists of 3 members.}
#'   \item{satisfaction}{Self-reported satisfaction level of the individual after the discussion.}
#'   \item{talk}{Amount of talk or utterances by the individual during the discussion.}
#'   \item{performance}{Evaluated performance of the group discussion (group-level variable).}
#'   \item{skill}{Evaluated conversation skill of the individual.}
#'   \item{condition}{Experimental condition indicator (e.g., 0 = control, 1 = treatment).}
#' }
#'
#' @source
#' Simulated dummy data created by the package author.
#'
#' @keywords datasets
"discussion"


#' Stations Data
#'
#' @description
#' A dataset representing distances or other relational metrics between stations.
#'
#' @format A data frame or matrix containing relational data between stations.
#'
#' @keywords datasets
"stations"
