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
#' \doi{10.2333/jbhmk.27.12}
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

#' Debate Simulation Data
#'
#' @description
#' A simulated dataset containing individual satisfaction scores and related
#' variables from a 3-person group debate experiment. This dummy data was
#' created by the package author to demonstrate hierarchical (multilevel) modeling.
#'
#' @format A data frame with 300 rows and 6 columns:
#' \describe{
#'   \item{group}{Group identifier. Each group consists of 3 members.}
#'   \item{sat}{Self-reported satisfaction level of the individual after the debate.}
#'   \item{talk}{Amount of talk or utterances by the individual during the debate.}
#'   \item{perf}{Evaluated performance of the group debate (group-level variable).}
#'   \item{skill}{Evaluated conversation skill of the individual.}
#'   \item{cond}{Experimental condition indicator (e.g., 0 = control, 1 = treatment).}
#' }
#'
#' @source
#' Simulated dummy data created by the package author.
#'
#' @keywords datasets
"debate"

#' Social Skills Training Data
#'
#' @description
#' A small dataset containing repeated outcome measurements from a social
#' skills training study. The outcome was measured at four time points for
#' participants assigned to a control or intervention condition.
#'
#' @format A data frame with 12 rows and 8 columns:
#' \describe{
#'   \item{ID}{Participant identifier.}
#'   \item{time1}{Outcome score at the first measurement occasion.}
#'   \item{time2}{Outcome score at the second measurement occasion.}
#'   \item{time3}{Outcome score at the third measurement occasion.}
#'   \item{time4}{Outcome score at the fourth measurement occasion.}
#'   \item{a}{Training condition indicator: 0 = control group, 1 = intervention group.}
#'   \item{b}{Qualitative baseline skill classification.}
#'   \item{age}{Participant age.}
#' }
#'
#' @source
#' Simulated dummy data created by the package author.
#'
#' @keywords datasets
"training"

