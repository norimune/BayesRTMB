test_that("message progress reports only 20 percent steps", {
  out <- capture.output({
    meter <- .rtmb_progress_meter(5, "message", label = "generated quantities")
    for (i in 1:5) meter$advance(1)
    meter$finish()
  })

  expect_identical(
    out,
    paste0("generated quantities: ", c(20, 40, 60, 80, 100), "%")
  )
})

test_that("progress modes resolve to the message implementation", {
  expect_identical(.rtmb_resolve_progress("auto"), "message")
  expect_identical(.rtmb_resolve_progress("bar"), "message")
  expect_identical(.rtmb_resolve_progress("message"), "message")
  expect_identical(.rtmb_resolve_progress("none"), "none")
})

test_that("parallel progress files retain all unread messages", {
  progress_file <- tempfile()
  .rtmb_write_progress_file(progress_file, "chain 1 started...")
  .rtmb_write_progress_file(progress_file, "chain 1: iter 200 warmup")
  .rtmb_write_progress_file(progress_file, "chain1 done")

  line_counts <- integer(1)
  out <- capture.output({
    line_counts <- .rtmb_report_progress_files(progress_file, line_counts)
  })

  expect_identical(out, c(
    "chain 1 started...",
    "chain 1: iter 200 warmup",
    "chain1 done"
  ))
  expect_identical(line_counts, 3L)

  out_again <- capture.output({
    line_counts <- .rtmb_report_progress_files(progress_file, line_counts)
  })
  expect_identical(out_again, character(0))
  expect_identical(line_counts, 3L)
})

test_that("VB progress arguments match MCMC-style controls", {
  mdl <- rtmb_lm(mpg ~ wt, data = mtcars)

  expect_true("progress" %in% names(formals(mdl$variational)))
  expect_true("progress" %in% names(formals(VB_Fit$public_methods$transformed_draws)))
  expect_true("progress" %in% names(formals(VB_Fit$public_methods$generated_quantities)))
})
