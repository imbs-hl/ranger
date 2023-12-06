
.onAttach = function(libname, pkgname) {
  if (!interactive()) {
    return()
  }
  
  threads_option <- getOption("ranger.num.threads")
  threads_env <- Sys.getenv("R_RANGER_NUM_THREADS")
  if (!is.null(threads_option)) {
    thread_string <- paste(threads_option, "threads as set by options(ranger.num.threads = N). Can be overwritten with num.threads.")
  } else if (threads_env != "") {
    thread_string <- paste(threads_env, "threads as set by environment variable R_RANGER_NUM_THREADS. Can be overwritten with num.threads.")
  } else {
    thread_string <- "2 threads (default). Change with num.threads in ranger() and predict(), options(ranger.num.threads = N) or environment variable R_RANGER_NUM_THREADS."
  }
  
  packageStartupMessage(paste("ranger", packageVersion("ranger"), "using", thread_string))
}
