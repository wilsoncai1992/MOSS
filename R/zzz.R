.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "MOSS v", utils::packageDescription("MOSS")$Version,
    ": Model One-Step Survival"
  ))
}
