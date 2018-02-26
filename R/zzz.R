.onAttach <- function(...) {
    packageStartupMessage(paste0(
        "MOSS v", utils::packageDescription("MOSS")$Version,
        ": model one-step survival"
    ))
}
