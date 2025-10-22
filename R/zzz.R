.onAttach <- function(libname, pkgname) {
  
  pkg_version <- utils::packageVersion("faers")
  
  welcome_msg <- c(
    paste0("Welcome to 'faers' package!"),
    paste0("================================================================"),
    paste0("You are using faers version ", pkg_version),
    paste0(""),
    paste0("Project home : https://github.com/WangLabCSU/faers"),
    paste0("================================================================")
  )
  
  packageStartupMessage(paste(welcome_msg, collapse = "\n"))
  
  if (interactive()) {
    tryCatch({
      current_version <- package_version(pkg_version)
    }, error = function(e) {
    })
  }
}
