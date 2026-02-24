#' Display the MsiDataset Logo
#'
#' This function displays the logo and version information for the
#' msidataset package.
#'
#' @return None. This function is called for its side effect of printing
#' the logo and messages.
#' @author Xiaoyu Xu \email{XIAOYU021@@e.ntu.edu.sg}
#' @export
#' @examples
#' msidataset_logo()

msidataset_logo <- function() {
  
  message("Thank you for using msidataset!")
  message(
    "Version ",
    msidataset_version,
    " (",
    update_date,
    ")"
  )
  message("More information: msidataset.tidymass.org")
  
  cat(
    c(
      "  __  __  _____ _____   _____        _            _        ",
      " |  \\/  |/ ____|_   _| |  __ \\      | |          | |       ",
      " | \\  / | (___   | |   | |  | | __ _| |_ __ _ ___| |_ ___  ",
      " | |\\/| |\\___ \\  | |   | |  | |/ _` | __/ _` / __| __/ __| ",
      " | |  | |____) |_| |_  | |__| | (_| | || (_| \\__ \\ |_\\__ \\ ",
      " |_|  |_|_____/|_____| |_____/ \\__,_|\\__\\__,_|___/\\__|___/ ",
      "                                                           ",
      "                                                           "
    ),
    sep = "\n"
  )
}