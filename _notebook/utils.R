# Utils R
#
# It defines some "util" functions.
#


#' Fix the "levels" of a dataframe.
#'
#' It fixes the "levels" issue of a "test" data set, to avoid the issue of "type of predictors in new data do not match".
#'
#' @param reference The "reference" data frame.
#' @param to_fixed The "target" data frame whose "levels" needs to be fixed.
#'
#' @return Fixed dataframe.
#'
fixLevels <- function(reference, to_fixed) {
    #
    for (col in names(to_fixed)) {
        #
        # Only update if it has "fewer" levels
        if ( col %in% names(reference) ) {
            if(length(levels(to_fixed[ ,col ])) <= length(levels(reference[ ,col ])) ) {
                levels(to_fixed[ ,col ]) = levels(reference[ ,col ])
            } else {
                levels(reference[ ,col ]) = levels(to_fixed[ ,col ])
            }
        }
    }

    return(list(ref = reference, fixed = to_fixed))
}


#' Parse a list of "string elements" from a given string.
#'
#' For each of the "string elements", it need to help trim the "space" before and after.
#'
#' @param orig_str The given string.
#'
#' @return THe list of "string elements".
#'
str_to_list <- function(orig_str) {
    #
    if (orig_str == '') {
        return(c())
    }

    #
    element_list = as.list(strsplit(orig_str, ",")[[1]])
    element_list = lapply(element_list, trimws)     # Remove "space"

    return(element_list)
}