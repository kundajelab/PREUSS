# Helper R
#
# It includes several helper functions.

#' Prep the dataframe.
#'
#' Several tasks to do:
#'
#' - Drop some columns, by \code{columns_to_drop}
#' - Correct some columns, by \code{columns_to_correct}
#' - Filter out records whose "predict" column has "NA" as value
#'
#' @param orig_data The "original" data frame.
#' @param column_to_predict The "column name" of "predict" column.
#' @param columns_to_keep The list of "column names" to be kept. If `columns_to_keep` is set, ignore `columns_to_drop`. 
#' @param columns_to_drop The list of "column names" to be dropped.
#' @param column_to_predict The list of "column names" to be corrected.
#'
#' @return Processed dataframe.
#'
prepData <- function(orig_data, column_to_predict="editing_level", 
                     columns_to_drop=NULL, columns_to_keep=NULL, 
                     columns_to_correct=c()) {
    # Drop the columns
    if (!is.null(columns_to_keep)) {
        # Also keep "predict column"
        columns_keep = c(column_to_predict, columns_to_keep)
        new_data = orig_data[ , (names(orig_data) %in% columns_keep)]
    } else if (!is.null(columns_to_drop)) {
        new_data = orig_data[ , !(names(orig_data) %in% columns_to_drop)]
    }
    
    # --- Criteria to help drop "invalid" rows
    # Remove records whose "editing level" is "NA"
    new_data = new_data[ is.na(new_data[ , column_to_predict])==FALSE, ]
    
    # --- Criteria to help drop "invalid" columns
    # Drop any columns that are composed of 'NA' entirely 
    new_data <- new_data[, colSums(is.na(new_data)) < nrow(new_data)]
    # Drop any columns with '0' variance 
    # remove_cols <- nearZeroVar(new_data, names = TRUE, 
    #                            freqCut = 2, uniqueCut = 0.0001)
    # all_cols <- names(new_data)
    # new_data <- new_data[ , setdiff(all_cols, remove_cols)]
    
    # Correct the columns
    # Convert "None" to "NA" for some columns which are supposed to be "number"
    # Avoid the "NA" error
    # data <- data.frame(lapply(data, function(x) {
    #             gsub("None", "NA", x)
    #         }))
    for (col in columns_to_correct) {
        new_data[ , col] <- suppressWarnings(as.numeric(as.character(new_data[ , col])))
    }
    
    # If a column is "not" numeric / factor, try to convert it to numeric
    for (col in names(new_data)) {
        if (typeof(new_data[ , col]) != 'integer') {
            new_data[ , col] <- suppressWarnings(as.numeric(as.character(new_data[ , col])))
        }
    }
    
    # Impute missing values in predictor data using proximity from randomForest.
    x = new_data[ ,2:ncol(new_data) ]
    y = new_data[ ,1 ]
    #
    tryCatch({
            new_data = rfImpute(x, y, iter = 5, ntree = 300)
        },
        warning = function(w) {},
        error = function(e) {
            print(e)
            
            # Change the "first" column name to "y"
            local_data = new_data
            colnames(local_data)[1] <- "y"
            new_data <<- local_data
        },
        finally = {}
    ) # END tryCatch

    #
    return(new_data)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



