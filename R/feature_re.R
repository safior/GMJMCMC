#' @export print.feature.ts
# Print correlation feature
print.feature.re <- function(re){
    if (is.list(re)){
        if (re[[2]]=="") {
            fString <- paste(re[[1]])
        }
        else{
            # If the feature contains a grouping structure for nlme
            fString <- paste(re[[1]], re[[2]], sep = ", ")
        }
        return(fString)
    }
    return(re)
}