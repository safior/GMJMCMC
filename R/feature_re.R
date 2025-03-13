#' @export print.feature.ts
print.feature.re <- function(re){
    #print(re)
    if (is.list(re)){
        #fString <- ""
        #for (i in 1:length(re)){
        #    fString <- paste(fString, re[[i]], sep=", ")
        #}
        fString <- paste(re[[1]], re[[2]], sep = ", ")
        return(fString)
    }
    return(re)
}