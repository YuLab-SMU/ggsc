.check_aes_nm <- function(object, data, prefix='label', aes.character = "colour"){
    if (prefix %in% colnames(data)){
        default_mapping <- aes(!!rlang::sym(prefix))
    }else{
        cli::cli_warn(c("'ident' is not in the {.cls {class(object)}}.",
                       "You can set `mapping = aes(color = AnnotationColumnName)` manually."))   
        default_mapping <- aes(NULL)
    }
    names(default_mapping) <- aes.character
    
    return(default_mapping)

}


.add_aes_ <- function(x, prefix = 'label', aes.character = 'x'){
    if (!aes.character %in% names(x)){
        new.aes <- aes(!!rlang::sym(prefix))
        names(new.aes) <- aes.character
        x <- modifyList(new.aes, x)
    }
    return(x)
}

.add_aes <- function(x, prefix = c('label', 'value'), aes.character = c('x', 'y')){
    for (i in seq(length(prefix))){
        x <- .add_aes_(x, prefix = prefix[i], aes.character = aes.character[i])
    }
    return(x)
}



