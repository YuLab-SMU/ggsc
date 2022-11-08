
sc_feature <- function(object, features, 
                    dims=c(1,2), reduction=NULL, 
                    cells=NULL, slot = "data") {

    if (is.null(reduction)) {
        reduction <- SeuratObject::DefaultDimReduc(object)
    }
    if (is.null(cells)) {
        cells <- colnames(object)
    }

    dims <- paste0(SeuratObject::Key(object = object[[reduction]]), dims)
    d <- SeuratObject::FetchData(object, vars = c(dims, "ident", 
        features), cells = cells, slot = slot)
 
    d2 <- tidyr::pivot_longer(d, 4:ncol(d))
    d2$name <- factor(d2$name, features)

    ggplot(d2, aes_string(dims[1], dims[2])) + 
        scattermore::geom_scattermore(aes(color=value), 
                    pointsize=2, pixels = c(700, 700)) +
        facet_wrap(~name, ncol=3) + 
        #scale_color_gradient2(low='blue', mid='grey', high='red') + 
        theme_bw() +
        theme(axis.text = element_blank(), 
            axis.title = element_blank(), 
            axis.ticks = element_blank(), 
            panel.grid = element_blank()
        )
}

## features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
## sc_feature(pbmc, features)
