#' Load colors for tumor types plotting
#'
#' @return Loaded R-objects containing the colors for the three domains.
#' @export
#'
loadColors <- function() {

    hematoCol <<- c("#FF5543", "#ff7477", "#ffafcc", "#D4ABEF", "#abc4ff", "#c1d3fe", "#95b8d1",
                   "#50BEC9", "#9bf6ff","#7DEF73", "#caffbf","#fcf5c7",
                   "#ffee93", "#f6bc66", "#F28041"
                    )

    neuroCol <<- c(  "#abc4ff",#"#c1d3fe",
                     "#95b8d1",
                                 "#50BEC9", "#9bf6ff","#7DEF73", "#caffbf",
                                 "#ffee93", "#f6bc66", "#F28041", "#FF452C", "#ff7477", "#ffafcc",
                                 "#D4ABEF"
                    )

    solidCol <<- c("#FF452C", "#ff7477", "#f08080", "#f4978e", #"#f8ad9d",

                   "#ffafcc",  "#ffc6ff", "#FC9FFF", "#D4ABEF","#a7bed3","#95b8d1", "#bde0fe",

                   "#abc4ff","#809bce",  "#a2d2ff",  #"#e1dbd6",#

                   "#68b6ef", "#65cbe9", "#50BEC9",

                   "#43BFA2", "#48EFAF", "#4FF246", "#adf7b6", "#caffbf", "#eeffc0",

                   "#fcf5c7", "#ffee93", "#FFE566","#EDC500",

                   "#f6bc66", "#F28041"

    )

}
