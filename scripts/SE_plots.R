library(ggplot2)
library(ggpubr)

read_input <- function(file.path){

    df <- as.data.frame(read.table(file.path))
    colnames(df) <- c('intron_name', 'SE', 'DRIP')
    df

}


plot_SE_vs_DRIP_scatter <- function(df){

    ggplot(df, aes(x=DRIP, y=SE)) + geom_point(alpha=0.4, color='dodgerblue') +
            theme_pubr() + labs(x='DRIPc Signal', y='Splicing Efficiency')

}