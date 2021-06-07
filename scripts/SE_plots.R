library(ggplot2)
library(ggpubr)


range01 <- function(x){
    (x-min(x))/(max(x)-min(x))
}

read_input <- function(file.path){

    df <- as.data.frame(read.table(file.path))
    colnames(df) <- c('chr', 'start', 'stop', 'intron_name', 'SE', 'strand', 'DRIP')
    df <- subset(df, SE != 'NA')  # remove introns with no SE
    df$SE <- as.numeric(df$SE)
    df[df == '.'] <- 0
    df$DRIP <- as.numeric(df$DRIP)
    df

}


plot_SE_vs_DRIP_scatter <- function(df){

    a <- ggplot(df, aes(x=DRIP, y=SE)) + geom_point(alpha=0.4, color='dodgerblue') +
            theme_pubr() + labs(x='DRIPc Signal', y='Splicing Efficiency')
    b <- ggplot(df, aes(x=SE)) + geom_density() + theme_pubr()
    ggarrange(a, b, nrow=2, ncol=1)

}


main <- function(){

    input.file <- as.character(snakemake@input)
    output.file <- as.character(snakemake@output)

    df <- read_input(input.file)
    plt <- plot_SE_vs_DRIP_scatter(df)
    ggsave(output.file, plt, dpi=500)


}

if (! interactive()){
    main()
}