# reach file coming in has intron portion and then DRIP score
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

read_input_file <- function(file.path){

    file.path <- as.character(file.path)
    message(file.path)
    split <- unlist(strsplit(file.path, "\\."))
    calc_type = split[2]
    pol_type = split[3]
    strand = split[4]
    df <- as.data.frame(read.table(file.path))
    colnames(df) <- c('chr', 'start', 'stop', 'intron_name', 'score', 'strand', 'DRIP')
    df$calc_type <- calc_type
    df$pol_type <- pol_type
    df

}

combine_input_files <- function(df.list){

    do.call('rbind', df.list)

}

drip_vs_no_drip_plots <- function(big.df){

    drip.df <- subset(big.df, DRIP > 0)
    drip.df$drip_coverage <- 'coverage'
    no_drip.df <- subset(big.df, DRIP == 0)
    no_drip.df$drip_coverage <- 'no_coverage'
    labeled.df <- rbind(drip.df, no_drip.df)
    a <- ggplot(subset(labeled.df, calc_type=='ADARediting'), aes(x=score, fill=drip_coverage)) + 
        geom_density(alpha=0.7) + facet_wrap(~pol_type) + theme_pubr() +
        scale_color_brewer(palette = "Dark2")

    b <- ggplot(subset(labeled.df, calc_type=='ADARediting'), aes(y=score, x=drip_coverage, fill=drip_coverage)) + 
        geom_boxplot() + facet_wrap(~pol_type) + theme_pubr() +
        scale_fill_brewer(palette = "Dark2") + stat_compare_means(method='t.test') +
        theme(legend.position = "none")

    c <- ggplot(subset(labeled.df, calc_type=='structureScore'), aes(x=score, fill=drip_coverage)) +
        geom_density(alpha=0.7) + theme_pubr() + scale_fill_brewer(palette = "Dark2")
    ggarrange(a, b, labels=c('A', 'B'), nrow=2, ncol=1)

    d <- ggplot(subset(labeled.df, calc_type=='structureScore'), aes(y=score, x=drip_coverage, fill=drip_coverage)) +
        geom_boxplot() + theme_pubr() + scale_fill_brewer(palette = "Dark2") + stat_compare_means(method='t.test') +
        theme(legend.position = "none")

    ggarrange(

        ggarrange(a, b, nrow=1, ncol=2, widths = c(2, 1), heights = c(1, 1)),
        ggarrange(c, d, nrow=1, ncol=2, widths = c(2, 1), heights = c(1, 1)),
        nrow=2, ncol=1, labels=c('A', 'B')


    )


}


adar_editing_vs_drip_plots <- function(big.df){

    adar.df <- subset(big.df, calc_type=='ADARediting')
    a <- ggplot(adar.df, aes(x=DRIP, y=score, color=pol_type)) + 
        geom_point(alpha=0.4) + geom_smooth(method='lm') + 
        scale_color_brewer(palette = "Dark2") + theme_pubr() +
        labs(y='Proportion edited bases')
    b <- ggplot(adar.df, aes(y=DRIP, x=pol_type)) + geom_boxplot() +
        theme_pubr() + stat_compare_means(method='t.test')
    c <- ggplot(adar.df, aes(y=score, x=pol_type)) + geom_boxplot() +
        theme_pubr() + stat_compare_means(method='t.test')
    ggarrange(a,                                                 # First row with scatter plot
          ggarrange(b, c, ncol = 2, labels = c("B", "C")), # Second row with box and dot plots
          nrow = 2, 
          labels = "A"                                        # Labels of the scatter plot
        ) 

}

SE_vs_drip_plots <- function(big.df){

    se.df <- subset(big.df, calc_type=='structureScore')
    a <- ggplot(se.df, aes(x=DRIP, y=score)) + 
        geom_point(alpha=0.4, color='dodgerblue') + geom_smooth(method='lm') +
        theme_pubr() + labs(y='Splicing Efficiency')
    a
}

ADAR_total_vs_drip_plots <- function(big.df){

    adar_total.df <- subset(big.df, calc_type=='ADAReditingTotalEdits')
    print(head(adar_total.df))
    a <- ggplot(adar_total.df, aes(x=DRIP, y=score)) + 
        geom_point(alpha=0.4, color='dodgerblue') + geom_smooth(method='lm') +
        theme_pubr() + labs(y='Total ADAR edits per intron')
    a
}

plot_big_df <- function(big.df, output.path){


    SE_plots <- SE_vs_drip_plots(big.df)
    ADAR_plots <- adar_editing_vs_drip_plots(big.df)
    ADAR_total <- ADAR_total_vs_drip_plots(big.df)
    cov_no_cov <- drip_vs_no_drip_plots(big.df)
    pdf(output.path, width=14 ,height=14)
    print(SE_plots)
    print(ADAR_total)
    print(cov_no_cov)
    print(ADAR_plots)
    dev.off()

}


main <- function(){

    input.files <- snakemake@input
    df.list <- list()
    for (i in 1:length(input.files)){
        df.list[[i]] <- read_input_file(input.files[[i]]) 
    }
    big.df <- combine_input_files(df.list)
    save.image('test.plot.Rdata')
    message(typeof(big.df))
    plot_big_df(big.df, as.character(snakemake@output))

}


if (! interactive() ){

    main()

}


