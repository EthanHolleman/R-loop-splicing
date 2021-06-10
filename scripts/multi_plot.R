# reach file coming in has intron portion and then DRIP score
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyr)
library(plyr)
library(dplyr)
library(ggridges)


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

# Calculate 95% confidence interval for a vector of data
# returns the margin of error
confidence_interval <- function(x){

    t.score <- qt(p=0.05/2, df=length(x)-1, lower.tail=F)
    se <- sd(x) / sqrt(length(x))
    margin_error <- t.score * se

    margin_error

}


read_gene_lengths <- function(genes.path){
    gene.df <- as.data.frame(read.table(genes.path))

    colnames(gene.df) <- c('gene_chr', 'gene_start', 'gene_end', 'gene_name', 'gene_score', 'gene_strand')

    gene.df <- separate(gene.df , 4, c('gene_name', 'ensemble_id'), sep = '_')
    gene.df$gene_length <- (gene.df$gene_end - gene.df$gene_start)

    gene.df

}


read_input_file <- function(file.path){
    # Read an individual input bed file. Should be a bunch of introns in
    # bed 6 plus an extra row. Bed6 score row will be the value for
    # whatever calculation was done (specified by file wildcards)
    # and last column should be DRIP peak value that mapped to the
    # intron (may be zero). Additional information is infered from
    # the wildcard values in the filepath so ordering of these
    # wildcards in the snakemake pipeline is important and
    # should be maintained / checked before changing.

    file.path <- as.character(file.path)
    split <- unlist(strsplit(file.path, "\\."))

    # get wildcard values present in filepaths
    calc_type = split[2]
    pol_type = split[3]
    strand = split[4]
    intron_inclusion = split[5]

    # read bed as a dataframe 
    df <- as.data.frame(read.table(file.path))

    # rename included columns 
    colnames(df) <- c('chr', 'start', 'stop', 'intron_name', 'score', 'strand', 'DRIP')

    # Add wildcard values as columns for later identification and plotting
    df$calc_type <- calc_type
    df$pol_type <- pol_type
    df$intron_inclusion = intron_inclusion
    df <- separate(df, 4, c('gene_name', 'ensemble_id', 'dist_from_start', 'dist_from_end'), sep = '_')
    df$dist_from_start <- as.numeric(df$dist_from_start)
    df$dist_from_end <- as.numeric(df$dist_from_end)
    df <- df %>% add_count(gene_name, name='gene_intron_count')

    df

}


intron_editing_and_SE_vs_intron_location <- function(big.df){
    # given the dataframe from a single input file, label
    # the introns as either first in gene, last in gene,
    # or intragene
    big.df.all_introns <- subset(big.df, intron_inclusion=='labledAllIntrons')
    big.df.all_introns.wt<- subset(big.df.all_introns, pol_type=='WT')

    # convert distance from start to factor in correct order
    # big.df.all_introns$dist_from_start <- factor( big.df.all_introns.wt$dist_from_start, 
    #                                                 levels=sort(unique( big.df.all_introns.wt$dist_from_start))
    #                                             )
    big.df.all_introns.wt$dist_from_start <- as.factor(big.df.all_introns.wt$dist_from_start)

                                        
    intron_location_vs_score <- ggplot(big.df.all_introns.wt, aes(x=dist_from_start, y=score, fill=calc_type)) +
                                geom_boxplot() + facet_wrap(~calc_type, scales = "free", nrow=3) +
                                theme_pubr() + scale_fill_brewer(palette='Dark2') + 
                                theme(legend.position = "none")

    hist_introns_per_dist <- ggplot(big.df.all_introns.wt, aes(x=dist_from_start, fill=calc_type)) +
                            geom_bar() + scale_fill_brewer(palette='Dark2') + 
                                theme(legend.position = "none") + facet_wrap(~calc_type, scales = "free", nrow=3) +
                                theme_pubr()
    
    avg_num_introns_per_gene <- ggplot(big.df.all_introns.wt, aes(y=gene_intron_count, x=calc_type, fill=calc_type)) +
                                  geom_boxplot() + scale_fill_brewer(palette='Dark2') + theme_pubr() +
                                  geom_hline(yintercept=7.8,  linetype="dashed", size=1.5)  # average number introns per gene 

    ggarrange(intron_location_vs_score,
              ggarrange(hist_introns_per_dist, avg_num_introns_per_gene, nrow=1, ncol=2, labels=c('B', 'C')),
              labels=c('A'), nrow=2, ncol=1)
    


    # intron_location_vs_score <- ggplot(big.df.all_introns.wt, aes(x=reorder(dist_from_start, sort(as.numeric(dist_from_start)), y=DRIP, fill=calc_type))) +
    #                             geom_boxplot() + facet_wrap(~calc_type, scales = "free", nrow=3) +
    #                             theme_pubr() + scale_fill_brewer(palette='Dark2') + 
    #                             theme(legend.position = "none")
    

}



compare_gene_lengths <- function(gene.df, big.df){

    big.df.all_introns <- subset(big.df, intron_inclusion=='labledAllIntrons')
    big.df.all_introns.wt<- subset(big.df.all_introns, pol_type=='WT')

    big.df.merge <- merge(big.df.all_introns.wt, gene.df, by="gene_name")
    big.df.merge$intron_density <- big.df.merge$gene_intron_count / big.df.merge$gene_length

    intron_density_comparison <- ggplot(big.df.merge, aes(y=intron_density, x=calc_type, fill=calc_type)) +
                                geom_boxplot() + scale_fill_brewer(palette='Dark2') + theme_pubr()
    
    gene_length_comparison <- ggplot(big.df.merge, aes(x=gene_length, y=calc_type, fill=calc_type)) +
                              geom_density_ridges(alpha=0.7) + theme_pubr() + 
                              scale_fill_brewer(palette='Dark2') +
                              theme(legend.position = "none")

    ggarrange(intron_density_comparison, gene_length_comparison, nrow=1, ncol=2)
                
}


remove_intron_by_proximity_to_gene_ends <- function(df.merged, dist_bp=0, prop_len=0){

    included_rows <- c()
    for (i in 1:nrow(df.merged)){
        row <- df.merged[i, ]
        if (prop_len > 0){
            cut_off <- prop_len * row$gene_length
        }
        else{
            cut_off <- dist_bp
        }

        if (row$start - row$gene_start > cut_off & row$gene_end - row$stop > cut_off){
            included_rows <- c(included_rows, i)
        }

    }

    df.merged[included_rows, ]


}


compare_metrics_drip_no_drip_introns_intron_order <- function(big.df){

    big.df.all_introns <- subset(big.df, intron_inclusion=='labledAllIntrons')  # use all introns
    big.df.all_introns.wt<- subset(big.df.all_introns, pol_type=='WT')  # use wildtype results
    big.df.all_introns.wt$drip_coverage <- as.factor(big.df.all_introns.wt[, 'DRIP'] > 0)
    # label drip coverage 

   # compare splicing efficieny across genes
    big.df.all_introns.wt <- subset(big.df.all_introns.wt, calc_type=='splicingEfficiency')

    # label first and last introns
    # require having at least 2 introns ? at least 3 introns?
    # remove genes without at least 2 introns
    big.df.2_introns <- subset(big.df.all_introns.wt, gene_intron_count >= 2)
    big.df.2_introns$total_introns <- nrow(big.df.2_introns)

    boxplot_compare <- function(df, i){
        means <- aggregate(score ~  drip_coverage, df, mean)
        boxplot <- ggplot(df, aes(y=score, x=drip_coverage, fill=drip_coverage)) +
                                geom_boxplot() + theme_pubr() + scale_fill_brewer(palette='Dark2') +
                                theme(legend.position = "none") + stat_compare_means(method='t.test') +
                                stat_summary(fun=mean, colour="black", show.legend=FALSE) + 
                                geom_text(data = means, aes(label = score, y = score + 0.08)) +
                                labs(title=paste(i, 'introns from flanks excluded'))
        change_introns_count <- data.frame(introns=c('All', 'Trimmed'), 
                                           number_introns=c(unique(df$total_introns), nrow(df))
                                          )
        print(change_introns_count)
        barplot <- ggplot(change_introns_count, aes(x=introns, y=number_introns, fill=introns)) +
                    geom_bar(stat='identity') + theme_pubr() + scale_fill_brewer(palette='Dark2') +
                    theme(legend.position = "none")
        ggarrange(boxplot, barplot)
    }

    boxplots <- list()
    for (i in 0:10){
        df <- subset(big.df.2_introns, as.numeric(dist_from_start) > i & as.numeric(dist_from_end) > i)
        boxplots[[i+1]] <- boxplot_compare(df, i)
    }

    ggarrange(plotlist=boxplots)

}


introns_with_and_without_drip_coverage <- function(big.df){

    # Note slow calculations only for ADAR editing
    big.df.all_introns <- subset(big.df, intron_inclusion=='labledAllIntrons')  # use all introns
    big.df.all_introns$drip_coverage <- as.factor(big.df.all_introns[, 'DRIP'] > 0)

    ggplot(big.df.all_introns, aes(x=drip_coverage, fill=calc_type)) +
            geom_bar(position=position_dodge(),  color="black") + theme_pubr() + 
            facet_wrap(~pol_type) + scale_fill_brewer(palette="Dark2")
            labs(x='DRIPc coverage', y='Number introns')
    
}


compare_metrics_drip_no_drip_introns_intron_proximity <- function(big.df, genes.df){

    big.df.all_introns <- subset(big.df, intron_inclusion=='labledAllIntrons')  # use all introns
    big.df.all_introns.wt<- subset(big.df.all_introns, pol_type=='WT')  # use wildtype results
    big.df.all_introns.wt$drip_coverage <- as.factor(big.df.all_introns.wt[, 'DRIP'] > 0)
    big.df.all_introns.wt <- merge(big.df.all_introns.wt, genes.df, by="gene_name")
    # label drip coverage 

   # compare splicing efficieny across genes
    big.df.all_introns.wt <- subset(big.df.all_introns.wt, calc_type=='splicingEfficiency')

    # label first and last introns
    # require having at least 2 introns ? at least 3 introns?
    # remove genes without at least 2 introns
    big.df.2_introns <- subset(big.df.all_introns.wt, gene_intron_count >= 2)
    big.df.2_introns$total_introns <- nrow(big.df.2_introns)
    print(colnames(big.df.2_introns))

    boxplot_compare <- function(df, i){
        means <- aggregate(score ~  drip_coverage, df, mean)
        boxplot <- ggplot(df, aes(y=score, x=drip_coverage, fill=drip_coverage)) +
                                geom_boxplot() + theme_pubr() + scale_fill_brewer(palette='Dark2') +
                                theme(legend.position = "none") + stat_compare_means(method='t.test') +
                                stat_summary(fun=mean, colour="black", show.legend=FALSE) + 
                                geom_text(data = means, aes(label = score, y = score + 0.08)) +
                                labs(title=paste(i*100, '% introns from flanks excluded'))
        change_introns_count <- data.frame(introns=c('All', 'Trimmed'), 
                                           number_introns=c(unique(df$total_introns), nrow(df))
                                          )
        print(change_introns_count)
        barplot <- ggplot(change_introns_count, aes(x=introns, y=number_introns, fill=introns)) +
                    geom_bar(stat='identity') + theme_pubr() + scale_fill_brewer(palette='Dark2') +
                    theme(legend.position = "none")
        ggarrange(boxplot, barplot)          
    }

    boxplots <- list()
    k <- 1
    for (i in seq(0, 1, 0.1)){
        df <- remove_intron_by_proximity_to_gene_ends(big.df.2_introns, prop_len=i)
        if (nrow(df) < 1){
            break
        }
        boxplots[[k]] <- boxplot_compare(df, i)
        k <- k + 1
    }
    print(length(boxplots))
    ggarrange(plotlist=boxplots)

}


drip_peak_score_by_intron_se <- function(big.df){

    big.df.all_introns <- subset(big.df, intron_inclusion=='labledAllIntrons')  # use all introns
    big.df.all_introns.se <- subset(big.df.all_introns, calc_type=='splicingEfficiency')
    big.df.all_introns.se$dist_from_start <- as.factor(big.df.all_introns.se$dist_from_start)

    boxplot <- ggplot(big.df.all_introns.se, aes(x=dist_from_start, y=DRIP)) +
            geom_boxplot(fill='dodgerblue') + theme_pubr() +
            theme(axis.text.x = element_text(angle = 90)) +
            labs(title='Mean DRIPc peak signal and 95% CI by intron dist to gene start in introns',
                x='Introns to start of gene',
                 y='Mean DRIPc signal'
                 )
    

    # calculate mean and margin of error for DRIP signal for each intron.
    # Group introns by distance from start of gene in 
    all_introns.agg.mean <- aggregate(big.df.all_introns.se$DRIP, 
                                      list(big.df.all_introns.se$dist_from_start),
                                      mean)
    all_introns.agg.margin_error <- aggregate(big.df.all_introns.se$DRIP, 
                                      list(big.df.all_introns.se$dist_from_start),
                                      confidence_interval)
    all_introns.mean.margin_error <- merge(all_introns.agg.mean, all_introns.agg.margin_error, by='Group.1')

    colnames(all_introns.mean.margin_error) <- c('dist_from_start', 'mean', 'margin_error')
    all_introns.mean.margin_error <- all_introns.mean.margin_error[complete.cases(all_introns.mean.margin_error), ]
    all_introns.mean.margin_error$mean <- as.numeric(all_introns.mean.margin_error$mean)
    all_introns.mean.margin_error$margin_error <- as.numeric(all_introns.mean.margin_error$margin_error)

    # Remove anything with margin error > 100 this is due to low sample size
     all_introns.mean.margin_error <- subset( all_introns.mean.margin_error, margin_error < 100)

    lineplot <- ggplot(all_introns.mean.margin_error, aes(x=dist_from_start, y=mean)) +
                geom_errorbar(aes(ymin=mean-margin_error, ymax=mean+margin_error), color='dodgerblue') +
                geom_point(color='black') + theme_pubr() +
                theme(axis.text.x = element_text(angle = 90)) + 
                labs(title='Mean DRIPc peak signal and 95% CI by intron dist to gene start in introns',
                     x='Introns to start of gene',
                     y='Mean DRIPc signal')

    ggarrange(boxplot, lineplot, nrow=2, ncol=1)

}


combine_input_files <- function(df.list){

    do.call('rbind', df.list)

}



plot_big_df <- function(big.df, genes.df, output.path){

    a <- intron_editing_and_SE_vs_intron_location(big.df)
    a_1 <- introns_with_and_without_drip_coverage(big.df)
    b <- compare_gene_lengths(genes.df, big.df)
    c <- compare_metrics_drip_no_drip_introns_intron_order(big.df)
    e <- compare_metrics_drip_no_drip_introns_intron_proximity(big.df, genes.df)
    d <- drip_peak_score_by_intron_se(big.df)

    pdf(output.path, width=14 ,height=14)

    print(a)
    print(a_1)
    print(b)
    print(c)
    print(d)
    print(e)

    dev.off()

}


main <- function(){
    print(snakemake@input)
    input.files <- snakemake@input$data
    print(input.files)
    gene.file <- unlist(snakemake@input$genes)
    gene.df <- read_gene_lengths(gene.file)
    save.image('plot_image.RData')
    df.list <- list()
    for (i in 1:length(input.files)){
        df.list[[i]] <- read_input_file(input.files[[i]]) 
    }
    big.df <- combine_input_files(df.list)
    plot_big_df(big.df, gene.df, as.character(snakemake@output))
    

}


if (! interactive() ){

    main()

}


