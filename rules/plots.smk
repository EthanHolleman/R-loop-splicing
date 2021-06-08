

plot_input_files = expand(
            'output/map_drip/map_drip.{calc_type}.{pol_type}.{strand}.{intron_inclusion}.bed',
            strand=['fwd', 'rev'],
            pol_type=['WT'],
            calc_type=['splicingEfficiency', 'ADAReditingTotalEdits'],
            allow_missing=True
            ) + expand(
            'output/map_drip/map_drip.{calc_type}.{pol_type}.{strand}.{intron_inclusion}.bed',
            strand=['fwd', 'rev'],
            pol_type=['WT', 'slow'],
            calc_type=['ADARediting'],
            allow_missing=True
        )

plot_input_files = expand(
    plot_input_files, intron_inclusion=['noFirstLast', 'labledAllIntrons']
    )



rule make_plots:
    conda:
        '../envs/R.yml'
    input:
        lambda wildcards: plot_input_files
    output:
        'output/plots/plot.pdf'
    script:'../scripts/multi_plot.R'

