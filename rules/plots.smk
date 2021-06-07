

plot_input_files = expand(
            'output/map_drip/map_drip.{calc_type}.{pol_type}.{strand}.bed',
            strand=['fwd', 'rev'],
            pol_type=['WT'],
            calc_type=['structureScore', 'ADAReditingTotalEdits']
            ) + expand(
            'output/map_drip/map_drip.{calc_type}.{pol_type}.{strand}.bed',
            strand=['fwd', 'rev'],
            pol_type=['WT', 'slow'],
            calc_type=['ADARediting']
        )


rule make_plots:
    conda:
        '../envs/R.yml'
    input:
        lambda wildcards: plot_input_files
    output:
        'output/plots/plot.pdf'
    script:'../scripts/multi_plot.R'

