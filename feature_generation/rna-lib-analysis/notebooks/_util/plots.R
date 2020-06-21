# Plots R
#
# It defines some "generic plot" functions.
#

wuplot_boxplot <- function(dataset, xcol, ycol,
                           boxplot_options = list(),
                           fonts = list(),
                           compare_means = list(),
                           exports = list()) {
    # Defaults
    default_boxplot_options = list(
        title = FALSE,
        subtitle = FALSE,
        caption = FALSE,
        tag = FALSE,
        xlab = FALSE,
        ylab = FALSE,
        
        # legend
        legend = "none",
        
        palette = "jco",
        add = "jitter",
        
        # 
        bxp.errorbar = TRUE,
        
        # Outlier
        outlier.shape = NA,
        
        # 
        ggtheme = theme_pubr()
    )
    default_fonts = list(
        title = list(size = 12, color = "black", face = "plain"),
        subtitle = list(size = 12, color = "black", face = "plain"),
        caption = list(size = 12, color = "black", face = "plain"),
        
        xlab = list(size = 12, color = "black", face = "plain"),
        ylab = list(size = 12, color = "black", face = "plain"),
        xy.text = list(size = 12, color = "black", face = "plain"),
        
        legend.title = list(size = 12, color = "black", face = "plain"),
        legend.text = list(size = 12, color = "black", face = "plain")
    )
    default_compare_means = list(
        # Options - t.test, wilcox.test, anova, kruskal.test
        # Default - wilcox.test
        method = 't.test',
        
        comparisons = FALSE,
        
        # Label
        # - p: the p-value
        # - p.format: the formatted p-value
        # - p.signif: the significance level.
        # Options - p.signif, ..p.signif.., p.format, ..p.format.., p, ..p..
        label = "p.signif",
        size = 4
    )
    default_exports = list(
        #
        # PDF settings
        res = 72,
        width = 7, # 504 px
        height = 7, # 504 px
        
        # PNG settings
        # res = 300,
        # width = 2000,
        # height = 2000,
        
        # 
        filename = NULL
    )
    
    # Update Options if needed
    if(length(boxplot_options) != 0) {
        updated_boxplot_options = list.merge(default_boxplot_options, boxplot_options)
    } else {
        updated_boxplot_options = default_boxplot_options
    }
    updated_fonts = list.merge(default_fonts, fonts)
    updated_compare_means = list.merge(default_compare_means, compare_means)
    updated_exports = list.merge(default_exports, exports)
    
    # Plot
    plot.boxplot <- ggboxplot(dataset.filtered,
                              y = ycol,
                              x = xcol,
                              color = xcol,
                              
                              #
                              title = updated_boxplot_options[['title']],
                              subtitle = updated_boxplot_options[['subtitle']],
                              caption = updated_boxplot_options[['caption']],
                              tag = updated_boxplot_options[['tag']],
                              xlab = updated_boxplot_options[['xlab']],
                              ylab = updated_boxplot_options[['ylab']],
                              
                              # Color Palette
                              palette = updated_boxplot_options[['palette']],
                              add = updated_boxplot_options[['add']],
                              
                              # 
                              bxp.errorbar = updated_boxplot_options[['bxp.errorbar']],
                              
                              # 
                              outlier.shape = updated_boxplot_options[['outlier.shape']],
                              
                              #
                              ggtheme = updated_boxplot_options[['ggtheme']]
    ) +
        # Add a baseline
        geom_hline(yintercept = mean(dataset[, ycol]), linetype = 2)
    
    # Fonts
    for (font_type in names(updated_fonts)) {
        plot.boxplot <- plot.boxplot + font(font_type, size = updated_fonts[[font_type]][['size']],
                                            color = updated_fonts[[font_type]][['color']], face = updated_fonts[[font_type]][['face']])
    }
    
    # Pairwise Comparisons
    plot.boxplot <- plot.boxplot + stat_compare_means(
        method = updated_compare_means[['method']],
        
        # Add comparisons
        comparisons = updated_compare_means[['comparisons']],
        
        #
        label = updated_compare_means[['label']],
        size = updated_compare_means[['size']]
    )
    
    #
    print(plot.boxplot)
    
    if ( !is.null(updated_exports[['filename']]) ) {
        ggarrange(plot.boxplot, ncol = 1) %>%
            ggexport(
                res = updated_exports[['res']],
                width = updated_exports[['width']],
                height = updated_exports[['height']],
                filename = updated_exports[['filename']]
            )
    }
}


