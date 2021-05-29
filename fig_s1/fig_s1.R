library(grid)
library(dplyr)
library(latex2exp)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

plot_fontsize = 12

df = read.csv('../gillespie_simulations/gillespie_data.csv')

df %>%
    dplyr::filter(ka == 1 & ks == 1) %>%
    dplyr::filter(pprod %in% c(100, 500, 1000)) %>%
    dplyr::filter(pdeg %in% c(0.5, 5, 50)) -> df_smol

p = ggplot(aes(x = Time, y = Fraction.On, color = factor(A.0.)),
           group = factor(A.0.),
           data =
               df_smol)
p = p + scale_color_viridis("A(0)",
                            discrete = TRUE,
                            begin = 0,
                            end = 0.9)
p = p + geom_line(size = 1)
p = p + facet_grid(pprod ~ pdeg, labeller = label_parsed)
p = p + coord_cartesian(xlim = c(0, 10), ylim = c(0, 1))
p = p + scale_x_continuous(breaks = c(0, 5, 10))
p = p + scale_y_continuous(breaks = c(0, 0.5, 1))
p = p + xlab("Time (days)")
p = p + ylab("Fraction with Active Promoter")
p = p + theme_bw(base_size = plot_fontsize)
p = p + theme(text = element_text(colour = "black"))

p


# everything below from https://stackoverflow.com/questions/36941197/overall-label-for-facets

# Labels
labelR = "CR Production Rate (AU/day)"
labelT = "CR Degradation Rate (1/day)"

# Get the ggplot grob
z <- ggplotGrob(p)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(z$layout, grepl("strip-r", name), select = t:r)
posT <- subset(z$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips,
# and a new row on top of current top strips
width <- z$widths[max(posR$r)]    # width of current right strips
height <- z$heights[min(posT$t)]  # height of current top strips

z <- gtable_add_cols(z, width, max(posR$r))
z <- gtable_add_rows(z, height, min(posT$t) - 1)

# Construct the new strip grobs
stripR <-
    gTree(name = "Strip_right", children = gList(rectGrob(gp = gpar(
        col = NA, fill = "grey85"
    )),
    textGrob(
        labelR,
        rot = -90,
        gp = gpar(fontsize = plot_fontsize, col = "grey10")
    )))

stripT <-
    gTree(name = "Strip_top", children = gList(rectGrob(gp = gpar(
        col = NA, fill = "grey85"
    )),
    textGrob(
        labelT, gp = gpar(fontsize = plot_fontsize, col = "grey10")
    )))

# Position the grobs in the gtable
z <-
    gtable_add_grob(
        z,
        stripR,
        t = min(posR$t) + 1,
        l = max(posR$r) + 1,
        b = max(posR$b) + 1,
        name = "strip-right"
    )
z <-
    gtable_add_grob(
        z,
        stripT,
        t = min(posT$t),
        l = min(posT$l),
        r = max(posT$r),
        name = "strip-top"
    )

# Add small gaps between strips
z <- gtable_add_cols(z, unit(1 / 5, "line"), max(posR$r))
z <- gtable_add_rows(z, unit(1 / 5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(z)

z

ggsave("fig_s1a.pdf",
       z,
       height = 4,
       width = 6)



df = read.csv('../gillespie_simulations/gillespie_single_cell_trajectories.csv')

df %>%
    dplyr::mutate(gene_scaled = 130 * gene) -> df

p = ggplot(data = df,
           aes(x = time))
p = p + geom_step(aes(y = gene_scaled, color = "Gene"),
                  direction = "hv",
                  size = 1)
p = p + geom_line(aes(y = prot, color = "Protein"), size = 0.5)
p = p + geom_line(aes(y = sstate), size = 0.5, color = "black")
p = p + facet_grid(ks ~ pdeg, labeller = label_parsed)
p = p + coord_cartesian(xlim = c(0, 10), ylim = c(0, 150))
p = p + scale_y_continuous(breaks = c(0, 50, 100, 150))
p = p + scale_x_continuous(breaks = c(0, 5, 10))
p = p + xlab("Time (days)")
p = p + ylab("Gene/Protein Level (AU)")
p = p + theme_bw(base_size = plot_fontsize)
p

# everything below from https://stackoverflow.com/questions/36941197/overall-label-for-facets

# Labels
labelR = "Silencing Rate (1/day)"
labelT = "CR Degradation Rate (1/day)"

# Get the ggplot grob
z <- ggplotGrob(p)

# Get the positions of the strips in the gtable: t = top, l = left, ...
posR <- subset(z$layout, grepl("strip-r", name), select = t:r)
posT <- subset(z$layout, grepl("strip-t", name), select = t:r)

# Add a new column to the right of current right strips,
# and a new row on top of current top strips
width <- z$widths[max(posR$r)]    # width of current right strips
height <- z$heights[min(posT$t)]  # height of current top strips

z <- gtable_add_cols(z, width, max(posR$r))
z <- gtable_add_rows(z, height, min(posT$t) - 1)

# Construct the new strip grobs
stripR <- gTree(name = "Strip_right", children = gList(rectGrob(gp = gpar(
    col = NA, fill = "grey85"
)),
textGrob(
    labelR,
    rot = -90,
    gp = gpar(fontsize = plot_fontsize, col = "grey10")
)))

stripT <- gTree(name = "Strip_top", children = gList(rectGrob(gp = gpar(
    col = NA, fill = "grey85"
)),
textGrob(
    labelT, gp = gpar(fontsize = plot_fontsize, col = "grey10")
)))

# Position the grobs in the gtable
z <-
    gtable_add_grob(
        z,
        stripR,
        t = min(posR$t) + 1,
        l = max(posR$r) + 1,
        b = max(posR$b) + 1,
        name = "strip-right"
    )
z <-
    gtable_add_grob(
        z,
        stripT,
        t = min(posT$t),
        l = min(posT$l),
        r = max(posT$r),
        name = "strip-top"
    )

# Add small gaps between strips
z <- gtable_add_cols(z, unit(1 / 5, "line"), max(posR$r))
z <- gtable_add_rows(z, unit(1 / 5, "line"), min(posT$t))

# Draw it
grid.newpage()
grid.draw(z)

z

ggsave("fig_s1b.pdf",
       z,
       height = 4,
       width = 6)

# read in data
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
df = read.csv('../gillespie_simulations/gillespie_protein_data.csv')
df %>%
    dplyr::mutate(prot_scaled = prot * pprod/pdeg) %>%
    dplyr::arrange(ks) -> df_new
df_new %>%
    dplyr::filter(pprod == 5000) %>%
    dplyr::filter(pdeg >= 50) -> phi_dhi
df_new %>%
    dplyr::filter(pprod == 5000) %>%
    dplyr::filter(pdeg <= 5) -> phi_dlo

# handy functions
fix_pdeg = function(df, pdval) {
    df %>%
        dplyr::filter(pdeg == pdval) -> df_new
    return(df_new)
}
fix_pdeg_hi = function(df) { return(fix_pdeg(df, 100)) }
fix_pdeg_lo = function(df) { return(fix_pdeg(df, 1)) }
fix_ks = function(df) { 
    df %>%
        dplyr::filter(ks == 1.0) -> df_new
    return(df_new)
}

# high production rate, high degradation rate
ksplot_phi_dhi = function() {
    df_smol = fix_pdeg_hi(phi_dhi)
    p = ggplot()
    p = p + geom_density(
        data = df_smol,
        aes(
            x = prot_scaled,
            color = factor(ks),
            fill = factor(ks)
        ),
        size = 1,
        alpha = 0
    )
    kslab = TeX("k_s")
    p = p + scale_color_viridis_d(kslab, end = 0.9)
    p = p + scale_fill_viridis_d(kslab, end = 0.9)
    p = p + scale_y_continuous(breaks = c(0, 0.05, 0.1), limits = c(0, 0.11))
    p = p + scale_x_continuous(breaks = c(0, 35, 70), limits = c(0, 70))
    p = p + xlab("Final CR Protein Level (AU)")
    p = p + ylab("Probability")
    p = p + theme_cowplot()
    p = p + guides(fill = guide_legend(override.aes = list(alpha = 0.75)))
    p = p + geom_vline(xintercept = 25, linetype="dashed", color="black", size = 1.25)
    
    df_smol %>%
        dplyr::group_by(ks) %>%
        dplyr::summarise(
            n_prots = n(),
            n_thresh = sum(prot > 0.5),
            p_thresh = 100 * n_thresh / n_prots,
            l_thresh = sprintf("%0.1f%%", p_thresh)
        ) -> df_lab
    pons = df_lab$p_thresh
    lvals = df_lab$l_thresh
    ycoords = c(0.0105, 0.0283, 0.055)
    cvals = viridis(3, begin = 0, end = 0.9)
    annotation_df = data.frame(cval = rev(cvals),
                               yval = ycoords,
                               lval = rev(lvals))
    p1 = p + geom_text(
        data = annotation_df,
        aes(
            y = yval,
            label = lval),
        size = 4.5,
        x = 50,
        show.legend = FALSE
    )
    return(p1)
}
pdplot_phi_dhi = function() {
    df_smol = fix_ks(phi_dhi)
    p = ggplot()
    p = p + geom_density(
        data = df_smol,
        aes(x = prot_scaled,
            color = factor(pdeg),
            fill = factor(pdeg)),
        size = 1,
        alpha = 0
    )
    p = p + scale_color_viridis_d("Deg. Rate", end = 0.9)
    p = p + scale_fill_viridis_d("Deg. Rate", end = 0.9)
    p = p + scale_y_continuous(breaks = c(0, 0.09, 0.18), limits = c(0, 0.18))
    p = p + scale_x_continuous(breaks = c(0, 75, 150), limits = c(0, 150))
    p = p + xlab("Final CR Protein Level (AU)")
    p = p + ylab("Probability")
    p = p + theme_cowplot()
    p = p + guides(fill = guide_legend(override.aes = list(alpha = 0.75)))
    # p = p + ggtitle("Prod = 5000, ks = 1")

    
    df_smol %>%
        dplyr::group_by(pdeg) %>%
        dplyr::summarise(
            n_prots = n(),
            n_thresh = sum(prot > 0.5),
            p_thresh = 100 * n_thresh / n_prots,
            l_thresh = sprintf("%0.1f%%", p_thresh)
        ) -> df_lab
    pons = df_lab$p_thresh
    lvals = df_lab$l_thresh
    xcoords = c(105, 55, 20)
    ycoords = c(0.025, 0.035, 0.07)
    cvals = viridis(3, begin = 0, end = 0.9)
    annotation_df = data.frame(cval = rev(cvals),
                               xval = xcoords,
                               yval = ycoords,
                               lval = lvals)
    p1 = p + geom_text(
        data = annotation_df,
        aes(x = xval,
            y = yval,
            label = lval),
        size = 4.5,
        show.legend = FALSE
    )
    return(p1)
}    

# high production rate, low degradation rate
ksplot_phi_dlo = function() {
    df_smol = fix_pdeg_lo(phi_dlo)
    p = ggplot()
    p = p + geom_density(
        data = df_smol,
        aes(x = prot_scaled,
            color = factor(ks),
            fill = factor(ks)),
        size = 1,
        alpha = 0
    )
    kslab = TeX("k_s")
    p = p + scale_color_viridis_d(kslab, end = 0.9)
    p = p + scale_fill_viridis_d(kslab, end = 0.9)
    p = p + scale_x_continuous(limits = c(0, 6000), breaks = c(0, 3000, 6000))
    p = p + scale_y_continuous(limits = c(0, 0.004), breaks = c(0, 0.002, 0.004))
    p = p + xlab("Final CR Protein Level (AU)")
    p = p + ylab("Probability")
    p = p + theme_cowplot()
    p = p + guides(fill = guide_legend(override.aes = list(alpha = 0.75)))
    p = p + geom_vline(xintercept = 2500, linetype="dashed", color="black", size = 1.25)
    # p = p + ggtitle("Prod = 5000, Deg = 1")

    
    df_smol %>%
        dplyr::group_by(ks) %>%
        dplyr::summarise(
            n_prots = n(),
            n_thresh = sum(prot > 0.5),
            p_thresh = 100 * n_thresh / n_prots,
            l_thresh = sprintf("%0.1f%%", p_thresh)
        ) -> df_lab
    pons = df_lab$p_thresh
    lvals = df_lab$l_thresh
    ycoords = c(0.0008, 0.0005, 0.0033)
    xcoords = c(750, 3500, 5000)
    cvals = viridis(3, begin = 0, end = 0.9)
    annotation_df = data.frame(cval = rev(cvals),
                               xval = xcoords,
                               yval = ycoords,
                               lval = rev(lvals))
    p1 = p + geom_text(
        data = annotation_df,
        aes(y = yval,
            x = xval,
            label = lval),
        size = 4.5,
        show.legend = FALSE
    )
    return(p1)
}   
pdplot_phi_dlo = function() {
    df_smol = fix_ks(phi_dlo)
    p = ggplot()
    p = p + geom_density(
        data = df_smol,
        aes(x = prot_scaled,
            color = factor(pdeg),
            fill = factor(pdeg)),
        size = 1,
        alpha = 0
    )
    p = p + scale_color_viridis_d("Deg. Rate", end = 0.9)
    p = p + scale_fill_viridis_d("Deg. Rate", end = 0.9)
    p = p + scale_x_continuous(limits = c(0, 1.1e4), breaks = c(0, 5000, 10000))
    p = p + scale_y_continuous(limits = c(0, 0.002), breaks = c(0, 0.001, 0.002))
    p = p + xlab("Final CR Protein Level (AU)")
    p = p + ylab("Probability")
    p = p + theme_cowplot()
    p = p + guides(fill = guide_legend(override.aes = list(alpha = 0.75)))

    
    df_smol %>%
        dplyr::group_by(pdeg) %>%
        dplyr::summarise(
            n_prots = n(),
            n_thresh = sum(prot > 0.5),
            p_thresh = 100 * n_thresh / n_prots,
            l_thresh = sprintf("%0.1f%%", p_thresh)
        ) -> df_lab
    pons = df_lab$p_thresh
    lvals = df_lab$l_thresh
    ycoords = c(0.000275, 0.000375, 0.0017)
    xcoords = c(7500, 3500, 1250)
    cvals = viridis(3, begin = 0, end = 0.9)
    annotation_df = data.frame(cval = rev(cvals),
                               xval = xcoords,
                               yval = ycoords,
                               lval = lvals)
    p1 = p + geom_text(
        data = annotation_df,
        aes(y = yval,
            x = xval,
            label = lval),
        size = 4.5,
        show.legend = FALSE
    )
    return(p1)
}    

p3 = ksplot_phi_dlo()
p4 = pdplot_phi_dlo()

ggsave("fig_s1c.pdf", p3, height = 3, width = 6)
ggsave("fig_s1d.pdf", p4, height = 3, width = 6)
