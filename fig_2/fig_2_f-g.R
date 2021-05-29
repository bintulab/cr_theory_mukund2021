library(grid)
library(gridExtra)
library(latex2exp)
library(patchwork)

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

p1 = ksplot_phi_dhi()
p2 = pdplot_phi_dhi()
p3 = ksplot_phi_dlo()
p4 = pdplot_phi_dlo()

ggsave("fig_2f.pdf", p1, height = 3, width = 6)
ggsave("fig_2g.pdf", p2, height = 3, width = 6)