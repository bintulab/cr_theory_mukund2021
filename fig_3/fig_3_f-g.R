library(patchwork)
library(ggExtra)
library(grid)
library(gridExtra)
library(latex2exp)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

df = read.csv('../loop_simulations/3f_response_peaks.csv')
df %>%
    dplyr::mutate(adx = case_when(ad <= 0 ~ 0,
                                  TRUE ~ ad)) -> newdf
p1 = ggplot(data = newdf)
p1 = p1 + geom_raster(aes(x = ka, y = ks, fill = pk))
p1 = p1 + scale_fill_viridis()
p1 = p1 + scale_fill_viridis(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0))
p1 = p1 + coord_cartesian(xlim = c(0, 20),
                          ylim = c(0, 20),
                          expand = 0)
p1 = p1 + scale_x_continuous(breaks = c(0, 10, 20))
p1 = p1 + scale_y_continuous(breaks=c(0, 10, 20))
p1 = p1 + labs(x = TeX("$k_p$ (day$^{-1})"),
               y = TeX("$k_s$ (day$^{-1})"),
               fill = "Peak")
p1 = p1 + theme_cowplot()
p1 = p1 + theme(axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black"),
                aspect.ratio = 1)
p1

ggsave("fig_3f.pdf",
       p1,
       dpi = 600,
       height = 3,
       width = 4)


df = read.csv('../loop_simulations/20_response_peaks.csv')
df %>%
    dplyr::mutate(adx = case_when(ad <= 0 ~ 0,
                                  TRUE ~ ad)) -> newdf
newdf %>%
    dplyr::filter(ki == 0.5) -> p2df
p2 = ggplot(data = p2df)
p2 = p2 + geom_point(aes(x = ka, y = ad, color = ks, group = ks), size = 1)
p2 = p2 + scale_y_log10(breaks = c(1, 10, 100))
p2 = p2 + scale_x_continuous(limits=c(0, 20.1), breaks=c(0, 10, 20))
p2 = p2 + scale_color_viridis(limits=c(0, 20.1), breaks = c(0, 10, 20))
p2 = p2 + labs(x = TeX("$k_p$ (day$^{-1})"),
               y = "Response Half-Life (days)",
               color = TeX("$k_s$ (day$^{-1})"))
p2 = p2 + theme_cowplot()
p2


ggsave("fig_3g.pdf",
       p2,
       dpi = 600,
       height = 3,
       width = 4)

