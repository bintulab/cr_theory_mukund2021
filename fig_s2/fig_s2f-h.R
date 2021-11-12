library(patchwork)
library(ggExtra)
library(grid)
library(gridExtra)
library(latex2exp)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

df = read.csv('../loop_simulations/20_response_peaks.csv')

df %>%
    dplyr::mutate(adx = case_when(ad <= 0 ~ 0,
                                  TRUE ~ ad)) -> newdf

newdf %>%
    dplyr::filter(ki == 0.5) -> p4df
p4 = ggplot(data = p4df)
p4 = p4 + geom_point(aes(x = pk, y = adx, color = ks, group = ks), size = 1)
p4 = p4 + coord_cartesian(xlim=c(0, 1.0))
p4 = p4 + scale_y_log10(limits = c(1, 200), breaks = c(1, 10, 100))
p4 = p4 + scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0))
p4 = p4 + scale_color_viridis(limits=c(0, 20.1), breaks=c(0, 10, 20))
p4 = p4 + theme_cowplot()
p4 = p4 + labs(x = "Peak",
               y = "Response Half-Life (days)",
               color = TeX("$k_s$ (day$^{-1})"))
p4

# p4 = p4 + annotate(GeomPoint, x = 0.169, y = 2.864, shape = 4, size = 3, color = "red", stroke = 2)
ggsave("fig_s2f.pdf",
       p4,
       dpi = 600,
       height = 3,
       width = 4)


newdf %>%
    dplyr::filter(ki == 0.5) -> p4df
p4 = ggplot(data = p4df)
p4 = p4 + geom_line(aes(x = pk, y = adx, color = ka, group = ka), size = 1)
p4 = p4 + coord_cartesian(xlim=c(0, 1.0))
p4 = p4 + scale_y_log10(limits = c(1, 200), breaks = c(1, 10, 100))
p4 = p4 + scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0))
p4 = p4 + scale_color_viridis(limits=c(0, 20.1), breaks=c(0, 10, 20))
p4 = p4 + theme_cowplot()
p4 = p4 + labs(x = "Peak",
               y = "Response Half-Life (days)",
               color = TeX("$k_p$ (day$^{-1})"))
p4
ggsave("fig_s2g.pdf",
       p4,
       dpi = 600,
       height = 3,
       width = 4)



newdf %>%
    dplyr::filter(ka == 1) -> p3df
p3 = ggplot(data = p3df)
p3 = p3 + geom_smooth(aes(x = ks, y = pk2, color = ki, group = ki), size = 1, se=F)
p3 = p3 + scale_color_viridis(limits = c(0, 2.1), breaks=c(0, 1, 2))
p3 = p3 + scale_x_continuous(limits=c(0, 20), breaks=c(0, 10, 20))
# p3 = p3 + scale_y_continuous(limits = c(0, 0.5), breaks=c(0, 0.25, 0.5))
p3 = p3 + scale_y_log10(breaks = c(0.000001, 0.001, 1),
                        labels = c(TeX("10^{-6}"), 
                                   TeX("10^{-3}"), 
                                   TeX("1")))
p3 = p3 + coord_cartesian(ylim=c(0.000001, 1))
p3 = p3 + labs(x = TeX("$k_s$ (day$^{-1})"),
               y = "Second Response Peak",
               color = TeX("$k_i$ (day$^{-1})"))
p3 = p3 + theme_cowplot()
p3
ggsave("fig_s2h.pdf",
       p3,
       dpi = 600,
       height = 3,
       width = 4)