library(grid)
library(gridExtra)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

df = read.csv('./recruitment_graded_response.csv')
p = ggplot(data=df, aes(x=Release.Duration, y=Recruitment.Duration))
p = p + geom_raster(aes(fill=Fraction))
p = p + scale_fill_viridis(limits=c(0,1), breaks=c(0, 0.5, 1))
p = p + facet_grid(. ~ Population)
p = p + coord_fixed(ratio=1, xlim=c(0, 5), ylim=c(0, 5), expand = FALSE)
p = p + xlab("Release Duration") + ylab("Recruitment Duration")
# p = p + theme_cowplot(font_family="Arial")
p = p + theme_bw(base_size=14)
p = p + theme(panel.spacing = unit(1, "lines"))
p

ggsave("../fig_s2b.pdf",
       p,
       height = 4,
       width = 8)