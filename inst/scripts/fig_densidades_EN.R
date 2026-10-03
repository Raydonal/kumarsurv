## English-label version of fig_densidades.R.
source("config_paths.R")
if (!exists("DIR_FIGURAS_EN")) stop("DIR_FIGURAS_EN not set")
suppressMessages({library(ggplot2); library(patchwork)})
dk <- function(y,a,k) a*k*y^(a-1)*(1-y^a)^(k-1)
qk <- function(u,a,k) (1-(1-u)^(1/k))^(1/a)
combos <- list(c(2,2), c(0.75,0.75), c(2,0.75), c(0.75,2), c(1,1))
labs <- c("(2, 2)","(0.75, 0.75)","(2, 0.75)","(0.75, 2)","(1, 1)")
cols <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#666666")
yy <- seq(0.001,0.999,length.out=400); uu <- seq(0.001,0.999,length.out=400)
dfd <- do.call(rbind, Map(function(c,l) data.frame(x=yy, val=dk(yy,c[1],c[2]), par=l), combos, labs))
dfq <- do.call(rbind, Map(function(c,l) data.frame(x=uu, val=qk(uu,c[1],c[2]), par=l), combos, labs))
dfd$par <- factor(dfd$par, levels=labs); dfq$par <- factor(dfq$par, levels=labs)

base <- theme_minimal(base_size=11, base_family="DejaVu Sans") +
  theme(panel.grid.minor=element_blank(),
        legend.position="bottom", legend.title=element_text(size=10),
        plot.title=element_text(size=11, hjust=0))

p1 <- ggplot(dfd, aes(x,val,colour=par)) + geom_line(linewidth=0.8) +
  scale_colour_manual(values=cols, name=expression((alpha*","~kappa)),
    labels=c(expression((2*","~2)),expression((0.75*","~0.75)),
             expression((2*","~0.75)),expression((0.75*","~2)),expression((1*","~1)))) +
  coord_cartesian(ylim=c(0,3)) +
  labs(x="y", y=expression(f(y*"|"*alpha*","*kappa)), title="(a) Density") + base

p2 <- ggplot(dfq, aes(x,val,colour=par)) + geom_line(linewidth=0.8) +
  scale_colour_manual(values=cols, name=expression((alpha*","~kappa)),
    labels=c(expression((2*","~2)),expression((0.75*","~0.75)),
             expression((2*","~0.75)),expression((0.75*","~2)),expression((1*","~1)))) +
  labs(x="u", y=expression(F^{-1}*(u*"|"*alpha*","*kappa)), title="(b) Quantile function") + base

g <- (p1 | p2) + plot_layout(guides="collect") & theme(legend.position="bottom")
ggsave(file.path(DIR_FIGURAS_EN, "Fig_densidades.pdf"), g, width=9, height=3.6, device=cairo_pdf)
cat("Fig_densidades (EN) ok\n")
