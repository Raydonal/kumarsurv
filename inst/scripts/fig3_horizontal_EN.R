## English-label version of fig3_horizontal.R. Writes Fig3.pdf into DIR_FIGURAS_EN.
source("config_paths.R")  # resolves dados/, resultados/, figuras/
if (!exists("DIR_FIGURAS_EN")) stop("DIR_FIGURAS_EN not set")
d<-read.csv("dda_platina.csv")
yr<-as.integer(substr(d$epiweek,1,4)); wk<-as.integer(substr(d$epiweek,7,8)); tnum<-yr+(wk-1)/52
cairo_pdf(file.path(DIR_FIGURAS_EN, "Fig3.pdf"), width=9, height=3.4, family="DejaVu Sans")
par(mar=c(3.4,3.9,0.8,1.0), mgp=c(2.3,0.7,0), cex.axis=0.9, las=1)
plot(tnum, d$cases, type="h", col="grey45", lwd=1.1, xlab="Year",
     ylab="Weekly ADD records", xaxs="i", ylim=c(0,max(d$cases)*1.02))
dev.off()
cat("Fig3 (EN) ok\n")
