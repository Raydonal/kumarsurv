source("config_paths.R")  # Rodada 33: resolve dados/, resultados/, figuras/
d<-read.csv("dda_platina.csv")
yr<-as.integer(substr(d$epiweek,1,4)); wk<-as.integer(substr(d$epiweek,7,8)); tnum<-yr+(wk-1)/52
cairo_pdf("Fig3.pdf", width=9, height=3.4, family="DejaVu Sans")
par(mar=c(3.4,3.9,0.8,1.0), mgp=c(2.3,0.7,0), cex.axis=0.9, las=1)
plot(tnum, d$cases, type="h", col="grey45", lwd=1.1, xlab="Ano",
     ylab="Registros semanais de DDA", xaxs="i", ylim=c(0,max(d$cases)*1.02))
dev.off()
cat("Fig3 horizontal ok\n")
