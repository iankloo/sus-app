#--functions
sus_converter <- function(d, odds, evens){
  x <- sum(d[odds] - 1, 5 - d[evens]) * 2.5
  return(x)
}

log.lik.tnorm<-function(mu, sig, dat){
  if ((length(dat) == 0)) {
    return(0)
  } else {
    sum(dtnorm(dat, mean = mu, sd = sig, lower = 0, upper = 100, log=T))
  }
}

lik2<-function(parvec, dat){
  cp<-c(mean=parvec[1], s.d.=parvec[2], gamma1= parvec[3])
  dp<-cp2dp(cp,family="SN")
  LB<-psn(0,dp=dp)
  UB<-psn(100,dp=dp)
  sum(dsn(dat,dp[1], dp[2], dp[3], log = T)-log(UB-LB))
}

prior<-function(mu,muprior,sdprior=20){
  varprior<-sdprior^2
  convert<-estBetaParams(muprior,sdprior)
  alpha<-convert$alpha
  beta<-convert$beta
  dbeta(mu/100,alpha,beta,log=TRUE)
}

estBetaParams <- function(mu.scaled, var.scaled) {
  mu<-mu.scaled/100
  var<-var.scaled/100^2
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

genPlot <- function(bootdata, conf_int, num_boots, raw_dat, title) {
  
  #Cairo(width = 10, height = 8, file="bootplot.png", type="png", pointsize=12,
  #      bg = "white", canvas = "white", units = "in", dpi = 300)
  
  h<-hist(bootdata, prob = F, xlim=c(0,100))
  h.y.min <- par('usr')[3]
  h.y.max <- par('usr')[4]
  
  h<-hist(bootdata, main = title, cex.main = 1, 
          xlab = 'SUS Score', cex.lab = 0.9, col = 'gray90', prob = F,
          cex.axis = 0.8, xaxt = 'n', yaxt = 'n',
          sub = paste0("SUS respondents = ", length(raw_dat), ", bootstrap sample size = ", num_boots),
          cex.sub = 0.75, xlim=c(0,100), ylim=c(h.y.min*5,h.y.max)) # the ylim argument creates space for the scales 
  axis(1, xaxp=c(0, 100, 10), las=1, cex.axis = 0.8)
  axis(2, yaxp=c(0, round_any(h.y.max,100,f=floor), 5), las=0, cex.axis = 0.8)
  
  # add in 95% CI bounds, mean, min, and max
  min.SUS<-min(raw_dat)
  max.SUS<-max(raw_dat)
  abline(v = min.SUS, col = 'darkgray', lty=2)
  abline(v = max.SUS, col = 'darkgray', lty=2)
  abline(v = conf_int, col = 'red', lty=2)
  abline(v = mean(bootdata), col = 'blue', lty=2)
  rect(0, 0.55*h.y.min+h.y.max, 100, -1.3*h.y.min+h.y.max, density = NULL, col = "white", border = NA, lty = 1, lwd = 1)
  text(conf_int, c(h.y.max,h.y.max), labels=c(paste(round(conf_int[1],1)),paste(round(conf_int[2],1))),pos=c(2,4),cex=0.7,col=c('red','red')) 
  text(conf_int, c(h.y.max-0.72*h.y.min,h.y.max-0.72*h.y.min), labels=c('95% CI LB','95% CI UB'),pos=c(2,4),cex=0.7,col=c('red','red')) 
  text(mean(bootdata), h.y.max, labels=paste(round(mean(bootdata),1)),cex=0.7,col='blue') 
  text(mean(bootdata), h.y.max-0.75*h.y.min, labels='Mean',cex=0.7,col='blue')
  text(c(min.SUS,max.SUS), c(h.y.max+2.25*h.y.min,h.y.max+2.25*h.y.min), labels=c(paste(round(min.SUS,1)),paste(round(max.SUS,1))),pos=c(2,4),cex=0.7,col=c('grey50','grey50')) 
  text(c(min.SUS,max.SUS), c(h.y.max+1.5*h.y.min,h.y.max+1.5*h.y.min), labels=c('Min','Max'),pos=c(2,4),cex=0.7,col=c('grey50','grey50'))
  
  # Scales to interpret SUS Scores 
  rect(0, 1.3*h.y.min, 50, 0.1*h.y.min, density = NULL, col = "indianred1", border = NULL, lty = 1, lwd = 1)
  rect(50, 1.3*h.y.min, 70, 0.1*h.y.min, density = NULL, col = "lightyellow", border = NULL, lty = 1, lwd = 1)
  rect(70, 1.3*h.y.min, 100, 0.1*h.y.min, density = NULL, col = "darkseagreen2", border = NULL, lty = 1, lwd = 1)
  text(c(25,60,85),c(0.75*h.y.min,0.75*h.y.min,0.75*h.y.min),labels=c('Unacceptable','Marginal','Acceptable'),cex=0.7)
  
  rect(0, 2.5*h.y.min, 51.7, 1.3*h.y.min, density = NULL, col = "indianred1", border = NULL, lty = 1, lwd = 1)
  rect(51.7, 2.5*h.y.min, 62.6, 1.3*h.y.min, density = NULL, col = MixColor("indianred1", "lightyellow", 0.5), border = NULL, lty = 1, lwd = 1)
  rect(62.6, 2.5*h.y.min, 72.5, 1.3*h.y.min, density = NULL, col = "lightyellow", border = NULL, lty = 1, lwd = 1)
  rect(72.5, 2.5*h.y.min, 78.8, 1.3*h.y.min, density = NULL, col = MixColor("darkseagreen2", "lightyellow", 0.5), border = NULL, lty = 1, lwd = 1)
  rect(78.8, 2.5*h.y.min, 100, 1.3*h.y.min, density = NULL, col = "darkseagreen2", border = NULL, lty = 1, lwd = 1)
  text(c(25.85,57.15,67.55,75.65,89.4),c(1.9*h.y.min),labels=c('F','D','C','B','A'),cex=0.7, adj=0.5)
  
  rect(0, 3.7*h.y.min, 100, 2.5*h.y.min, density = NULL, col = "white", border = NULL, lty = 1, lwd = 1)
  rect(15.58, 3.7*h.y.min, 25.02, 2.5*h.y.min, density = NULL, col = "indianred1", border = NA, lty = 1, lwd = 1)
  rect(32.79, 3.7*h.y.min, 38.61, 2.5*h.y.min, density = NULL, col = MixColor("indianred1", "lightyellow", 0.5), border = NA, lty = 1, lwd = 1)
  rect(49.04, 3.7*h.y.min, 52.76, 2.5*h.y.min, density = NULL, col = "lightyellow", border = NA, lty = 1, lwd = 1)
  rect(70.18, 3.7*h.y.min, 72.62, 2.5*h.y.min, density = NULL, col = MixColor("darkseagreen2", "lightyellow", 0.5), border = NA, lty = 1, lwd = 1)
  rect(84.30, 3.7*h.y.min, 86.70, 2.5*h.y.min, density = NULL, col = "darkseagreen2", border = NA, lty = 1, lwd = 1)
  text(c(20.3,35.7,50.9,71.4,85.5),c(3.1*h.y.min,3.1*h.y.min,3.1*h.y.min,3.1*h.y.min,3.1*h.y.min),labels=c('Awful','Poor','OK','Good','Excellent'),cex=0.7)
  rect(0, 3.7*h.y.min, 100, 2.5*h.y.min, density = NULL, col = NULL, border = "black", lty = 1, lwd = 1)
  
  rect(0, 4.9*h.y.min, 100, 3.7*h.y.min, density = NULL, col = "white", border = NULL, lty = 1, lwd = 1)
  text(seq(from=20, to=90, by=10),rep(4.3*h.y.min, 8),labels=c(0.01,0.02,0.06,0.13,0.29,0.56,0.88,0.998),cex=0.6,font=1)
  text(-0.5,4.35*h.y.min,labels="Score Percentiles",cex=0.6,font=1, pos=4)
  
  rect(0, 6.1*h.y.min, 100, 4.9*h.y.min, density = NULL, col = NA, border = NA, lty = 1, lwd = 1)
  rect(min.SUS-1, 6*h.y.min, max.SUS+1, 5*h.y.min, density = NULL, col = "azure1", border = NA, lty = 1, lwd = 1)
  freqs<-count(raw_dat)
  text(freqs[,1],rep(5.5*h.y.min, length(freqs[,1])),labels=freqs[,2],cex=0.6,font=1)
  text(-0.5,5.6*h.y.min,labels="Frequencies",cex=0.6,font=1, pos=4)
  
  grid.text("Acceptability ranges (Bangor et al., 2008, p. 592)\nLetter grades (Sauro & Lewis, 2016, p. 204)\nAdjective ratings (Bangor et al., 2009, p. 118)\nScore percentiles (Sauro & Lewis, 2016, p. 203)",
            x = unit(0.015, "npc"), y = unit(0.042, "npc"),
            just="left", gp=gpar(fontsize=6, col="grey50"))
  
  #dev.off()
  
}

genPlotBayes <- function(bootdata, conf_int, raw_dat, title) {
  
  #Cairo(width = 10, height = 8, file="bootplot.png", type="png", pointsize=12,
  #      bg = "white", canvas = "white", units = "in", dpi = 300)
  
  h<-hist(bootdata, prob = F, xlim=c(0,100))
  h.y.min <- par('usr')[3]
  h.y.max <- par('usr')[4]
  
  h<-hist(bootdata, main = title, cex.main = 1, 
          xlab = 'SUS Score', cex.lab = 0.9, col = 'gray90', prob = F,
          cex.axis = 0.8, xaxt = 'n', yaxt = 'n',
          sub = paste0("SUS respondents = ", length(raw_dat)),
          cex.sub = 0.75, xlim=c(0,100), ylim=c(h.y.min*5,h.y.max)) # the ylim argument creates space for the scales 
  axis(1, xaxp=c(0, 100, 10), las=1, cex.axis = 0.8)
  axis(2, yaxp=c(0, round_any(h.y.max,100,f=floor), 5), las=0, cex.axis = 0.8)
  
  # add in 95% CI bounds, mean, min, and max
  min.SUS<-min(raw_dat)
  max.SUS<-max(raw_dat)
  abline(v = min.SUS, col = 'darkgray', lty=2)
  abline(v = max.SUS, col = 'darkgray', lty=2)
  abline(v = conf_int, col = 'red', lty=2)
  abline(v = mean(bootdata), col = 'blue', lty=2)
  rect(0, 0.55*h.y.min+h.y.max, 100, -1.3*h.y.min+h.y.max, density = NULL, col = "white", border = NA, lty = 1, lwd = 1)
  text(conf_int, c(h.y.max,h.y.max), labels=c(paste(round(conf_int[1],1)),paste(round(conf_int[2],1))),pos=c(2,4),cex=0.7,col=c('red','red')) 
  text(conf_int, c(h.y.max-0.72*h.y.min,h.y.max-0.72*h.y.min), labels=c('95% CI LB','95% CI UB'),pos=c(2,4),cex=0.7,col=c('red','red')) 
  text(mean(bootdata), h.y.max, labels=paste(round(mean(bootdata),1)),cex=0.7,col='blue') 
  text(mean(bootdata), h.y.max-0.75*h.y.min, labels='Mean',cex=0.7,col='blue')
  text(c(min.SUS,max.SUS), c(h.y.max+2.25*h.y.min,h.y.max+2.25*h.y.min), labels=c(paste(round(min.SUS,1)),paste(round(max.SUS,1))),pos=c(2,4),cex=0.7,col=c('grey50','grey50')) 
  text(c(min.SUS,max.SUS), c(h.y.max+1.5*h.y.min,h.y.max+1.5*h.y.min), labels=c('Min','Max'),pos=c(2,4),cex=0.7,col=c('grey50','grey50'))
  
  # Scales to interpret SUS Scores 
  rect(0, 1.3*h.y.min, 50, 0.1*h.y.min, density = NULL, col = "indianred1", border = NULL, lty = 1, lwd = 1)
  rect(50, 1.3*h.y.min, 70, 0.1*h.y.min, density = NULL, col = "lightyellow", border = NULL, lty = 1, lwd = 1)
  rect(70, 1.3*h.y.min, 100, 0.1*h.y.min, density = NULL, col = "darkseagreen2", border = NULL, lty = 1, lwd = 1)
  text(c(25,60,85),c(0.75*h.y.min,0.75*h.y.min,0.75*h.y.min),labels=c('Unacceptable','Marginal','Acceptable'),cex=0.7)
  
  rect(0, 2.5*h.y.min, 51.7, 1.3*h.y.min, density = NULL, col = "indianred1", border = NULL, lty = 1, lwd = 1)
  rect(51.7, 2.5*h.y.min, 62.6, 1.3*h.y.min, density = NULL, col = MixColor("indianred1", "lightyellow", 0.5), border = NULL, lty = 1, lwd = 1)
  rect(62.6, 2.5*h.y.min, 72.5, 1.3*h.y.min, density = NULL, col = "lightyellow", border = NULL, lty = 1, lwd = 1)
  rect(72.5, 2.5*h.y.min, 78.8, 1.3*h.y.min, density = NULL, col = MixColor("darkseagreen2", "lightyellow", 0.5), border = NULL, lty = 1, lwd = 1)
  rect(78.8, 2.5*h.y.min, 100, 1.3*h.y.min, density = NULL, col = "darkseagreen2", border = NULL, lty = 1, lwd = 1)
  text(c(25.85,57.15,67.55,75.65,89.4),c(1.9*h.y.min),labels=c('F','D','C','B','A'),cex=0.7, adj=0.5)
  
  rect(0, 3.7*h.y.min, 100, 2.5*h.y.min, density = NULL, col = "white", border = NULL, lty = 1, lwd = 1)
  rect(15.58, 3.7*h.y.min, 25.02, 2.5*h.y.min, density = NULL, col = "indianred1", border = NA, lty = 1, lwd = 1)
  rect(32.79, 3.7*h.y.min, 38.61, 2.5*h.y.min, density = NULL, col = MixColor("indianred1", "lightyellow", 0.5), border = NA, lty = 1, lwd = 1)
  rect(49.04, 3.7*h.y.min, 52.76, 2.5*h.y.min, density = NULL, col = "lightyellow", border = NA, lty = 1, lwd = 1)
  rect(70.18, 3.7*h.y.min, 72.62, 2.5*h.y.min, density = NULL, col = MixColor("darkseagreen2", "lightyellow", 0.5), border = NA, lty = 1, lwd = 1)
  rect(84.30, 3.7*h.y.min, 86.70, 2.5*h.y.min, density = NULL, col = "darkseagreen2", border = NA, lty = 1, lwd = 1)
  text(c(20.3,35.7,50.9,71.4,85.5),c(3.1*h.y.min,3.1*h.y.min,3.1*h.y.min,3.1*h.y.min,3.1*h.y.min),labels=c('Awful','Poor','OK','Good','Excellent'),cex=0.7)
  rect(0, 3.7*h.y.min, 100, 2.5*h.y.min, density = NULL, col = NULL, border = "black", lty = 1, lwd = 1)
  
  rect(0, 4.9*h.y.min, 100, 3.7*h.y.min, density = NULL, col = "white", border = NULL, lty = 1, lwd = 1)
  text(seq(from=20, to=90, by=10),rep(4.3*h.y.min, 8),labels=c(0.01,0.02,0.06,0.13,0.29,0.56,0.88,0.998),cex=0.6,font=1)
  text(-0.5,4.35*h.y.min,labels="Score Percentiles",cex=0.6,font=1, pos=4)
  
  rect(0, 6.1*h.y.min, 100, 4.9*h.y.min, density = NULL, col = NA, border = NA, lty = 1, lwd = 1)
  rect(min.SUS-1, 6*h.y.min, max.SUS+1, 5*h.y.min, density = NULL, col = "azure1", border = NA, lty = 1, lwd = 1)
  freqs<-count(raw_dat)
  text(freqs[,1],rep(5.5*h.y.min, length(freqs[,1])),labels=freqs[,2],cex=0.6,font=1)
  text(-0.5,5.6*h.y.min,labels="Frequencies",cex=0.6,font=1, pos=4)
  
  grid.text("Acceptability ranges (Bangor et al., 2008, p. 592)\nLetter grades (Sauro & Lewis, 2016, p. 204)\nAdjective ratings (Bangor et al., 2009, p. 118)\nScore percentiles (Sauro & Lewis, 2016, p. 203)",
            x = unit(0.015, "npc"), y = unit(0.042, "npc"),
            just="left", gp=gpar(fontsize=6, col="grey50"))
  
  #dev.off()
  
}

genScatter <- function(dat){
  
  scatter.option <- 2 # 1 uses dots for the scatter plots; 2 uses counts (as text)
  
  dat_sub <- dat[colnames(dat)[grep('X[0-9]', colnames(dat))]]
  colnames(dat_sub) <- paste0('Q',1:length(colnames(dat_sub)))
  
  odds <- seq(1, ncol(dat_sub), by = 2)
  evens <- seq(2, ncol(dat_sub), by = 2)
  
  dat_sub <- dat_sub[,c(odds, evens)]
  
  dat_sub <- rbind(dat_sub, rep(0,10), rep(6,10))
  real_data_points <- nrow(dat_sub) - 2

  cols = rev(brewer.pal(11, "RdBu"))   # reverses the pallete to go from blue (positive correlation) to white to red (negative correlation)cols = brewer.pal(11, "PuOr")
  pal = colorRampPalette(cols)
  cor_colors = data.frame(correlation = seq(-1,1,0.01), 
                          correlation_color = pal(201)[1:201])  # assigns a color for each r correlation value
  cor_colors$correlation_color = as.character(cor_colors$correlation_color)
  
  # function to report spearman correlations (with significance testing) in the upper triangular matrix
  panel.cor <- function(x, y, digits=2, cex.cor)
  {
    x<-x[1:real_data_points]
    y<-y[1:real_data_points]
    par(usr = c(0, 1, 0, 1))
    u <- par('usr') 
    names(u) <- c("xleft", "xright", "ybottom", "ytop")
    r <- cor(x, y, method="spearman", use="complete.obs")
    op <- options(warn = (-1))
    test <- cor.test(x,y, method="spearman")
    options(op)
    bgcolor = cor_colors[2+(-r+1)*100,2] # converts correlation into a specific color
    bgcolor = adjustcolor(bgcolor,alpha.f = 0.5)
    do.call(rect, c(col = bgcolor, as.list(u))) # colors the correlation box
    
    if (test$p.value> 0.05){
      text(0.5,0.5,"Insignificant",cex=1.1)
    } else{
      text(0.5, 0.75, paste("rho =",round(r,2)),cex=1.1) # prints correlatoin coefficient
      text(.5, .25, paste("p =",formatC(test$p.value, format = "e", digits = 1)),cex=1.1)  
      abline(h = 0.5, lty = 2) # draws a line between correlatoin coefficient and p value
    }
  }
  
  # function to create scatterplots (with lowess) in the lower triangular matrix
  panel.smooth<-function (x, y, col = "black", bg = NA, pch = 19, cex = 1.5, # 0.75 lower size
                          col.smooth = "blue", span = 2/3, iter = 3, ...) {
    x<-x[1:real_data_points]
    y<-y[1:real_data_points]
    par(usr = c(0, 6, 0, 6)) # test lines to set the x and y axes of the scatterplots from 1 to 5
    # using adjustcolor increases how dark the points are based on the number of occurences
    points(x, y, pch = pch, col = adjustcolor(col,alpha.f = 0.2), bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
      lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), lwd=2.5, col = col.smooth, ...)
  }
  
  # function to create scatterplots (with linear models) in the lower triangular matrix
  reg <- function(x, y, col, lwd) {
    clip(min(x),max(x),min(y),max(y))
    abline(lm(y~x), col=col, lwd=lwd)
  }
  panel.lm =  function (x, y, col = "black", bg = NA, pch = 19, 
                        cex = 1, col.lm = "blue", span = 2/3, iter = 3, ...)  {
    x<-x[1:real_data_points]
    y<-y[1:real_data_points]
    par(usr = c(0, 6, 0, 6)) 
    if (scatter.option == 1) {
      # Option 1: Use dots to represent the data, applying adjustcolor to increase how dark the points are based on the number of occurences
      points(x, y, pch = pch, col = adjustcolor(col,alpha.f = 0.2), bg = bg, cex = cex)  
    } else {
      # Option 2: Use text to show the count of occurrences
      freqs<-count(data.frame(x,y))
      text(freqs[,1],freqs[,2],labels=freqs[,3],cex=0.9,font=2)
    } # end if
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) reg(x[ok], y[ok], adjustcolor(col.lm,alpha.f = 0.35), 1.5)
  }
  
  # function to create histograms along the main diagonal
  panel.hist.user.defined <- function(x, ...)
  {
    x<-x[1:real_data_points]
    counts<-rep(0,5)
    for(i in 1:5) {
      counts[i]<-length(which(x == i))
    }
    par(usr = c(0,6,0,1.5))
    y <- counts; y <- counts/max(counts)
    rect(seq(from = 0.5, 4.5, by = 1), 0, seq(from = 1.5, 5.5, by = 1), y, 
         col=c("indianred1","pink","gray90","darkseagreen1","darkseagreen2"), ...)
  }
  
  
  
  
  # Cairo(width = 10, height = 8, file="TestScatter3.png", type="png", pointsize=12,
  #       bg = "white", canvas = "white", units = "in", dpi = 300)
  
  pairs(dat_sub, lower.panel=panel.lm, upper.panel=panel.cor, 
        diag.panel=panel.hist.user.defined, cex.labels=1.5, font.labels=2, gap=0.5, main = "Scatterplot Matrix of SUS Responses")
  
  grid.lines(c(0.5,0.5), c(0.06,.91), gp=gpar(col="gray",lwd=2))
  grid.lines(c(0.05,0.95), c(0.4841,0.4841), gp=gpar(col="gray",lwd=2))
  
  bump_right <- 0.025
  grid.text("Legends", x = unit(0.122+bump_right, "npc"), y = unit(0.02, "npc"),
            just="right", gp=gpar(fontsize=10, col="black", fontface="bold"))
  
  if (scatter.option == 1) {
    grid.text("1. Scatterplot Freq:  Low", x = unit(0.13+bump_right, "npc"), y = unit(0.02, "npc"),
              just="left", gp=gpar(fontsize=10, col="black", fontface="italic"))
    grid.circle(x=0.295+bump_right, y=0.02, r=0.004, gp=gpar(fill="gray90"))
    grid.circle(x=0.305+bump_right, y=0.02, r=0.004, gp=gpar(fill="gray70"))
    grid.circle(x=0.315+bump_right, y=0.02, r=0.004, gp=gpar(fill="gray50"))
    grid.circle(x=0.325+bump_right, y=0.02, r=0.004, gp=gpar(fill="gray30"))
    grid.circle(x=0.335+bump_right, y=0.02, r=0.004, gp=gpar(fill="gray10"))
    grid.text("High,", x = unit(0.345+bump_right, "npc"), y = unit(0.02, "npc"),
              just="left", gp=gpar(fontsize=10, col="black", fontface="italic"))
  } else {
    grid.text("1. Scatterplot Freq:  # equates to count,", x = unit(0.13+bump_right, "npc"), y = unit(0.02, "npc"),
              just="left", gp=gpar(fontsize=10, col="black", fontface="italic"))
  } # end if
  
  
  grid.text("2. Hist Response:  Disagree", x = unit(0.39+bump_right, "npc"), y = unit(0.02, "npc"),
            just="left", gp=gpar(fontsize=10, col="black", fontface="italic"))
  grid.rect(x=0.575+bump_right, y=0.02, width=0.008, height=0.008, gp=gpar(fill="indianred1"))
  grid.rect(x=0.586+bump_right, y=0.02, width=0.008, height=0.008, gp=gpar(fill="pink"))
  grid.rect(x=0.597+bump_right, y=0.02, width=0.008, height=0.008, gp=gpar(fill="gray90"))
  grid.rect(x=0.608+bump_right, y=0.02, width=0.008, height=0.008, gp=gpar(fill="darkseagreen1"))
  grid.rect(x=0.619+bump_right, y=0.02, width=0.008, height=0.008, gp=gpar(fill="darkseagreen2"))
  grid.text("Agree,", x = unit(0.63+bump_right, "npc"), y = unit(0.02, "npc"),
            just="left", gp=gpar(fontsize=10, col="black", fontface="italic"))
  
  grid.text("3. Correlation:  -1", x = unit(0.685+bump_right, "npc"), y = unit(0.02, "npc"),
            just="left", gp=gpar(fontsize=10, col="black", fontface="italic"))
  grid.rect(x=0.805+bump_right, y=0.02, width=0.008, height=0.008, gp=gpar(fill=adjustcolor(cor_colors[201,2],alpha.f = 0.5)))
  grid.rect(x=0.816+bump_right, y=0.02, width=0.008, height=0.008, gp=gpar(fill=adjustcolor(cor_colors[150,2],alpha.f = 0.5)))
  grid.rect(x=0.827+bump_right, y=0.02, width=0.008, height=0.008, gp=gpar(fill=adjustcolor(cor_colors[100,2],alpha.f = 0.5)))
  grid.rect(x=0.838+bump_right, y=0.02, width=0.008, height=0.008, gp=gpar(fill=adjustcolor(cor_colors[50,2],alpha.f = 0.5)))
  grid.rect(x=0.849+bump_right, y=0.02, width=0.008, height=0.008, gp=gpar(fill=adjustcolor(cor_colors[1,2],alpha.f = 0.5)))
  grid.text("+1", x = unit(0.86+bump_right, "npc"), y = unit(0.02, "npc"),
            just="left", gp=gpar(fontsize=10, col="black", fontface="italic"))
  
  grid.text("Agreement = Evidence of Usability", x = unit(0.014, "npc"), y = unit(0.69, "npc"),
            rot = 90, just="center", gp=gpar(fontsize=10, col="darkgreen"))
  grid.text("Disagreement = Evidence of Usability", x = unit(0.014, "npc"), y = unit(0.275, "npc"),
            rot = 90, just="center", gp=gpar(fontsize=10, col="indianred1"))
  
  grid.lines(x = unit(c(0.025, 0.025), "npc"), y = unit(c(0.49, 0.895), "npc"), 
             arrow = arrow(angle=90, length=unit(0.05, "inches"), ends="both", type="closed"), gp=gpar(col="darkgreen"))
  grid.lines(x = unit(c(0.025, 0.025), "npc"), y = unit(c(0.07, 0.48), "npc"), 
             arrow = arrow(angle=90, length=unit(0.05, "inches"), ends="both", type="closed"), gp=gpar(col="indianred1"))  
  
  grid.lines(x = unit(c(0.056, 0.497), "npc"), y = unit(c(0.933, 0.933), "npc"), 
             arrow = arrow(angle=90, length=unit(0.05, "inches"), ends="both", type="closed"), gp=gpar(col="darkgreen"))
  grid.lines(x = unit(c(0.503, 0.944), "npc"), y = unit(c(0.933, 0.933), "npc"), 
             arrow = arrow(angle=90, length=unit(0.05, "inches"), ends="both", type="closed"), gp=gpar(col="indianred1"))
  
  #dev.off()
}

error_checker <- function(dat){
  
  odds <- seq(1, ncol(dat), by = 2)
  evens <- seq(2, ncol(dat), by = 2)
  
  dat <- dat[, 2:ncol(dat)]
  
  find_sl <- function(d){
    (max(d) - min(d)) < 2
  }
  
  sl <- apply(dat, 1, find_sl)
  
  return(sl)
}




