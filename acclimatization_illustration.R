x <- seq(295, 298, length.out = 10) + rnorm(10, 0, 0.2) ##Tair
y <- seq(294, 296, length.out = 10) + rnorm(10, 0, 0.5) ##Topt

xmin <- min(c(x,y)); xmax <- max(c(x,y)) 
x1 <- seq(xmin-0.05, xmax+0.05, length.out=1000)

dxy1 = data.frame(x=x, y=y)
dxy = dxy1[order(dxy1$x),]
txy = data.frame(x=x1, y=x1)

dist <- NA
for(i in 1:nrow(dxy)){
  dist[i] <- min(apply(txy, 1, function(x) 
    (sqrt((x[1] - dxy[i,1])^2 + (x[2] - dxy[i,2])^2)) ))
}

main_title1 <- expression(T[air]~"vs"~ T[opt])
xlab1 <- expression(T[air]~"(°C)")
ylab1 <- expression(T[opt]~"(°C)")

par(mfrow=c(1,2), oma = c(0.5,1, 0.2, 0.2))
plot(dxy, xlim=c(xmin, xmax), ylim=c(xmin, xmax), pch = 19, cex=2, col= 'royalblue4',
     main = main_title1, xlab = xlab1, ylab = ylab1,
     cex.main=2, cex.lab=1, cex.axis=1)
abline(a=0, b=1, lty=2, lwd=3)

acc <- c(mean(dist), cor(dxy$x, dist, method='kendall'))
mod <- summary(lm(dist ~ dxy$x))
rsq <- mod$r.squared
main_title2 <- bquote("Kendall"~tau~"="~ .(rsq))
xlab2 <- expression(T[air]~"(°C)")
ylab2 <- expression(Delta~"(°C)")
plot(dxy$x, dist, ty='o', col= 'royalblue4', pch = 19, cex=2,
     main = main_title2, xlab = xlab2, ylab = ylab2,
     cex.main=2, cex.lab=1, cex.axis=1)
abline(a = mod$coefficients[1], b= mod$coefficients[2], lty=1, lwd=2)
abline(a = acc[1], b=0, lty=2, lwd=3, col='red')

ap1 <- ggplot(dxy, aes(x = x, y = y))+
  geom_point(size=4, colour= "blue")+
  geom_abline(linetype=5, colour='black', size = 2)+
  xlim(xmin, xmax) +  ylim(xmin, xmax) +
  xlab(xlab1) + ylab(ylab1) + 
  ggtitle(main_title1)+
  theme_bw() +
  theme(
    panel.grid = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text( size=12, colour='black',vjust=0.70),
    axis.text.y = element_text(size=12, colour='black'),
    axis.title.y = element_text(size=14, colour='black'),
    plot.title = element_text(size=16, colour = 'black', hjust = 0.5)
  )

p2df <- data.frame(x = dxy$x, y = dist)

ap2 <- ggplot(p2df, aes(x = x, y = y))+
  geom_point(size=4, colour= "brown1")+
  geom_hline(yintercept =  mean(p2df$y), linetype=4, colour='blue', size = 1)+
  geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2],
              size=1, colour = "brown1")+
  xlab(xlab2) + ylab(ylab2) + 
  ggtitle(main_title2)+
  theme_bw() +
  theme(
    panel.grid = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text( size=12, colour='black',vjust=0.70),
    axis.text.y = element_text(size=12, colour='black'),
    axis.title.y = element_text(size=14, colour='black'),
    plot.title = element_text(size=16, colour = 'black', hjust = 0.5)
  )
ap2

library(gridExtra)
grid.arrange(ap1, ap2, ncol=2)


##Acclimation illustration##
x <- seq(295, 298, length.out = 10) + rnorm(10, 0, 0.4) ##Tair
y <- seq(295, 298, length.out = 10) + rnorm(10, 0, 0.05) ##Topt

xmin <- min(c(x,y)); xmax <- max(c(x,y)) 
x1 <- seq(xmin-0.05, xmax+0.05, length.out=1000)

dxy1 = data.frame(x=x, y=y)
dxy = dxy1[order(dxy1$x),]
txy = data.frame(x=x1, y=x1)

dist <- NA
for(i in 1:nrow(dxy)){
  dist[i] <- min(apply(txy, 1, function(x) 
    (sqrt((x[1] - dxy[i,1])^2 + (x[2] - dxy[i,2])^2)) ))
}

main_title1 <- expression(T[air]~"vs"~ T[opt])
xlab1 <- expression(T[air]~"(°C)")
ylab1 <- expression(T[opt]~"(°C)")

par(mfrow=c(1,2), oma = c(0.5,1, 0.2, 0.2))
plot(dxy, xlim=c(xmin, xmax), ylim=c(xmin, xmax), pch = 19, cex=2, col= 'royalblue4',
     main = main_title1, xlab = xlab1, ylab = ylab1,
     cex.main=2, cex.lab=1, cex.axis=1)
abline(a=0, b=1, lty=2, lwd=3)

acc <- c(mean(dist), cor(dxy$x, dist, method='kendall'))
mod <- summary(lm(dist ~ dxy$x))
rsq <- acc[2]
main_title2 <- bquote("Kendall"~tau~"="~ .(rsq))
xlab2 <- expression(T[air]~"(°C)")
ylab2 <- expression(Delta~"(°C)")
plot(dxy$x, dist, ty='o', col= 'royalblue4', pch = 19, cex=2,
     main = main_title2, xlab = xlab2, ylab = ylab2,
     cex.main=2, cex.lab=1, cex.axis=1)
abline(a = mod$coefficients[1], b= mod$coefficients[2], lty=1, lwd=2)
abline(a = acc[1], b=0, lty=2, lwd=3, col='red')

ap1 <- ggplot(dxy, aes(x = x, y = y))+
  geom_point(size=4, colour= "blue")+
  geom_abline(linetype=5, colour='black', size = 2)+
  xlim(xmin, xmax) +  ylim(xmin, xmax) +
  xlab(xlab1) + ylab(ylab1) + 
  ggtitle(main_title1)+
  theme_bw() +
  theme(
    panel.grid = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text( size=12, colour='black',vjust=0.70),
    axis.text.y = element_text(size=12, colour='black'),
    axis.title.y = element_text(size=14, colour='black'),
    plot.title = element_text(size=16, colour = 'black', hjust = 0.5)
  )

p2df <- data.frame(x = dxy$x, y = dist)

ap2 <- ggplot(p2df, aes(x = x, y = y))+
  geom_point(size=4, colour= "brown1")+
  geom_hline(yintercept =  mean(p2df$y), linetype=4, colour='blue', size = 1)+
  geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2],
              size=1, colour = "brown1")+
  xlab(xlab2) + ylab(ylab2) + 
  ggtitle(main_title2)+
  theme_bw() +
  theme(
    panel.grid = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text( size=12, colour='black',vjust=0.70),
    axis.text.y = element_text(size=12, colour='black'),
    axis.title.y = element_text(size=14, colour='black'),
    plot.title = element_text(size=16, colour = 'black', hjust = 0.5)
  )
ap2

library(gridExtra)
grid.arrange(ap1, ap2, ncol=2)

