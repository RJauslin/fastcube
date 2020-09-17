
rm(list = ls())


library(Matrix)
library(ggplot2)
library(sampling)
library(viridis)

N <- 500

cat1 <- rep(1:50,each = N/50)
cat2 <- rep(1:100,each = N/100)
cat3 <- rep(1:125,each = N/125)
cat1 <- cat1[sample.int(N,N)]
cat2 <- cat2[sample.int(N,N)]
cat3 <- cat3[sample.int(N,N)]

Xcat <-as.matrix(data.frame(cat1 = cat1,
                            cat2 = cat2,
                            cat3 = cat3))


p <- 30
X <- matrix(rnorm(N*p),ncol = 30)


Xcat_tmp <- disjMatrix(Xcat)
Xred <- as.matrix(cbind(X,Xcat_tmp))

pik <- rep(100/N,N)
system.time(s <- ffphase(X,pik,Xcat))

system.time(s <- fcube(X,pik,Xcat))

system.time(s <- fcube(Xred,pik))
system.time(s <- BalancedSampling::cube(pik,as.matrix(Xred)))
system.time(s <- ReducedSamplecube(Xred,pik))
sum(s)


Xred <- as(Xred,"sparseMatrix")
dat_moran <- summary(Xred)
dat_moran$type <- rep("Tore",nrow(dat_moran))

theme_wave <- function(...) {
  theme_minimal() +
    theme(
      text = element_text(family="sans",color = "black",size = 9),
      panel.spacing = unit(2, "lines"),
      # title
      plot.title = element_text(hjust = 0.5,size = 9),
      # axes
      axis.line=element_blank(),
      axis.ticks=element_blank(),
      # legend
      legend.position="bottom",
      legend.title = element_text(size = 9,vjust = +1.0),
      legend.key.size = unit(0.3, "cm"),
      legend.key.width = unit(0.7,"cm") ,
      # background colors
      panel.background=element_blank(),
      panel.border=element_rect(colour = "black",fill = "transparent"),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      # keep edge black facet_wrap
      # strip.background = element_rect(fill="white"),
      strip.text =element_text(color = "black",size = 8)
    )
}

p_2 <-ggplot(data = dat_moran) +
  geom_tile(aes(x =j,y = i, fill = as.factor(x)))+
  scale_y_reverse()+
  scale_fill_grey(start = 0.01,end = 0.5)+
  # scale_fill_viridis_d()+
  coord_fixed()+
  labs(fill = "$a_{kl}$")+
  ggtitle("Stratification matrix")+
  theme_wave()+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none")
p_2

