
rm(list = ls())


library(Matrix)
library(ggplot2)
library(sampling)
library(viridis)

N <- 1000
pik <- rep(0.5,N)

Xcat <- as.matrix(rep(1:500,each = N/500))
Xcat_tmp <- disjMatrix(Xcat)
p <- 5
X <- cbind(pik,matrix(rnorm(N*p),ncol = p))
# X <- as.matrix(pik)
A <- cbind(X,Xcat_tmp)/pik

Xred <- as.matrix(cbind(X,Xcat_tmp))

system.time(s <- balcat(X,Xcat,pik))
round(s,9)
# s <- landingLP(X,s,pik,Xcat)

t(A)%*%s

t(A)%*%pik
# system.time(s <- fcube(X,pik,Xcat))
#
# system.time(s <- fcube(Xred,pik))
# system.time(s <- BalancedSampling::cube(pik,as.matrix(Xred)))
# system.time(s <- ReducedSamplecube(Xred,pik))
# sum(s)


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

