library(tidyverse)
library(ggh4x)


ddm2w <- function(d_v, d_h, d_b, thres, nDT, tIn, bias, vd, hd, bd, vecOut) {
  T <- 5.2
  dt <- 0.008
  lt <- as.integer(T / dt)
  
  vec_tsecBF <- rep(1, lt)
  vec_tprimBF <- rep(1, lt)
  vec_tBonus <- rep(1, lt)
  aux <- abs(as.integer(tIn / dt))## turns the atribute entering second to 0. 
  
  if (tIn > 0) {
    vec_tsecBF[1:aux] <- 0
  } else if (tIn < 0) {
    vec_tprimBF[1:aux] <- 0
  }
  
  allX <- vector("list", length(vecOut))
  
  for (i in seq_along(vecOut)) {
    vecOut[i] <- T
    X <- numeric()
    X <- c(X, bias)
    flag <- 0
    cont <- 0
    
    while (flag == 0 && cont < lt) {
      noise <- rnorm(1, 0, sqrt(dt))
      newX <- X[length(X)] + (d_v * vd * vec_tprimBF[cont + 1] + d_h * hd * vec_tsecBF[cont + 1] + d_b * bd * vec_tBonus[cont + 1]) * dt + noise
      X <- c(X, newX)
      
      if (newX > 1) {
        flag <- 1
        vecOut[i] <- nDT + cont * dt
      } else if (newX < -1) {
        flag <- 1
        vecOut[i] <- -nDT - cont * dt
      }
      cont <- cont + 1
    }
    allX[[i]] <- X
  }
  
  return(list(vecOut = vecOut, X = allX,s_bf=aux))
}

# Example usage in R
d_v <- 0.5# primary bf drift
d_h <- 0.3# scondary bf drift
d_b <- 0.1# bonus drift
thres <- 1.0
nDT <- 0.1
tIn <- 0.8
bias <- 0
vd <- 1.0
hd <- -2.0
bd <- 0
vecOut <- numeric(10) # Example output vector
set.seed(123)

result <- ddm2w(d_v, d_h, d_b, thres, nDT, tIn, bias, vd, hd, bd, vecOut)

# Access the results
vecOut_result <- result$vecOut
X_result <- result$X
aux<-result$s_bf
# Print the results
#print(vecOut_result)
#print(X_result)


viz<-X_result%>%map_dfr(.,~{
  tibble(val=.x)%>%
    mutate(x=1:length(.x))%>%
    mutate(Bfs=x<aux)
},.id = "d")

x_breaks=c(0,50,100,150,200,250,300,350,400,450,500)
x_labels=x_breaks*0.008*1000


colors<-RColorBrewer::brewer.pal(8,"Dark2")[c(1,3,8)]

viz%>%mutate(val=case_when(
  val>=1 ~ 1,
  val<=-1~ -1,
  T~val
)
)%>%
  ggplot(aes(x=x,y=val,group=d,color=Bfs))+
  geom_line(size=1)+
  scale_color_manual(name="",values=c("#1B9E77","#7570B3"), breaks=(c(F,T)),  labels = c("\u03BB  second", "\u03BB  first"))+
  scale_x_continuous(name="Time [ms]",breaks=x_breaks,labels=x_labels)+
  scale_y_continuous(name="Decision variable [au]")+
  # annotate("text", x = 200, y = 1.2, label = "Threshold decide for ingroup",size=10)+
  #annotate("text", x = 200, y = -1.2, label = "Threshold decide for outgroup",size=10)+
  geom_hline(aes(yintercept = 1.1),size=8,color="grey90")+
  geom_hline(aes(yintercept = -1.1),size=8,color="grey90")+
  geom_vline(aes(xintercept = 100),size=2,color="darkgrey")+
  theme_minimal(18)+
  guides( y = "axis_truncated")+
  annotate("segment", x = 10, xend = 80, y = 0, yend = 0.8,
           arrow = arrow(type = "closed", length = unit(0.2, "inches")),
           color = "black",size=10)+
  annotate("segment", x = 10, xend = 80, y = 0, yend = 0.8,
           arrow = arrow(type = "closed", length = unit(0.2, "inches")),
           color = "#7570B3",size=3)+
  annotate("segment", x = 150, xend = 250, y = 0.5, yend = 0.4,
           arrow = arrow(type = "closed", length = unit(0.2, "inches")),
           color = "black",size=10)+
  annotate("segment", x = 150, xend = 250, y = 0.5, yend = 0.4,
           arrow = arrow(type = "closed", length = unit(0.2, "inches")),
           color = "#1B9E77",size=3)+
  geom_segment(aes(x = -30, y = 0, xend = 0, yend = 0), color = "purple", size = 1) +
  annotate("text", x = 107, y = 1.1, label = "\u03B8", color = "darkgrey", size = 10) +
  coord_cartesian(xlim=c(0,310))+
  theme(
    legend.position = c(0.89,0.5),
    legend.text = element_text(size=19),
    #legend.title = element_blank(),
    #legend.background = element_rect(colour = "black")
  )->ddm_fig
ddm_fig


ggsave(ddm_fig,filename=here::here("figures","ddm.png"),height = 7,width=7)
#   theme(
#     axis.line.x = element_line(color = "black"),
#     axis.line.y = element_line(color = "black"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks.length = unit(0.25, "cm"),
#     axis.ticks.x = element_line(color = "black"),
#     axis.ticks.y = element_line(color = "black"),
#     axis.ticks.margin = unit(0.5, "cm")
#   ) 




