
# UpSetR ------------------------------------------------------------------

library(UpSetR)
library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)


# Upregulated Data preparing ----------------------------------------------------------

# u2h
# change logFC column name to logFC.2h and same for others

u2h<-read.csv("Data/Upregulated.in.2h.pathogen.less.3.zero.csv")
u4h<<-read.csv("Data/Upregulated.in.4h.pathogen.less.3.zero.csv")
u8h<<-read.csv("Data/Upregulated.in.8h.pathogen.less.3.zero.csv")
u16h<<-read.csv("Data/Upregulated.in.16h.pathogen.less.3.zero.csv")
u24h<<-read.csv("Data/Upregulated.in.24h.pathogen.less.3.zero.csv")

u2h$logFC.2h[u2h$logFC.2h>0]<-1
u4h$logFC.4h[u4h$logFC.4h>0]<-1
u8h$logFC.8h[u8h$logFC.8h>0]<-1
u16h$logFC.16h[u16h$logFC.16h>0]<-1
u24h$logFC.24h[u24h$logFC.24h>0]<-1
u<-rbind.fill(u2h,u4h,u8h,u16h,u24h)
u[is.na(u[1:6])]<-0

U <- u %>% pivot_longer(logFC.2h:logFC.24h # (2:6)
                         , names_to = "timepoint", values_to = "sign") %>%
  dplyr::filter(sign == 1) %>%
  pivot_wider(names_from = timepoint, values_from = sign, values_fill = 0) %>% t()


write.csv(U,file = "Upregulated-ED-transposes.pathogen.csv",row.names = TRUE)


####################################################################################################### Read Upregulated data

df <-t(read.csv("Upregulated-ED-transposes.pathogen.csv")) %>% as.data.frame()
colnames(df)<- df[1,]
df <- df[-1,] %>% mutate_if(is.character, as.numeric)
upset(df)


####################################################################################################### Setting colors

main_bar_col <- c("violetred4")
sets_bar_col <- c("turquoise4")
matrix_col <- c("slateblue4")
shade_col <- c("wheat4")



####################################################################################################### Setting Set Variables

mb_ratio1 <- c(0.55,0.45)


####################################################################################################### Upregulated Plot

upset(df,
      mb.ratio = mb_ratio1,
      mainbar.y.label = "Interaction of upregulated genes",
      sets.x.label = "Number of genes",
      order.by = "freq",
      # show.numbers = TRUE,
      point.size = 2,
      line.size = 1,
      main.bar.color = main_bar_col,
      sets.bar.color = sets_bar_col,
      matrix.color = matrix_col,
      shade.color = shade_col )



# Downregulated Data preparing ----------------------------------------------------------

d2h<-read.csv("Data/Downregulated.in.2h.pathogen.less.3.zero.csv")
d4h<<-read.csv("Data/Downregulated.in.4h.pathogen.less.3.zero.csv")
d8h<<-read.csv("Data/Downregulated.in.8h.pathogen.less.3.zero.csv")
d16h<<-read.csv("Data/Downregulated.in.16h.pathogen.less.3.zero.csv")
d24h<<-read.csv("Data/Downregulated.in.24h.pathogen.less.3.zero.csv")



d2h$logFC.2h[d2h$logFC.2h<0]<-1
d4h$logFC.4h[d4h$logFC.4h<0]<-1
d8h$logFC.8h[d8h$logFC.8h<0]<-1
d16h$logFC.16h[d16h$logFC.16h<0]<-1
d24h$logFC.24h[d24h$logFC.24h<0]<-1
d<-rbind.fill(d2h,d4h,d8h,d16h,d24h)
d[is.na(d[1:6])]<-0


D <- d %>% pivot_longer(logFC.2h:logFC.24h, names_to = "timepoint", values_to = "sign") %>%
  dplyr::filter(sign == 1) %>%
  pivot_wider(names_from = timepoint, values_from = sign, values_fill = 0) %>% t()


write.csv(D,file = "Downregulated-ED-transposes.pathogen.csv",row.names = TRUE)


####################################################################################################### Read Downregulated data

d <-t(read.csv("Downregulated-ED-transposes.pathogen.csv")) %>% as.data.frame()
colnames(d)<- d[1,]
d <- d[-1,] %>% mutate_if(is.character, as.numeric)
upset(d)


####################################################################################################### Setting colors

 main_bar_col <- c("blue4")
 sets_bar_col <- c("coral1")
 matrix_col <- c("forestgreen")
 shade_col <- c("wheat4")



####################################################################################################### Setting Set Variables

mb_ratio1 <- c(0.55,0.45)

####################################################################################################### Downregulated Plot

upset(d,
      mb.ratio = mb_ratio1,
      mainbar.y.label = "Interaction of downregulated genes",
      sets.x.label = "Number of genes",
      order.by = "freq",
      # show.numbers = TRUE,
      point.size = 2,
      line.size = 1,
      main.bar.color = main_bar_col,
      sets.bar.color = sets_bar_col,
      matrix.color = matrix_col,
      shade.color = shade_col )

