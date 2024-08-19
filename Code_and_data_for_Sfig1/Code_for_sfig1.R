library(ggplot2)

# Panel B GAM models

## Left : mus models
df <- read.table("data_for_timemodels.txt")

df_mus <- df[df$species=="mus",]
pred.gam<-read.table("data_pred.gam_mus_all.txt")

weightstage_mus=(ggplot(pred.gam) + 
       geom_line(aes(x = weight, y = predicted)) +          # slope
       geom_ribbon(aes(x = weight, ymin = predicted - std.error, ymax = predicted + std.error), 
                   fill = "lightgrey", alpha = 0.5) +  # error band
       geom_point(data = df_mus,                      # adding the predicted data (scaled values)
                  aes(x = weight, y = stage),color="black",alpha=0.1) + 
       labs(x = "Weight", y = "stage", 
            title = "Weight and stage model MOUSE") + 
       theme_minimal()
)
ggsave(file=paste0("FigS1_PanelB_left.pdf"),weightstage_mus,width = 17,
       height = 15,units="cm")



## Right : ham models
df_ham <- df[df$species=="ham",]
pred.gam<-read.table("data_pred.gam_ham_all.txt")


weightstage_ham=(ggplot(pred.gam) + 
                   geom_line(aes(x = weight, y = predicted)) +          # slope
                   geom_ribbon(aes(x = weight, ymin = predicted - std.error, ymax = predicted + std.error), 
                               fill = "lightgrey", alpha = 0.5) +  # error band
                   geom_point(data = df_ham,                      # adding the predicted data (scaled values)
                              aes(x = weight, y = stage),color="black",alpha=0.1) + 
                   labs(x = "Weight", y = "stage", 
                        title = "Weight and stage model HAMSTER") + 
                   theme_minimal()
)
ggsave(file=paste0("FigS1_PanelB_right.pdf"),weightstage_ham,width = 17,
       height = 15,units="cm")


# Panel C Numerical ages

## Left : mus models

pred.gam.cusp<-read.table("data_cusp_mus_all.txt")
pred.gam.rnaseq<- read.table("data_rnaseq_mus_all.txt")
  
weighttime_mus=(ggplot(pred.gam.cusp) + 
         geom_line(aes(x = weight, y = fit)) +         
         geom_ribbon(aes(x = weight, ymin = lower1, ymax = upper1), 
                     fill = "lightgrey", alpha = 0.5) +  # error band
         
         geom_point(data = pred.gam.cusp,                      # adding the predicted data (scaled values)
                    aes(x = weight, y = fit),color="dark grey",alpha=0.8) + 
         
         geom_point(data = pred.gam.rnaseq,                      # adding the predicted data (scaled values)
                    aes(x = weight, y = fit),color="dark red",alpha=0.8) + 
         geom_errorbar(data = pred.gam.rnaseq,                      # adding the predicted data (scaled values)
                       aes(x = weight, ymin = lower1, ymax = upper1), colour="black", width=.1) +
         geom_point(data = pred.gam.rnaseq,                      # adding the dpc 
                    aes(x = weight, y = stage),color="dark orange",alpha=0.8)+
         labs(x = "Weight", y = "stage", 
              title = "Weight and stage model MOUSE") + 
         theme_minimal()
  )
 
ggsave(file=paste0("FigS1_PanelC_left.pdf"),weighttime_mus,width = 17,
       height = 15,units="cm")


## Right : hamster models

pred.gam.cusp<-read.table("data_cusp_ham_all.txt")
pred.gam.rnaseq<- read.table("data_rnaseq_ham_all.txt")

weighttime_mus=(ggplot(pred.gam.cusp) + 
                  geom_line(aes(x = weight, y = fit)) +         
                  geom_ribbon(aes(x = weight, ymin = lower1, ymax = upper1), 
                              fill = "lightgrey", alpha = 0.5) +  # error band
                  
                  geom_point(data = pred.gam.cusp,                      # adding the predicted data (scaled values)
                             aes(x = weight, y = fit),color="dark grey",alpha=0.8) + 
                  
                  geom_point(data = pred.gam.rnaseq,                      # adding the predicted data (scaled values)
                             aes(x = weight, y = fit),color="dark red",alpha=0.8) + 
                  geom_errorbar(data = pred.gam.rnaseq,                      # adding the predicted data (scaled values)
                                aes(x = weight, ymin = lower1, ymax = upper1), colour="black", width=.1) +
                  geom_point(data = pred.gam.rnaseq,                      # adding the dpc 
                             aes(x = weight, y = stage),color="dark orange",alpha=0.8)+
                  labs(x = "Weight", y = "stage", 
                       title = "Weight and stage model HAMSTER") + 
                  theme_minimal()
)

ggsave(file=paste0("FigS1_PanelC_right.pdf"),weighttime_mus,width = 17,
       height = 15,units="cm")

