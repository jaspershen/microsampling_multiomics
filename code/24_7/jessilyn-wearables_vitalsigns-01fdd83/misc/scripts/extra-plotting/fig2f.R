library("dplyr")
library("tidyr")
library("ggplot2")
source("ggplot-theme.R")

load("res.cca.Rda")
df.res = data.frame(res.cca) %>%
  gather(group, value, Electrolytes:Hematologic) %>%
  group_by(group) %>%
  summarise(mean=mean(value), sd=sd(value), pval=pnorm(0, mean(value), sd(value))) %>%
  arrange(desc(mean))
df.res = data.frame(df.res)
#df.res[order(levels(df.res$group)),]$mean = res.cca.means[order(res.cca.means$test),]$value
df.res = df.res[-5,]
df.res$group = c("Hepatic","Hematologic","Metabolic","Cardiovascular","Immune","Electrolytes")
df.res$group = factor(df.res$group, levels = df.res$group)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(3)

p=ggplot(df.res, aes(x=group, y=mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, color=cols[1]) +
  theme(legend.title = element_blank()) +
  geom_point(size=7, color=cols[1]) +
  weartals_theme + 
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold")) + 
  ylim(0,0.5) +
  labs(x = NULL, y = NULL) #"Correlation Coefficient")
p
geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.5)
ggsave(paste0("fig2f.png"),p,width=5.5,height=4)
