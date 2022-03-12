library("dplyr")
library("ggplot2")
library("tidyr")

load("figure5c.Rda")

weartals_theme = theme_classic() + theme(text = element_text(size=18), panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

lambda.choice <- "lasso.min"

fig.2c.plot <- gather(res$fig.2c.corr.coefs, variable, value, -test)
fig.2c.plot = fig.2c.plot[fig.2c.plot$variable!="experiment",]
fig.2c.plot = fig.2c.plot %>% group_by(test,variable) %>% summarise(mean = mean(value), sd = sd(value))

tmp = fig.2c.plot[fig.2c.plot$variable == "personal.mean",]
lvls = as.character(tmp$test[order(-tmp$mean)])

fig.2c.plot$test = factor(fig.2c.plot$test, levels = lvls)
#^ Ran out of time, but I can simplify this later, which will probably rid the error.
fig.2c <- fig.2c.plot

ggplot(fig.2c[fig.2c$variable %in% c("rf.pers","personal.mean"),], aes(x=test, y=mean, color = variable)) +
  geom_errorbar(size = 0.8, aes(ymin=mean-sd, ymax=mean+sd), width=.8, position=position_dodge(width=0.7)) +
  geom_point(size = 3, position=position_dodge(width=0.7)) + #, aes(shape=variable)
  weartals_theme +
  ylim(0,1) +
  scale_color_manual(values=gg_color_hue(5)[c(3,1,2,5,4)]) +
  labs(x = "Lab tests",y = expression(paste("Sqrt of % Variance Explained")))
