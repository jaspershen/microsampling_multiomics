library("dplyr")
library("tidyr")
library("ggplot2")
source("ggplot-theme.R") # just to make things look nice

load("figure4b.Rda")

# Some preprocessing
res = res.iPOP[,c(1,7,9)]
colnames(res)[2] = "personal mean (iPOP)"

toplot.iPOP = gather(res, model, value, -test) %>%
  group_by(test, model) %>%
  summarise(mean = mean(value), sd = sd(value), pval = mean(value<0))

toplot.iPOP = toplot.iPOP[order(-toplot.iPOP$mean),]
toplot.iPOP = toplot.iPOP[order(toplot.iPOP$model,decreasing = TRUE),] # order by the personal mean
toplot.iPOP$test = factor(as.character(toplot.iPOP$test), levels = unique(as.character(toplot.iPOP$test)))

toplot.30k = group_by(res.30k,test,model) %>%
  summarise(mean = mean(value), sd = sd(value), pval = mean(value<0) )
toplot = rbind(toplot.30k,toplot.iPOP)

# Order by 
toplot = toplot[order(-toplot$mean),]
toplot = toplot[order(toplot$model),] # order by the personal mean
toplot$test = factor(as.character(toplot$test), levels = unique(as.character(toplot$test)))

# Plot
pp = ggplot(toplot, aes(test, mean, group = model, color = model)) +
  ylab(expression(sqrt("Variance explained"))) +
  xlab("Lab test") +
  geom_point(size = 3, position=position_dodge(width=0.7)) +
  geom_errorbar(size = 0.8, aes(ymin=mean-sd, ymax=mean+sd), width=.8, position=position_dodge(width=0.7)) +
  weartals_theme + 
  theme(text = element_text(size=14))

print(pp)

# Compare vitals and vitals + personal
t.test(toplot[toplot$model=="vitals" & toplot$test == "MCHC",]$mean,toplot[toplot$model=="vitals + personal mean and slope" & toplot$test == "MCHC",]$mean)
