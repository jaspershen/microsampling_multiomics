###############
#  Figure 5A #
###############
# Visits vs R^2
generate5A = function(clin,dataset = "30k",min_visits=10,cap=200){
  if (dataset == "iPOP"){
    identifier = "iPOP_ID"
    corDf.tmp = iPOPcorDf[!is.na(iPOPcorDf[[clin]]),]
  }
  else{
    identifier = "ANON_ID"
    corDf.tmp = corDf[!is.na(corDf[[clin]]),]
  }
  
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["Pulse"]]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["Temp"]]),]
  
  # Here we select people with the largest number of observations
  toppat = table(corDf.tmp[[identifier]])
  toppat = toppat[toppat > min_visits]
  toppat = names(toppat)
  toppat = toppat[1:min(cap,length(toppat))]
  
  dd = corDf.tmp[corDf.tmp[[identifier]] %in% toppat ,c(identifier,"Temp","Pulse",clin,"Clin_Result_Date")]
  
  # Compute R for individual models
  res = c()
  for (pat in toppat){
    frm = paste0(clin," ~ Pulse + Pulse^2 + Temp + Temp^2")
    
    # Crossval
    allpreds = c()
    for (i in 1:10){
      ind.data = dd[dd[[identifier]] == pat,]
      nsmpl = nrow(ind.data)
      train.idx = sample(nsmpl)[1:floor(0.8 * nsmpl)]
      model = lm(frm, data = ind.data[train.idx,])
      preds = predict(model,ind.data[-train.idx,])
      
      allpreds = rbind(allpreds, cbind(ind.data[-train.idx,clin], preds))
    }
    
    err = sqrt(1 - var(allpreds[,1] - allpreds[,2])/var(allpreds[,1] - mean(allpreds[,1])))
    vis = sum(dd[[identifier]] == pat)
    span = max(as.Date(ind.data$Clin_Result_Date)) - min(as.Date(ind.data$Clin_Result_Date))
    res = rbind(res, c(vis, span, err))
  }
  dres = data.frame(span = res[,2]/365, r = res[,3])
  ggplot(dres, aes(span, r)) + 
    weartals_theme + theme(text = element_text(size=25)) +
    geom_point(size=2) + 
    geom_smooth(method="lm", formula = y ~ x, size=1)
  ggsave(paste0("plots/Figure-5A-",clin,"-",dataset,".png"),width = 9,height = 6,units = "in")
}
generate5A("CHOL","30k",cap=100,min_visits = 10)

#####################################
#   Figure 5A (Timecourse from 2D)  #
#####################################

#Run after running individual time course from 2C (reading in various files for wear by timespans) that produced pct_var, corr_coeffs, & num_Records files   #
weartals_theme = theme_bw() + theme(text = element_text(size=18), panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
# read in each of the corr_coeffs from the different time windows
# had to manually add back in headers in the with demog files that I ran on scg: 20180327_corr_coeffs_AllData.csv , 20180327_corr_coeffs_MonthPrior.csv, 20180327_corr_coeffs_TwoWeekPrior.csv, 20180327_corr_coeffs_WeekPrior.csv,20180327_corr_coeffs_DayPrior.csv, 20180327_corr_coeffs_ThreeDayPrior.csv, 
# save as pdf 4x12.5"
# data <-read.table("/Users/jessilyn/Desktop/framework_timecourse/with_resting_bugfix_and_demographics/20180420_corr_coeffs_AllData_demog.csv",
#                   header=TRUE,sep=',',stringsAsFactors=FALSE)
fig.2c.df <-read.csv("/Users/jessilyn/Desktop/framework_timecourse/with_restingbugfix_demog_pctdev/20180507_AllData_pct_var.csv",
                     header=TRUE,sep=',',stringsAsFactors=FALSE)
fig.2c.plot <- melt(fig.2c.df)
fig.2c.plot[,3][is.nan(fig.2c.plot[,3])] <- 0 #replace % var explained of NaN w/ 0
fig.2c$test <- fig.2c.plot[order(-fig.2c.plot[,3]),] # reorder by LM Vitals
fig.2c$test = factor(fig.2c$test, levels = order(-fig.2c.plot[,3]))

# Plot the % var explained
ggplot(fig.2c, aes(x=test, y=value, color = variable)) + geom_point(size = 5, aes(shape=variable, color=variable)) +
  weartals_theme +
  ylim(0,1) +
  scale_shape_discrete(breaks=c("vitals", "lasso", "rf"),
                       labels=c("LM vitals", "LASSO", "RF")) +
  scale_color_discrete(breaks=c("vitals", "lasso", "rf"),
                       labels=c("LM vitals", "LASSO", "RF")) +
  labs(x = "Lab tests",y = expression(paste("Sqrt of % Variance Explained"))) 


### THE PART BELOW IS UNDER CONSTRUCTION
# Look at differences between DayPrior and AllData
data.allData <-read.table("/Users/jessilyn/Desktop/framework_timecourse/with_restingbugfix_demog_pctdev/20180503_pct_var_Dayprior.csv",
                          header=TRUE,sep=',',stringsAsFactors=FALSE)
data.dayPrior <-read.table("/Users/jessilyn/Desktop/framework_timecourse/with_restingbugfix_demog_pctdev/20180507_AllData_pct_var.csv",
                           header=TRUE,sep=',',stringsAsFactors=FALSE)
#data.dayPrior$model<- gsub("vitals", "vitals.ipop", data$model)
#data$r_squared <- pmax(data$r_squared, 0)
#data.nodemog$r_squared <- pmax(data$r_squared, 0)
#df <- merge(data.dayPrior, data.allData, by = c("model", "test"))
#colnames(df)[3:4] <- c("DayPrior", "AllData")
df <- data.dayPrior[,3:4] - data.allData[,3:4]
df2 <- melt(df)
df2<- df2[order(-df2$value),]
df2<-cbind(data.dayPrior[,1], df2)
colnames(df2)<-c("test", "model", "pct.var")
ggplot(df2, aes(test,pct.var, color = model)) + geom_point(size = 5, aes(shape=model, color=model)) +
  weartals_theme + 
  ylim(-0.5,0.5) +
  scale_shape_discrete(breaks=c("rf", "lasso"),
                       labels=c("RF", "LASSO")) +
  scale_color_discrete(breaks=c("rf", "lasso"),
                       labels=c("RF", "LASSO")) +
  labs(x = "Lab tests",y = expression(atop("Increase in Corr Coeff", paste("by using Day Prior vs All Watch Data"))))

### THE PART BELOW IS UNDER CONSTRUCTION
# make the case that if we could do the RF etc on all data (dont need to be individualized models) and we could combine the individualized models we could do an awesome job at preciting the clinical labs.
# can we create one more layer of mixed effects models in the iPOP analysis here?

# Look at differences between with and without demog
data.nodemog <-read.table("/Users/jessilyn/Desktop/framework_timecourse/with_resting_bugfix_no_demographics/20180420_corr_coeffs_TwoWeekPrior.csv",
                          header=TRUE,sep=',',stringsAsFactors=FALSE)
data <-read.table("/Users/jessilyn/Desktop/framework_timecourse/with_resting_bugfix_and_demographics/20180420_corr_coeffs_TwoWeekPrior_demog.csv",
                  header=TRUE,sep=',',stringsAsFactors=FALSE)
data$model<- gsub("vitals", "vitals.ipop", data$model)
#data$r_squared <- pmax(data$r_squared, 0)
#data.nodemog$r_squared <- pmax(data$r_squared, 0)
df <- merge(data, data.nodemog, by = c("model", "test"))
colnames(df)[3:4] <- c("r_w_demog", "r_no_demog")
df$delta <- df[,3] - df[,4]
df<- df[order(-df$delta),]

ggplot(df, aes(test,delta, color = model)) + geom_point(size = 5, aes(shape=model, color=model)) +
  weartals_theme + 
  ylim(-0.5,0.5) +
  scale_shape_discrete(breaks=c("all-rf", "lasso-rf", "all-lm", "lasso-lm", "vitals"),
                       labels=c("RF all variables", "RF + LASSO", "LM all variables", "LM + LASSO", "LM vitals")) +
  scale_color_discrete(breaks=c("all-rf", "lasso-rf", "all-lm", "lasso-lm", "vitals"),
                       labels=c("RF all variables", "RF + LASSO", "LM all variables", "LM + LASSO", "LM vitals")) +
  labs(x = "Lab tests",y = expression(atop("Increase in Corr Coeff", paste("by adjusting for demographics"))))

###############
#  Figure 5B #
###############
getEvents = function(dres, codes)
{
  dres$date = as.Date(dres$date)
  dres$change = c(0, dres$slope[2:length(dres$slope)] - dres$slope[2:length(dres$slope) - 1])
  dres$change = abs(dres$change)
  dres$time_diff = c(0, dres$date[2:length(dres$slope)] - dres$date[2:length(dres$slope) - 1])
  newpat = c(TRUE,dres$identifier[2:length(dres$slope)] != dres$identifier[2:length(dres$slope) - 1])
  dres = dres[!newpat,]
  dres = dres[order(-dres$change),]
  
  res.events = c()
  
  for(i in 1:10){
    events = codes[as.character(dres$identifier[i]) == codes$ANON_ID,]
    events$dist = abs(events$date - dres$date[i])
    events = events[order(events$dist),]
    res.events = rbind(res.events, events[1:10,])
  }
  res.events
}


# Temporal evolution of the estimate of the mean
# Method: Get last 50 (param: window) observations and compute adjasted rsquared of the lm
generate5C = function(clin,vit,dataset = "30k",window=50,filter = NULL,col = 1)
{
  if (dataset == "iPOP"){
    identifier = "iPOP_ID"
    corDf.tmp = iPOPcorDf[!is.na(iPOPcorDf[[clin]]),]
  }
  else{
    identifier = "ANON_ID"
    corDf.tmp = corDf[!is.na(corDf[[clin]]),]
  }
  
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[[vit]]),]
  if (!is.null(filter)){
    rows = corDf.tmp[[identifier]] %in% filter
    corDf.tmp = corDf.tmp[rows,]
  }
  
  # Here we select people with the largest number of observations
  toppat = table(corDf.tmp[[identifier]])
  toppat = names(sort(-toppat))[1:length(filter)]
  
  dd = corDf.tmp[corDf.tmp[[identifier]] %in% toppat, c(identifier,"Clin_Result_Date",vit,clin)]
  
  dates = c()
  slopes = c()
  rsquared = c()
  ids = c()
  vtrue = c()
  vpred = c()
  vitals = c("Pulse","Temp","Systolic","Diastolic","Respiration")
  
  for (pat in toppat){
    dd.pat = dd[dd[[identifier]] == pat,]
    for (i in (window+1):(nrow(dd.pat)-1) ){
      model = lm(formula(paste(clin,"~",paste(vitals,sep=" + "))), data = dd.pat[(i-window):i, ])
      slopes = c(slopes, model$coefficients[1])
      rsquared = c(rsquared, sqrt(summary(model)$r.squared))
      dates = c(dates, dd.pat$Clin_Result_Date[i])
      vtrue = c(vtrue, dd.pat[i+1,clin])
      vpred = c(vpred, predict(model, newdata = dd.pat[i+1,]) )
      ids = c(ids, pat)
    }
  }
  
  dres = data.frame(date = as.Date(as.POSIXct(dates)), slope = slopes, identifier = ids, rsquared = rsquared, true = vtrue,
                    predicted = vpred)
  dres$time = as.numeric(difftime(dres$date,min(dres$date),units="days")) / 365.0
  
  dres.pred = gather(dres, values, measurement, true:predicted, factor_key=TRUE)
  
  plt_pred = ggplot(dres.pred, aes(time, measurement, group = values, color = values)) + 
    weartals_theme + theme(text = element_text(size=25)) +
    geom_line(size=1.3) +
    geom_point(size=2) 
  print(plt_pred)
  
  plt_slope = ggplot(dres, aes(time, slope, group = identifier, color = identifier)) + 
    weartals_theme + theme(text = element_text(size=25)) +
    geom_line(size=1.3) +
    geom_point(size=2) 
  write.table(dres, file=paste0("data/Figure-5C-",clin,"-",vit,"-",pat,".csv"),sep = ',')
  
  cols = gg_color_hue(3)
  
  
  plt_rsq = ggplot(dres, aes(time, rsquared, group = identifier, color = identifier)) + 
    weartals_theme + theme(text = element_text(size=25)) +
    ylab(expression(sqrt("Variance explained")))
  
  if (dataset != "iPOP"){
    plt_rsq = plt_rsq + geom_point(size=4, colour=cols[col]) +
      geom_line(size=2, colour=cols[col]) 
  }
  else{
    plt_rsq = plt_rsq + geom_point(size=4) +
      geom_line(size=1.3) + ylim(c(0,0.6)) + theme(legend.position = "none")
  }
  
  events = NULL
  if (dataset != "iPOP")
    events = getEvents(dres, codes)
  
  if (dataset == "iPOP"){
    ggsave(paste0("plots/Figure-5C-",clin,'-',vit,"-",window,"-",dataset,"-slopes.png"), 
           plot = plt_slope, width = 6, height = 6,units = "in")
    ggsave(paste0("plots/Figure-5C-",clin,'-',vit,"-",window,"-",dataset,"-rsqured.png"),
           plot=plt_rsq,width = 6, height = 6,units = "in")
    ggsave(paste0("plots/Figure-5C-",clin,'-',vit,"-",window,"-",dataset,"-predictions.png"),
           plot=plt_pred,width = 6, height = 6,units = "in")
  }
  list(dres = dres, plt_slope = plt_slope, plt_rqs = plt_rsq, plt_pred = plt_pred, events = events)
}

## Load codes (initial)
codes = read.csv(paste0(dir,"initial_MI.csv"),stringsAsFactors=FALSE,header = FALSE)
colnames(codes) = c("ANON_ID","ICD_DATE","ICD9","ICD10","DX_NAME")
codes$date = as.Date(as.POSIXct(codes$ICD_DATE,format="%d-%b-%Y")) + 1 # POSIX no time mapped by Date to the previous day so add 1

# Find Patients with enough visits
getTop10Patients = function(){
  pats = unique(codes$ANON_ID)
  corDf.tmp = corDf[corDf$ANON_ID %in% pats,]
  counts = aggregate(corDf.tmp$ANON_ID, by=list(corDf.tmp$ANON_ID), FUN=length)
  counts = counts[order(-counts$x),]
  
  npats = 10
  pats = counts[counts[,2] > 50,1][1:10] # top 5 with at least 50 visits
  pats
}
#pats = getTop10Patients()

generate5Cevents = function(pats,col,dataset="30k"){
  dres = generate5C("HCT","Pulse",dataset,
                    window = 30,filter=pats,col)
  
  if(dataset == "30k"){
    codes_pats = codes[codes$ANON_ID %in% pats,]
    ids.tmp = c("!",codes_pats$ANON_ID)
    first = c(ids.tmp[-1] != ids.tmp[-length(ids.tmp)])
    last = c(first[-1],TRUE)
    
    # correct for D-148 (last event instead of first)
    if (pats == "D-148"){
      first[] = FALSE
      first[6] = TRUE
    }
    codes_pats$time = as.numeric(codes_pats$date - min(dres$dres$date))/365
    codes_pats = codes_pats[first,]
  }
  else{#Mike
    times = c(7,43,52,64)
    first = rep(FALSE, length(dres$dres$time))
    first[times] = TRUE
    codes_pats = list(time = as.numeric(dres$dres$time[times]), ICD10 = NULL)
    print(dres$dres$date[times])
  }
  filename = paste0("plots/Figure-5C-rsqured-",pats,".png")
  
  plt_cur = dres$plt_rqs + theme(legend.position="none")
  for (evid in 1:sum(first)){
    plt_cur = plt_cur + geom_vline(xintercept = codes_pats$time[evid],color="red",size=1) +
      geom_text(aes_q(x=codes_pats$time[evid], label=paste0("\n",codes_pats$ICD10[evid]),
                      y=max(dres$dres$rsquared) - sd(dres$dres$rsquared)/2),
                colour="black", angle=90, text=element_text())
  }
  
  ## Other events (TODO: very manual for now...)
  print(pats)
  if (pats == "D-148"){
    eid = 28
    plt_cur = plt_cur + geom_vline(xintercept = dres$dres$time[eid],color="blue",size=1)
    print(dres$dres$date[eid])
  }
  if (pats == "D-145"){
    for (eid in c(12, 17)){
      plt_cur = plt_cur + geom_vline(xintercept = dres$dres$time[eid],color="blue",size=1)
      print(dres$dres$date[eid])
    }
  }
  if (pats == "PD-6145"){
    for (eid in c(25)){
      plt_cur = plt_cur + geom_vline(xintercept = dres$dres$time[eid],color="blue",size=1)
      print(dres$dres$date[eid])
    }
  }
  
  print(plt_cur)
  ggsave(paste0(filename),
         plot=plt_cur,width = 6,height = 6,units = "in")
  filename = paste0("plots/Figure-5C-pred-",pats,".png")
  ggsave(paste0(filename),
         plot=dres$plt_pred,width = 6,height = 6,units = "in")
  
  dres
}

pats = c("D-145","D-148","PD-6145") #, "N-3683")
for (i in 1:length(pats))
  generate5Cevents(pats[i],i)

# Models depending on the number of observations
# Method:
# * Get only patients with at least nobs + neval observations
# * Start from nstart observations and go to nobs observations for training
generate5D = function(pat = "1636-69-001",lab.test = "HCT", nobs = 25, neval = 3, nstart = 6, shift = 2){
  ## 1636-69-001
  
  ## How much data we need for the estimates
  corDf.tmp = iPOPcorDf[!is.na(iPOPcorDf[[lab.test]]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["Pulse"]]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["Temp"]]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["systolic"]]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["diastolic"]]),]
  tbl = table(corDf.tmp$iPOP_ID)
  patids = names(tbl[tbl >= nobs + neval])
  wall = which(corDf.tmp$iPOP_ID == patids[1])
  half2 = wall[ceiling(length(wall)/2):length(wall)]
  corDf.tmp[half2,]$iPOP_ID = "1636-69-001-late"
  patids = c("1636-69-001-late", patids)
  patids = patids[1:2]
  
  # Here we select people with the largest number of observations
  vts = c("Pulse","Temp","systolic","diastolic") #,"Respiration")
  
  res = matrix(1, length(patids), nobs)
  for (i in nstart:nobs){
    for (j in 1:length(patids)){
      pat = patids[j]
      d = corDf.tmp[corDf.tmp$iPOP_ID %in% pat, c(lab.test,vts)]
      ntest = nrow(d) - shift
      
      train = shift + (ntest - neval + 1 - i):(ntest-neval)
      test = shift + (ntest - neval + 1):ntest
      
      model = lm(paste0(lab.test," ~ ."),data=d[train,])
      
      sm = summary(model)
      pp = predict(model, newdata=d[ntest,])
      err = mean((pp - d[(ntest - neval + 1):ntest,lab.test])**2) / mean((mean(d[train,lab.test]) - d[test,lab.test])**2)
      res[j,i] = err
    }
  }
  mse = colMeans(res[,1:nobs])
  mse[mse > 1] = 1
  df = data.frame(observations = 1:nobs, RPVE = sqrt(1 - mse))
  #plot(df,ylab="RPVE of HCT prediction",xlab="observation used for the model")
  
  plt = ggplot(df, aes(observations, RPVE)) + 
    geom_point(size=3) +
    weartals_theme + theme(text = element_text(size=25)) +
    theme(legend.position="none") +
    scale_x_continuous(breaks = pretty(df$observations, n = 13)) 
  print(plt)
  ggsave(paste0("plots/Figure-5D-",pat,"-",lab.test,".png"),width = 6, height = 6)
  write.table(df, file=paste0("data/Figure-5D-",pat,"-",lab.test,".csv"))
  df
}
res = generate5D("1636-69-001", "HCT")


# An alternative version of the one above, now with number of days monitored
generate5D_time = function(pat = "1636-69-001",lab.test = "HCT", nobs = 25, neval = 3, nstart = 6, shift = 2){
  ## 1636-69-001
  
  ## How much data we need for the estimates
  corDf.tmp = iPOPcorDf[!is.na(iPOPcorDf[[lab.test]]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["Pulse"]]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["Temp"]]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["systolic"]]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["diastolic"]]),]
  tbl = table(corDf.tmp$iPOP_ID)
  patids = names(tbl[tbl >= nobs + neval])
  wall = which(corDf.tmp$iPOP_ID == patids[1])
  #half2 = wall[ceiling(length(wall)/2):length(wall)]
  #corDf.tmp[half2,]$iPOP_ID = "1636-69-001-late"
  #patids = c("1636-69-001-late", patids)
  patids = patids[1]
  
  # Here we select people with the largest number of observations
  vts = c("Pulse","Temp","systolic","diastolic") #,"Respiration")
  
  res = c()
  d = corDf.tmp[corDf.tmp$iPOP_ID == patids, c(lab.test, vts, "Clin_Result_Date")]
  pdate = as.POSIXct(d[["Clin_Result_Date"]])
  d["time"] = (as.numeric(pdate))/(60*60*24)
  d["time"] = d["time"] - min(d["time"])

  
  windows = 2:20*50
  res = list()
  
  for (window in windows){
    res[[as.character(window)]] = c()
    for (i in 2:(nrow(d)-neval)){
      test = i:(i+neval)#nrow(d)
      
      train = 1:(i-1)
      enddate = d[["time"]][i]
      startdate = enddate - window
      
      train_data = d[train,]
      train_data = train_data[train_data[["time"]]>=startdate,]
      
      if (nrow(train_data) < length(vts) + 2)
        next
      if (nrow(train_data) == i-1)
        next
      
      model = lm(paste0(lab.test," ~ ."),data=train_data[,c(lab.test, vts)])
        
      sm = summary(model)
      pp = predict(model, newdata=d[test,c(lab.test, vts)])
      err = mean((pp - d[test,lab.test])**2) / mean((mean(d[train,lab.test]) - d[test,lab.test])**2)
      res[[as.character(window)]] = c(res[[as.character(window)]], err)
    }
  }
  mse = c()
  for (window in windows){
    mse = c(mse, mean(res[[as.character(window)]]))
  }
  mse[mse>1] = 1
  
  #mse = colMeans(res[,1:nobs])
  #mse[mse > 1] = 1
  df = data.frame(days_monitored = windows, RPVE = sqrt(1 - mse))
  #plot(df,ylab="RPVE of HCT prediction",xlab="observation used for the model")
  
  plt = ggplot(df, aes(days_monitored, RPVE)) + 
    geom_point(size=3) +
    weartals_theme + theme(text = element_text(size=25)) +
    theme(legend.position="none") +
    scale_x_continuous(breaks = pretty(df$days_monitored, n = 13)) 
  print(plt)
  ggsave(paste0("plots/Figure-5D-",pat,"-",lab.test,".png"),width = 6, height = 6)
  write.table(df, file=paste0("data/Figure-5D-",pat,"-",lab.test,".csv"))
  df
}
res = generate5D_time("1636-69-001", "HCT")

# corDf.tmp = corDf.tmp[corDf.tmp[["iPOP_ID"]] == "1636-69-001",]
# org = as.numeric(as.POSIXct(corDf.tmp[["Clin_Result_Date"]]),origin = "1970-01-01")
# dist = (org - lag(org))/(60*60*24)
# boxplot(dist)

###############
#   Figure 5 #
###############
# run after reading in and cleaning data and running Figure 2D section to get top.names

#weartals_theme = theme_bw() + theme(text = element_text(size=18), panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

#corr.coefs <-read.table("../SECURE_data/20180322_ranked_models_test_lm.csv",row.names=1, sep=",")
corr.coefs <-read.table("../SECURE_data/20180403_ranked_models_ipop_lm.csv",row.names=1, sep=",")
top.names<-rownames(corr.coefs) # names of lab tests from either the 30k or the iPOP simple bivariate models
top.names<-top.names[top.names %in% names(wear)] # only keep the lab names that are also present in the iPOP data

## Univariate Mixed-effect: True vs predicted 
# !! Only patients with at least min_visits = 20

min_visits = 20
mm.corr.coefs <- c()
lr.corr.coefs <- c()
id.corr.coefs <- c()
clin.idx <- c()
#for (i in 1:length(top.names)){
for (i in 1:4){
  clin = top.names[i]
  patients = sort(table(corDf[!is.na(corDf[[clin]]),]$ANON_ID))
  labs.vitals.tmp = corDf[corDf$ANON_ID %in% names(patients[patients > min_visits]),]
  labs.vitals.tmp$ANON_ID = factor(labs.vitals.tmp$ANON_ID)
  
  nn = nrow(labs.vitals.tmp)
  smp = sample(nn)
  test = smp[(1+floor(nn*0.9)):nn]
  train = smp[1:floor(nn*0.9)]
  
  frm = paste0(clin," ~ Pulse + Temp + (Pulse + Temp|ANON_ID)")
  print(frm)
  if (nrow(labs.vitals.tmp[train,]) && length(unique(labs.vitals.tmp[train,]$ANON_ID)) > 1){
    clin.idx <-c(clin.idx, clin)
    mm = lmer(frm, data = labs.vitals.tmp[train,])
    cf = coef(mm)
    vit = "Pulse"
    #qq = qplot(cf$ANON_ID[vit], geom="histogram")  + weartals_theme + xlab(paste0(top8[i]," ~ ",vit)) + ylab("count")
    #print(qq) 
    #, vp = viewport(layout.pos.row = matchidx$row,
    #                         layout.pos.col = matchidx$col))
    tt = labs.vitals.tmp[test,clin]
    
    # Evaluate LR model
    frm = paste0(clin," ~ Pulse + Temp")
    m0 = lm(frm, labs.vitals.tmp[train,])
    pp = predict(m0, newdata = labs.vitals.tmp[test,])
    #plot(pp, tt)
    lr.corr.coefs <- c(lr.corr.coefs, cor(pp,tt,use = "na.or.complete")) # corr coef of LR model
    
    # Evaluate MM model
    pp = predict(mm, newdata = labs.vitals.tmp[test,])
    #plot(pp, tt)
    mm.corr.coefs<-c(mm.corr.coefs, cor(pp,tt,use = "na.or.complete")) # corr coef of MM model
    
    # Evaluate LR model with ID
    frm = paste0(clin," ~ ANON_ID")
    m0 = lm(frm, labs.vitals.tmp[train,])
    pp = predict(m0, newdata = labs.vitals.tmp[test,])
    #plot(pp, tt)
    id.corr.coefs <- c(id.corr.coefs, cor(pp,tt,use = "na.or.complete")) # corr coef of LR model only with patient ID
  }
}
indiv.corr.coefs <- cbind(lr.corr.coefs, mm.corr.coefs, id.corr.coefs)
rownames(indiv.corr.coefs) <- clin.idx
write.table(indiv.corr.coefs, "../SECURE_data/20180329_indiv_30k_corr_coeffs.csv",row.names=TRUE,col.names=TRUE, sep=",")
data <-read.table("../SECURE_data/20180330/20180329_indiv_30k_corr_coeffs.csv",
                  header=TRUE,sep=',',stringsAsFactors=FALSE)
d <- melt(data, id.vars="X")
ggplot(d, aes(x=X, y=value, col=variable, shape=variable))+
  geom_point(cex=2.5) + 
  weartals_theme



# Build a RF HCT model for "1636-69-001", plot R^2 over time
generate5E = function(){
  
  wear.tmp = wear
  
  patients = unique(wear.tmp$iPOP_ID)
  #onehot = t(simplify2array(lapply(wear.tmp$iPOP_ID, FUN = function(x){x == patients})))
  #colnames(onehot) = paste("p",patients,sep="")
  
  wear.variables <- unlist(read.table(paste0(dir,"FinalLasso_153WearableFactors.csv"), stringsAsFactors = FALSE))[1:4] # the table of model features we want to work with
  #wear.tmp = cbind(wear.tmp, onehot)
  
  wear.variables = c("HCT","iPOP_ID","Clin_Result_Date",wear.variables)
  
  #change gender and ethnicity to dummy variables
  gender <- data.frame(model.matrix( ~ Gender - 1, data=wear.tmp))
  ethn <- data.frame(model.matrix( ~ Ethn - 1, data=wear.tmp))
  
  #remove the least populated gender and ethnicity (NCOL-1)
  cache <- names(gender)[which(sapply(gender,sum)==max(sapply(gender,sum)))]
  gender <- data.frame(cache=gender[which(sapply(gender,sum)==max(sapply(gender,sum)))])
  ethn <- ethn[,-which(sapply(ethn,sum)==min(sapply(ethn,sum)))]
  
  #store names as demo.variables
  demo.variables <- c("AgeIn2016", names(gender), names(ethn))
  
  #cbind new gender and ethnicity variables to "data"
  wear.tmp <- cbind(wear.tmp,gender,ethn)
  wear.tmp = wear.tmp[,c(wear.variables,demo.variables)]
  wear.tmp = na.omit(wear.tmp)
  
  test.patient = patients[1]
  test = wear.tmp$iPOP_ID == test.patient
  
  
  onehot = model.matrix(as.formula(paste0("~ iPOP_ID - 1")), wear.tmp)
  colnames(onehot) = paste0("ID",1:ncol(onehot))
  wear.tmp = cbind(wear.tmp, onehot)
  wear.tmp = wear.tmp[,-2]
  
  dates = c()
  rpve = c()
  ntest = 10
  for (i in 1:(sum(test) - ntest - 1) ){
    trainobs = which(test)[(0):(i)]
    testobs = which(test)[(i+1):(i+ntest)]
    
    train = rep(TRUE, nrow(wear.tmp))
    train[test] = FALSE
    train[trainobs] = TRUE
    
    wear.dates = wear.tmp[,c(2)]
    wear.data = wear.tmp[,-c(2)]
    
    model = randomForest(wear.data[train,-1],wear.data[train,1])
    
    tt = wear.data[testobs,1]
    pp = predict(model,newdata = wear.tmp[testobs,-1])
    rpve = c(rpve, sqrt(max(1 - mean((tt - pp)**2) / var(wear.data[,1]), 0)))
    dates = c(dates, wear.dates[testobs[ntest]])
  }
  dres.pred = data.frame(time = (as.Date(dates) - min(as.Date(dates)))/365.0, rpve = rpve, dates = as.Date(dates))
  plt_pred = ggplot(dres.pred, aes(time, rpve)) + 
    weartals_theme + theme(text = element_text(size=25)) +
    geom_line(size=1.3,color="red") +
    ylim(0, 1) +
    geom_point(size=3,color="red")
  ggsave("plots/Figure-5E.png",plot = plt_pred,width = 6, height = 6)
  list(plt_pred=plt_pred,  dres.pred=dres.pred)
}

# visits = aggregate(iPOPcorDf$iPOP_ID,by=list(iPOPcorDf$iPOP_ID),length)
# visits = visits[order(-visits$x),]
# pats = visits[1:1,1]

# FIGURE 6
# Starts at 2015-04-21 ends at 2017-06-14
generate6 = function(){
  dres = generate5C("HCT","Pulse","iPOP",
                    window = 10,filter="1636-69-001",1)
  dres$plt_rqs
  
  # Starts at 2014-12-07 ends at 2015-06-17
  figure45.varyingr2 = generate5E() 
  figure45.varyingr2$plt_pred
  
  # COMBINE
  dres.combined = dres$dres[,c("date","rsquared")]
  dres.combined$rpve = sqrt(dres.combined$rsquared)
  dres.combined = dres.combined[,c("date","rpve")]
  dres.combined$type="cVS"
  dres.combined.wear = figure45.varyingr2$dres.pred[,c("dates","rpve")]
  colnames(dres.combined.wear)[1] = "date"
  dres.combined.wear$type = "wVS"
  dres = rbind(dres.combined,dres.combined.wear)
  dres$time = (as.Date(dres$date) - min(as.Date(dres$date)))/365.0
  
  plt = ggplot(dres, aes(time, rpve, group=type, colour=type)) + 
    weartals_theme + theme(text = element_text(size=25)) +
    ylim(0, 1) +
    geom_line(size=1.3) +
    geom_point(size=3) +
    scale_color_manual(values=gg_color_hue(2)[c(2,1)]) +
    ylab("Multiple correlation\ncoefficient (R)") +
    xlab("Time (Years)")
  ggsave("plots/Figure-6.png",plot = plt,width = 8, height = 6)
  list(plt,dres)
}
generate6()
