#%%%%%%%%%%%%%%%%%%%%
# F2 Earthworms
# Data inspection
# Exploratory Data Analysis
#%%%%%%%%%%%%%%%%%%%%


# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/GatherSource/F2_EW_MakeLikeFile.R")
source("Data/GatherSource/Functions/panelutils.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Data processing ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Subsetting 

## Explanatory variables are the soil properties of the plots
df.exp <- data[,c("clay","sand","pH","mc","cn")]

df.groups <- data[,c("age_class","samcam","field.ID", "hole")]

df.climate <- data[,c("ats1", "ata1", "atb1","hum1", "rad1", "prec1")]

## Response variables for Abundances
dt.rsp.abn <- as.data.table(data[,c("anc", "ancad", "endo", "endad", "N")])

#Juveniles
dt.rsp.abn[, anc.juv:=anc-ancad]
dt.rsp.abn[, endo.juv:=endo-endad]
dt.rsp.abn[,juv:= anc.juv + endo.juv]


## Response variables for biomass
dt.rsp.bms <- as.data.table(data[,c("anc.bm", "ancad.bm", "endo.bm", "endad.bm", "N.bm")])

# Juveniles
dt.rsp.bms[, anc.juv.bm:=anc.bm-ancad.bm]
dt.rsp.bms[, endo.juv.bm:=endo.bm-endad.bm]
dt.rsp.bms[,juv.bm:= anc.juv.bm + endo.juv.bm]


## Response variables for biodiversity
dt.rsp.bdv <- as.data.table(data[,c("SR", "H", "J")])

# Put together
dt.rsp1 <- cbind(dt.rsp.abn, dt.rsp.bms, dt.rsp.bdv)

data = cbind(data, dt.rsp1[,c("anc.juv", "endo.juv", "juv", "anc.juv.bm", "endo.juv.bm", "juv.bm"), with=F])

## Variable selection lists 
time <- c("age_class", "age", "samcam")
climate <- c("ats1", "rad1", "hum1", "prec1")
soil <- c("mc", "pH", "cn", "clay", "sand", "silt")

abn.sfg <- c("anc", "endo", "epi")
abn.lifecycle <- c("ancad", "endad", "juv", "N")
bms.sfg <- c("anc.bm", "endo.bm", "epi.bm")
bms.lifecycle <- c("ancad.bm", "endad.bm", "juv.bm", "N.bm")

bdv <- c("SR", "H")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Inspect climate data in th course of the investigation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Climate variiables vs date
df.m.date <- melt(cbind(date=data$date, location=data$location, data[,climate]), id.vars=c("date", "location"))
ggplot(df.m.date, aes(x=date, y=value, group=location)) + 
  geom_point() + geom_line(aes(col=location)) + facet_wrap(~variable, ncol=2, scales="free")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Plot response vs explanatory variables ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Plot defaults 
default.ggplot <- list(geom_point(),
                       geom_smooth(),
                       geom_smooth(method=lm, fullrange=T, col="yellow"),
                       geom_smooth(method=lm, fullrange=T, formula = y ~ x + I(x^2), col="red"),
                       facet_wrap(~variable, ncol=2, scales="free"))

## Abundance ####
### Abundance Soil functional groups vs climate ####
df.melt <- melt(cbind(data[,abn.sfg], data[,climate]), id.vars=c("anc", "endo", "epi"))
ggplot(df.melt, aes(x=value, y=anc, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=endo, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=epi, group=variable)) + default.ggplot

df.melt <- melt(cbind(data[,abn.lifecycle], data[,climate]), id.vars=c("ancad", "endad", "juv", "N"))
ggplot(df.melt, aes(x=value, y=ancad, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=endad, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=juv, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=N, group=variable)) + default.ggplot

### Abundance vs soil variables ####
df.melt <- melt(cbind(data[,abn.sfg], data[,soil]), id.vars=c("anc", "endo", "epi"))
ggplot(df.melt, aes(x=value, y=anc, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=endo, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=epi, group=variable)) + default.ggplot

df.melt <- melt(cbind(data[,abn.lifecycle], data[,soil]), id.vars=c("ancad", "endad", "juv", "N"))
ggplot(df.melt, aes(x=value, y=ancad, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=endad, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=juv, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=N, group=variable)) + default.ggplot

### Abundance vs time ####
df.melt <- melt(cbind(data[,abn.sfg], data[,time], date=data$date), id.vars=c("age_class", "date", "age", "samcam"))
ggplot(df.melt, aes(x=age_class, y=value, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=date, y=value, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=age, y=value, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=samcam, y=value, group=variable)) + default.ggplot

df.melt <- melt(cbind(data[,abn.lifecycle], data[,time], date=data$date), id.vars=c("age_class", "date", "age", "samcam"))
ggplot(df.melt, aes(x=age_class, y=value, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=date, y=value, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=age, y=value, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=samcam, y=value, group=variable)) + default.ggplot


## Biomass ####
### Biomass vs. climate ####
df.melt <- melt(cbind(data[,bms.sfg], data[,climate]), id.vars=c("anc.bm", "endo.bm", "epi.bm"))
ggplot(df.melt, aes(x=value, y=anc.bm, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=endo.bm, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=epi.bm, group=variable)) + default.ggplot

df.melt <- melt(cbind(data[,bms.lifecycle], data[,climate]), id.vars=c("ancad.bm", "endad.bm", "juv.bm", "N.bm"))
ggplot(df.melt, aes(x=value, y=ancad.bm, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=endad.bm, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=juv.bm, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=N.bm, group=variable)) + default.ggplot

### Biomass vs. soil ####
df.melt <- melt(cbind(data[,bms.sfg], data[,soil]), id.vars=c("anc.bm", "endo.bm", "epi.bm"))
ggplot(df.melt, aes(x=value, y=anc.bm, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=endo.bm, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=epi.bm, group=variable)) + default.ggplot

df.melt <- melt(cbind(data[,bms.lifecycle], data[,soil]), id.vars=c("ancad.bm", "endad.bm", "juv.bm", "N.bm"))
ggplot(df.melt, aes(x=value, y=ancad.bm, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=endad.bm, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=juv.bm, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=N.bm, group=variable)) + default.ggplot

### Biomass vs. time ####
df.melt <- melt(cbind(data[,bms.sfg], data[,time], date=data$date), id.vars=c("age_class", "date", "age", "samcam"))
ggplot(df.melt, aes(x=age_class, y=value, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=date, y=value, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=age, y=value, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=samcam, y=value, group=variable)) + default.ggplot

df.melt <- melt(cbind(data[,bms.lifecycle], data[,time], date=data$date), id.vars=c("age_class", "date", "age", "samcam"))
ggplot(df.melt, aes(x=age_class, y=value, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=date, y=value, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=age, y=value, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=samcam, y=value, group=variable)) + default.ggplot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




















str(data)

matplot( data[,climate], data$ancad)

x11()
Mydotplot(data[,soil])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


View(data2)

### Biodiversity vs.t climate ####
df.melt <- melt(cbind(data[,bdv], data[,climate]), id.vars=c("SR", "H"))
ggplot(df.melt, aes(x=value, y=SR, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=H, group=variable)) + default.ggplot

### biodiversity vs. soil ####
df.melt <- melt(cbind(data[,bdv], data[,soil]), id.vars=c("SR", "H"))
ggplot(df.melt, aes(x=value, y=SR, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=value, y=H, group=variable)) + default.ggplot

### biodiversity vs. time ####
df.melt <- melt(cbind(data[,bdv], data[,time], date=data$date), id.vars=c("age_class", "date", "age", "samcam"))
ggplot(df.melt, aes(x=age_class, y=value, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=date, y=value, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=age, y=value, group=variable)) + default.ggplot
ggplot(df.melt, aes(x=samcam, y=value, group=variable)) + default.ggplot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Data summary ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p <- ncol(dt.rsp1)

df.rsp1 <- as.data.frame(dt.rsp1)

for(i in 1:(p-1)) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  y <- df.rsp1[,i]
  plot(y)
  boxplot(y)
  hist(y, main="")
  plot(y^2)
  title(names(df.rsp1)[i],outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Multicollinearity between explanatory variables #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pairs(data[,c(time, climate, soil)], 
      panel=panel.smooth, 
      upper.panel=panel.cor, method="pearson",
      diag.panel=panel.hist, 
      main = "Bivariate Plots with histogram and smooth curves")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Species occurrence, spatial abundance ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spatial <- c("dgO", "dgN", "ALO", "ARO", "ACA", "ACH", "OCY", "OLA", "LTR", "LCA", "LRU")

df.melt <- melt(data[,spatial], id.vars=c("dgO", "dgN"))

library(mapdata)
germany <- fortify(map("worldHires", "Germany",fill = TRUE, plot = FALSE))

ggplot(df.melt, aes(dgO, dgN))  + 
  geom_point(aes(size=value*100, col="red")) +
  facet_wrap(~ variable, ncol=3, scales="free") +
  geom_map(data=germany, map=germany, aes(x=long, y=lat,map_id=region), fill="white", alpha=0.1, col="black")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






