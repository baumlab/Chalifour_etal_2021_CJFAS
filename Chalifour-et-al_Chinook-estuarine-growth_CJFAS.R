## Code for Chalifour L, Scott DC, MacDuffee M, Stark S, Dower JF, Beacham TD, Martin TG, and Baum JK. Chinook salmon exhibit long-term rearing and early marine growth in the Fraser River, B.C., a large urban estuary. 

library(here)
##load data
chin<- read.csv("data/chin.csv") #all Chinook catch 2016 and 2017 including stock

##summarize Harrison otoliths extracted by month for Table 1
table(chin$oto.extracted[chin$month=="March"], chin$Habitat[chin$month=="March"])
table(chin$oto.extracted[chin$month=="April"], chin$Habitat[chin$month=="April"])
table(chin$oto.extracted[chin$month=="May"], chin$Habitat[chin$month=="May"])
table(chin$oto.extracted[chin$month=="June"], chin$Habitat[chin$month=="June"])
table(chin$oto.extracted[chin$month=="July"], chin$Habitat[chin$month=="July"])

##prep data
chin<- chin[chin$hatchery !=("Y"),] #remove 18 fin-clipped fish (16 from 2016, 2 from 2017)
chin<- chin[chin$hatchery !=c("T"),] #remove 5 additional fish that had potential thermal marks on otoliths (2016)
oto<- read.csv("data/otolith.measurements.csv") #otolith measurements 2016 wild Harrison Chinook

#re-order habitat and month factors
chin$Habitat<- factor(chin$Habitat, levels = c("Marsh", "Sand flat", "Eelgrass"))
chin$month<- factor(chin$month, levels = c("March", "April", "May", "June", "July", "October"))

library(plyr)  #combine otolith and 2016 Chinook catch data 
oto <- merge(chin, oto, by = "fishid")  

####################### Back-calculation of Forklength at entry using Linear Regression against otolith widths (better relationship than radius for our data). Principals derived from Francis (1990). Back-calculation of fish length: a critical review. J. Fish Biol. 36:883 -902. Code for back-calculations adapted from sample code provided by Cameron Freshwater. #######################

### Proportional Calculation #1  
#Calculation using proportional hypothesis #1 - Scale Proportional Hypothesis (See Francis 1990)  
#SPH assumes that if the scale were n% larger/smaller when the fish was caught than the avg scale for that size of fish, the scale would remain n% larger/smaller than normal throughout  life. This principal can be applied to otolith radius or width as well.  

# if f(L) is mean scale radius for fish of length L and a and b are regression coefficients, >f(L) = a + bL or 
# >f(Li) = (Si/Sc) x f(Lc) 
# -> mean scale radius for fish of length L at time i = scale radius at time i / scale radius at time of capture x mean scale radius for fish of length L at time of capture. 

#Convert oto widths to mm from um  
oto.width<- oto$avg.w/1000  
EntryWidth<- oto$avg.fw/1000
#PART 1: SPH
### 1) Calculate otolith width based on fork length
OW.lm<- lm(oto.width ~ Forklength.mm, data=oto) 
### 2) Calculate mean OW for each fish based on regression
OW.PopCapture <- (coef(OW.lm)[2]*oto$Forklength.mm) + coef(OW.lm)[1] 
### 3) Calculate ratio of observed to expected
OW.Ratio<- oto.width/OW.PopCapture
### 4) Adjust observed measured otolith estuarine entry width by ratio to calculate expected
OW.ExpEntry <- EntryWidth/OW.Ratio
### 5) Estimate entry length based on proportion
OW.EntryFL <- (OW.ExpEntry - coef(OW.lm)[1])/coef(OW.lm)[2]  

### Proportional Calculation #2    
#Calculation using proportional hypothesis #2 - Body Proportional Hypothesis (See Francis 1990)  
#BPH assumes that if the body length were n% larger/smaller when the fish was caught than the avg length of fish for that size of scale, the fish would remain n% larger/smaller than normal throughout life.  

#Linear forms of equation: 
#  if g(s) is mean body length for fish with scale radius S and c and d are regression coefficients, >g(S) = c + dS or Li = [g(Si)/(gS)c]/Lc
# -> Length at time i = (mean length at scale radius i / mean length at scale radius capture) / mean length at capture   
# >in linear regression this also becomes Li = [(c + dSi)/(c + dSc)]*Lc 
# -> so Length at time i = [(intercept + scale radius(i) x slope)/(intercept + scaleradius(capture) x slope)] x Length at capture

### 1) Calculate FL based on OW
FL.lm<- lm(Forklength.mm ~ oto.width, data=oto)
### 2) Calculate mean FL for each fish based on regression
FL.PopCapture <- (coef(FL.lm)[2]*oto.width) + coef(FL.lm)[1]
### 3) Calculate ratio of observed to expected
FL.Ratio <- oto$Forklength.mm/FL.PopCapture
### 4) Calculate mean length for individual with observed entry width
FL.EntryLength <- (coef(FL.lm)[2]*EntryWidth) + coef(FL.lm)[1] 
### 5) Adjust by ratio
FL.EntryFL <- FL.Ratio * FL.EntryLength  

## Calculating uncertainty by comparing/averaging both methods (See Francis 1990)  
#  Take mean of both methods to reduce error - this is final entry fork length used in all future analyses  
oto$EntryFL <- (FL.EntryFL + OW.EntryFL)/2 

### Grab slope and intercept from Width SPH model
summary(OW.lm)
W_slope<- coef(OW.lm)[2]
W_intercept<- coef(OW.lm)[1] 

### Grab slope and intercept from Width BPH model
summary(FL.lm)
FL_slope<- coef(FL.lm)[2]
FL_intercept<- coef(FL.lm)[1]  

### Summary of two methods
summary(lm(FL.EntryFL~OW.EntryFL))

##### Plot exploration of methodology - Figure S1   
#save to tiff
tiff("final.figures/FigS1.tiff", width=7.2, height = 8, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 600)
#low res figure for initial submission
tiff("final.figures/FigS1_lowres.tiff", width=7.2, height = 8, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 100)

par(mfrow = c(2,2), oma = c(0.5,0.5,0.5,0.5), mar = c(3,2.4,0,0), mgp = c(1.5,0.5,0), bty="l", cex=1, cex.axis=0.85)
#Plot SPH (blue)
plot(oto.width~Forklength.mm, ylab = "Otolith width at capture (mm)", xlab = "Size at capture (mm)", data=oto, xlim = c(30, 115), ylim = c(0.5, 1.8), pch = c(1,0,2)[as.factor(oto$Habitat)], col = adjustcolor(c( "blue", "orange","green")[as.factor(oto$Habitat)]));abline(a=W_intercept, b=W_slope,lwd=1, col="blue")
text(110, 1.8,"A")

#Plot BPH (red)
par(mar=c(3,2.4,0,0))
plot(Forklength.mm~oto.width, ylab = "Size at capture (mm)", xlab = "Otolith width at capture (mm)", data=oto, xlim = c(0.5, 1.8), ylim = c(30, 115), pch = c(1,0,2)[as.factor(oto$Habitat)], col = adjustcolor(c( "blue", "orange","green")[as.factor(oto$Habitat)]));abline(a=FL_intercept, b=FL_slope, lwd=1, col = "red")
text(1.7, 115,"B")

#Plot forklength vs forklength for two methods
par(mar = c(3,2.4,0.5,0))
plot(FL.EntryFL~OW.EntryFL, xlim = c(0,80), ylim = c(0,85), ylab = "Entry size predicted by BPH (mm)", xlab = "Entry size predicted by SPH (mm)", pch = c(1,0,2)[as.factor(oto$Habitat)], col = adjustcolor(c( "blue", "orange","green")[as.factor(oto$Habitat)])); abline(lm(FL.EntryFL~OW.EntryFL))
legend(0,80, pch = c(1,0,2), y.intersp = 1, col = c( "blue", "orange","green"), c( "Marsh", "Sand flat","Eelgrass"), bty = "n", cex=0.8, pt.cex = 1.2)
text(80, 85,"C")

#Plot MEAN OF 2 METHODS (i.e. final entry fork length used); Size at entry vs size at capture
par(mar=c(3,2.4,0.5,0))
plot(EntryFL~Forklength.mm, ylab = "Size at entry (mm)", xlab = "Size at capture (mm)", ylim = c(15, 85), xlim = c(15, 110), data = oto, pch = c(1,0,2)[as.factor(oto$Habitat)], col = adjustcolor(c( "blue", "orange","green")[as.factor(oto$Habitat)]))
text(110, 85,"D")

dev.off()
#############
hist(oto$EntryFL, xlab = "Forklength at entry (mm)", ylim = c(0,40), xlim = c(10,80))

### Calculate daily growth in Forklength by dividing growth since entry by days in estuary
oto$Avg.day.FL<- (oto$Forklength.mm-oto$EntryFL)/oto$ndays.me 
summary(oto$Avg.day.FL)  #min 0.2837, median 0.5571, mean 0.57445, max 0.9211
sd(oto$Avg.day.FL)  #sd 0.129

### Size at capture and size at entry 
summary(oto$Forklength.mm) #min: 40.00, median: 62.00, mean: 62.04, max: 110.00 mm
summary(oto$EntryFL) #min: 18.53, median: 35.79, mean: 37.69, max: 67.57 mm
summary(oto$mass) #min: 0.513, median: 2.275, mean: 2.64, max: 14.822 g

#strong relationship between mass and fork.length -- use to derive mass at entry using forklength at entry
lmm <- lm(log(mass) ~ Forklength.mm, data=oto)
lengthvalues<- seq(30,120,1)
mass.exp<- exp(predict(lmm, list(Forklength.mm=lengthvalues)))

plot(mass ~ Forklength.mm, data = oto, main = "Mass to Forklength relationship", ylab = "Mass at capture (g)", xlab = "Forklength at capture (mm)", pch=16);text(85,.5,c("Adj. R-sq: 0.91"), cex=0.8);legend("topleft", pch = 16, col = c(1,3,2), c("Marsh", "Eelgrass", "Sand flat"), bty = "n"); lines(lengthvalues,mass.exp)

#Use intercept and slope of mass-fork length relationship to derive mass at entry from fork length at entry
oto$mass.entry<-0
for(i in 1:length(oto$EntryFL)){
  oto$mass.entry[i]<- exp(-2.269057+0.049116*oto$EntryFL[i])
}
#plot(mass~mass.entry, data=oto)

summary(oto$mass.entry) #min: 0.2569, median: 0.5997, mean: 0.7558, max: 2.8573 g

##################################################################################################################################
## Note Figure 1 (map) was created in QGIS; Figure 2 (otolith images) was created in Photoshop


## Figure 3: Chinook salmon emigration patterns in the Fraser River estuary - CPUE summarized by month for all Chinook salmon in 2016 and 2017; stock composition (% Harrison) by month pooled for both years.
summary(chin$stock)

library(reshape)
outmigration<- subset(chin, select=c("month", "Year"))
stock<- subset(chin, select=c("Year","month", "stock"))
stock$freq<- rep(1)
stock.table<-cast(stock, stock~month + Year)
write.csv(stock.table, "data/stock.table.raw.csv", row.names = FALSE)

#manually convert single fish caught in October to 3 so that it will be visible on plot
chin3<- chin
fish<- subset(chin, month=="October")
chin3<- rbind(chin3, fish,fish)
chin3$month<- factor(chin3$month, c("March", "April", "May", "June", "July", "August", "September", "October"))
chin3$freq<- 1

### note add sampling effort
effort<- read.csv("data/effort_year_month.csv")
effort$Month<- factor(effort$Month, c("March", "April", "May", "June", "July", "August", "September", "October"))

#add effort and summarise catch by abundance
chin4<- subset(chin3, select=c("Year", "month","freq"))
chin4<- ddply(chin4, .drop=FALSE, .(Year, month), summarise, abundance=sum(freq)) #use .drop=FALSE to retain 0 value combos
effort$month_year<- paste(effort$Month, effort$Year, sep="_") #create month-year identifier
chin4$month_year<- paste(chin4$month, chin4$Year, sep="_")
chin4$effort<- effort[match(chin4$month_year, effort$month_year),2] #add correct effort values for each month to chin4
chin4$CPUE<- chin4$abundance/chin4$effort

#stock composition as a proportion of total catch
chin.stock<- chin
levels(chin.stock$stock)<- c("Other","Other","Harrison","Other","Other","Other","Other","Other", "Other", "Other","Other","Other","Other", "Other","Other","Harrison") #re-name stocks to labels we will be using ("Harrison" and "W_Chilliwack" both labelled "Harrison")
chin.stock<- chin.stock[chin.stock$month != "October",]
chin.stock$month<- factor(chin.stock$month)
chin.stock$stock<- factor(chin.stock$stock, levels = c("Harrison", "Other")) #re-order levels
chin.stock.table<- table(chin.stock$stock, chin.stock$month)
chin.stock.table<- prop.table(chin.stock.table,2)
chin.stock.table


####### FINAL FIGURE 
o<- order(chin4$month) #re-order catch data frame so that bars are ordered by month & year
ForBars<- chin4[o,]

#save to tiff
tiff("final.figures/Figure3.tiff", width=7.2, height = 4, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 600)
#low res figure for initial submission
#tiff("final.figures/Figure3_lowres.tiff", width=7.2, height = 4, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 100)
par(mfrow=c(1,2), mgp = c(2.2,0.6,0), mar = c(3.5,3.2,2.5,0), bty = "l", las=2, xpd=NA, cex =  1, cex.axis=0.85);
b<-barplot(ForBars$CPUE, space=c(rep(c(1,0),8)), col=c("white", "black"), ylim = c(0,13), ylab = "Chinook CPUE 2016-2017", xaxt = "n", xlim=c(1,24))
axis(1, at = seq(2,23,3), c("March", "April", "May", "June", "July", "Aug.", "Sept.", "Oct."))
text(22, 14.3, "A")
legend(13,14, fill=c("white", "black"), legend=c("2016", "2017"),
       bty="n", cex=0.8)

par(mar = c(3.5,5.5,3,0.5), xpd=TRUE, las=2)
barplot(chin.stock.table, xlab = "", border=NA, col=c("grey30","darkgrey"), ylab = "Proportion of GSI samples")
text(5.5, 1.15, "B")
legend(1.5,1.2, pch = c(15), legend=c("Harrison","Other stocks"), col=c("grey30","darkgrey"), bty="n", x.intersp = 0.5, y.intersp = 1, text.width = 0.3, xpd=NA, cex=0.8, pt.cex=1, trace=TRUE)

dev.off()
#############################

## Figure 4: Portfolio of estuarine entry timing and residency period for Harrison Chinook. Estuarine entry point was determined using chemical analysis of sagittal otoliths via LA-ICP-MS. Daily growth and number of days after entry were subsequently determined using visual measurements of the otoliths

#Add column for estuarine entry date, based on the date caught and the number of days in the estuary derived from otolith analyses
oto$EntryDay<- oto$J.date-oto$ndays.me

# Assess relationship between residency period and entry date, removing outlier fish that were caught early in the season at very small size (Set cut-off at 45 mm fork length, as this size and below were never caught outside of the marsh and are assumed to be below the minimum body size required for successful ocean emigration, therefore they would have had longer residency if we had not caught them). Smallest sandflat fish is 54, smallest eelgrass is 56mm, so by setting min to 45mm we are being conservative in estimating a minimum size threshold to handle saltier water in outer habitats. This clarifies relationship below, that residency period is longest for fish that enter earliest.   

out<- oto[oto$Forklength.mm<46,] # 6 fish, entering in March (1) and April (5) at an entering fork length of between 33 and 40 mm. 
oto_res<- subset(oto, Forklength.mm > 45)

##Assess relationship between residency and entry date
#all otos
lm_res<- lm(ndays.me~EntryDay, data=oto_res)  #R2 = 0.5473, p=6.07E-16, slope -0.64 day residency for every day later of entry

#habitat as fixed effect - examined with and without eelgrass, and found in both cases that the apparent difference between the sandflat and marsh fish is not significant
lm_res_hab<- lm(ndays.me~EntryDay + Habitat, data = oto_res) #R2 = 0.6094, p < 2.2E-16, HabitatSand flat Estimate = -4.4, Pr(>|t|) = 0.17055

#sand flat fish only
lm_res_sf<- lm(ndays.me~EntryDay, data = oto_res[oto_res$Habitat=="Sand flat",]) #R2 = 0.6998, p=1.514E-05, slope -1.03 day residency for every day later of entry

#marsh fish only
lm_res_ma<- lm(ndays.me~EntryDay, data = oto_res[oto_res$Habitat=="Marsh",]) #R2 = 0.619, p=5.638E-14, slope -0.69 day residency for every day later of entry

#Use identify function to help examine patterns in the two clusters of marsh fish -- looked at variety of factors and determined catch date (J.date) has an apparent pattern
plot(ndays.me~EntryDay, data=oto_res[oto_res$Habitat=="Marsh",], bty="l")
identify(oto_res$ndays.me~oto_res$EntryDay, labels = as.character(oto_res$J.date))

#group fish by catch date to separate two clusters
ma_early<- subset(oto_res, J.date<145)
ma_late<- subset(oto_res, J.date>144)

#Test to see if groups represent clusters -- they do.
plot(ndays.me~EntryDay, data=ma_early, ylim = c(0, 90), xlim = c(60, 140), col="black")
par(new=T)
plot(ndays.me~EntryDay, data=ma_late, yaxt="n", xaxt="n", ylim = c(0, 90), xlim = c(60, 140),col="blue")

#Assess relationship
lm_entry_ma_catchday<- lm(ndays.me~EntryDay + J.date, data = oto_res[oto_res$Habitat=="Marsh",]) #R2 = 1, p<2.2E-16, slope -1 day residency for every day later of entry and +1 day residency for every day later of capture

ma_early$group<- "early"
ma_late$group<- "late"
oto_res2<- rbind(ma_early, ma_late)

#compare groups by residency period
catch_res_anova<- aov(ndays.me~group, data=oto_res2)
TukeyHSD(catch_res_anova) #late-early comparison is not significant (p=0.299)

#compare groups by entry timing
catch_ent_anova<- aov(EntryDay~group, data=oto_res2)
TukeyHSD(catch_ent_anova) #late-early comparison is significant (diff = 18.96, p=2.1E-06)
#### Based on above, stick to plotting overall residency vs entry timing relationship, with habitats highlighted using symbols.

#Bin Julian date of entry into months (note 2016 leap year but no anomalies along edges)
oto$EntryMonth<- NA 
for (i in 1:nrow(oto)) {
  if (oto$EntryDay[i] < 60) {
    oto$EntryMonth[i]<- "February"
  } else if (oto$EntryDay[i] > 59 & oto$EntryDay[i] < 91) {
    oto$EntryMonth[i]<- "March"
  } else if (oto$EntryDay[i] > 90 & oto$EntryDay[i] < 121) {
    oto$EntryMonth[i]<- "April"
  } else if (oto$EntryDay[i] > 120 & oto$EntryDay[i] < 152) {
    oto$EntryMonth[i]<- "May"
  } else {
    oto$EntryMonth[i]<- "NA"
  }
}
oto$EntryMonth<- factor(oto$EntryMonth, levels = c("February", "March", "April", "May"))

##### plot Fig. 4 residency, growth
#save to tiff
tiff("final.figures/Figure4.tiff", width=7.2, height = 8, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 600)
#low res figure for initial submission
#tiff("final.figures/Figure4_lowres.tiff", width=7.2, height = 8, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 100)

par(mfrow = c(2,2), mgp = c(1.75,0.5,0), oma = c(0.5,0.5,0.5,0.5), mar = c(3,3,0,0), bty = "l", xpd=FALSE, cex=1, cex.axis=0.85); 
hist(oto$EntryDay, ylim = c(0,21), xlim = c(40,160), main = "", xlab = "Estuarine entry (Julian day)", ylab = "Frequency (number of fish)");
abline(v=mean(oto$EntryDay), col="blue", lwd=1.5); 
abline(v=median(oto$EntryDay), col="darkblue", lty = 2, lwd = 1.5);
text(155, 21, "A");

par(mar = c(3,1.5,0,0));
hist((oto$ndays.me), xlab = "Estuarine residency (days)", ylab = "", yaxt = "n", ylim = c(0,21), xlim = c(0,100), main = "");
abline(v=mean(oto$ndays.me), col="blue", lwd=1.5); 
abline(v=median(oto$ndays.me), col="darkblue", lty = 2, lwd=1.5); 
axis(2, at = c(0,5,10,15,20), labels = FALSE);
text(98, 21, "B");
legend(60,20, lty = c(1,2), lwd = 1.5, col = c("blue", "darkblue"), c("Mean", "Median"), bty="n", cex = 0.8, x.intersp = 0.5)

par(mar = c(3,3,0,0), bty="l");
plot(ndays.me~EntryDay, data=oto[oto$Forklength.mm>45,], ylim = c(0,120), xlim = c(40,150),ylab = "Estuarine residency (days)", xlab = "Estuarine entry (Julian day)", pch = c(16,15,17)[as.factor(oto$Habitat)], col = adjustcolor(c( "blue", "orange","green")[as.factor(oto$Habitat)], alpha.f = 0.8));
abline(lm_res, lwd=1.5); 
legend(40,128, pch = c(16,15,17), col = c("blue", "orange", "green"), c("Marsh", "Sand flat", "Eelgrass"), bty = "n", cex = 0.8, pt.cex = 1.2)
text(145, 121,"C")

par(mar = c(3,1.5,0,0), bty="l");
plot(ndays.me~EntryMonth, data=oto[oto$Forklength.mm>45,], ylim = c(0,120), ylab = "", xlab = "Estuarine entry (month)", yaxt = "n");
axis(2, at = c(0,20,40,60,80,100,120), labels = FALSE);
text(4.4, 121,"D")
dev.off()
#######################################

## Figure 5: Otolith edge microchemistry in relationship to habitat
#examine difference between otolith edge Sr by habitat using ANOVA
Sr_aov<-aov(Sr.edge~Habitat, data=oto)
TukeyHSD(Sr_aov) #only Sand flat - Marsh comparison is significant (p = 0.0001482)

#linear regressions of Sr to residency for each habitat - only significant for marsh
marSr<- lm(Sr.edge~ndays.me, data=oto[oto$Habitat==c("Marsh"),]) #Rsq 0.1818, p 0.0003205
sandSr<-lm(Sr.edge~ndays.me, data=oto[oto$Habitat==c("Sand flat"),]) #Rsq 0.1643, p 0.0952
eelSr<- lm(Sr.edge~ndays.me, data=oto[oto$Habitat==c("Eelgrass"),]) #Rsq 0.1608, p 0.4307
#overall relationship between Sr increasing as ndays.me increases is only consistent pattern
Sr<- lm(Sr.edge~ndays.me, data=oto) #Rsq 0.2378, p 9.475E-7

#######
#Fig 5 2-panel plot

#save to tiff
tiff("final.figures/Figure5.tiff", width=7.2, height = 4, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 600)
#low res figure for initial submission
tiff("final.figures/Figure5_lowres.tiff", width=7.2, height = 4, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 100)

par(mfrow = c(1,2), oma = c(0.5,0.5,1,1), mar = c(3,3,0,0), mgp = c(1.8,0.5,0), bty = "l", cex=1, cex.axis=0.85)

#plot edge Sr:Ca by habitat
plot(Sr.edge~Habitat, data=oto, ylab = c("Otolith edge Sr:Ca (mmol/mol)"), xlab = c("Habitat of capture"), ylim = c(0,4));text(3.3, 4,"A")

#plot edge Sr:Ca over residency
par(mar = c(3,1,0,0))
plot(Sr.edge~ndays.me, data=oto,  ylab = "", yaxt="n", xlab = c("Estuarine residency (days)"), ylim = c(0,4), xlim = c(0, 100), pch = c(16,15,17)[as.factor(oto$Habitat)], col = adjustcolor(c( "blue", "orange","green")[as.factor(oto$Habitat)], alpha.f = 0.8), bty = "l"); axis(2, at = c(0,1,2,3,4), labels = FALSE); abline(Sr, lty=1, lwd = 1.5, col="black");text(98, 4,"B")

dev.off()
#######################################

## Figure 6: Size of Harrison Chinook at capture in relationship to habitat; panel A: boxplot showing 2016 and 2017 Harrison Chinook catch by habitat
harrison<- chin[which(chin$stock%in%c("Harrison", "W_Chilliwack")),]

#examine difference between mean size in each habitat using ANOVA
harr_aov<-aov(Forklength.mm~Habitat, data=harrison)
tukey<- TukeyHSD(harr_aov) #all relationships significant
tukey$Habitat[,"p adj"] #Sand flat-Marsh: 7.696543e-11     Eelgrass-Marsh: 7.693690e-11 Eelgrass-Sand flat:2.745533e-06 

#add otolith analyzed identifier to harrison catch
oto_list<- data.frame(oto = "Y", fishid = oto$fishid)
run<- data.frame(oto = "N", fishid = rep_len("NA", length.out = 504-91))
oto_l<- rbind(oto_list,run)
harrison$oto<- NA
harrison$oto<- oto_l[match(harrison$fishid,oto_l$fishid),1]

#linear regressions of forklength for each habitat
mar<- lm(Forklength.mm~J.date, data=harrison[harrison$Habitat==c("Marsh"),]) #R2 = 0.4646, p = 2.2E-16
sand<-lm(Forklength.mm~J.date, data=harrison[harrison$Habitat==c("Sand flat"),]) #R2 = 0.3809, p = 6.069E-9
eel<- lm(Forklength.mm~J.date, data=harrison[harrison$Habitat==c("Eelgrass"),]) #R2 = 0.1517, p = 1.617E-4

#######
#Fig 6 2-panel plot

#save to tiff
tiff("final.figures/Figure6.tiff", width=7.2, height = 4, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 600)
#low res figure for initial submission
tiff("final.figures/Figure6_lowres.tiff", width=7.2, height = 4, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 100)

par(mfrow = c(1,2), oma = c(0.5,0.5,1,1), mar = c(3,3,0,0), mgp = c(1.8,0.5,0), bty = "l", cex=1, cex.axis=0.85)
plot(Forklength.mm~Habitat, data=harrison, ylab = c("Fork length at capture (mm)"), xlab = c("Habitat of capture"), ylim = c(20,141));text(3.3, 140,"A")

#plot both years size by habitat
par(mar = c(3,1,0,0))
plot(Forklength.mm~J.date, data=harrison,  ylab = "", yaxt="n", xlab = c("Date of capture (Julian day)"), ylim = c(20,141), xlim = c(80,201), pch = c(1,0,2)[as.factor(harrison$Habitat)], col = adjustcolor(c( "blue", "orange","green")[as.factor(harrison$Habitat)], alpha.f = 0.8), bty = "l"); points(Forklength.mm~J.date, data=harrison[harrison$oto == "Y",], col = adjustcolor(c("blue", "orange","green")[as.factor(harrison$Habitat)], alpha.f = 0.8), pch = c(16,15,17)[as.factor(harrison$Habitat)]); axis(2, at = c(20,40,60,80,100,120,140), labels = FALSE); legend(80,147, pch = c(1,0,2,16), y.intersp = 1, col = c( "blue", "orange","green", "black"), c( "Marsh", "Sand flat","Eelgrass", "Otolith analyzed"), bty = "n", cex=0.8, pt.cex = 1.2);abline(eel, lty=1, lwd = 1.5, col="green");abline(mar, lty=2, lwd = 1.5, col="blue");abline(sand, lty=6, lwd = 1.5, col="orange");text(200, 140,"B") 

dev.off()

#########################################

#Figure 7 Estuarine growth as a function of time and body size

### A. growth over time
summary(oto$J.date)
date<- seq(98, 192, length.out = 91)
#exponential relationship is stronger for Jday
expgrday<- lm(log(Avg.day.FL)~J.date, data=oto) #Adjusted R-squared:  0.1605,p-value: 8.32e-05 
pred.day<- exp(predict(expgrday, list(J.date=date)))

### B. growth by body size
gr_FL<- lm(Avg.day.FL~Forklength.mm, data=oto) #R2=0.343, p=1.075E-09

######## Total estuarine growth over size at entry - i.e. proportional growth for season or body size increase over the season... Fork length at capture - fork length at entry / fork length at entry; expressed as %
oto$total<- NA
oto$total<- ((oto$Forklength.mm - oto$EntryFL)/oto$EntryFL)*100

### C. Total growth by residency time - exponential
tot_res<- lm(total~ndays.me, data=oto) #R: 0.7023, p=2.2E-16
tot_res_exp<- lm(log(total)~ndays.me, data=oto) #R: 0.782, p=2.2E-16
res<- seq(5.25, 89.16, length.out = 91)
respred<- exp(predict(tot_res_exp, list(ndays.me = res)))

### D. Total estuarine growth by entry size
entry<- seq(18.5,68,length.out = 91)
sizeentry<- lm(log(total)~EntryFL, data=oto) #R = 0.4294, p=1.103 E-12
sizepred<- exp(predict(sizeentry, list(EntryFL = entry)))

#########################################
#save to tiff
tiff("final.figures/Figure7.tiff", width=7.2, height = 8, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 600)
#low res figure for initial submission
tiff("final.figures/Figure7_lowres.tiff", width=7.2, height = 8, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 100)

par(mfrow = c(2,2), oma = c(0.5,0.5,0.5,0.5), mar = c(3,2.4,0,0), mgp = c(1.6,0.5,0), bty="l", cex=1, cex.axis=0.85)
plot(Avg.day.FL~J.date, data=oto, ylab="Mean daily growth (mm)", xlab = "Date of capture (Julian day)", ylim = c(0.25,1.05), yaxt = "n", pch = c(16,15,17)[as.factor(oto$Habitat)], col = adjustcolor(c( "blue", "orange","green")[as.factor(oto$Habitat)], alpha.f = 0.8)); 
lines(date, pred.day, lwd=1.5); 
text(190,1.05, "A"); 
axis(2, at = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), labels = c("",0.4,"",0.6,"",0.8,"","1.0"));
legend(100,1.1, pch = c(16,15,17), y.intersp = 1, col = c( "blue", "orange","green"), c( "Marsh", "Sand flat","Eelgrass"), bty = "n", cex=0.8, pt.cex = 1.2) 

par(mar=c(3,2.3,0,0.5))
plot(Avg.day.FL~Forklength.mm, ylab="", yaxt = "n", xlab="Fork length at capture (mm)", ylim = c(0.25,1.05), pch = c(16,15,17)[as.factor(oto$Habitat)], col = adjustcolor(c( "blue", "orange","green")[as.factor(oto$Habitat)], alpha.f = 0.8), data=oto); 
abline(gr_FL, lwd=1.5);
axis(2, at = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), labels = FALSE);
text(107, 1.05,"B")

par(mar = c(3,2.4,0.5,0))
plot(total~ndays.me, data=oto,  ylab = "Total growth (%)", xlab="Estuarine residency (days)", ylim = c(0,265), xlim = c(5,90), pch = c(16,15,17)[as.factor(oto$Habitat)], col = adjustcolor(c( "blue", "orange","green")[as.factor(oto$Habitat)], alpha.f = 0.8)); 
axis(1, at = c(10,20,30,40,50,60,70,80,90), labels = FALSE);
lines(res, respred, lwd=1.5);
text(88, 260,"C")

par(mar=c(3,2.3,0.5,0.5))
plot(total~EntryFL, data=oto, xlab = "Fork length at entry (mm)", ylab = "", yaxt = "n", ylim = c(0,265), pch = c(16,15,17)[as.factor(oto$Habitat)], col = adjustcolor(c( "blue", "orange","green")[as.factor(oto$Habitat)], alpha.f = 0.8)); 
lines(entry, sizepred, lwd=1.5);
axis(1, at = c(20,30,40,50,60,70), labels = FALSE);
axis(2, at = c(0,50,100,150,200,250), labels = FALSE);
text(65, 260,"D")

dev.off()
########################################################################


############### SUPPLEMENTAL FIGURES AND ANALYSES ###########################################
#note that Figure S1 is above in methodology. Table S1 is a subset of field measurements for water quality parameters in the estuary in 2016 and 2017 and was adapted from Table S1 of Chalifour et al (2019) Habitat use by juvenile salmon, other migratory fish, and resident fish species underscores the importance of estuarine habitat mosaics. Marine Ecology Progress Series 625: 145–162


###############Figure S2 - Evaluation of otolith microchemistry in relation to estuarine residency and habitat

### EVALUATION OF OTOLITH MICROCHEMISTRY IN RELATION TO ESTUARINE RESIDENCY AND HABITAT

tiff("final.figures/FigureS2.tiff", width=7.2, height = 8, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 600)
#low res
tiff("final.figures/FigureS2_lowres.tiff", width=7.2, height = 8, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 100)

par(mfrow = c(2,2), oma = c(0.5,0.5,1,1), mar = c(3,3,0,1), mgp = c(1.5,0.5,0), bty = "l", cex=1, cex.axis=0.85)

#Miller et al 2010 experiment documented brackish zone as Sr:Ca values above 1.55 mmol/mol (but some measured below) AND Ba:Ca values above 2umol/mol (range 0.35 - 1.31 umol/mol). If Ba:Ca is above but Sr:Ca is below, this is supposed to indicate freshwater, and if Ba:Ca is below but Sr:Ca is above, this indicates full marine (a few eelgrass fish meet this). Note that we also have eelgrass and sand flat fish in the "fresh" category, which seems highly unlikely. The system that Miller looked at likely had high natural levels of these isotopes relative to the Harrison and Fraser.

#Otolith microchemistry at otolith edge (just prior to capture)
plot((Ba.edge*10)~Sr.edge, data=oto, pch = c(1,0,2)[as.factor(oto$Habitat)], col=c("blue", "orange","green")[as.factor(oto$Habitat)], ylim = c(0, 44.5), xlim = c(0.5,3.5), ylab="Edge Ba:Ca (umol/mol)", xlab="Edge Sr:Ca (mmol/mol)"); 
abline(v=1.55); 
abline(h=2); 
text(1, 30, "Fresh", cex=0.7); 
text(2.6, 20, "Brackish", cex=0.7); 
text(3, 0.55, "Marine", cex=0.7); 
text(3.3, 44.5, "A"); 
legend(1.7, 48, pch = c(1,0,2), col=c("blue", "orange", "green"), c("Marsh", "Sand flat", "Eelgrass"), cex=0.8, pt.cex=1, bty="n")

#Otolith microchemistry at estuarine entry
plot((Ba.me*10)~Sr.me, data=oto, pch = c(1,0,2)[as.factor(oto$Habitat)], col=c("blue", "orange","green")[as.factor(oto$Habitat)], bty="l", ylab="Estuarine entry Ba:Ca (umol/mol)", xlab="Estuarine entry Sr:Ca (mmol/mol)", ylim = c(0,44.5), xlim = c(0.5,3.5)); 
abline(v=1.55); 
abline(h=2); 
text(1, 30, "Fresh", cex=0.7); 
text(2.6, 20, "Brackish", cex=0.7); 
text(3, 0.55, "Marine", cex=0.7); 
text(3.3, 44.5, "B")

#Only fish that were caught in the marsh had Ba:Ca concentrations at the otolith edge that were above 1.1 umol/mol Ca, reaffirming that they were in a less saline environment up until capture. The highest edge Ba:Ca levels were approaching 2 um, suggesting that the marsh has a high level of 138 Ba and can mimic freshwater conditions. These fish also had Sr:Ca levels approaching or above 1 mmol/mol, which affirms the estuarine environment and supports that the brackish/marine signature of the estuary did imprint on the otolith, and we are not just seeing a remnant mainstem signature. 

#Sr:Ca at the otolith edge increases, and Ba:Ca at the otolith edge decreases with increasing residency in the estuary. This suggests that tidal fluctuations in the estuary and an ontogenetic shift outward towards more saline habitats influence otolith microchemistry over time as the fish reside in this estuarine environment. 

#otolith edge Ba:Ca vs residency
plot((Ba.edge*10)~ndays.me, data=oto, ylab="Edge Ba:Ca (umol/mol)", xlab="Estuarine residency (days)", ylim = c(0,44.5), xlim = c(0, 100), pch = c(1,0,2)[as.factor(oto$Habitat)], col=c("blue", "orange","green")[as.factor(oto$Habitat)]); text(98, 44.5, "C")

#otolith edge Sr:Ca vs residency
plot(Sr.edge~ndays.me, data=oto,  ylab = "Edge Sr:Ca (mmol/mol)", xlab = c("Estuarine residency (days)"), ylim = c(0,4.45), xlim = c(0, 100), pch = c(1,0,2)[as.factor(oto$Habitat)], col = (c( "blue", "orange","green")[as.factor(oto$Habitat)]), bty = "l");text(98, 4.45,"D")

dev.off()
#########

#Further exploration and validation
babies<-subset(oto, ndays.me<14) #only 3 fish had lower than 14 days residency before capture, and of these one had Sr:Ca over 1.5 mmol/mol - strongly suggesting estuarine (caught in March after 13.6 days). All were caught in the marsh at forklengths 40-44 mm.Two caught in March, one in April. The one with the lowest Sr:Ca (0.8668 mmol/mol) was only in the estuary for an estimated 5.3 days, and had the highest Ba:Ca, so could very possibly have a partilally influenced mainstem signature. Essentially we saw an inflection right near the edge, and also caught the fish in the marsh, but this could indicate that the mainstem causes an inflection. Interestingly though, the Ba:Ca at estuarine entry (me) was less than 1 umol:mol as well, which would be very low for the mainstem, suggesting that it was the entry into brackish water that caused the mutual uptake of Ba with Sr towards the edge.
plot(Ba.edge~Sr.edge, data = babies, col=c("blue", "yellow","green")[as.factor(babies$Habitat)])
plot(Ba.me~Sr.me, data=babies)

############################################

###############Figure S3 - Relationship between mean daily growth and estuarine entry timing

tiff("final.figures/FigureS3.tiff", width=7.2, height = 4, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 600)
#low res
tiff("final.figures/FigureS3_lowres.tiff", width=7.2, height = 4, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 100)

par(mfrow = c(1,2), oma = c(0.5,0.5,0.5,0.5), mar = c(3,2.4,0,0), mgp = c(1.6,0.5,0), bty="l", cex=1, cex.axis=0.85)

plot(Avg.day.FL~EntryDay, data=oto, ylim = c(0.25,1.05), ylab="Mean daily growth (mm)", xlab="Estuarine entry (Julian day)", pch = c(1,0,2)[as.factor(oto$Habitat)], col = c("blue", "orange","green")[as.factor(oto$Habitat)]); 
text(140,1.05, "A");
axis(2, at = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), labels = c("",0.4,"",0.6,"",0.8,"","1.0"));
legend(40,1.1, pch = c(1,0,2), y.intersp = 1, col = c( "blue", "orange","green"), c( "Marsh", "Sand flat","Eelgrass"), bty = "n", cex=0.8, pt.cex = 1.2) 

par(mar=c(3,2.3,0,0.5))
plot(Avg.day.FL~ndays.me, ylim = c(0.25,1.05), xlim = c(0,100), ylab="", xlab="Estuarine residency (days)", pch = c(1,0,2)[as.factor(oto$Habitat)], col = c( "blue", "orange","green")[as.factor(oto$Habitat)], data=oto); 
text(98,1.05, "B");
axis(2, at = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), labels = c("","","","","","","",""))

dev.off()

############################################

###############Figure S4 - Allometric patterns of growth by estuarine entry
oto$Avg.day.prop<- oto$Avg.day/oto$Forklength.mm
gr_prop<- lm(Avg.day.prop~J.date, data=oto) #R2=0.090, p=0.003816 
gr_FL_prop<- lm(Avg.day.prop~Forklength.mm, data=oto) #R2=0.2195, p=2.82E-06


tiff("final.figures/FigureS4.tiff", width=7.2, height = 8, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 600)
tiff("final.figures/FigureS4_lowres.tiff", width=7.2, height = 8, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 100)

par(mfrow = c(2,2), oma = c(0.5,0.5,0.5,0.5), mar = c(3,2.4,0,0), mgp = c(1.6,0.5,0), bty="l", cex=1, cex.axis=0.85)

plot(Avg.day.FL~J.date, data=oto, ylim = c(0.25,1.05), ylab="Mean daily growth (mm)", xlab="", pch = c(1,0,2)[as.factor(oto$Habitat)], col = c("blue", "orange","green")[as.factor(oto$Habitat)]); 
lines(date, pred.day);
text(190,1.05, "A");
axis(2, at = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), labels = c("",0.4,"",0.6,"",0.8,"","1.0"));
legend(98,1.1, pch = c(1,0,2), y.intersp = 1, col = c( "blue", "orange","green"), c( "Marsh", "Sand flat","Eelgrass"), bty = "n", cex=0.8, pt.cex = 1.2) 

par(mar=c(3,2.3,0,0.5))
plot(Avg.day.FL~Forklength.mm, ylim = c(0.25,1.05), ylab="", xlab="", pch = c(1,0,2)[as.factor(oto$Habitat)], col = c( "blue", "orange","green")[as.factor(oto$Habitat)], data=oto);
text(107,1.05, "B"); abline(gr_FL) 

par(mar=c(3,2.4,0,0))
plot(Avg.day.prop~J.date, data=oto, ylim = c(0.025,0.105), ylab = "Proportional daily growth", xlab="Date of capture (Julian day)", pch = c(1,0,2)[as.factor(oto$Habitat)], col = c( "blue", "orange","green")[as.factor(oto$Habitat)]); text(190,0.105, "C"); 
abline(gr_prop)

par(mar=c(3,2.3,0,0.5))
plot(Avg.day.prop~Forklength.mm, ylim = c(0.025,0.105), xlab="Fork length at capture (mm)", ylab="",data=oto, pch = c(1,0,2)[as.factor(oto$Habitat)], col = c( "blue", "orange","green")[as.factor(oto$Habitat)]); text(107,0.105, "D"); 
abline(gr_FL_prop)

dev.off()
############################################

###############Figure S5 - Average daily width of otolith increments during the freshwater, early estuarine entry, and late estuary just prior to capture

## create new data frame that compiles mean daily growth for freshwater, early estuary, and late estuary

fishid<- oto$fishid

early<- oto[,c("fishid","EntryDay", "Avg.day.early", "EntryFL", "Avg.day.FL")]
colnames(early)<- c("fishid","EntryDay","Avg.day", "EntryFL", "Avg.day.FL")
early$fishid<-paste(fishid,"early", sep = "_")
early<- na.omit(early)
early$stage<- "early"

late<- oto[,c("fishid","EntryDay", "Avg.day.late", "EntryFL", "Avg.day.FL")]
colnames(late)<- c("fishid","EntryDay","Avg.day","EntryFL", "Avg.day.FL")
late$fishid<-paste(fishid,"late", sep = "_")
late<- na.omit(late)
late$stage<- "late"

fresh<- oto[,c("fishid","EntryDay", "Avg.day.fresh", "EntryFL", "Avg.day.FL")]
colnames(fresh)<- c("fishid","EntryDay","Avg.day","EntryFL", "Avg.day.FL")
fresh$fishid<-paste(fishid,"fresh", sep = "_")
fresh<- na.omit(fresh)
fresh$stage<- "fresh"

growth<-rbind(early, late, fresh) #complete data frame, now with mean daily growth in single column, labeled by factor of growth stage in second column

gr_anova<- aov(Avg.day~stage, data=growth)
TukeyHSD(gr_anova) #All comparisons are significant (<0.001)

gr_anova_fl<- aov(Avg.day.FL~stage, data=growth)
TukeyHSD(gr_anova_fl) #No comparisons are significant when converted to FL


growth$stage<- factor(growth$stage, levels = c("fresh", "early", "late"))

tiff("final.figures/FigureS5.tiff", width=6, height = 4, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 600)
tiff("final.figures/FigureS5_lowres.tiff", width=6, height = 4, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 100)

par(oma = c(0.5,0.5,0.5,0.5), mar = c(3,2.4,0,0), mgp = c(1.6,0.5,0), bty="l", cex=1, cex.axis=0.85)

plot(Avg.day~EntryDay, data=growth, xlim = c(40,150), col=c("lightblue", "blue", "darkblue")[growth$stage], pch=c(15,16,17)[growth$stage], ylab = "Mean daily growth (um)", xlab = "Estuarine entry date", bty="l"); 
legend(37,5.6, title = "Growth period", c("Freshwater", "Early estuary", "Late estuary"), col = c("lightblue", "blue", "darkblue"), pch = c(15,16,17), cex=0.8, pt.cex = 1.2, bty="n")

dev.off()
############################################

###############Table S2: Individual juvenile Harrison Chinook that experienced a decline of >5% in mean daily growth over time (Early:Fresh or Late:Early). 

#EVALUATION OF DECREASED DAILY GROWTH OVER TIME IN INDIVIDUAL FISH
## look at differential in mean daily growth for otoliths where we could measure freshwater (7-14 days prior to entry), early estuary (7-14 days after entry) and late estuary dailies (7-14 days prior to capture)
oto$early.fresh<- NA
oto$early.fresh<- oto$Avg.day.early/oto$Avg.day.fresh
hist(oto$early.fresh) #majority between 0 - 20% increase, but 5 below 0 indicating possible decrease

oto$late.early<- NA
oto$late.early<- oto$Avg.day.late/oto$Avg.day.early
hist(oto$late.early) #majority ~ 30% increase, but several below 0 and one considerably lower

## A value near 0 does not necessarily indicate starvation, but simply consistency, especially given that these are averages across 7-14 daily increments. Set threshold of "true decrease" at a negative differential of at least 5%, to allow for some margin of error.  

## lower growth measured in early estuarine entry vs fresh - 1 fish
lowfresh<- subset(oto, early.fresh<0.95)
## lower growth measured in late estuarine entry vs early - 4 fish
lowest<- subset(oto, late.early<0.95) 

## combine and view
low<- rbind(lowfresh, lowest)

low$Site<- factor(low$Site)
plot(EntryFL~Forklength.mm, data=low, pch = c(1:3)[as.factor(low$Site)]); legend(48, 48, pch = c(1:3), c("M2", "M4", "SF3"), bty = "n") #4 marsh and one sand flat fish; one entered at 34 mm and the rest between 44-47 mm. Three of the marsh fish were caught at site M2 (Westham marsh) and one at M4 (Brunswick Point), which are both marsh channels that are more exposed - possibly experienced greater tidal and temperature fluctuations. 


plot(Avg.day.FL~ndays.me, data=low); abline(h=0.575); abline(h=0.575+0.129, lty=2); abline(h=0.575-0.129, lty=2) #note that two are above mean daily growth for population, and one below. Other two are within 1 SD of the mean. 

#What stands out is that growth is variable for these fish over their estuarine residency. One fish came in small, and had poor growth throughout, which was the individual with high Sr:Ca edge. Interestingly, the fish that had low fresh-early growth recovered and had normal late-early growth. 

############# save data frame low for Table S2 #######################
low<- subset(low, select = c(Site, Sr.me, Sr.edge, EntryDay, J.date, ndays.me, EntryFL, Forklength.mm, Avg.day.FL, Avg.day, Avg.day.fresh, Avg.day.early, Avg.day.late, early.fresh, late.early))
write.csv(low, "data/TableS2.csv", row.names = FALSE)

