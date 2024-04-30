################################################################################
#################### Load libraries  ##############################
################################################################################

library(dplyr)
library(mice) #Imputation
library(metafor) # forest plot
library(Hmisc)
library(rms)
library(survival)
library(survminer)
library(zoo)
library(ggplot2)
library(tidymodels)
library(gtsummary)
library(dcurves)
library(ggsci)
library(cowplot)
library(survminer)
library(caret)
################################################################################
#################### Functions  ##############################
################################################################################

calibration.plot <- function( cox.model, year ){
  
  cal <- calibrate(cox.model, method = "boot", u= year*365, m=10, B=1000, bw=FALSE)
  pred <- as.data.frame(1 - cal[,'pred'])
  colnames(pred) <- c('x')
  
  cal.base.correct <- as.data.frame(1 - cal[,'calibrated.corrected'])
  colnames(cal.base.correct) <- c('Value')
  cal.base.correct$Type <- 'clinical'
  
  cal.base.observed <- as.data.frame(1 - cal[,'calibrated'])
  colnames(cal.base.observed) <- c('Value')
  cal.base.observed$Type <- 'observed'
  
  cal_plot_data <- rbind(cbind(pred, cal.base.correct), cbind(pred,cal.base.observed))
  
  cal.plot <- ggplot(cal_plot_data, aes(x = x, y = Value, color = Type)) +
    geom_line(linewidth = 2) +
    scale_x_continuous(name =paste("Predicted ", as.character(year) ,"-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7), limits=c(0, 0.7)) +
    scale_y_continuous(name=paste("Observed ", as.character(year) ,"-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7), limits=c(0, 0.7)) +
    labs(color = '', title = "Calibration curve") +
    scale_color_manual(values = c("clinical" = "#374e55", "observed" = "#df8f44"),labels = c("Internal validation", "Original")) +
    geom_abline(aes(intercept = 0, slope = 1, color='black'), guide = "none")	+
    theme_classic() +
    theme(plot.subtitle = element_text(size = 14),
          plot.title = element_text(size = 16, hjust = 0.5, face="bold", margin=margin(b = 20, unit = "pt")),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          panel.grid.major.y = element_line(),
          panel.grid.minor.y = element_line()) +
    theme(legend.position="bottom") 
  return(cal.plot)
}

calibration.values <- function( cox.model, year ){
  
  cal <- calibrate(cox.model, method = "boot", u= year*365, m=10, B=1000, bw=FALSE)
  output <- as.data.frame(1 - cal[,'pred'])
  colnames(output) <- c('x')
  output$y <- 1 - cal[,'calibrated.corrected']
  
  return(output)
}


validation.values <- function( cox.model, uhr1000.data, year ){
  val <- validate(cox.model, method="boot", B=1000, bw=FALSE)
  
  val.outcomes <- val.surv(cox.model, uhr1000.data, Surv(uhr1000.data$transdays, uhr1000.data$transtat), u = year*365)
  val.data <- data.frame(1 - val.outcomes$actual, 1-val.outcomes$p)
  names(val.data) <- c('observed','predicted')
  cal.large <- mean(val.data$observed) - mean(val.data$predicted)
  c.val = 0.5 + val[1,5]/2
  cal.slope = val[3,5]
  return(data.frame(c.val = c.val, 
              cal.large = cal.large,
              cal.slope = cal.slope))
}

val.calibration.slope <- function(cox.model, cox.late, uhr1000.data, year){
  # Calculate PI for both models first
  x.early <- model.matrix(cox.model) %*% coef(cox.model)
  xavg.early <- sum(coef(cox.model)*cox.model$means) 
  PI.early <- x.early - xavg.early # centered PI in development dataset (for discrimination later)
  
  # Calculate PI in validation set
  x.val.early <- model.matrix(cox.late) %*% coef(cox.model) # determine Xb in validation data 
  # by using coefficients of development data 
  PI.early.val <- x.val.early - xavg.early # center PI by using mean of PI of development data 
  
  # Regress out PI from model and calculate calibration slope in validation set
  uhr1000.data$PI.early.val <- PI.early.val
  fit.val.early <- cph(Surv(transdays, transtat) ~ PI.early.val,data = uhr1000.data, x=TRUE, y=TRUE, surv = TRUE, time.inc = year*365)
  #cat("Calibration slope for model developed with data from 1995 - 2016")
  return(fit.val.early$coefficients)
}

validation.plot <- function( cox.model, year, number, uhr1000.data){

  val.outcomes <- val.surv(cox.model, uhr1000.data, Surv(uhr1000.data$transdays, uhr1000.data$transtat), u = year*365)
  
  val.data <- data.frame(1 - val.outcomes$actual, 1-val.outcomes$p)
  colnames(val.data) <- c('x','Value')
  val.data$Type <- number
  
  return(val.data)
}

################################################################################
#################### Load UHR 1000 tables  ##############################
################################################################################

uhr1000 <- read.csv("path/to/dataset/", header = TRUE, sep=",")

# Calculate follow-up in months
uhr1000$transmonths <- (as.yearmon(as.character(uhr1000$enddate), format = "%d/%m/%Y") - as.yearmon(as.character(uhr1000$startdate), format = "%d/%m/%Y"))*12

################################################################################
#################### Convert and prepare variables  ##############################
################################################################################

# Convert GAF scores to SOFAS scores
gaf <- read.csv(paste(dirname(rstudioapi::getSourceEditorContext()$path),"/GAF_SOFAS.csv", sep = ""), header = TRUE, sep=",") # Load conversion table
names(gaf) <- c("GAF","SOFAS")
gaf$SOFAS <- round(gaf$SOFAS)

uhr1000$gaf_sofas <- uhr1000$sofas_currscore_0

for(index in c(1:nrow(gaf))){
  uhr1000[uhr1000$gafrat_0 == gaf[index,"GAF"] & !is.na(uhr1000$gafrat_0 == gaf[index,"GAF"]) & is.na(uhr1000$sofas_currscore_0), "gaf_sofas"] <- gaf[index,"SOFAS"]
}

# Convert Hamilton Depression scores to MADRS scores
hamd <- read.csv(paste(dirname(rstudioapi::getSourceEditorContext()$path),"/MADRS_HAMD.csv", sep = ""), header = TRUE, sep=",") #Load conversion table
names(hamd) <- c("MADRS","HAMD")

uhr1000$hamdep_madrs <- uhr1000$madrs_tot_0

for(index in c(1:nrow(hamd))){
  uhr1000[uhr1000$hamdep_tot_0 == hamd[index,"HAMD"] & !is.na(uhr1000$hamdep_tot_0 == hamd[index,"HAMD"]) & is.na(uhr1000$madrs), "hamdep_madrs"] <- hamd[index,"MADRS"]
}

# Calculate UHR category based on dominating UHR group
uhr1000$uhr_cat <- uhr1000$uhr_group

uhr1000[(uhr1000$uhr_group == 3 | uhr1000$uhr_group == 5  | uhr1000$uhr_group == 7) & !is.na(uhr1000$uhr_group), c('uhr_cat')] <- 1
uhr1000[uhr1000$uhr_group == 6 & !is.na(uhr1000$uhr_group), c('uhr_cat')] <- 2
uhr1000[uhr1000$uhr_group == 4 & !is.na(uhr1000$uhr_group), c('uhr_cat')] <- 3

# Code gender variable as 0 male and 1 female
uhr1000$gender <- uhr1000$gender - 1  
 
# Set categorical variables and reference
uhr1000$gender <- as.factor(uhr1000$gender)
uhr1000$group <- as.factor(uhr1000$group)
uhr1000$uhr_cat <- as.factor(uhr1000$uhr_cat)
uhr1000$uhr_cat <- relevel(uhr1000$uhr_cat, ref = 1)

# Nelson aalen estimator
uhr1000$nelsonaalen <- nelsonaalen(uhr1000, transdays, transtat)

# Extract baseline year
uhr1000$year <- format(as.Date(uhr1000$assessment_date_0, format="%d/%m/%Y"),"%Y")

# Assign subjects to groups for temporal validation
uhr1000$year.group <- with(uhr1000, ifelse(year  < 2008 , 1, 
                                           ifelse(year > 2016, 3, 2)))

# Save dataset prior to imputation
uhr1000.org <- uhr1000

################################################################################
#################### Select variables  ##############################
################################################################################

uhr1000 <- uhr1000 %>% dplyr::select(year.group, site, caarms_DS_sev_0, caarms_PA_sev_0, caarms_UTC_sev_0, bprs_tot_0, sans_tot_0, gaf_sofas, assessment_age, timeSxService, uhr_cat, transtat, transdays, transmonths, group, gender, study, nelsonaalen)

# Log transform time to service
uhr1000[uhr1000$timeSxService == 0 & !is.na(uhr1000$timeSxService), "timeSxService"] <- 1
uhr1000$timeSxService <- log(uhr1000$timeSxService)

# List missing data percentage
p_missing <- unlist(lapply(uhr1000, function(x) sum(is.na(x))))/nrow(uhr1000)
sort(p_missing[p_missing > 0], decreasing = TRUE)

# Exclude invdividuals with no follow-up time information
uhr1000 <- uhr1000[!is.na(uhr1000$transdays),]

################################################################################
#################### Imputation  ##############################
################################################################################

#### Multiple Imputation ####
imp <- mice(uhr1000, maxit=0)

# Extract predictorMatrix and methods of imputation 
predM <- imp$predictorMatrix
meth <- imp$method

# Setting values of variables to 0 in the predictor matrix that don't need to be imputed
#predM[, c("assessment_age")] <- 0
predM[, c("transtat")] <- 0
#predM[, c("gender")] <- 0
predM[, c("transdays")] <- 0
predM[, c("transmonths")] <- 0
predM[, c("year.group")] <- 0
predM[, c("site")] <- 0
predM[c("transdays"),] <- 0
predM[c("transmonths"),] <- 0
predM[c("year.group"),] <- 0
predM[c("site"),] <- 0
predM[, c("group")] <- 0

# Ordered categorical variables 
poly <- c("uhr_cat")

# Turn their methods matrix into the specified imputation models
meth[poly] <- "polyreg"

# Impute training set
uhr1000.imp <- mice(uhr1000, m = 50, maxit = 30, 
                  predictorMatrix = predM, 
                  method = meth, print =  FALSE, seed = 1235)

# Check logged events
#uhr1000.imp$loggedEvents

#Plot convergence
#plot(uhr1000.imp, layout = c(5, 5))

# Show density plots
#densityplot(uhr1000.imp)

data.comp <- complete(uhr1000.imp, action = "long", include = FALSE)

################################################################################
#################### Rubin's Rule  ##############################
################################################################################

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ caarms_DS_sev_0 + caarms_PA_sev_0 + caarms_UTC_sev_0 + bprs_tot_0 + 
                                 sans_tot_0 + gaf_sofas + group + timeSxService + assessment_age + gender + uhr_cat))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

################################################################################
#################### Forest Plot  ##############################
################################################################################

plot.data <- summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

plot.data <- plot.data[c(1:12),]

plot.data <- plot.data[,c("term", "estimate", "2.5 %", "97.5 %", "p.value")]

names(plot.data) <- c("Predictor", "HR", "lower", "upper", "P.value")
plot.data$Predictor <- as.character(plot.data$Predictor)

plot.data$P.value <- sprintf("%.3f", round(plot.data$P.value,3))

plot.data[nrow(plot.data)+1,] <- c(0, as.double(1.0), as.double(1.0), as.double(1.0), as.double(sprintf("%.3f",1.000)))

plot.data$Predictor <- c("CAARMS Disorganized Speech, severity",
                         "CAARMS Perceptual Abnormalities, severity",
                         "CAARMS Unusual Thought Content, severity",
                         "BPRS total score",
                         "SANS total score",
                         "SOFAS score",
                         "Received non-standard intervention treatment",
                         "Time to UHR service (log transformed)",
                         "Age at baseline",
                         "Female",
                         "UHR intake group - APS",
                         "UHR intake group - Trait",
                         "UHR intake group - BLIPS (ref)")

# colors to be used in the plot
colp <- "#6b58a6"
coll <- "#a7a9ac"

k <- nrow(plot.data)

forest(plot.data$HR, ci.lb=plot.data$lower, ci.ub = plot.data$upper, slab=plot.data$Predictor, annotate=TRUE, cex = 1, 
       cex.main = 1.5, cex.lab = 1, cex.axis = 1, cex.sub = 3, 
       efac = 0, 
       xlim = c(-2,3),
       at = c(0,1,2),
       alim = c(0, 2), 
       xlab="Hazard ratio (95% CI)", 
       psize = plot.data$HR*0.5,
       level = 95,
       refline = 1) # ilab = plot.data$P.value, ilab.xpos = c(5)
text(c(-2,1.8), c(14.5, 14.5), pos=4, 
     c("", "Hazard ratio (95% CI)"), font=2)

# redraw the CI lines and points in the chosen color
segments(plot.data$lower, k:1, plot.data$upper, k:1, col="black", lwd=1.5)
points(plot.data$HR, k:1, pch=18, cex = plot.data$HR*2,  col="white")
points(plot.data$HR, k:1, cex = plot.data$HR*2, pch=18, col="black")

#### Sensitivity analysis with study sites ####

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ caarms_DS_sev_0 + factor(study)))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ caarms_UTC_sev_0 +  factor(study)))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ caarms_PA_sev_0 +  factor(study)))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ gender +  factor(study)))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ bprs_tot_0 +  factor(study)))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ sans_tot_0 +  factor(study)))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ gaf_sofas +  factor(study)))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ assessment_age +  factor(study)))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ timeSxService +  factor(study)))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ group +  factor(study)))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ uhr_cat +  factor(study)))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ caarms_DS_sev_0 + caarms_UTC_sev_0 + sans_tot_0 + uhr_cat +  factor(study)))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

################################################################################
#################### Internal validation ######################
################################################################################

# Single imputation 
imp.nr <- sample(1:50,1, replace=F) 
#imp.nr = 19 # select specific imputation to reproduce results from paper
uhr1000.single <- subset(data.comp, .imp == imp.nr )

# Preperation for cox model fit
dd <- datadist(uhr1000.single)
options(datadist='dd')

# Set year for validation
year = 10

#### Create model with 3 CPP ####
uhr1000.single$study <- factor(uhr1000.single$study)

uhr.cox <- cph(Surv(transdays, transtat) ~ caarms_DS_sev_0 +  caarms_UTC_sev_0 + timeSxService + sans_tot_0 + gaf_sofas + uhr_cat , data = uhr1000.single, x=TRUE, y=TRUE, surv = TRUE, time.inc = year*365)

uhr.cox.pred <- survest(uhr.cox, uhr1000.single, times = year*365)
risk.uhr <- 1 - uhr.cox.pred$surv

# Output original model
uhr.cox

# Internal validation using bootstrapping
val.uhr <- validate(uhr.cox, method="boot", B=1000, bw=FALSE)
val.uhr

# Calibration plot
cal.uhr <- calibration.plot(uhr.cox, year)
cal.uhr

# Calculate PPV and NPV for different thresholds
uhr1000.boot <- bootstraps(uhr1000.single,
                           times = 1000,
                           strata = transtat,
                           apparent = FALSE)

pred.table.clin <- array(rep(NaN, 1000*12*4), c(1000, 12, 4))

for(x in c(1:1000)){
  
  data.b <- as.data.frame(uhr1000.boot[[x,1]])
  pred <- survest(uhr.cox, data.b, times = year*365)
  
  predicted.y <- pred$surv
  gt <- as.numeric(data.b$transtat == 0)
  i <- 1
  for (th in c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6)){
    predicted.th <- as.numeric((1-predicted.y) < th)
    conf_matrix<-confusionMatrix(factor(predicted.th),factor(gt), positive = "1")
    pred.table.clin[x,i,1] <- conf_matrix$byClass[1]
    pred.table.clin[x,i,2] <- conf_matrix$byClass[2]
    pred.table.clin[x,i,3] <- conf_matrix$byClass[3]
    pred.table.clin[x,i,4] <- conf_matrix$byClass[4]  
    i <- i + 1
  }
}

pred.Mean <- apply(pred.table.clin, c(2,3), mean, na.rm = TRUE)
pred.SD <- apply(pred.table.clin, c(2,3), sd, na.rm = TRUE)

pred.Mean <- round(pred.Mean*100, digits = 2)
pred.SD <- round(pred.SD*100, digits = 2)

# Decision curves
uhr1000.decision <-
  uhr1000.single %>%
  mutate(
    pr_failure18 =
      1 - summary(survfit(uhr.cox, newdata = uhr1000.single), times = year*365)$surv[1, ]
  )

dca(Surv(transdays, transtat) ~ pr_failure18,
    data = uhr1000.decision,
    time = year*365,
    thresholds = seq(0.05, 0.7, 0.01),
    label = list(pr_failure18 = "Prediction Model")
) %>%
  plot(smooth = TRUE)

dca(Surv(transdays, transtat) ~ pr_failure18,
    data = uhr1000.decision,
    time = year*365,
    thresholds = seq(0.05, 0.7, 0.01),
    label = list(pr_failure18 = "Prediction Model")
) %>%
  net_intervention_avoided() %>%
  plot(smooth = TRUE)

tmp <- dca(Surv(transdays, transtat) ~ pr_failure18,
           data = uhr1000.decision,
           time = year*365,
           thresholds = seq(0.05, 0.7, 0.01),
           label = list(pr_failure18 = "Prediction Model")
) %>%
  net_intervention_avoided()

 ggplot(tmp$dca, aes(x = threshold, y = net_intervention_avoided, color = label)) +
  stat_summary(fun.y = 'mean',  geom = 'line', linewidth = 1.5) +
   geom_line(linewidth = 1.5) +
  #stat_summary(fun.data = 'mean_sdl', alpha = 0.2) +
  scale_x_continuous(name = "Threshold Probability", breaks = c(0.1,0.2,0.3,0.4, 0.5,0.6,0.7), limits=c(0.1, 0.7)) +
  scale_y_continuous(name="Net reduction in intervention", breaks = c(0.00,0.25,0.50,0.75,1), limits=c(-0.2, 1)) +
  scale_color_manual(values = c("Treat All" = "#AEAEAE", "Treat None" = "#374e55", "Prediction Model" = "#df8f44"),
                     labels = c("Treat All", "Treat None", "UHR 1000+ model")) +
  labs(colour = "") + 
  theme_classic() +
  theme(plot.subtitle = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face="bold", margin=margin(b = 20, unit = "pt")),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line()) +
  theme(legend.position="bottom") 

 dca_plot <- ggplot(tmp$dca, aes(x = threshold, y = net_benefit, color = label)) +
  #stat_summary(fun.y = 'mean',  geom = 'line', linewidth = 1.5) +
  geom_smooth(method = "loess", se = FALSE , linewidth = 1.5) +
  #geom_line(linewidth = 1.5) +
  #stat_summary(fun.data = 'mean_sdl', alpha = 0.2) +
   scale_x_continuous(name = "Threshold Probability", breaks = c(0.1,0.2,0.3,0.4, 0.5,0.6,0.7), limits=c(0.1, 0.7)) +
   scale_y_continuous(name="Net benefit", breaks = c(0.00,0.05, 0.10, 0.15, 0.20), limits=c(-0.05, 0.20)) +
  scale_color_manual(values = c("Treat All" = "#AEAEAE", "Treat None" = "#374e55", "Prediction Model" = "#df8f44"),
                     labels = c("Treat All", "Treat None", "UHR 1000+ model")) +
  labs(colour = "", 
       title = "Decision curve analysis") + 
  theme_classic() +
  theme(plot.subtitle = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face="bold", margin=margin(b = 20, unit = "pt")),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line()) +
  theme(legend.position="bottom") 

ggarrange(cal.uhr3, dca_plot, ncol= 2)



################################################################################
#################### Temporal validation  #######################
################################################################################

set.seed(4826)

val.early <- vector()
val.all <- vector()

cal.large.early <- vector()
cal.large.all <- vector()

c.early <- vector()
c.all <- vector()

mape.early <- vector()
mape.all <- vector()

plot.early <- data.frame()
plot.all <- data.frame()

cal.data.early <- data.frame()
cal.data.all <- data.frame()

cal.early <- data.frame()
cal.all <- data.frame()

interv.avoid.early <- data.frame()
net.benefit.early <- data.frame()

interv.avoid.all <- data.frame()
net.benefit.all <- data.frame()

for(imp.nr in c(1:50)){
  uhr1000.single <- subset(data.comp, .imp == imp.nr )
  # Only use subjects from Melbourne
  uhr1000.single <- subset(uhr1000.single, site == 1 )
  # Splitting up by years
  uhr1000.early <- subset(uhr1000.single, year.group == 1 )
  uhr1000.all <- subset(uhr1000.single, year.group < 3)
  uhr1000.late <- subset(uhr1000.single, year.group == 3)
  # Splitting up by studies
  #uhr1000.early <- subset(uhr1000.single, study < 6 )
  #uhr1000.all <- subset(uhr1000.single, study < 7 | study == 11 | study == 13 | study == 8 )
  #uhr1000.late <- subset(uhr1000.single, study == 7 | study == 10 | study == 9 | study == 12)
  
  # Develop models with different baseline years
  dd <- datadist(uhr1000.early)
  options(datadist='dd')
  cox.early <- cph(Surv(transdays, transtat) ~ caarms_DS_sev_0 +  caarms_UTC_sev_0 + timeSxService + sans_tot_0 + gaf_sofas + uhr_cat, data = uhr1000.early, x=TRUE, y=TRUE, surv = TRUE, time.inc = 2*365)
  if(imp.nr == 1){
    tmp <- calibration.values(cox.early, 2)
    cal.data.early <- rbind(cal.data.early, tmp)
  }
  tmp <- validation.values( cox.early, uhr1000.early, 2 )
  cal.early <- rbind(cal.early, tmp)
  
  # Calculate c-index for both models in validation set
  surv.obj <- with(uhr1000.late,Surv(transdays,transtat)) 
  estimates <- survest(cox.early,newdata=uhr1000.late,times=2*365)$surv
  c.early[imp.nr] <- rcorr.cens(x=estimates,S=surv.obj)[1]
  
  # Calculate calibration curve and calibration-in-the-large
  uhr1000.late.double <- rbind(uhr1000.late, uhr1000.late[c(1:200),]) # Enlarge dataset with same samples to trick val.surv into getting a better hare fit
  tmp <- validation.plot(cox.early, 2, imp.nr, uhr1000.late.double)
  tmp <- tmp[1:nrow(uhr1000.late),] # Get original sample again
  mape.early[imp.nr] <- mean(abs(tmp$Value - tmp$x))
  cal.large.early[imp.nr] <- mean(tmp$x) - mean(tmp$Value)
  plot.early <- rbind(plot.early, tmp)
  
  # Decision curves internal
  uhr1000.early <-
    uhr1000.early %>%
    mutate(
      trans_pred =
        1 - summary(survfit(cox.early, newdata = uhr1000.early), times = 2*365)$surv[1, ]
    )
  
  tmp <- dca(Surv(transdays, transtat) ~ trans_pred,
             data = uhr1000.early,
             time = 2*365,
             thresholds = seq(0.05, 0.4, 0.01),
             label = list(trans_pred = "Prediction Model")
  ) %>%
    net_intervention_avoided()
  
  tmp.interv.early <-   cbind(tmp$dca$net_intervention_avoided, tmp$dca$threshold, tmp$dca$label)
  tmp.net.benefit.early <-  cbind(tmp$dca$net_benefit, tmp$dca$threshold, tmp$dca$label)
  
  # Decision curves temporal
  uhr1000.late <-
    uhr1000.late %>%
    mutate(
      trans_pred =
        1 - summary(survfit(cox.early, newdata = uhr1000.late), times = 2*365)$surv[1, ]
    )
  
  tmp <- dca(Surv(transdays, transtat) ~ trans_pred,
             data = uhr1000.late,
             time = 2*365,
             thresholds = seq(0.05, 0.4, 0.01),
             label = list(trans_pred = "Prediction Model")
  ) %>%
    net_intervention_avoided()
  
  tmp.interv.early <-   rbind(tmp.interv.early, cbind(tmp$dca$net_intervention_avoided[73:108], tmp$dca$threshold[73:108],rep(4,36)))
  tmp.net.benefit.early <-  rbind(tmp.net.benefit.early, cbind(tmp$dca$net_benefit[73:108], tmp$dca$threshold[73:108], rep(4,36)))
  
  interv.avoid.early <- rbind(interv.avoid.early, tmp.interv.early)
  net.benefit.early <- rbind(net.benefit.early, tmp.net.benefit.early)
  
  
  dd <- datadist(uhr1000.all)
  options(datadist='dd')
  cox.all <- cph(Surv(transdays, transtat) ~ caarms_DS_sev_0 +  caarms_UTC_sev_0 + timeSxService + sans_tot_0 + gaf_sofas + uhr_cat, data = uhr1000.all, x=TRUE, y=TRUE, surv = TRUE, time.inc = 2*365)
  if(imp.nr == 1){
    tmp <- calibration.values(cox.all, 2)
    cal.data.all <- rbind(cal.data.all, tmp)
  }
  tmp <- validation.values( cox.all, uhr1000.all, 2 )
  cal.all <- rbind(cal.all, tmp)
  
  # Calculate c-index for both models in validation set
  surv.obj <- with(uhr1000.late,Surv(transdays,transtat)) 
  estimates <- survest(cox.all,newdata=uhr1000.late,times=2*365)$surv
  c.all[imp.nr] <- rcorr.cens(x=estimates,S=surv.obj)[1]
  
  # Calculate calibration curve and calibration-in-the-large
  tmp <- validation.plot(cox.all, 2, imp.nr, uhr1000.late)
  cal.large.all[imp.nr] <- mean(tmp$x) - mean(tmp$Value)
  mape.all[imp.nr] <- mean(abs(tmp$Value - tmp$x))
  plot.all <- rbind(plot.all, tmp)  
  
  # Decision curves internal
  uhr1000.all <-
    uhr1000.all %>%
    mutate(
      trans_pred =
        1 - summary(survfit(cox.all, newdata = uhr1000.all), times = 2*365)$surv[1, ]
    )
  
  tmp <- dca(Surv(transdays, transtat) ~ trans_pred,
             data = uhr1000.all,
             time = 2*365,
             thresholds = seq(0.05, 0.4, 0.01),
             label = list(trans_pred = "Prediction Model")
  ) %>%
    net_intervention_avoided()
  
  tmp.interv.all<-   cbind(tmp$dca$net_intervention_avoided, tmp$dca$threshold, tmp$dca$label)
  tmp.net.benefit.all<-  cbind(tmp$dca$net_benefit, tmp$dca$threshold, tmp$dca$label)
  
  # Decision curves
  uhr1000.late <-
    uhr1000.late %>%
    mutate(
      trans_pred =
        1 - summary(survfit(cox.all, newdata = uhr1000.late), times = 2*365)$surv[1, ]
    )
  
  tmp <- dca(Surv(transdays, transtat) ~ trans_pred,
             data = uhr1000.late,
             time = 2*365,
             thresholds = seq(0.05, 0.4, 0.01),
             label = list(trans_pred = "Prediction Model")
  ) %>%
    net_intervention_avoided()
  
  tmp.interv.all<-   rbind(tmp.interv.all, cbind(tmp$dca$net_intervention_avoided[73:108], tmp$dca$threshold[73:108],rep(4,36)))
  tmp.net.benefit.all<-  rbind(tmp.net.benefit.all, cbind(tmp$dca$net_benefit[73:108], tmp$dca$threshold[73:108], rep(4,36)))
  
  interv.avoid.all<- rbind(interv.avoid.all, tmp.interv.all)
  net.benefit.all<- rbind(net.benefit.all, tmp.net.benefit.all)
  
  # Create Cox model with same predictors in validation set
  dd <- datadist(uhr1000.late)
  options(datadist='dd')
  cox.late <- cph(Surv(transdays, transtat) ~ caarms_DS_sev_0 +  caarms_UTC_sev_0 + timeSxService + sans_tot_0 + gaf_sofas + uhr_cat , data = uhr1000.late, x=TRUE, y=TRUE, surv = TRUE, time.inc = 2*365)
  
  #### Calculate calibration curve and c-index in validation set ####
  
  val.early[imp.nr] <- val.calibration.slope(cox.early, cox.late, uhr1000.late, 2)
  val.all[imp.nr] <- val.calibration.slope(cox.all, cox.late, uhr1000.late, 2)
}

cal.large <- as.data.frame(cbind(cal.large.early, cal.large.all))
names(cal.large) <- c('early','all')
cal.large$Type <- 'Calibration-in-the-large (Temporal)'

cal.slope <- as.data.frame(cbind(val.early,  val.all))
names(cal.slope) <- c('early','all')
cal.slope$Type <- 'Slope (Temporal)'

c.index <- as.data.frame(cbind(c.early, c.all))
names(c.index) <- c('early','all')
c.index$Type <- 'Harrell C index (Temporal)'

val.results <- rbind(c.index, cal.slope, cal.large)

internal.c <- data.frame(cbind(cal.early[,1], cal.all[,1], rep('Harrell C index (Internal)', imp.nr)))
names(internal.c) <- names(val.results)
val.results <- rbind( val.results, internal.c)

internal.large <-  data.frame(cbind(cal.early[,2], cal.all[,2], rep('Calibration-in-the-large (Internal)', imp.nr)))
names(internal.large) <- names(val.results)
val.results <- rbind( val.results, internal.large)

internal.slope <- data.frame( cbind(cal.early[,3], cal.all[,3], rep('Slope (Internal)', imp.nr)))
names(internal.slope) <- names(val.results)
val.results <- rbind( val.results, internal.slope)

val.results[, c(1:2)] <- sapply(val.results[, c(1:2)], as.numeric)
val.results$Type <- as.factor(val.results$Type)

table_validation <- val.results %>% 
  select(early, all, Type) %>% # keep only columns of interest
  tbl_summary(
    by = "Type",
    statistic = all_continuous() ~ "{mean}",        # stats and format for continuous columns
    digits = all_continuous() ~ 3,  
    label  = list(                                              # display labels for column names
      early ~ "1995 - 2007",
      all  ~ "1995 - 2016"  
    ),
    missing = "no", # don't list missing data separately;
  ) %>%
  add_ci(pattern = "{stat} ({ci})",
         style_fun = everything() ~ purrr::partial(style_number, digits = 3))%>% # add column with total number of non-missing observations
  modify_header(label = "Variable") %>% # update the column header
  bold_labels()

table_validation

table_validation %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = ".../validation_table_years.docx")

cal.early.plot <- ggplot() + 
  geom_smooth(aes(x=cal.data.early$x, y=cal.data.early$y, linetype = "Internal"),method="loess", se=FALSE, linewidth = 1.2, color = "black") +
  geom_smooth(aes(x=plot.early$Value, y=plot.early$x, linetype = "Temporal (2017 - 2020)"),method="loess", se=FALSE, linewidth = 1.2, color = "black") +
  #geom_line(aes(x=plot.early$Value, y = plot.early$lwl), color = "#ea801c", linetype = "dashed") +
  #geom_line(aes(x=plot.early$Value, y = plot.early$upl), color = "#ea801c", linetype = "dashed") +
  scale_x_continuous(name =paste("Predicted 2-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5), limits=c(0, 0.5)) +
  scale_y_continuous(name=paste("Observed 2-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5), limits=c(0, 0.5)) +
  geom_abline(aes(intercept = 0, slope = 1), color="black")	+
  labs(linetype = "",
       title = "1995 - 2007 (PACE 400)",
       subtitle = "A) Calibration curves")+ 
  scale_linetype_manual(values = c( "Internal" = "dotted" , "Temporal (2017 - 2020)" = "twodash"),
                        labels = c("1995 - 2007", "2017 - 2020")) +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face="bold", margin=margin(b = 10, unit = "pt")),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major.y = element_line(),
        legend.key.width = unit(2, "line"),
        panel.grid.minor.y = element_line()) +
  theme(legend.position="bottom") 

cal.all.plot <- ggplot() + 
  geom_smooth(aes(x=cal.data.all$x, y=cal.data.all$y, linetype = "Internal"),method="loess", se=FALSE, linewidth = 1.2, color = "black") +
  geom_smooth(aes(x=plot.all$Value, y=plot.all$x, linetype = "Temporal (2017 - 2020)"),method="loess", se=FALSE, linewidth = 1.2, color = "black") +
  #geom_line(aes(x=plot.all$Value, y = plot.all$lwl), color = "#ea801c", linetype = "dashed") +
  #geom_line(aes(x=plot.all$Value, y = plot.all$upl), color = "#ea801c", linetype = "dashed") +
  scale_x_continuous(name =paste("Predicted 2-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5), limits=c(0, 0.5)) +
  scale_y_continuous(name=paste("Observed 2-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5), limits=c(0, 0.5)) +
  geom_abline(aes(intercept = 0, slope = 1), color="black")	+
  labs(linetype = "",
       title = "1995 - 2016",
       subtitle = "A) Calibration curves")+  
  scale_linetype_manual(values = c( "Internal" = "dotted" , "Temporal (2017 - 2020)" = "twodash"),
                     labels = c("1995 - 2016", "2017 - 2020")) +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face="bold", margin=margin(b = 10, unit = "pt")),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major.y = element_line(),
        legend.key.width = unit(2, "line"),
        panel.grid.minor.y = element_line()) +
  theme(legend.position="bottom") 

net.benefit.early$V3 <- as.factor(net.benefit.early$V3)
net.benefit.early <- subset(net.benefit.early, V3!='1')

dc.early.plot <- ggplot(net.benefit.early, aes(x = V2, y = V1, linetype = V3)) +
  geom_smooth(method = "loess", se = FALSE , linewidth = 1.2, color = "black") +
  #stat_summary(fun.y = 'mean',  geom = 'line', linewidth = 1.5) +
  #(fun.data = 'mean_sdl', alpha = 0.2) +
  theme_classic() +
  scale_x_continuous(name = "Threshold Probability", breaks = c(0.1,0.2,0.3,0.4), limits=c(0.05, 0.4)) +
  scale_y_continuous(name="Net Benefit", breaks = c(0.00,0.05,0.1, 0.10), limits=c(-0.01, 0.12)) +
  scale_linetype_manual(values = c( "2" = "solid", "3" = "dotted" , "4" = "twodash"),
                        labels = c("Treat None", "1995 - 2007", "2017 - 2020")) +
  labs(linetype = "",
       subtitle = "B) Decision curve analysis")+ 
  theme(plot.subtitle = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face="bold", margin=margin(b = 10, unit = "pt")),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major.y = element_line(),
        legend.key.width = unit(2, "line"),
        panel.grid.minor.y = element_line()) +
  theme(legend.position="bottom") 

net.benefit.all$V3 <- as.factor(net.benefit.all$V3)
net.benefit.all <- subset(net.benefit.all, V3!='1')

dc.all.plot <- ggplot(net.benefit.all, aes(x = V2, y = V1, linetype = V3)) +
  geom_smooth(method = "loess", se = FALSE , linewidth = 1.2, color = "black") +
  #stat_summary(fun.y = 'mean',  geom = 'line', linewidth = 1.5) +
  #(fun.data = 'mean_sdl', alpha = 0.2) +
  theme_classic() +
  scale_x_continuous(name = "Threshold Probability", breaks = c(0.1,0.2,0.3,0.4), limits=c(0.05, 0.4)) +
  scale_y_continuous(name="Net Benefit", breaks = c(0.00,0.05,0.1, 0.10), limits=c(-0.01, 0.12)) +
  scale_linetype_manual(values = c( "2" = "solid", "3" = "dotted" , "4" = "twodash"),
                     labels = c("Treat None", "1995 - 2016", "2017 - 2020")) +
  labs(linetype = "",
       subtitle = "B) Decision curve analysis")+ 
  theme(plot.subtitle = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5, face="bold", margin=margin(b = 10, unit = "pt")),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid.major.y = element_line(),
        legend.key.width = unit(2, "line"),
        panel.grid.minor.y = element_line()) +
  theme(legend.position="bottom") 


plot_grid(cal.early.plot, cal.all.plot, dc.early.plot, dc.all.plot, ncol = 2, nrow = 2)

################################################################################
#################### Tables  ##############################
################################################################################

#### UHR1000 descriptive table ####

table_desc <- uhr1000.org %>% 
  select(marital_status, children, education_highest, accom_current, occup_current, Immigrated, Immigrated_parents, assessment_age, gender, timeSxService, uhr_group, transtat, group, caarms_DS_sev_0, caarms_PA_sev_0, caarms_UTC_sev_0, bprs_tot_0, sans_tot_0, gaf_sofas) %>% # keep only columns of interest
  mutate(
    gender = factor(gender, labels = c("Male", "Female")) ,
    uhr_group = factor(uhr_group, labels = c("BLIPS","Attenuated Psychosis","Atten. Psychosis+BLIPS","Vulnerability","BLIPS+Vulnerability","Vulnerability+Atten. Pscychosis","BLIPS+Atten. Psychosis+Vulnerability")) ,
    marital_status = factor(marital_status, labels = c("Married/De facto", "Not Married/De facto")),
    education_highest = factor(education_highest, labels = c("primary completed, secondary ongoing","secondary completed","TAFE/diploma/certificate unfinished/finished","undergrad unfinished","undergrad finished","post-grad unfinished","post grad finished")),
    accom_current = factor(accom_current, labels = c("None","crisis accomodation","boarding/rented room/foster","rented house/flat","own house/flat","house/flat with family of origin","Institution","other")),
    occup_current = factor(occup_current, labels = c("unemployed","full-time paid employed","part-time paid employed/casual/apprenticeship","student secondary","student post-secondary","caregiver")),
    Immigrated = factor(Immigrated, labels = c("Yes", "No")),
    Immigrated_parents = factor(Immigrated_parents, labels = c("Yes", "No")),
    
    group = factor(group, labels = c("Standard intervention treatment", "Non-standard intervention treatment")),
    transtat = factor(transtat, labels = c("Not transitioned", "Transitioned"))
  ) %>% 
  tbl_summary(  
    by = transtat,
    statistic = list(all_continuous() ~ "{mean} ({sd})",        # stats and format for continuous columns
                     all_categorical() ~ "{n} / {N} ({p}%)"),   # stats and format for categorical columns
    digits = all_continuous() ~ 1,                              # rounding for continuous columns
    type   = list(c(gender, group, transtat) ~ "dichotomous",
                  c(assessment_age, timeSxService, caarms_DS_sev_0, caarms_PA_sev_0, caarms_UTC_sev_0, bprs_tot_0, sans_tot_0, gaf_sofas) ~ "continuous",
                  c(uhr_group) ~ "categorical"),
    value = list(gender ~ "Female",
                 transtat ~ "Transitioned",
                 group ~ "Non-standard intervention treatment"),
    label  = list(                                              # display labels for column names
      assessment_age   ~ "Age (years)",                           
      gender ~ "Gender",
      timeSxService  ~ "Time between first symptoms and acceptance/treatment at UHR/PACE service",    
      uhr_group ~ "UHR syndrome",
      group ~ "Received non-standard intervention treatment as part of trial",
      caarms_DS_sev_0 ~ "CAARMS Disorganized Speech, Severity",
      caarms_PA_sev_0 ~ "CAARMS Perceptual Abnormalities, Severity", 
      caarms_UTC_sev_0 ~ "CAARMS Unusual Thought Content, Severity",
      bprs_tot_0 ~  "BPRS total score",
      sans_tot_0 ~ "SANS total score", 
      gaf_sofas ~ "GAF-SOFAS score",
      transtat ~ "Transition status",
      marital_status    ~ "Marital Status",
      children      ~ "Number of children",
      education_highest  ~ "Educational level",
      accom_current    ~ "Current accomodation",
      occup_current      ~ "Current occupation",
      Immigrated  ~ "Participant immigrated from another country",     
      Immigrated_parents  ~ "One or both parents of participant immigrated from another country"
    ),
    missing = "no", # don't list missing data separately;
  ) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "Variable") %>% # update the column header
  bold_labels()

table_desc %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = ".../descriptive_table.docx")

#### Demographic table ####

table_demog <- uhr1000.org %>% 
  select(assessment_age, gender, timeSxService, group, uhr_group, marital_status, children, education_highest, accom_current, occup_current, Immigrated, Immigrated_parents) %>% # keep only columns of interest
  mutate(
    gender = factor(gender, labels = c("Male", "Female")) ,
    uhr_group = factor(uhr_group, labels = c("BLIPS","Attenuated Psychosis","Atten. Psychosis+BLIPS","Vulnerability","BLIPS+Vulnerability","Vulnerability+Atten. Pscychosis","BLIPS+Atten. Psychosis+Vulnerability")) ,
    marital_status = factor(marital_status, labels = c("Married/De facto", "Not Married/De facto")),
    education_highest = factor(education_highest, labels = c("primary completed, secondary ongoing","secondary completed","TAFE/diploma/certificate unfinished/finished","undergrad unfinished","undergrad finished","post-grad unfinished","post grad finished")),
    accom_current = factor(accom_current, labels = c("None","crisis accomodation","boarding/rented room/foster","rented house/flat","own house/flat","house/flat with family of origin","Institution","other")),
    occup_current = factor(occup_current, labels = c("unemployed","full-time paid employed","part-time paid employed/casual/apprenticeship","student secondary","student post-secondary","caregiver")),
    Immigrated = factor(Immigrated, labels = c("Yes", "No")),
    Immigrated_parents = factor(Immigrated_parents, labels = c("Yes", "No")),
    group = factor(group, labels = c("Refused/Control", "Intervention treatment"))
  ) %>% 
  tbl_summary(     
    statistic = list(all_continuous() ~ "{mean} ({sd})",        # stats and format for continuous columns
                     all_categorical() ~ "{n} / {N} ({p}%)"),   # stats and format for categorical columns
    digits = all_continuous() ~ 1,                              # rounding for continuous columns
    type   = list(c(gender, marital_status, Immigrated, Immigrated_parents) ~ "dichotomous",
                  c(assessment_age, timeSxService) ~ "continuous",
                  c(uhr_group, education_highest, accom_current, occup_current, children) ~ "categorical"),
    value = list(gender ~ "Female",
                 marital_status ~ "Married/De facto",
                 Immigrated ~ "Yes",
                 Immigrated_parents ~ "Yes",
                 group ~ "Intervention treatment"),
    label  = list(                                              # display labels for column names
      assessment_age   ~ "Age (years)",                           
      gender ~ "Gender",
      timeSxService  ~ "Time between first symptoms and acceptance/treatment at UHR/PACE service",    
      uhr_group ~ "UHR category",
      marital_status    ~ "Marital Status",
      children      ~ "Number of children",
      education_highest  ~ "Educational level",
      accom_current    ~ "Current accomodation",
      occup_current      ~ "Current occupation",
      Immigrated  ~ "Participant immigrated from another country",     
      Immigrated_parents  ~ "One or both parents of participant immigrated from another country",
      group ~ "Part of cohort trial/Received placebo treatment vs. Received intervention treatment"
    ),
    missing = "no", # don't list missing data separately;
  ) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "Variable") %>% # update the column header
  bold_labels()

table_demog %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = ".../demog_table.docx")

#### UHR table ####

table_uhr <- uhr1000.org %>% 
  select(transtat, assessment_age, timeSxService, uhr_group) %>% # keep only columns of interest
  mutate(transtat = factor(transtat, labels = c("No transition", "Transition to FEP")),
         uhr_group = factor(uhr_group, labels = c("BLIPS","Attenuated Psychosis","Atten. Psychosis+BLIPS","Vulnerability","BLIPS+Vulnerability","Vulnerability+Atten. Pscychosis","BLIPS+Atten. Psychosis+Vulnerability"))) %>% 
  tbl_summary(
    by = uhr_group,
    statistic = all_continuous() ~ "{mean} ({sd})",        # stats and format for continuous columns
    digits = all_continuous() ~ 1,  
    type   = list(c(transtat) ~ "dichotomous",
                  c(assessment_age, timeSxService) ~ "continuous"),# rounding for continuous columns
    value = list(transtat ~ "Transition to FEP"),
    label  = list(                                              # display labels for column names
      transtat ~ "Transition status",
      assessment_age   ~ "Age (years)",                           
      timeSxService  ~ "Time between first symptoms and acceptance/treatment at UHR/PACE service"  
    ),
    missing = "no", # don't list missing data separately;
  ) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "Variable") %>% # update the column header
  bold_labels()

table_uhr %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = ".../uhr_table.docx")


#### Baseline table by study ####

table_bl <- uhr1000.org %>% 
  select(assessment_age, gender, timeSxService, group, uhr_cat, transtat, study, caarms_UTC_sev_0, caarms_NBI_sev_0, caarms_PA_sev_0, caarms_DS_sev_0, bprs_tot_0, bprs_ps_0, sans_tot_0, sans_af_0, sans_al_0, sans_av_0, sans_an_0, sans_at_0, sofas_currscore_0, gafrat_0, gf_social_score_0, gf_role_score_0) %>% # keep only columns of interest
  mutate(transtat = factor(transtat, labels = c("No transition", "Transition to FEP")),
         gender = factor(gender, labels = c("Male", "Female")) ,
         uhr_cat = factor(uhr_cat, labels = c("BLIPS","Attenuated Psychosis","Vulnerability")) ,
         group = factor(group, labels = c("Refused/Control", "Intervention treatment")),
         study = factor(study, labels = c("Prediction Study","First Intervention Study","Stress Cortisol Study","Lithium/Monitoring/Ris-Aus-9","Self","Step","SANE","Speak","PROSCAN","Neurapro","Charms","EU-GEI"))
  ) %>% 
  tbl_summary(
    by = study,
    statistic = all_continuous() ~ "{mean} ({sd})",        # stats and format for continuous columns
    digits = all_continuous() ~ 1,  
    type   = list(c(transtat, gender, group) ~ "dichotomous",
                  c(uhr_cat) ~ "categorical",
                  c(assessment_age, caarms_UTC_sev_0, caarms_NBI_sev_0, caarms_PA_sev_0, caarms_DS_sev_0, bprs_tot_0, bprs_ps_0, sans_tot_0, sans_af_0, sans_al_0, sans_av_0, sans_an_0, sans_at_0, sofas_currscore_0, gafrat_0, gf_social_score_0, gf_role_score_0) ~ "continuous"),# rounding for continuous columns
    value = list(transtat ~ "Transition to FEP",
                 gender ~ "Female",
                 group ~ "Intervention treatment"),
    label  = list( 
      assessment_age   ~ "Age (years)",                           
      gender ~ "Female",
      timeSxService  ~ "Time between first symptoms and acceptance/treatment at UHR/PACE service",    
      uhr_cat ~ "UHR category",
      group ~ "Part of cohort trial/Received placebo treatment vs. Received intervention treatment",
      transtat ~ "Transition status",
      caarms_UTC_sev_0 ~ "Unusual Thought Content, Severity Baseline",
      caarms_NBI_sev_0 ~ "Non-Bizarre Ideas, Severity Baseline",
      caarms_PA_sev_0 ~ "Perceptual Abnormalities, Severity Baseline",
      caarms_DS_sev_0 ~ "Disorganized Speech, Severity Baseline",
      bprs_tot_0 ~ "BPRS Total",
      bprs_ps_0 ~ "BPRS Psychotic Subscale",
      sans_tot_0 ~ "SANS Total",
      sans_af_0 ~ "SANS Affective flattening or blunting",
      sans_al_0 ~ "SANS Alogia",
      sans_av_0 ~ "SANS Avolition-Apathy",
      sans_an_0 ~ "SANS Anhedonia-Asociality",
      sans_at_0 ~ "SANS Attention",
      sofas_currscore_0 ~ "SOFAS (Social and Occupational Functioning Assessment Scale) current score",
      gafrat_0 ~ "GAF score (current)",
      gf_social_score_0 ~ "Global Functioning: Social Scale Current level (past month) Baseline",
      gf_role_score_0 ~ "Global Functioning: Role Scale Current level (past month) Baseline"
    ),
    missing = "no", # don't list missing data separately;
  ) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "Variable") %>% # update the column header
  bold_labels()

table_bl %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = ".../baseline_study_table.docx")


################################################################################
#################### UHR subgroup plot  ##############################
################################################################################

uhr1000.sub <- subset(uhr1000.org, !(is.na(uhr_cat)))
uhr1000.sub$transyears <- round(as.numeric(uhr1000.sub$transdays)/365.25, digits = 1)

fit <- survfit(Surv(transyears, transtat) ~ uhr_cat, data = uhr1000.sub)

cumhaz_plot <- ggsurvplot(fit,
                          conf.int = FALSE,
                          pval = TRUE,
                          pval.method = TRUE,
                          pval.coord = c(12.5,0.1),
                          pval.method.coord = c(12.5,0.15),
                          pval.method.size = 5,
                          xlab = "Time in years",   # customize X axis label.
                          ylab = "Risk of transition to psychotic disorder",
                          legend.labs = c("BLIPS", "APS or APS+Trait", "Trait"), 
                          risk.table.col = "strata", # Change risk table color by groups
                          ggtheme = theme_bw(), # Change ggplot2 theme
                          palette = c("#374e55", "#df8f44", "#00a1d5"),
                          fun = "event",
                          font.x = c(20),
                          font.y = c(20),
                          font.tickslab = c(15),
                          font.legend = c(15))

cumhaz_plot