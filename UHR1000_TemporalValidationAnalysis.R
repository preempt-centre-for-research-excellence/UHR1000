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
  
  cal.plot <- ggplot(cal_plot_data, aes(x = x, y = Value, color = Type, linetype = Type)) +
    geom_line(linewidth = 2) +
    scale_x_continuous(name =paste("Predicted ", as.character(year) ,"-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7), limits=c(0, 0.7)) +
    scale_y_continuous(name=paste("Observed ", as.character(year) ,"-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7), limits=c(0, 0.7)) +
    theme_bw() +
    labs(color = '', linetype = '') +
    scale_color_manual(values = c("clinical" = "#225ea8", "observed" = "#a1dab4"),
                       labels = c("Optimism corrected", "Original")) +
    scale_linetype_manual(values = c("clinical" = "solid", "observed" = "dotted"), labels = c("Optimism corrected", "Original")) +
    geom_abline(aes(intercept = 0, slope = 1, color='black'), guide = "none")	+
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 14),
          axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 13)) +
    theme(legend.position="bottom")
  return(cal.plot)
}

instability.plot <- function(imp.data, cox.formula, risk.orig, year, B ){
  
  set.seed(123)
  uhr1000.boot <- bootstraps(imp.data,
                             times = B,
                             strata = transtat,
                             apparent = FALSE)
  risk = data.frame()
  for(x in c(1:B)){
    data.b <- as.data.frame(uhr1000.boot[[x,1]])
    # Run Cox model
    cox.bl <- cph(cox.formula, data = data.b, x=TRUE, y=TRUE, surv = TRUE, time.inc = year*12)
    pred <- survest(cox.bl, imp.data, times = year*12)
    if(x == 1){
      risk <- 1 - pred$surv
      insta.plot <- plot(x=risk.orig, y=1 - pred$surv, 
           xlab = "Estimated risk from developed model", ylab = "Estimated risk from bootstrap models",
           pch = 19, frame = TRUE, xlim = c(0,1), ylim = c(0,1), bg = "gray", col = "gray", cex = 0.5, cex.axis = 1.25, cex.lab=1.5 ) 
    }else if(x %% 10 == 1){
      points(x=risk.orig, y=1 - pred$surv,  bg = "gray", col = "gray", cex = 0.5, pch = 19)
      risk <- cbind(risk, 1 - pred$surv)
    }else{
      risk <- cbind(risk, 1 - pred$surv)
    }
  }
  grid(nx = NA, ny = 5, lty = 2, lwd = 0.5, col = "gray")
  abline(coef = c(0,1), lwd = 2)
  risk.lower <- vector()
  risk.upper <- vector()
  risk.mape <- vector()
  for(ind in c(1:length(risk.orig))){
    tmp <- quantile(risk[ind,], probs = c(0.025, 0.975))
    risk.lower[ind] <- tmp[1]
    risk.upper[ind] <- tmp[2]
    risk.mape[ind] <- mean(abs(risk[ind,] - risk.orig[ind]))
  }
  risk.all <- as.data.frame(cbind(risk.orig, risk.lower, risk.upper))
  risk.all <- risk.all[order(risk.all$risk.orig),]
  lw.lower <- loess(risk.lower ~ risk.orig, data = risk.all, span = 0.5)
  lw.upper <- loess(risk.upper ~ risk.orig, data = risk.all, span = 0.3)
  lines(risk.all$risk.orig, lw.lower$fitted, lwd = 2, lty = 2)
  lines(risk.all$risk.orig, lw.upper$fitted, lwd = 2, lty = 2)
  
  return(risk.mape)
}

val.calibration.slope <- function(cox.model, cox.late, uhr1000.late){
  # Calculate PI for both models first
  x.early <- model.matrix(cox.model) %*% coef(cox.model)
  xavg.early <- sum(coef(cox.model)*cox.model$means) 
  PI.early <- x.early - xavg.early # centered PI in development dataset (for discrimination later)
  
  # Calculate PI in validation set
  x.val.early <- model.matrix(cox.late) %*% coef(cox.model) # determine Xb in validation data 
  # by using coefficients of development data 
  PI.early.val <- x.val.early - xavg.early # center PI by using mean of PI of development data 
  
  # Regress out PI from model and calculate calibration slope in validation set
  uhr1000.late$PI.early.val <- PI.early.val
  fit.val.early <- cph(Surv(transdays, transtat) ~ PI.early.val,data = uhr1000.late, x=TRUE, y=TRUE, surv = TRUE, time.inc = 2*365)
  #cat("Calibration slope for model developed with data from 1995 - 2016")
  return(fit.val.early$coefficients)
}

validation.plot <- function( cox.model, year, number, uhr1000.late){

  val.outcomes <- val.surv(cox.model, uhr1000.late, Surv(uhr1000.late$transdays, uhr1000.late$transtat), u = year*365)
  
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

# Extract just Melbourne sample
uhr1000 <- subset(uhr1000, site == 1)

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

################################################################################
#################### Tables  ##############################
################################################################################

#### UHR1000 descriptive table ####

table_desc <- uhr1000 %>% 
  select(assessment_age, gender, timeSxService, uhr_group, transtat, group, caarms_DS_sev_0, caarms_PA_sev_0, caarms_UTC_sev_0, caarms_NBI_sev_0, bprs_tot_0, sans_tot_0, gaf_sofas) %>% # keep only columns of interest
  mutate(
    gender = factor(gender, labels = c("Male", "Female")) ,
    uhr_group = factor(uhr_group, labels = c("BLIPS","Attenuated Psychosis","Atten. Psychosis+BLIPS","Vulnerability","BLIPS+Vulnerability","Vulnerability+Atten. Pscychosis","BLIPS+Atten. Psychosis+Vulnerability")) ,
    group = factor(group, labels = c("Standard intervention treatment", "Non-standard intervention treatment")),
    transtat = factor(transtat, labels = c("Not transitioned", "Transitioned"))
  ) %>% 
  tbl_summary(  
    by = transtat,
    statistic = list(all_continuous() ~ "{mean} ({sd})",        # stats and format for continuous columns
                     all_categorical() ~ "{n} / {N} ({p}%)"),   # stats and format for categorical columns
    digits = all_continuous() ~ 1,                              # rounding for continuous columns
    type   = list(c(gender, group) ~ "dichotomous",
                  c(assessment_age, timeSxService, caarms_DS_sev_0, caarms_PA_sev_0, caarms_UTC_sev_0, caarms_NBI_sev_0, bprs_tot_0, sans_tot_0, gaf_sofas) ~ "continuous",
                  c(uhr_group) ~ "categorical"),
    value = list(gender ~ "Female",
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
      caarms_NBI_sev_0 ~ "CAARMS Non-Bizarre Ideas, Severity",
      bprs_tot_0 ~  "BPRS total score",
      sans_tot_0 ~ "SANS total score", 
      gaf_sofas ~ "GAF-SOFAS score"
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

table_demog <- uhr1000 %>% 
  select(assessment_age, gender, timeSxService, group, uhr_group, marital_status, children, education_highest, accom_current, occup_current, Immigrated, Immigrated_parents) %>% # keep only columns of interest
  mutate(
    gender = factor(gender, labels = c("Male", "Female")) ,
    uhr_group = factor(uhr_group, labels = c("BLIPS","Attenuated Psychosis","Atten. Psychosis+BLIPS","Vulnerability","BLIPS+Vulnerability","Vulnerability+Atten. Pscychosis","BLIPS+Atten. Psychosis+Vulnerability")) ,
    marital_status = factor(marital_status, labels = c("Married/De facto", "Not Married/De facto")),
    education_highest = factor(education_highest, labels = c("primary completed, secondary ongoing","secondary completed","TAFE/diploma/certificate unfinished/finished","undergrad unfinished","undergrad finished","post-grad unfinished","post grad finished")),
    accom_current = factor(accom_current, labels = c("None","crisis accomodation","boarding/rented room/foster","rented house/flat","own house/flat","house/flat with family of origin","other")),
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

table_uhr <- uhr1000 %>% 
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
  flextable::save_as_docx(path = ".../UHR1000/uhr_table.docx")


#### Baseline table by study ####

table_bl <- uhr1000 %>% 
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
  flextable::save_as_docx(path = "U:/PostDoc/Research/Publications/Own/UHR1000/baseline_study_table.docx")


################################################################################
#################### UHR subgroup plot  ##############################
################################################################################

uhr1000.sub <- subset(uhr1000, !(is.na(uhr_cat)))
uhr1000.sub$transyears <- round(as.numeric(uhr1000.sub$transdays)/365.25, digits = 1)

fit <- survfit(Surv(transyears, transtat) ~ uhr_cat, data = uhr1000.sub)

cumhaz_plot <- ggsurvplot(fit,
                          conf.int = TRUE,
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
                          palette = c("#225ea8", "#41b6c4", "#a1dab4"),
                          fun = "event",
                          font.x = c(20),
                          font.y = c(20),
                          font.tickslab = c(15),
                          font.legend = c(15))

cumhaz_plot

################################################################################
#################### Select variables  ##############################
################################################################################

uhr1000 <- uhr1000 %>% dplyr::select(caarms_DS_sev_0, caarms_PA_sev_0, caarms_UTC_sev_0, caarms_NBI_sev_0, bprs_tot_0, sans_tot_0, gaf_sofas, assessment_age, timeSxService, uhr_cat, transtat, transdays, transmonths, group, study, nelsonaalen)

# Log transform time to service
uhr1000[uhr1000$timeSxService == 0 & !is.na(uhr1000$timeSxService), "timeSxService"] <- 1
uhr1000$timeSxService <- log(uhr1000$timeSxService)

# List missing data percentage
p_missing <- unlist(lapply(uhr1000, function(x) sum(is.na(x))))/nrow(uhr1000)
sort(p_missing[p_missing > 0], decreasing = TRUE)

# Exclude invdividuals with no follow-up time information
uhr1000 <- uhr1000[!is.na(uhr1000$transdays),]

# Save dataset prior to imputation
uhr1000.org <- uhr1000

# Limit to 2-year probability (right censor all individuals with a follow-up longer than 2 years)
uhr1000$transtat <- uhr1000$transtat * (uhr1000$transday < 2*365.25)
uhr1000[uhr1000$transdays > 2*365.25 & !is.na(uhr1000$transdays), 'transdays'] <- 2*365 

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
predM[, c("transdays")] <- 0
predM[, c("transmonths")] <- 0
predM[c("transdays"),] <- 0
predM[c("transmonths"),] <- 0
predM[, c("group")] <- 0

# Ordered categorical variables 
poly <- c("uhr_cat")

# Turn their methods matrix into the specified imputation models
meth[poly] <- "polyreg"

# Impute training set
uhr1000.imp <- mice(uhr1000, m = 50, maxit = 30, 
                  predictorMatrix = predM, 
                  method = meth, print =  FALSE, seed = 123)

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

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ caarms_DS_sev_0 + caarms_PA_sev_0 + caarms_UTC_sev_0 + caarms_NBI_sev_0 + bprs_tot_0 + 
                                 sans_tot_0 + gaf_sofas + group + timeSxService + assessment_age + uhr_cat))
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

plot.data$Predictor <- c("CAARMS Disorganized Speech, Severity",
                         "CAARMS Perceptual Abnormalities, Severity",
                         "CAARMS Unusual Thought Content, Severity",
                         "CAARMS Non-Bizarre Ideas, Severity",
                         "BPRS total score",
                         "SANS total score",
                         "GAF-SOFAS score",
                         "Received non-standard intervention treatment",
                         "Time to UHR service (log transformed)",
                         "Age at baseline",
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
text(c(-2,2.2), c(14.5, 14.5), pos=4, 
     c("Variable", "Hazard ratio (95% CI)"), font=2)

# redraw the CI lines and points in the chosen color
segments(plot.data$lower, k:1, plot.data$upper, k:1, col=colp, lwd=1.5)
points(plot.data$HR, k:1, pch=18, cex = plot.data$HR*2,  col="white")
points(plot.data$HR, k:1, cex = plot.data$HR*2, pch=18, col=colp)

#### Sensitivity analysis with study sites ####

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ caarms_DS_sev_0 + factor(study)))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ caarms_UTC_sev_0 +  factor(study)))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ caarms_PA_sev_0 +  factor(study)))
summary(pool(models), conf.int = 0.95, exponentiate = TRUE)

models <- with(uhr1000.imp,coxph(Surv(transdays, transtat) ~ caarms_NBI_sev_0 +  factor(study)))
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

#### Create model with 3 CPP ####
uhr.3 <- cph(Surv(transdays, transtat) ~ caarms_DS_sev_0 + caarms_UTC_sev_0 + sans_tot_0, data = uhr1000.single, x=TRUE, y=TRUE, surv = TRUE, time.inc = 2*365)
pred.uhr3 <- survest(uhr.3, uhr1000.single, times = 2*365)
risk.uhr3 <- 1 - pred.uhr3$surv

# Output original model
uhr.3

# Internal validation using bootstrapping
val.uhr3 <- validate(uhr.3, method="boot", B=1000, bw=FALSE)
val.uhr3

# Calibration plot
cal.uhr3 <- calibration.plot(uhr.3, 2)
cal.uhr3

# Instability plot
insta.uhr3 <- instability.plot(uhr1000.single, as.formula(Surv(transdays, transtat) ~ caarms_DS_sev_0 + caarms_UTC_sev_0 + sans_tot_0), risk.uhr3, 2, 1000)
mean(insta.uhr3)

################################################################################
#################### Temporal validation   ##############################
################################################################################

val.early <- vector()
val.middle <- vector()
val.all <- vector()

c.early <- vector()
c.middle <- vector()
c.all <- vector()

plot.early <- data.frame()
plot.middle <- data.frame()
plot.all <- data.frame()

for(imp.nr in c(1:50)){
  uhr1000.single <- subset(data.comp, .imp == imp.nr )
  uhr1000.early <- subset(uhr1000.single, study < 7 )
  uhr1000.all <- subset(uhr1000.single, study < 7 | study == 11 | study == 13 | study == 12 | study == 8  )
  uhr1000.middle <- subset(uhr1000.single, study == 5 | study == 6 | study == 11 | study == 13 | study == 12 | study == 8  )
  uhr1000.late <- subset(uhr1000.single, study == 7 | study == 10 | study == 9 )
  
  # Develop models with different baseline years
  dd <- datadist(uhr1000.early)
  options(datadist='dd')
  cox.early <- cph(Surv(transdays, transtat) ~ caarms_DS_sev_0 + caarms_UTC_sev_0 + sans_tot_0 , data = uhr1000.early, x=TRUE, y=TRUE, surv = TRUE, time.inc = 2*365)
  
  dd <- datadist(uhr1000.middle)
  options(datadist='dd')
  cox.middle <- cph(Surv(transdays, transtat) ~ caarms_DS_sev_0 + caarms_UTC_sev_0 + sans_tot_0 , data = uhr1000.middle, x=TRUE, y=TRUE, surv = TRUE, time.inc = 2*365)
  
  dd <- datadist(uhr1000.middle)
  options(datadist='dd')
  cox.all <- cph(Surv(transdays, transtat) ~ caarms_DS_sev_0 + caarms_UTC_sev_0 + sans_tot_0 , data = uhr1000.all, x=TRUE, y=TRUE, surv = TRUE, time.inc = 2*365)
  
  # Create Cox model with same predictors in validation set
  dd <- datadist(uhr1000.late)
  options(datadist='dd')
  cox.late <- cph(Surv(transdays, transtat) ~ caarms_DS_sev_0 + caarms_UTC_sev_0 + sans_tot_0 , data = uhr1000.late, x=TRUE, y=TRUE, surv = TRUE, time.inc = 2*365)
  
  #### Calculate calibration curve and c-index in validation set ####
  
  val.early[imp.nr] <- val.calibration.slope(cox.early, cox.late, uhr1000.late)
  val.middle[imp.nr] <- val.calibration.slope(cox.middle, cox.late, uhr1000.late)
  val.all[imp.nr] <- val.calibration.slope(cox.all, cox.late, uhr1000.late)

  # Calculate c-index for both models in validation set
  surv.obj <- with(uhr1000.late,Surv(transdays,transtat)) 
  
  estimates <- survest(cox.early,newdata=uhr1000.late,times=2*365)$surv
  c.early[imp.nr] <- rcorr.cens(x=estimates,S=surv.obj)
  
  estimates <- survest(cox.middle,newdata=uhr1000.late,times=2*365)$surv
  c.middle[imp.nr] <- rcorr.cens(x=estimates,S=surv.obj)
  
  estimates <- survest(cox.all,newdata=uhr1000.late,times=2*365)$surv
  c.all[imp.nr] <- rcorr.cens(x=estimates,S=surv.obj)
  
  # Calculate calibration curve
  tmp <- validation.plot(cox.early, 2, imp.nr, uhr1000.late)
  plot.early <- rbind(plot.early, tmp)
  
  tmp <- validation.plot(cox.middle, 2, imp.nr, uhr1000.late)
  plot.middle <- rbind(plot.middle, tmp)
  
  tmp <- validation.plot(cox.all, 2, imp.nr, uhr1000.late)
  plot.all <- rbind(plot.all, tmp)
  
}

plot.early$Type <- as.factor(plot.early$Type)
plot.middle$Type <- as.factor(plot.middle$Type)
plot.all$Type <- as.factor(plot.all$Type)

ggplot(plot.early, aes(x = x, y = Value, group = Type)) +
  geom_smooth(method="loess", se=FALSE, fullrange=FALSE, linewidth = 1, color = "#8c8c8c") +
  scale_x_continuous(name =paste("Predicted 2-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5), limits=c(0, 0.5)) +
  scale_y_continuous(name=paste("Observed 2-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5), limits=c(0, 0.5)) +
  geom_abline(aes(intercept = 0, slope = 1), color="black")	+
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13)) +
  theme(legend.position="bottom") 

ggplot(plot.middle, aes(x = x, y = Value, group = Type)) +
  geom_smooth(method="loess", se=FALSE, fullrange=FALSE, linewidth = 1, color =  "#8c8c8c") +
  scale_x_continuous(name =paste("Predicted 2-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5), limits=c(0, 0.5)) +
  scale_y_continuous(name=paste("Observed 2-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5), limits=c(0, 0.5)) +
  geom_abline(aes(intercept = 0, slope = 1), color="black")	+
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13)) +
  theme(legend.position="bottom") 

ggplot(plot.all, aes(x = x, y = Value, group = Type)) +
  geom_smooth(method="loess", se=FALSE, fullrange=FALSE, linewidth = 1, color =  "#8c8c8c") +
  scale_x_continuous(name =paste("Predicted 2-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5), limits=c(0, 0.5)) +
  scale_y_continuous(name=paste("Observed 2-year Probability of Transition to Psychosis",sep = ""), breaks = c(0.1,0.2,0.3,0.4,0.5), limits=c(0, 0.5)) +
  geom_abline(aes(intercept = 0, slope = 1), color="black")	+
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13)) +
  theme(legend.position="bottom") 
