# Hello! This will be the center for our analysis of VSRR data!

# Load packages
library(tidyverse)
library(naniar) #missingness visualization
library(skimr) 
library(labelled)
library(sjlabelled) #for labelling variables
library(tableone)
library(survival)
library(survminer)
library(longCatEDA)
library(rstatix) #pipe-freindly simple stats
library(PMCMRplus)
library(ggpubr)#allows p-values to be placed anywhere
library(cmprsk) #gives cumulative incidence function
library(JM) #prbc dataset example

getOption("continue")
options(continue = "  ")

### Load data -----------------------------------------------------------------

df <- valvesparingrootrepl_data_2021.06.06_2201 %>%
  dplyr::select(hx_bav, mrn, age, sex, hx_htn, hx_diabetes, hx_smoking, nyha, lvef_pre, lvidd_pre, 
         lvids_pre, diameter_annulus_pre, diameter_stj_pre, 
         diameter_sov_pre, asc_pre, ai_pre, ai_degree_pre) %>% 
  tibble()

### Data Cleaning (T1) --------------------------------------------------
  
gg_miss_var(df) #graph of missing variables
skim(df) #table showing missing n, IQR, range, mean, SD

#### HTN ================================================
df %>% 
  dplyr::select(hx_htn) %>% 
  group_by(hx_htn) %>% 
  summarise(count = n()) #seeing how many "2s" are in HTN column

help("select")

df %>% 
  filter(hx_htn == 2) #filtering to see if there are any "2s" in HTN column

df$hx_htn <- as.double(df$hx_htn) #redefining hx_htn to double

df <- df %>%
  mutate(hx_htn = if_else(hx_htn == 2, 0, hx_htn)) %>%
  print() #printing new variable

typeof(df$hx_htn) #checking on rtype of variable


barplot(table(df$hx_htn)) # barplot showing separation of htn hx pre-op

#### Diabetes ================================================

df %>% 
  dplyr::select(hx_diabetes) %>% 
  group_by(hx_diabetes) %>% 
  summarise(count = n()) #seeing how many "2s" are in diabetes column

df %>% 
  filter(hx_diabetes == 2) #filtering to see if there are any "2s" in diabetes column

df$hx_diabetes <- as.double(df$hx_diabetes) #redefining hx_diabetes to double

df <- df %>%
  mutate(hx_diabetes = if_else(hx_diabetes == 2, 0, hx_diabetes)) %>%
  print() #printing new variable

barplot(table(df$hx_diabetes)) # barplot showing separation of diabetes hx pre-op

####Smoking ================================================

df %>% 
  dplyr::select(hx_smoking) %>% 
  group_by(hx_smoking) %>% 
  summarise(count = n()) #seeing how many "2s" are in diabetes column

df %>% 
  filter(hx_smoking == 2) #filtering to see if there are any "2s" in diabetes column

df$hx_smoking <- as.double(df$hx_smoking) #redefining hx_diabetes to double

df <- df %>%
  mutate(hx_smoking = if_else(hx_smoking == 2, 0, hx_smoking)) %>%
  print() #printing new variable

barplot(table(df$hx_smoking)) # barplot showing separation of htn hx pre-op

#### Annulus ================================================

df <- df %>% 
  mutate(diameter_annulus_pre = if_else(diameter_annulus_pre<10, diameter_annulus_pre*10, diameter_annulus_pre)) 
  # mutating annulus diameter to standardize all to mm

hist(df$diameter_annulus_pre) #histogram to show change for annulus diameter

#### STJ ================================================

df <- df %>% 
  mutate(diameter_stj_pre = if_else(diameter_stj_pre<20, diameter_stj_pre*10, diameter_stj_pre))
  # mutating stj diameter to standardize all to mm

hist(df$diameter_stj_pre) #histogram to show change for stj diameter

#### SOV ================================================

df <- df %>% 
  mutate(diameter_sov_pre = if_else(diameter_sov_pre<15, diameter_sov_pre*10, diameter_sov_pre))
  # mutating sov diameter to standardize all to mm

hist(df$diameter_sov_pre) #histogram to show change for sov diameter


#### Asc. Aorta ================================================

hist(df$asc_pre) #histogram for ascending aorta length

df <- df %>% mutate(asc_pre = asc_pre/10) # creating a new variable by multiplying all ascending lengths by 10 to mm

df %>% filter(diameter_annulus_pre > 45) # filtering for the few (n < 5) patients w/ annulus diameter > 45

#### AI Degree ================================================

barplot(table(df$ai_degree_pre)) # barplot showing separation of ai degree pre-op

typeof(df$ai_degree_pre) #inspection test to see type of variable

df %>%
  dplyr::select(ai_degree_pre) %>% 
  group_by(ai_degree_pre) %>%
  summarise(count = n()) #pipe in order to investigate AI degree count

### Table 1 ------------------------------------------

df.bl_characters_all <- c("age", "sex", "hx_htn", "hx_diabetes", "hx_smoking", "nyha", "lvef_pre", "lvidd_pre", 
               "lvids_pre", "diameter_annulus_pre", "diameter_stj_pre", "diameter_sov_pre", "asc_pre", 
               "ai_pre", "ai_degree_pre") # list all baseline variables 

df.bl_characters_factors <- c("sex", "hx_htn", "hx_diabetes", "hx_smoking", "nyha", 
                              "ai_pre", "ai_degree_pre") #listing all variables that are NOT continuous as factors

df.stratified <- c("hx_bav") # creating our straifier, history of bicuspid AV

notnormal <- c("age", "lvef_pre","lvids_pre") #based on VSRR_Low_ef code, preop variables that are not normal

tab1 <- CreateTableOne(data = df, 
                       vars = df.bl_characters_all,
                       factorVars = df.bl_characters_factors,
                       strata = df.stratified,
                       includeNA=TRUE) #creating tab1 in order to then print it


print(tab1, 
      nonnormal = notnormal,
      varLabels = TRUE, #Chi-square for catagorical, one-way anova for normal cont, kruskal for nonnormal cont.
      dropEqual = TRUE,
      smd = TRUE) #does not show the outputed level (ex. COPD = Yes), smd needs to be in the print!!

print(tab1, 
      nonnormal = notnormal,
      varLabels = TRUE,
      dropEqual = TRUE,
      quote = TRUE, noSpaces = TRUE,
      smd = TRUE) #for export to Excel.

# Table 2 ------------------------------------------------------------

df2 <- valvesparingrootrepl_data_2021.06.06_2201 %>%
  dplyr::select(cpb, xc, circarrest, 
circarrest_time, archgraft, procasc, prochemi, procz1, procz2, proctotarch, 
procmitral, proccabg, leaflet_repair, cusprepair_reference___5, 
cusprepair_reference___2, cusprepair_reference___1, cusprepair_reference___3, 
cusprepair_reference___4, cusprepair_reference___6) %>%
  tibble()

skim(df2)

df.bl_characters_all <- c("cpb", "xc", "circarrest", 
                          "circarrest_time", "archgraft", "procasc", "prochemi", "procz1", "procz2", "proctotarch", 
                          "procmitral", "proccabg", "leaflet_repair", "cusprepair_reference___5", 
                          "cusprepair_reference___2", "cusprepair_reference___1", "cusprepair_reference___3", 
                          "cusprepair_reference___4", "cusprepair_reference___6") # list all intraop variables 

df.bl_characters_factors <- c( "circarrest", "procasc", "prochemi", "procz1", "procz2", "proctotarch", 
                               "procmitral", "proccabg", "leaflet_repair", "cusprepair_reference___5", 
                               "cusprepair_reference___2", "cusprepair_reference___1", "cusprepair_reference___3", 
                               "cusprepair_reference___4", "cusprepair_reference___6") #listing all variables that are NOT continuous as factors

df.stratified <- c("hx_bav") # creating our straifier, history of bicuspid AV

notnormal <- c("cpb", "xc", "circarrest_time") #based on VSRR_Low_ef code, preop variables that are not normal

tab2 <- CreateTableOne(data = valvesparingrootrepl_data_2021.06.06_2201, 
                       vars = df.bl_characters_all,
                       factorVars = df.bl_characters_factors,
                       strata = df.stratified,
                       includeNA=TRUE) #creating tab1 in order to then print it


print(tab2, 
      nonnormal = notnormal,
      varLabels = TRUE, #Chi-square for catagorical, one-way anova for normal cont, kruskal for nonnormal cont.
      dropEqual = TRUE,
      smd = TRUE) #does not show the outputed level (ex. COPD = Yes), smd needs to be in the print!!

print(tab2, 
      nonnormal = notnormal,
      varLabels = TRUE,
      dropEqual = TRUE,
      quote = TRUE, noSpaces = TRUE,
      smd = TRUE) #for export to Excel.



# Table 3 ------------------------------------------------------------

df3 <- valvesparingrootrepl_data_2021.06.06_2201 %>%
  dplyr::select(mrn, hx_bav, mort_all, reop, pg_dc, mg_dc, coaptheight_post, ai_dc, 
                ai_degree_dc, postop_neuro, reop_pacer, echotype_1yr,
                lvids_1yr, lvidd_1yr, lvef_1yr, pg_1yr, mg_1yr) %>% 
  tibble() 

skim(df3$lvidd_1yr) #table showing missing n, IQR, range, mean, SD


class(df3$lvids_1yr)
class(df3$lvidd_1yr)

df3$lvidd_1yr <- as.double(df3$lvidd_1yr)

df3 <- df3 %>%
  mutate(lvidd_1yr = if_else(lvidd_1yr == "4,8", 4.8, lvidd_1yr))

df3 <- df3 %>%
  mutate(echotype_1yr = if_else(echotype_1yr == 2, 0, echotype_1yr)) %>%
  print() #printing new variable





df3_characters_all <- c("mort_all", "reop", "pg_dc", "mg_dc", "coaptheight_post", "ai_dc", 
                          "ai_degree_dc", "postop_neurotype_2", "reop_pacer", "echotype_1yr",
                          "lvids_1yr", "lvidd_1yr", "lvef_1yr", "pg_1yr", "mg_1yr") # list all postop variables 

df3_characters_factors <- c("mort_all", "reop", "ai_dc", 
                              "ai_degree_dc", "postop_neurotype_2", "reop_pacer", 
                              "echotype_1yr") #listing all variables that are NOT continuous as factors

df3.stratified <- c("hx_bav") # creating our straifier, history of bicuspid AV

notnormal <- c("pg_dc", "mg_dc", "coaptheight_post", "lvids_1yr", "lvidd_1yr", "lvef_1yr", "pg_1yr", "mg_1yr") #based on VSRR_Low_ef code, preop variables that are not normal

tab3 <- CreateTableOne(data = df3, 
                       vars = df3_characters_all,
                       factorVars = df3_characters_factors,
                       strata = df3.stratified,
                       includeNA=TRUE) #creating tab1 in order to then print it


print(tab3, 
      nonnormal = notnormal,
      varLabels = TRUE, #Chi-square for catagorical, one-way anova for normal cont, kruskal for nonnormal cont.
      dropEqual = TRUE,
      smd = TRUE) #does not show the outputed level (ex. COPD = Yes), smd needs to be in the print!!

print(tab3, 
      nonnormal = notnormal,
      varLabels = TRUE,
      dropEqual = TRUE,
      quote = TRUE, noSpaces = TRUE,
      smd = TRUE) #for export to Excel.

### Figures -------------------------------------------------

#### Yale Graph --------------------------------------------- 

df.surv$dos <- as.Date.character(df.surv$dos, "%Y-%m-%d")

df.surv <- df.surv %>%
  mutate(potential_follow_up = (Sys.Date() - dos)/365.25)

typeof(Sys.Date())
class(Sys.Date())
typeof(df.surv$dos)
class(df.surv$dos)

typeof(follow_up_days)
        
ggplot(data = df.surv, mapping = aes(x = valvesparingrootrepl_data_2021.06.06_2201$followup_days/365.25, y = potential_follow_up*-1)) + geom_point() + geom_abline(intercept = 0.3, slope = -1)

ggplot(df.surv, aes(x = potential_follow_up, y = time_to_event_years, color = as.factor(mort_all))) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1) + 
  geom_abline(intercept = -1, slope = 1) + 
  geom_abline(intercept = -3, slope = 1) + 
  geom_abline(intercept = -5, slope = 1) 


#### KM Curves ----------------------------------------------

##### Mortality --------------------
df.surv <- valvesparingrootrepl_data_2021.06.06_2201 %>%
  dplyr::select(dos, mort_date, mort_all, alive_date, hx_bav, age, reop, reop_date,
                reop_aortic___1, reop_aortic___2, reop_aortic___3)
    #creating the survival variable

class(df.surv$mort_date)

df.surv$mort_date <- as.Date.character(df.surv$mort_date, "%Y-%m-%d") 
df.surv$dos <- as.Date.character(df.surv$dos, "%Y-%m-%d")
df.surv$alive_date <- as.Date.character(df.surv$alive_date, "%Y-%m-%d")
#redefining variables for analysis (converting characters to date)

df.surv <- df.surv %>%
  mutate(event_date = if_else(is.na(mort_date),alive_date,mort_date),
         time_to_event_days = (event_date - dos),
         time_to_event_years = time_to_event_days/365.25, 
         event = if_else(mort_all == "1", 1, 0, missing = NULL))
    #mutating variables

df.surv %>%
  group_by(event)%>%
  skim(time_to_event_years)

surv.all <- Surv(time = df.surv$time_to_event_years, event = df.surv$event)
surv.all.fit <- survfit(surv.all ~ 1)
print(surv.all.fit, print.rmean=TRUE)
summary(surv.all.fit)
plot(surv.all.fit, fun = "s") #overall 95% CI

surv.strat <- survfit(surv.all ~ df.surv$hx_bav)
print(surv.strat, print.rmean = TRUE)
summary(surv.strat,c(0,1,2,3,4,5,6,7,8, 9, 10))
plot(surv.strat, fun = "s", xlab = "Time (years)", ylab = "Survival",
     main = "Survival Curve")

##formal plot switch to this

ggsurvplot(surv.strat, data = df.surv,
           conf.int = TRUE,
           xlab = "Time in Years",risk.table = TRUE, break.time.by =1, xlim = c(0, 10), pval = TRUE, pval.method = TRUE,
           legend = c(.95, .5), legend.labs = c("TAV", "BAV"), title = "Freedom from Mortality")


help("ggsurvplot")

           break.time.by = 1,
           xlim = c(0, 12),  #change this to adjust x-axis
           legend.labs = c("≤45","46-55",">55"),
           #pval = "p = 0.30 (log-rank test)",
           pval = TRUE,
           pval.method = TRUE,
           risk.table = TRUE,
           risk.table.height = 0.25,
           palette = c("#c6dbef", "#6baed6","#08519c"))

##### Reoperation --------------------

class(df.surv$reop_date)

df.surv$reop_date <- as.Date.character(df.surv$reop_date, "%Y-%m-%d")

df.surv <- df.surv %>%
  mutate(reop_date_KM = if_else(is.na(reop_date),alive_date, reop_date),
         time_to_reop_days = (reop_date_KM - dos),
         time_to_reop_years = time_to_event_days/365.25, 
         event_reop = if_else(reop_aortic___1 | reop_aortic___2 | reop_aortic___3 == "1", 1, 0, missing = NULL))

surv.all_reop <- Surv(time = df.surv$time_to_reop_years, event = df.surv$reop)
surv.all.fit_reop <- survfit(surv.all_reop ~ 1)
print(surv.all.fit_reop, print.rmean=TRUE)
summary(surv.all.fit_reop)
plot(surv.all.fit_reop, fun = "s")

surv.strat_reop <- survfit(surv.all_reop ~ df.surv$hx_bav)
print(surv.strat_reop, print.rmean = TRUE)
summary(surv.strat_reop,c(0,1,2,3,4,5,6,7, 8, 9, 10, 11))
plot(surv.strat_reop, fun = "s", xlab = "Time (years)", ylab = "Freedom from Reoperation",
     main = "Reoperation Curve")

ggsurvplot(surv.strat_reop, data = df.surv,
           conf.int = TRUE,
           xlab = "Time in Years",risk.table = TRUE, break.time.by =1, xlim = c(0, 10), pval = TRUE, pval.method = TRUE,
           legend = c(.95, .5), legend.labs = c("TAV", "BAV"), title = "Freedom from Aortic Reoperation")

## specifically code for aortic reop!!!!



#### Freedom from AI
df.freedom_ai <- valvesparingrootrepl_data_2021.06.06_2201 %>%
  dplyr::select(mrn, dos, mort_date, mort_all, alive_date, hx_bav, age, reop, reop_date, 
                dc_echodate, ai_dc, ai_degree_dc, echodate_30day, ai_30day,	ai_degree_30day, 
                echodate_1yr, ai_1yr, ai_degree_1yr, echotype_1yr,
                echodate_2yr, ai_2yr, ai_degree_2yr, echotype_2yr,
                echodate_3yr, ai_3yr, ai_degree_3yr, echotype_3yr,
                echodate_4yr, ai_4yr, ai_degree_4yr, echotype_4yr,
                echodate_5yr, ai_5yr, ai_degree_5yr, echotype_5yr,
                echodate_6yr, ai_6yr, ai_degree_6yr, echotype_6yr,
                echodate_7yr, ai_7yr, ai_degree_7yr, echotype_7yr,
                echodate_8yr, ai_8yr, ai_degree_8yr, echotype_8yr,
                echodate_9yr, ai_9yr, ai_degree_9yr, echotype_9yr,
                echodate_10yr, ai_10yr, ai_degree_10yr, echotype_10yr,
                echodate_11yr, ai_11yr, ai_degree_11yr,
                echodate_12yr, ai_12yr, ai_degree_12yr)

typeof(df$dc_echodate)
class(df.surv$echodate_7yr)

df.freedom_ai$dc_echodate <- as.Date.character(df.freedom_ai$dc_echodate, "%Y-%m-%d")
df.freedom_ai$echodate_30day <- as.Date.character(df.freedom_ai$echodate_30day, "%Y-%m-%d")
df.freedom_ai$echodate_1yr <- as.Date.character(df.freedom_ai$echodate_1yr, "%Y-%m-%d")
df.freedom_ai$echodate_2yr <- as.Date.character(df.freedom_ai$echodate_2yr, "%Y-%m-%d")
df.freedom_ai$echodate_3yr <- as.Date.character(df.freedom_ai$echodate_3yr, "%Y-%m-%d")
df.freedom_ai$echodate_4yr <- as.Date.character(df.freedom_ai$echodate_4yr, "%Y-%m-%d")
df.freedom_ai$echodate_5yr <- as.Date.character(df.freedom_ai$echodate_5yr, "%Y-%m-%d")
df.freedom_ai$echodate_6yr <- as.Date.character(df.freedom_ai$echodate_6yr, "%Y-%m-%d")
df.freedom_ai$echodate_7yr <- as.Date.character(df.freedom_ai$echodate_7yr, "%Y-%m-%d")
df.freedom_ai$echodate_8yr <- as.Date.character(df.freedom_ai$echodate_8yr, "%Y-%m-%d")
df.freedom_ai$echodate_9yr <- as.Date.character(df.freedom_ai$echodate_9yr, "%Y-%m-%d")
df.freedom_ai$echodate_10yr <- as.Date.character(df.freedom_ai$echodate_10yr, "%Y-%m-%d")
df.freedom_ai$echodate_11yr <- as.Date.character(df.freedom_ai$echodate_11yr, "%Y-%m-%d")
df.freedom_ai$echodate_12yr <- as.Date.character(df.freedom_ai$echodate_12yr, "%Y-%m-%d")
df.freedom_ai$dos <- as.Date.character(df.freedom_ai$dos, "%Y-%m-%d")
df.freedom_ai$alive_date <- as.Date.character(df.freedom_ai$alive_date, "%Y-%m-%d")

df.freedom_ai <- df.freedom_ai %>%
  mutate(late_ai_new_days = if_else(ai_degree_10yr > 0, echodate_10yr - dos, alive_date - dos),
         late_ai_new_days = if_else(ai_degree_9yr > 0, echodate_9yr - dos, alive_date - dos),
         late_ai_new_days = if_else(ai_degree_8yr > 0, echodate_8yr - dos, alive_date - dos),
         late_ai_new_days = if_else(ai_degree_7yr > 0, echodate_7yr - dos, alive_date - dos),
         late_ai_new_days = if_else(ai_degree_6yr > 0, echodate_6yr - dos, alive_date - dos),
         late_ai_new_days = if_else(ai_degree_5yr > 0, echodate_5yr - dos, alive_date - dos),
         late_ai_new_days = if_else(ai_degree_4yr > 0, echodate_4yr - dos, alive_date - dos),
         late_ai_new_days = if_else(ai_degree_3yr > 0, echodate_3yr - dos,alive_date - dos),
         late_ai_new_days = if_else(ai_degree_2yr > 0, echodate_2yr - dos, alive_date - dos),
         late_ai_new_days = if_else(ai_degree_1yr > 0, echodate_1yr - dos, alive_date - dos),
         late_ai_new_days = if_else(ai_degree_30day > 0, echodate_30day - dos, alive_date - dos),
         late_ai_new_days = if_else(ai_degree_dc > 0, dc_echodate - dos, alive_date - dos),
          event_new_ai = if_else(late_ai_new_days != alive_date - dos, 1, 0), 
         late_ai_new_years = late_ai_new_days/365.25)

view(df.freedom_ai$late_ai_new_days)

surv.all_freedom_ai <- Surv(time = df.freedom_ai$late_ai_new_years, event = df.freedom_ai$event_new_ai)
surv.all.fit_freedom_ai <- survfit(surv.all_freedom_ai ~ 1)
print(surv.all.fit_freedom_ai, print.rmean=TRUE)
summary(surv.all.fit_freedom_ai)
plot(surv.all.fit_freedom_ai, fun = "s")




surv.strat_freedom_ai <- survfit(surv.all_freedom_ai ~ df.surv$hx_bav)
print(surv.strat_freedom_ai, print.rmean = TRUE)
summary(surv.strat_freedom_ai,c(0,1,2,3,4,5,6,7, 8))
plot(surv.strat_freedom_ai, xlab = "Time (years)", ylab = "Freedom from AI",
     main = "Freedom from AI", lty = 1:2)

ggsurvplot(surv.strat_freedom_ai, data = df.surv,
           conf.int = TRUE,
           xlab = "Time in Years",risk.table = TRUE, break.time.by =1, xlim = c(0, 10), pval = TRUE, pval.method = TRUE,
           legend = c(.95, .5), legend.labs = c("TAV", "BAV"), title = "Freedom from AI >0")


help(legend)

#### Bar Graph w/ AI Distribution ---------------------------

df.distribution <- df.freedom_ai %>%
  pivot_longer(c(ai_degree_dc, ai_degree_30day, ai_degree_1yr,
                ai_degree_2yr, ai_degree_3yr, ai_degree_4yr,
                ai_degree_5yr, ai_degree_6yr, ai_degree_7yr, 
                ai_degree_8yr, ai_degree_9yr, ai_degree_10yr,
                ai_degree_11yr, ai_degree_12yr), names_to = "AI_Follow_up",
                values_to = "Ai_degree") # reformat 

df.distribution <- df.freedom_ai %>%
  pivot_longer(c(ai_dc, ai_30day, ai_1yr,
                 ai_2yr, ai_3yr, ai_4yr,
                 ai_5yr, ai_6yr, ai_7yr, 
                 ai_8yr, ai_9yr, ai_10yr,
                 ai_11yr, ai_12yr), names_to = "AI_Follow_up_time",
               values_to = "Ai_yes_no")

df.distribution %>%
  dplyr::select(Ai_degree) %>%
  group_by(Ai_degree) %>%
  summarise(count = n())

df.distribution2 %>%
  dplyr::select(Ai_yes_no) %>%
  group_by(Ai_yes_no) %>%
  summarise(count = n())

df.distribution %>%
  dplyr::select(a) %>%
  group_by(a) %>%
  summarise(count = n())


IF ai_at_dc_any==0 & ai_dc_level = na THEN ai_dc_level == 0
Mutate(a = ai_discharge * 1
       a = if_else(a == 1, degree_discharge, a)
       
df.distribution <- df.distribution %>%
  mutate(a= Ai_degree*1,
         a= as.double(a),
         a= if_else(df.distribution2$Ai_yes_no==1 & Ai_degree==0, 0, a))

AI_follow_up_ordered <- factor(df.distribution$AI_Follow_up, levels=c("ai_degree_dc","ai_degree_30day", "ai_degree_1yr",
                        "ai_degree_2yr", "ai_degree_3yr", "ai_degree_4yr", "ai_degree_5yr", "ai_degree_6yr", "ai_degree_7yr",
                        "ai_degree_8yr", "ai_degree_9yr", "ai_degree_10yr", "ai_degree_11yr",
                        "ai_degree_12yr"), 
                        labels=c("DC","30 Day", "1 Yr", "2 Yr", "3 Yr",
                        "4 Yr", "5 Yr", "6 Yr", "7 Yr", "8 Yr", "9 Yr",
                        "10 Yr", "11 Yr","12 Yr")

df.distribution %>%
  dplyr::filter(!is.na(Ai_degree))%>% 
  ggplot(data = df.distribution, mapping = aes(fill = as.factor(a), x=AI_follow_up_ordered, y=a)) + 
  geom_bar(position = "stack", stat = "identity")
  
df.distribution$Ai_degree <- as.factor(df.distribution$Ai_degree)



help(is.na)

#### Bar Graph (redo)

df.freedom_ai$ai_degree_dc <- as.factor(df.freedom_ai$ai_degree_dc)
df.freedom_ai$ai_degree_30day <- as.factor(df.freedom_ai$ai_degree_30day)
df.freedom_ai$ai_degree_1yr <- as.factor(df.freedom_ai$ai_degree_1yr)
df.freedom_ai$ai_degree_2yr <- as.factor(df.freedom_ai$ai_degree_2yr)
df.freedom_ai$ai_degree_3yr <- as.factor(df.freedom_ai$ai_degree_3yr)

df.freedom_ai <- df.freedom_ai %>%
  mutate(ai_degree_dc = fct_recode(ai_degree_dc),
    ai_degree_dc = if_else(ai_dc==0 & ai_degree_dc!=0|1|2|3|4,0, fct_recode(ai_degree_dc)),
    ai_degree_30day = fct_recode(ai_degree_30day),
    ai_degree_30day = if_else(ai_30day==0 & ai_degree_30day!=0|1|2|3|4,0, ai_degree_30day),
    ai_degree_1yr = fct_recode(ai_degree_1yr),
    ai_degree_1yr = if_else(ai_1yr==0 & ai_degree_1yr!=0|1|2|3|4,0, ai_degree_1yr),
    ai_degree_2yr = fct_recode(ai_degree_2yr),
    ai_degree_2yr = if_else(ai_2yr==0 & ai_degree_2yr!=0|1|2|3|4,0, ai_degree_2yr),
    ai_degree_3yr = fct_recode(ai_degree_3yr),
    ai_degree_3yr = if_else(ai_3yr==0 & ai_degree_3yr!=0|1|2|3|4,0, ai_degree_3yr),
    ai_degree_4yr = fct_recode(ai_degree_4yr),
    ai_degree_4yr = if_else(ai_4yr==0 & ai_degree_4yr!=0|1|2|3|4,0, ai_degree_4yr),
    ai_degree_5yr = fct_recode(ai_degree_5yr),
    ai_degree_5yr = if_else(ai_5yr==0 & ai_degree_5yr!=0|1|2|3|4,0, ai_degree_5yr),
    ai_degree_6yr = fct_recode(ai_degree_6yr),
    ai_degree_6yr = if_else(ai_6yr==0 & ai_degree_6yr!=0|1|2|3|4,0, ai_degree_6yr),
    ai_degree_7yr = fct_recode(ai_degree_7yr),
    ai_degree_7yr = if_else(ai_7yr==0 & ai_degree_7yr!=0|1|2|3|4,0, ai_degree_7yr),
    ai_degree_8yr = fct_recode(ai_degree_8yr),
    ai_degree_8yr = if_else(ai_8yr==0 & ai_degree_8yr!=0|1|2|3|4,0, ai_degree_8yr),
    ai_degree_9yr = fct_recode(ai_degree_9yr),
    ai_degree_9yr = if_else(ai_9yr==0 & ai_degree_9yr!=0|1|2|3|4,0, ai_degree_9yr),
    ai_degree_10yr = fct_recode(ai_degree_10yr),
    ai_degree_10yr = if_else(ai_10yr==0 & ai_degree_10yr!=0|1|2|3|4,0, ai_degree_10yr),
    ai_degree_11yr = fct_recode(ai_degree_11yr),
    ai_degree_11yr = if_else(ai_11yr==0 & ai_degree_11yr!=0|1|2|3|4,0, ai_degree_11yr),
    ai_degree_12yr = fct_recode(ai_degree_12yr),
    ai_degree_12yr = if_else(ai_12yr==0 & ai_degree_12yr!=0|1|2|3|4,0, ai_degree_12yr)
    )

df.new <- df.freedom_ai %>%
  pivot_longer(c(ai_degree_dc, ai_degree_30day, ai_degree_1yr,
                 ai_degree_2yr, ai_degree_3yr, ai_degree_4yr,
                 ai_degree_5yr, ai_degree_6yr, ai_degree_7yr, 
                 ai_degree_8yr, ai_degree_9yr, ai_degree_10yr,
                 ai_degree_11yr, ai_degree_12yr), names_to = "AI_Follow_up",
               values_to = "Ai_degree") # reformat 

df.new$Ai_degree <- as.factor(df.new$Ai_degree)    

df.freedom_ai %>%
  dplyr::select(ai_degree_dc) %>%
  group_by(ai_degree_dc) %>%
  summarise(count = n())

df.new$Ai_degree <- as.factor(df.new$Ai_degree)

df.new %>%
  dplyr::filter(!is.na(Ai_degree))%>% 
  ggplot(data = df.new, mapping = aes(fill = as.factor(Ai_degree), x=AI_follow_up_ordered, y=Ai_degree)) + 
  geom_bar(position = "stack", stat = "identity")


