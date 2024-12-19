####################################################################################
## Part 4. Advanced topics
## Chpt 7. Structural equation modeling
####################################################################################

# Load libraries
library(tidyverse)     # Data manipulation 
library(survey)        # Complex survey analysis: Base-R approach
library(srvyr)         # Complex survey analysis: Tidyverse approach 
library(lavaan)        # Structural equation modeling (SEM)
library(lavaan.survey) # Complex survey analysis of SEMs 
library(mice)          # Multiple imputation (MI)
library(naniar)        # Summary statistics of missing data

# Import data
# setwd("~/CSA_Code")
dat <- haven::read_spss("kyrbs2020.sav")
dat

# Lower cases for convenience
names(dat) <- tolower(names(dat))

# Preprocess data
mydata <- dat %>% 
  mutate(
    female=ifelse(sex==1,0,1) %>% as.factor(), # male/female 
    agem=age_m/12, # age in years
    ses=e_ses, # socioeconomic status
    sch_type=factor(stype,labels=c("Co-ed school","All-boys school","All-girls school")), # type of schools
    sch_mdhg=factor(mh,labels=c("High school","Middle school")), # middle/high school
    achieve=e_s_rcrd, # grade 
    use_smart=ifelse(int_spwd==1&int_spwk==1,0,1) %>% as.factor(), # use smartphone(1) or not(0)
    int_spwd_tm=ifelse(is.na(int_spwd_tm),0,int_spwd_tm), # use smartphone on weekdays(1) or not(0)
    int_spwk=ifelse(is.na(int_spwk),0,int_spwk), # use smartphone on weekend(1) or not(0)
    spend_smart=(5*int_spwd_tm+2*int_spwk)/(7*60), # daily smartphone usage (hrs) 
    spend_smart=ifelse(use_smart==0,NA,spend_smart), # NA if not using smartphone 
    scale_gad=5-rowMeans(dat %>% select(starts_with("m_gad"))), # GAD (Generalized Anxiety Disorder)
    scale_spaddict=5-rowMeans(dat %>% select(starts_with("int_sp_ou"))), # smartphone addiction
    scale_spaddict=ifelse(use_smart==0,NA,scale_spaddict), # NA if not using smartphone 
    covid_suffer=e_covid19, # financial hardship due to covid 19
    exp_sad=ifelse(m_sad==2,1,0)%>% as.factor(), # experienced sadness/despair
    exercise=pa_tot-1, # weekly workout days 
    happy=6-pr_hd, # perceived happiness (reversed; 1[NOT]-5[HAPPY])
  ) %>% 
  select(female:happy,
         starts_with("m_gad"),
         starts_with("int_sp_ou"),
         cluster,strata,w,fpc)

# Complex survey analysis: Tidyverse approach　
cs_design <- mydata %>% 
  as_survey_design(ids=cluster,     # cluster 
                   strata=strata,   # strata
                   weights=w,       # weight 
                   fpc=fpc)         # FPC 

# Choose smartphone users only 
subpop <- mydata %>% filter(use_smart==1)
cs_subpop <- cs_design %>% filter(use_smart==1)

###########################################################
# Confirmatory Factor Analysis (CFA)
# Step 1: Set model formula 
mycfa <- "
# Lambda matrix
GAD =~ m_gad_1+m_gad_2+m_gad_3+m_gad_4+m_gad_5+m_gad_6+m_gad_7
ADDICT =~ int_sp_ou_1+int_sp_ou_2+int_sp_ou_3+int_sp_ou_4+int_sp_ou_5
ADDICT =~ int_sp_ou_6+int_sp_ou_7+int_sp_ou_8+int_sp_ou_9+int_sp_ou_10
# Phi matrix 
GAD ~~ ADDICT 
# Theta matrix (Based on modification indices)
int_sp_ou_1 ~~  int_sp_ou_2  
int_sp_ou_2 ~~  int_sp_ou_3
int_sp_ou_1 ~~ int_sp_ou_3
int_sp_ou_5 ~~ int_sp_ou_6
int_sp_ou_7 ~~ int_sp_ou_9
int_sp_ou_7 ~~ int_sp_ou_8
int_sp_ou_8 ~~ int_sp_ou_9
m_gad_4 ~~ m_gad_5
m_gad_5 ~~ m_gad_7
"
# Step 2: Estimation 
fit_cfa_n <- sem(mycfa, data=subpop, estimator="MLM")
# modindices(fit_cfa_n) %>% arrange(desc(mi)) %>% head()
# Step 3: Summary
summary(fit_cfa_n,fit=TRUE,standard=TRUE)  # model fits, standardized values

# Check out all possible indices 
fitMeasures(fit_cfa_n)

# Step 2.5: Survey-weighted estimation
fit_cfa_y <- lavaan.survey(fit_cfa_n,cs_subpop)
# Step 3: Summary
summary(fit_cfa_y,fit=TRUE,standard=TRUE)

# mycfa <- "
# # Lambda matrix
# GAD =~ m_gad_1+m_gad_2+m_gad_3+m_gad_4+m_gad_5+m_gad_6+m_gad_7
# ADDICT =~ int_sp_ou_1+int_sp_ou_2+int_sp_ou_3+int_sp_ou_4+int_sp_ou_5
# ADDICT =~ int_sp_ou_6+int_sp_ou_7+int_sp_ou_8+int_sp_ou_9+int_sp_ou_10
# # Phi matrix 
# GAD ~~ ADDICT 
# # Intercepts 
# int_sp_ou_1~1; int_sp_ou_2~1; int_sp_ou_3~1; int_sp_ou_4~1; int_sp_ou_5~1
# int_sp_ou_6~1; int_sp_ou_7~1; int_sp_ou_8~1; int_sp_ou_9~1; int_sp_ou_10~1
# # Theta matrix (Based on modification indices) 
# int_sp_ou_1 ~~  int_sp_ou_2  
# int_sp_ou_2 ~~  int_sp_ou_3
# int_sp_ou_1 ~~ int_sp_ou_3
# int_sp_ou_5 ~~ int_sp_ou_6
# int_sp_ou_7 ~~ int_sp_ou_9
# int_sp_ou_7 ~~ int_sp_ou_8
# int_sp_ou_8 ~~ int_sp_ou_9
# m_gad_4 ~~ m_gad_5
# m_gad_5 ~~ m_gad_7
# "

# Footnote: Calculate model fits
tibble(
  source=names(fitMeasures(fit_cfa_n)),
  cs_n=fitMeasures(fit_cfa_n),
  cs_y=fitMeasures(fit_cfa_y)
) %>% 
  filter(source=="chisq.scaled"|source=="df"|   
           source=="cfi"|source=="tli"|
           source=="rmsea"|source=="rmsea.ci.lower"|
           source=="rmsea.ci.upper"|source=="srmr") %>% 
  mutate(across(
    .cols=starts_with("cs_"),
    .fns=function(x){round(x,5)}
  ))

# Compare the results
cfa_compare <- inner_join(
  parameterEstimates(fit_cfa_n) %>% select(lhs,rhs,op,est,se,z),
  parameterEstimates(fit_cfa_y) %>% select(lhs,rhs,op,est,se,z),
  by=c("lhs","op","rhs")) %>% 
  drop_na(z.x) %>% 
  select(lhs,op,rhs,starts_with("est"),starts_with("se"),starts_with("z"))
cfa_compare

# Tidy the model outpus 
cfa_compare <- cfa_compare %>% 
  mutate(
    ratioF=est.y/est.x,
    conserv=ifelse(abs(z.x)>abs(z.y),TRUE,FALSE)
  )
quantile(cfa_compare$ratioF,prob=c(0,0.25,0.5,0.75,1))
count(cfa_compare, conserv)

# Path Modeling (PM)
# Step 1: Set model formula 
mypm <- "
# Control variables 
covid_suffer ~ achieve
happy ~ achieve
# X -> M 
covid_suffer ~ a*ses 
# X -> Y ; M -> Y 
happy ~ b*covid_suffer+c*ses 
# Intercept 
covid_suffer~1 
happy~1 
# Estimate IE/DE/TE
IE := a*b   # indirect effect
DE := c     # direct effect
TE := DE+IE # total effect
"
# Step 2: Estimation 
fit_pm_n <- sem(mypm,data=subpop,estimator="MLM")
# Step 3: Summary
summary(fit_pm_n) 

# Step 2.5: Survey-weighted estimation
fit_pm_y <- lavaan.survey(fit_pm_n,cs_subpop)
# Step 3: Summary
summary(fit_pm_y)

# User-defined function for model summary 
wrap_sem_function <- function(fit_n,fit_y){
  inner_join(parameterestimates(fit_n) %>% 
               select(lhs,rhs,op,est,se,z,pvalue),
             parameterestimates(fit_y) %>% 
               select(lhs,rhs,op,est,se,z,pvalue),
             by=c("lhs","op","rhs")) %>% 
    select(lhs,op,rhs,starts_with("est"),starts_with("se"),
           starts_with("z"),starts_with("pvalue")) %>% 
    drop_na(z.x) %>% 
    mutate(
      myreportx=str_c(format(round(est.x,4),nsmall=4),
                      cut(pvalue.x,c(-Inf,0.001,0.01,0.5,1),labels=c("***","**","**","")),
                      "\n(",
                      format(round(se.x,4),nsmall=4),
                      ")"),
      myreporty=str_c(format(round(est.y,4),nsmall=4),
                      cut(pvalue.y,c(-Inf,0.001,0.01,0.5,1),labels=c("***","**","**","")),
                      "\n(",
                      format(round(se.y,4),nsmall=4),
                      ")"),
      R2x=ifelse(op=="~~",format(round(1-est.x,4),nsmall=4),""),
      R2y=ifelse(op=="~~",format(round(1-est.y,4),nsmall=4),""),
      ratioF=est.y/est.x,
      conserv=ifelse(abs(z.x)>abs(z.y),TRUE,FALSE)
    ) %>% 
    select(lhs,op,rhs,myreportx,myreporty,R2x,R2y,ratioF,conserv) 
}
pm_compare <- wrap_sem_function(fit_pm_n,fit_pm_y)
write_excel_csv(pm_compare, "Table_Part4_Ch7_2_compare_PM_CS.csv")

# Structural Equation Modeling (SEM)
# Step 1: Set model formula 
mysem <- "
GAD =~ m_gad_1+m_gad_2+m_gad_3+m_gad_4+m_gad_5+m_gad_6+m_gad_7
ADDICT =~ int_sp_ou_1+int_sp_ou_2+int_sp_ou_3+int_sp_ou_4+int_sp_ou_5
ADDICT =~ int_sp_ou_6+int_sp_ou_7+int_sp_ou_8+int_sp_ou_9+int_sp_ou_10
# Theta matrix (Based on modification indices)
int_sp_ou_1 ~~  int_sp_ou_2  
int_sp_ou_2 ~~  int_sp_ou_3
int_sp_ou_1 ~~ int_sp_ou_3
int_sp_ou_5 ~~ int_sp_ou_6
int_sp_ou_7 ~~ int_sp_ou_9
int_sp_ou_7 ~~ int_sp_ou_8
int_sp_ou_8 ~~ int_sp_ou_9
m_gad_4 ~~ m_gad_5
m_gad_5 ~~ m_gad_7
# X -> Y1, Y2 
ADDICT ~ GAD+covid_suffer
spend_smart ~ achieve+covid_suffer
# Y1 <--> Y2  
ADDICT ~ p1*spend_smart
spend_smart ~ p2*ADDICT
# Test reciprocal relationships 
mypara := p1-p2 
"
# Step 2: Estimation 
fit_sem_n <- sem(mysem,data=mydata,estimator="MLM")
# Step 3: Summary
summary(fit_sem_n,fit=TRUE,standard=TRUE)

# Step 2.5: Survey-weighted estimation
fit_sem_y <- lavaan.survey(fit_sem_n,survey.design=cs_design)
# Step 3: Summary
summary(fit_sem_y,fit=TRUE,standard=TRUE)

# Footnote: Calculate model fits
tibble(
  source=names(fitMeasures(fit_sem_n)),
  cs_n=fitMeasures(fit_sem_n),
  cs_y=fitMeasures(fit_sem_y)
) %>% 
  filter(source=="chisq.scaled"|source=="df"|   
           source=="cfi"|source=="tli"|
           source=="rmsea"|source=="rmsea.ci.lower"|
           source=="rmsea.ci.upper"|source=="srmr") %>% 
  mutate(across(
    .cols=starts_with("cs_"),
    .fns=function(x){round(x,5)}
  ))

# Tidy the model outputs
sem_compare <- wrap_sem_function(fit_sem_n,fit_sem_y)
quantile(sem_compare$ratioF,c(0,0.25,0.50,0.75,1.00))
count(sem_compare, conserv)
sem_compare %>%
  filter(op=="~"|op==":=") %>% 
  write_excel_csv("Table_Part4_Ch7_4_compare_SEM_CS.csv")

###########################################################
# Missing data analysis
n_case_miss(mydata); n_case_complete(mydata)
miss_var_summary(mydata) # systematic missingness (use_smart==0)

# Choose smartphone users only
subpop <- mydata %>% filter(use_smart==1) # same with CFA's
miss_var_summary(subpop) # agem has missing values

# Step 1: Set model formula 
mysem2 <- str_replace(mysem,
                      "spend_smart ~ achieve",
                      "spend_smart ~ agem+achieve") # control variables
# Step 2: Estimation 
fit_sem_n2 <- sem(mysem2,data=subpop,estimator="MLM")
# Step 2.5: Survey-weighted estimation (Listwise deletion)
fit_sem_y2 <- lavaan.survey(fit_sem_n2,survey.design=cs_design)

# Multiple imputation (MI)
myimp <- mice(subpop,m=20,seed=1234,print=FALSE) # impute 20 times
# Combine m=20 data into one list
myimplist <- myimp %>%
  complete("all") %>%
  mitools::imputationList()
myimplist

# Complex survey analysis: Base-R approach 　
cs_design_imp <- svydesign(ids=~cluster,  # cluster 
                           strata=~strata,# strata
                           weights=~w,    # weight
                           fpc=~fpc,      # FPC
                           data=myimplist)

# Step 2.5: Survey-weighted estimation (Multiple imputation)
fit_sem_mi_y2 <- lavaan.survey(fit_sem_n2,survey.design=cs_design_imp)

# Step 3: Summary
summary(fit_sem_y2,fit=TRUE,standard=TRUE)
summary(fit_sem_mi_y2,fit=TRUE,standard=TRUE)

# Tidy the model outputs
sem_compare_no2 <- wrap_sem_function(fit_sem_n2,fit_sem_y2)
sem_compare_mi2 <- wrap_sem_function(fit_sem_n2,fit_sem_mi_y2)
left_join(sem_compare_no2 %>% select(lhs:myreporty) %>% 
            rename(myreporty_NO=myreporty),
          sem_compare_mi2 %>% select(lhs:myreporty) %>% 
            rename(myreporty_MI=myreporty)) %>% 
  filter(op=="~"|op==":=") %>% 
  write_excel_csv("Table_Part4_Ch7_5_compare_SEM_MI.csv")

quantile(sem_compare_no2$ratioF,c(0,0.25,0.50,0.75,1.00))
quantile(sem_compare_mi2$ratioF,c(0,0.25,0.50,0.75,1.00))

count(sem_compare_no2, conserv)
count(sem_compare_mi2, conserv)

###############################################
# Save the results
# save(fit_sem_n2,fit_sem_y2,fit_sem_mi_y2,
#      file="object_SEM_MI.RData")
# load("object_SEM_MI.RData")
