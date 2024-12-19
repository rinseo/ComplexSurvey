####################################################################################
## Part 4. Advanced topics
## Chpt 9. Propensity score analysis
####################################################################################

# Load libraries
library(tidyverse)  # Data manipulation 
library(survey)     # Complex survey analysis: Base-R approach
library(srvyr)      # Complex survey analysis: Tidyverse approach 
library(MatchIt)    # Propensity score matching

# Import data
# setwd("~/CSA_Code")
dat <- haven::read_spss("kyrbs2020.sav")

# Lower cases for convenience
names(dat) <- tolower(names(dat))

# Preprocess data
mydata <- dat %>% 
  mutate(
    achieve=e_s_rcrd, # dependent variable
    stype=factor(stype,labels=c("Co-ed school","All-boys school","All-girls school")), # type of schools
    sch_sep=ifelse(stype=="Co-ed school",0,1), # cause: 1(male/female only), 0(mixed)
    female=ifelse(sex==1,0,1) %>% as.factor(), # male/female  
    ses=e_ses,
    sch_mdhg=as.factor(mh),
    use_smart=ifelse(int_spwd==1&int_spwk==1,1,0) %>% as.factor(),
    scale_gad=rowMeans(dat %>% select(starts_with("m_gad"))),
    exercise=pa_tot-1
  ) %>% 
  select(achieve:exercise,
         cluster,strata,w,fpc)

# Scenario 1 (NO)
lm_naive <- lm(achieve~sch_sep,mydata) 
coef(lm_naive)
confint(lm_naive)

# Scenario 1 (CS)
## Tidyverse approach　
cs_design <- mydata %>% 
  as_survey_design(ids=cluster,     # cluster 
                   strata=strata,   # strata
                   weights=w,       # weight 
                   fpc=fpc)         # FPC 
lm_cs <- svyglm(achieve~sch_sep,
                design=cs_design) 
coef(lm_cs)
confint(lm_cs)

## ATT: DuGoff et al. (2014)
# Step 1: Propensity score matching
set.seed(1234)
m_out <- matchit(sch_sep~female+ses+sch_mdhg+scale_gad+exercise+w,
                 data=mydata,
                 method="nearest",
                 link="logit",
                 replace=FALSE)
# Return matched data 
m_data <- match.data(m_out)

# Step 2: Complex survey analysis
match_design_att <- m_data %>% 
  as_survey_design(ids=cluster, 
                   strata=strata,
                   weights=w,
                   fpc=fpc)

# Step 3: Estimate ATT
matchATT <- svyglm(achieve~sch_sep,
                     design=match_design_att) 
coef(matchATT)
confint(matchATT)

## ATE: Ridgeway et al., (2015)
# Step 1: Propensity score weighting
psmodel <- svyglm(sch_sep~female+ses+sch_mdhg+scale_gad+exercise,
                  design=cs_design,family=quasibinomial)

# Step 2: Complex survey analysis
# Calculate IPTW, final weights
cs_design2 <- cs_design %>% 
  mutate(
    psc=predict(psmodel, cs_design, type="response"),
    iptw=ifelse(sch_sep==1,1/psc,1/(1-psc)), # IPTW
    atewt=iptw*w  # combine IPTW and survey weights
  ) 

# Step 3: Estimate ATE
cs_design_ate <- cs_design2$variables %>% 
  as_survey_design(ids=cluster,  
                   strata=strata,
                   weights=atewt,
                   fpc=fpc) 
PSW_ATE <- svyglm(achieve~sch_sep,
                    design=cs_design_ate)
coef(PSW_ATE)
confint(PSW_ATE)

# Compare the effect sizes 
myfig <- list("Complex survey analysis: No"=lm_naive,
              "Complex survey analysis: Yes"=lm_cs,
              "Complex survey analysis: Yes,\nPropensity score matched"=matchATT,
              "Complex survey analysis: Yes,\nPropensity score weighted"=PSW_ATE) %>% 
  map_df(~c(coef(.)[2],confint(.)[2,]),.id="method") %>% 
  mutate(method=fct_reorder(method,row_number()),
         CS=ifelse(str_detect(method,"Yes"),"Yes","No"),
         PS=ifelse(str_detect(method,"Propensity score"),"Yes","No"))
myfig %>%
  ggplot(aes(x=method,y=sch_sep,lty=PS,shape=CS))+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=`2.5 %`,ymax=`97.5 %`),size=0.5,width=0.2)+
  geom_hline(yintercept=0,lty=3)+
  scale_linetype_manual(values=c(2,1))+
  theme_classic()+
  labs(x="\nMethod",y="Treatment effect size (Mean with 95% CI)",
       lty="Propensity score analysis",shape="Complex survey analysis")
ggsave("Figure_Part4_Ch9_1_compare_PSA_CS.png",width=15,height=13,units='cm')
