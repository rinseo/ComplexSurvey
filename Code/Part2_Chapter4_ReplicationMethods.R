####################################################################################
## Part 2. Analysis: Descriptive statistics 
## Chpt 4. Variance estimation using replication methods
####################################################################################

# Load libraries
library(tidyverse)  # Data manipulation 
library(survey)     # Complex survey analysis: Base-R approach
library(srvyr)      # Complex survey analysis: Tidyverse approach 

# Import data
# setwd("~/CSA_Code")
dat <- haven::read_spss("kyrbs2020.sav")

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
    happy=pr_hd, # perceived happiness (1[Happy]-5[Unhappy])
  ) %>% 
  select(female:happy,
         starts_with("m_gad"),
         starts_with("int_sp_ou"),
         cluster,strata,w,fpc)

###################################################
# BRR (balanced repeated replication)
# Begin with complex survey analysis using parametric approach
cs_design <- mydata %>% 
  as_survey_design(ids=cluster,     # cluster 
                   strata=strata,   # strata 
                   weights=w,       # weight 
                   fpc=fpc)         # FPC

# But an error pops up...
cs_design %>%
  as_survey_rep(type="BRR")

# ...because you can't split-half clusters with odd-number sizes
# Please select even-number sized clusters only (as an exercise)
even_cluster <- count(mydata, strata, cluster) %>% 
  group_by(strata) %>% 
  mutate(
    select=max(row_number())%%2 # left-overs
  ) %>% filter(select==0) %>% # only even numbers
  ungroup() %>% 
  select(cluster)
even_cluster
cs_design_even <- mydata %>% inner_join(even_cluster) %>% # choose even-numbered 
  as_survey_design(id=cluster,      # cluster 
                   strata=strata,   # strata 
                   weights=w)       # weight, without FPC
cs_BRR <- cs_design_even %>% as_survey_rep(type="BRR")
cs_BRR

# Compare with H-T estimation: Categorical variables 
cs_design_even %>% survey_count(sch_type,vartype='ci') %>% 
  data.frame() %>% mutate(range=n_upp-n_low)
cs_BRR %>% survey_count(sch_type,vartype='ci') %>% 
  data.frame() %>% mutate(range=n_upp-n_low)

# Compare with H-T estimation: Continuous variables 
cs_design_even %>% summarise(happy=survey_mean(happy,vartype='ci')) %>% 
  data.frame() %>% mutate(range=happy_upp-happy_low)
cs_BRR %>% summarise(happy=survey_mean(happy,vartype='ci')) %>% 
  data.frame() %>% mutate(range=happy_upp-happy_low)

# Repeat the procedure to all continuous variables 
HT_MCI <- cs_design_even %>% 
  select(female:happy) %>% 
  summarise(across(
    .cols=where(is.double),
    .fns=function(x){survey_mean(x,vartype='ci',
                                 na.rm=TRUE)}
  )) 
BRR_MCI <- cs_BRR %>% 
  select(female:happy) %>% 
  summarise(across(
    .cols=where(is.double),
    .fns=function(x){survey_mean(x,vartype='ci',
                                 na.rm=TRUE)}
  )) 
HT_BRR <- bind_rows(
  HT_MCI %>% mutate(approach="HT"),
  BRR_MCI %>% mutate(approach="BRR")
)
myfig <- HT_BRR %>% 
  pivot_longer(-approach) %>% 
  mutate(
    type="mean",
    type=ifelse(str_detect(name,"_low"),"ll",type),
    type=ifelse(str_detect(name,"_upp"),"ul",type),
    name=str_remove(name,"_low|_upp")
  ) %>% 
  pivot_wider(names_from="type",values_from="value") %>% 
  mutate(approach=fct_reorder(approach, row_number()))
myfig %>% 
  ggplot(aes(x=approach,y=mean,color=approach))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=ll,ymax=ul),width=0.2)+
  facet_wrap(~name,scales="free")+
  coord_flip()+
  scale_color_grey()+
  labs(x="",y="Mean with 95% CI")+
  theme_bw()+
  guides(color="none")
ggsave("Figure_Part2_Ch4_1_BRR_M_CI.png",width=25,height=22,units='cm')

# Fay's method: Handling with BRR's potential threats 
cs_FAY <- cs_design_even %>% 
  as_survey_rep(type="Fay", rho=0.30)  # Rho=0.30 (commonly used)
# Categorical variables
cs_FAY %>% survey_count(sch_type,vartype='ci')
# Continuous variables
cs_FAY %>% summarise(happy=survey_mean(happy,vartype='ci'))

## Base-R approach
even_mydata <- mydata %>% inner_join(even_cluster)
cs_design_classic <- svydesign(ids=~cluster,  # cluster 
                               strata=~strata,# strata 
                               weights=~w,    # weight 
                               data=even_mydata)
cs_BRR_classic <- as.svrepdesign(cs_design_classic,
                                 type="BRR")

# Categorical variable 
svytotal(~sch_type, cs_BRR_classic) %>% data.frame() # frequency with SE 
svytotal(~sch_type, cs_BRR_classic) %>% confint()    # frequency with 95% CI 
# Continuous variables
svymean(~happy, cs_BRR_classic) # mean 
svymean(~happy, cs_BRR_classic) %>% confint() # 95% CI 

########################################################
# JRR (Jackknife repeated replication)
# What is "Jackknife"?
set.seed(1234)
x <- rnorm(50,0,1)
mean(x); var(x)

# Exclude one and take the average of 49; repeat for 50 times
mean_x <- rep(NA,50)
for (i in 1:50){
  mean_x[i] <- mean(x[-i])
}
mean_x

mean(mean_x)
(50-1)*sum((mean_x - mean(mean_x))^2)

# JRR (JKn, if stratified sampled)
cs_JRR <- cs_design %>% 
  as_survey_rep(type="JKn")
cs_JRR

# Compare with H-T estimation: Categorical variables  
cs_design %>% survey_count(sch_type,vartype='ci') %>% 
  data.frame() %>% mutate(range=n_upp-n_low)
cs_JRR %>% survey_count(sch_type,vartype='ci') %>% 
  data.frame() %>% mutate(range=n_upp-n_low)

# Compare with H-T estimation: Continuous variables 
cs_design %>% summarise(happy=survey_mean(happy,vartype='ci')) %>% 
  data.frame() %>% mutate(range=happy_upp-happy_low)
cs_JRR %>% summarise(happy=survey_mean(happy,vartype='ci')) %>% 
  data.frame() %>% mutate(range=happy_upp-happy_low)

# Repeat the procedure to all continuous variables 
HT_MCI <- cs_design %>% 
  select(female:happy) %>% 
  summarise(across(
    .cols=where(is.double),
    .fns=function(x){survey_mean(x,vartype='ci',
                                 na.rm=TRUE)}
  )) 
JRR_MCI <- cs_JRR %>% 
  select(female:happy) %>% 
  summarise(across(
    .cols=where(is.double),
    .fns=function(x){survey_mean(x,vartype='ci',
                                 na.rm=TRUE)}
  )) 
HT_JRR <- bind_rows(
  HT_MCI %>% mutate(approach="HT"),
  JRR_MCI %>% mutate(approach="JRR")
) 
myfig <- HT_JRR %>% 
  pivot_longer(-approach) %>% 
  mutate(
    type="mean",
    type=ifelse(str_detect(name,"_low"),"ll",type),
    type=ifelse(str_detect(name,"_upp"),"ul",type),
    name=str_remove(name,"_low|_upp")
  ) %>% 
  pivot_wider(names_from="type",values_from="value") %>% 
  mutate(approach=fct_reorder(approach, row_number()))
myfig %>% 
  ggplot(aes(x=approach,y=mean,color=approach))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=ll,ymax=ul),width=0.2)+
  facet_wrap(~name,scales="free")+
  coord_flip()+
  scale_color_grey()+
  labs(x="",y="Mean with 95% CI")+
  theme_bw()+
  guides(color="none")
ggsave("Figure_Part2_Ch4_2_JRR_M_CI.png",width=25,height=22,units='cm')

## Base-R approach
cs_design_classic <- svydesign(ids=~cluster,  # cluster 
                               strata=~strata,# strata 
                               weights=~w,    # weight 
                               data=mydata)
cs_JRR_classic <- as.svrepdesign(cs_design_classic,
                                 type="JKn")

# Categorical variable 
svytotal(~sch_type, cs_JRR_classic) %>% data.frame() # frequency with SE 
svytotal(~sch_type, cs_JRR_classic) %>% confint()    # frequency with 95% CI 
# Continuous variables 
svymean(~happy, cs_JRR_classic)  # mean 
svymean(~happy, cs_JRR_classic) %>% confint()  # 95% CI 

###############################################
# Bootstrapping
# Resample for 100 times (increase>1,000 for practical use) 
set.seed(1234)
cs_boot <- cs_design %>% 
  as_survey_rep(type="bootstrap",
                replicates=100)

# Compare with H-T estimation: Categorical variables 
cs_design %>% survey_count(sch_type,vartype='ci') %>% 
  data.frame() %>% mutate(range=n_upp-n_low)
cs_boot %>% survey_count(sch_type,vartype='ci') %>% 
  data.frame() %>% mutate(range=n_upp-n_low)

# Compare with H-T estimation: Continuous variables 
cs_design %>% summarise(happy=survey_mean(happy,vartype='ci')) %>% 
  data.frame() %>% mutate(range=happy_upp-happy_low)
cs_boot %>% summarise(happy=survey_mean(happy,vartype='ci')) %>% 
  data.frame() %>% mutate(range=happy_upp-happy_low)

boot_MCI <- cs_boot %>% 
  select(female:happy) %>% 
  summarise(across(
    .cols=where(is.double),
    .fns=function(x){survey_mean(x,vartype='ci',
                                 na.rm=TRUE)}
  )) 
HT_boot <- bind_rows(
  HT_MCI %>% mutate(approach="HT"),
  JRR_MCI %>% mutate(approach="JRR"),
  boot_MCI %>% mutate(approach="Bootstrap")
) 
myfig <- HT_boot %>% 
  pivot_longer(-approach) %>% 
  mutate(
    type="mean",
    type=ifelse(str_detect(name,"_low"),"ll",type),
    type=ifelse(str_detect(name,"_upp"),"ul",type),
    name=str_remove(name,"_low|_upp")
  ) %>% 
  pivot_wider(names_from="type",values_from="value") %>% 
  mutate(approach=fct_reorder(approach, row_number()))
myfig %>% 
  ggplot(aes(x=approach,y=mean,color=approach))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=ll,ymax=ul),width=0.2)+
  facet_wrap(~name,scales="free")+
  coord_flip()+
  labs(x="",y="Mean with 95% CI")+
  theme_bw()+
  guides(color="none")
ggsave("Figure_Part2_Ch4_3_boot_M_CI.png",width=25,height=22,units='cm')

## Base-R approach 
set.seed(1234)
cs_boot_classic <- as.svrepdesign(cs_design_classic,
                                 type="bootstrap",
                                 replicates=100)

# Categorical variables
svytotal(~sch_type, cs_boot_classic) %>% data.frame() # frequency with SE 
svytotal(~sch_type, cs_boot_classic) %>% confint()    # frequency with 95% CI 
# Continuous variables
svymean(~happy, cs_boot_classic)  # mean 
svymean(~happy, cs_boot_classic) %>% confint()  # 95% CI 

###############################################
# Subpopulation analysis (female x stype)
# Categorical variables: Proportions
prop_sad <- list("HT"=cs_design,
                 "JRR"=cs_JRR, # takes a while
                 "Bootstrap"=cs_boot) %>% 
  map_df(~.x %>% group_by(female,sch_type) %>% 
           survey_count(exp_sad) %>%
           mutate(prop=n/sum(n)),.id="method")
prop_sad %>% arrange(female,sch_type,exp_sad) # point estimates are same

myfig <- prop_sad %>% 
  mutate(
    female=ifelse(female==0,"Boys","Girls"),
    exp_sad=ifelse(exp_sad==0,"Without sadness","With sadness"),
    subpop=str_c(sch_type,",\n", female),
    method=fct_reorder(method,row_number())
  )
myfig %>% ggplot(aes(x=subpop,y=prop,fill=method))+
  geom_bar(stat='identity',position=position_dodge(width=0.8),alpha=0.7)+
  labs(x="Subpopulation",y="Proportiona",fill="Method")+
  coord_cartesian(ylim=c(0.1,0.9))+
  theme_bw()+theme(legend.position="top")+
  facet_grid(~exp_sad)
ggsave("Figure_Part2_Ch4_4_compare_subpop_prop.png",width=20,height=15,units='cm')

# Continuous variables: Mean with 95% CI
mci_happy <- list("HT"=cs_design,
                  "JRR"=cs_JRR, # takes a while
                  "Bootstrap"=cs_boot) %>% 
  map_df(~.x %>% group_by(female,sch_type) %>% 
           summarise(mean=survey_mean(happy,vartype='ci')) %>%
           rename(ll=mean_low,ul=mean_upp),.id="method")
myfig <- mci_happy %>% 
  mutate(
    female=ifelse(female==0,"Boys","Girls"),
    subpop=str_c(sch_type,",\n", female),
    method=fct_reorder(method,row_number())
  )
myfig %>%
  ggplot(aes(x=subpop,y=mean,shape=method,color=method))+
  geom_point(size=2,position=position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=ll,ymax=ul),width=0.1,
                position=position_dodge(width=0.3))+
  labs(x="Subpopulation",y="Perceived happiness (Mean with 95% CI)",
       shape="Method",color="Method")+
  coord_cartesian(ylim=c(2.0,2.4))+
  theme_bw()+theme(legend.position="top")
ggsave("Figure_Part2_Ch4_5_compare_subpop_M_CI.png",width=20,height=15,units='cm')

## Base-R approach 
# Categorical variables
svytable(~exp_sad+sch_type+female,cs_JRR)
svytable(~exp_sad+sch_type+female,cs_boot)
# Continuous variables
svyby(~happy, by=~female+sch_type, 
      cs_JRR, svymean)
svyby(~happy, by=~female+sch_type, 
      cs_JRR, svymean) %>% confint() 
svyby(~happy, by=~female+sch_type, 
      cs_boot, svymean)
svyby(~happy, by=~female+sch_type, 
      cs_boot, svymean) %>% confint() 

###############################################
# Save results
# save(cs_BRR,cs_BRR_classic,BRR_MCI,
#      cs_JRR,cs_JRR_classic,JRR_MCI,
#      file="object_BRR_JRR.RData")
# load("object_BRR_JRR.RData")
