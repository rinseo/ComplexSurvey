####################################################################################
## Part 3. Analysis: Inferential statistics
## Chpt 5. Association tests
####################################################################################

# Load libraries
library(tidyverse)  # Data manipulation 
library(survey)     # Complex survey analysis: Base-R approach
library(srvyr)      # Complex survey analysis: Tidyverse approach 
library(jtools)     # Calculating correlation coefficient

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
    happy=pr_hd, # perceived happiness (1[Happy]-5[Unhappy])
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

# Take a subset of smartphone non-users 
cs_design_subset <- cs_design %>%  
  filter(use_smart==0) 
# Make frequency tables
mytable <- cs_design_subset %>% 
  survey_count(female,exp_sad)
mytable %>% group_by(female) %>% 
  mutate(prop=n/sum(n)) %>% 
  select(exp_sad,female,prop) %>% 
  pivot_wider(names_from="exp_sad",values_from="prop")

# Conventional chi-square tests 
xtabs(~female+exp_sad, mydata %>% filter(use_smart==0)) %>% 
  chisq.test() 
# Survey-weighted association tests: Wald statistics
svychisq(~female+exp_sad, cs_design_subset, statistic="Wald")
svychisq(~female+exp_sad, cs_design_subset, statistic="adjWald")
# Survey-weighted association tests: Rao-Scott statistics
svychisq(~female+exp_sad, cs_design_subset, statistic="Chisq")
svychisq(~female+exp_sad, cs_design_subset, statistic="F")
svychisq(~female+exp_sad, cs_design_subset, statistic="lincom")
# If an error occurs, try: install.packages("CompQuadForm")
svychisq(~female+exp_sad, cs_design_subset, statistic="saddlepoint")

tibble(
  source=c("NoCS","Wald","AdjWald","RS-1st","RS-2nd","lincom","saddlep"),
  pvalue=c(chisq.test(xtabs(~female+exp_sad, mydata %>% filter(use_smart==0)))$p.value, 
           c("Wald","adjWald","Chisq","F","lincom","saddlepoint") %>% 
             map_dbl(~svychisq(~female+exp_sad, cs_design_subset, statistic=.x)$p.value))
) %>% 
  arrange(pvalue) %>% data.frame() 

# T-tests
# One-sample t-test 
# Scenario 1 (NO)
t.test(mydata$happy,mu=2)
# Scenario 2 (CS)
svyttest(I(happy-2)~0, cs_design) 

# Paired sample t-test 
# Scenario 1 (NO)
mydata %>% summarise(m1=mean(m_gad_1),m2=mean(m_gad_2))
t.test(mydata$m_gad_1,mydata$m_gad_2,paired=TRUE)
# Scenario 2 (CS)
cs_design %>% 
  summarise(
    gad1=survey_mean(m_gad_1,vartype="ci"),
    gad2=survey_mean(m_gad_2,vartype="ci")
  )
svyttest(I(m_gad_1-m_gad_2)~0, cs_design) 

# Independent sample t-test 
# Scenario 1 (NO)
t.test(happy~sch_mdhg, mydata)

# Scenario 2 (CS)
svyttest(happy~sch_mdhg, cs_design)
cs_design %>% group_by(sch_mdhg) %>% 
  summarise(happy=survey_mean(happy, vartype='ci'))

########################################################
# Pearson's r 
# Visualize the relationship between happiness & financial 
myfig <- cs_design %>% 
  survey_count(happy,covid_suffer)
myfig
myfig %>% 
  ggplot(aes(x=covid_suffer, y=happy, alpha=n))+
  geom_tile()+
  scale_x_continuous(expand = expansion(mult = c(0, 0)))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)))+
  labs(x="Level of financial hardship due to COVID-19\n(1=High, 4=Low)",
       y="Perceived happiness\n(1=Happy, 5=Unhappy)")+
  theme_minimal()+
  theme(legend.position="none")
ggsave("Figure_Part3_Ch5_1_correlation_tileplot.png",width=10,height=10,units='cm')

# Scenario 1 (NO)
result_cor <- cor.test(~happy+covid_suffer, mydata)
result_cor$estimate     # r 
result_cor$statistic    # t-value 
result_cor$p.value      # p-value 

# Scenario 2 (CS)
svycor(~happy+covid_suffer, cs_design) 
result_svycor <- svycor(~happy+covid_suffer, cs_design, 
                        sig.stats=TRUE) # weights::wtd.cor()
# install.packages("weights") # Try if an error occurs
result_svycor$cors      # r 
result_svycor$t.values  # t-value 
result_svycor$p.values  # p-value 

########################################################
# Cronbach's alpha 
# Scenario 1 (NO)
psych::alpha(mydata %>% select(starts_with("m_gad")))
psych::alpha(mydata %>% select(starts_with("m_gad")))$total 
psych::alpha(mydata %>% filter(use_smart==1) %>% 
               select(starts_with("int_sp_ou")))$total

# Scenario 2 (CS)
svycralpha(~m_gad_1+m_gad_2+m_gad_3+m_gad_4+m_gad_5+m_gad_6+m_gad_6+m_gad_7,
           cs_design)
svycralpha(~int_sp_ou_1+int_sp_ou_2+int_sp_ou_3+int_sp_ou_4+int_sp_ou_5+
             int_sp_ou_6+int_sp_ou_7+int_sp_ou_8+int_sp_ou_9+int_sp_ou_10,
           cs_design %>% filter(use_smart==1))
