####################################################################################
## Part 4. Advanced topics
## Chpt 8. Multilevel modeling
####################################################################################

# Load libraries
library(tidyverse)   # Data manipulation 
library(lme4)        # Multilevel modeling (MLM)
library(WeMix)       # Survey-weighted MLM
library(broom.mixed) # Tidying MLM results

# Import data
# setwd("~/CSA_Code")
dat <- haven::read_spss("kyrbs2020.sav")
dat

# Lower cases for convenience
names(dat) <- tolower(names(dat))

# Preprocess data
mydata <- dat %>% 
  mutate(
    ses=e_ses, # socioeconomic status 
    sch_type=factor(stype,labels=c("Co-ed school","All-boys school","All-girls school")), # type of schools
    sch_mdhg=factor(mh,labels=c("High school","Middle school")), # middle/high school
    happy=6-pr_hd, # perceived happiness (reversed; 1[NOT]-5[HAPPY]) 
  ) %>% 
  select(ses,sch_type,sch_mdhg,happy,
         cluster,w)

# Re-scaling (Carle, 2009)
## Level-1 (student) weights
mydata <- mydata %>% 
  mutate(sqw=w^2) %>% 
  group_by(cluster) %>% 
  mutate(
    sumsqw=sum(sqw),
    sumw=sum(w),
    nj=length(w),
    aw1=w*nj/sumw,
    bw1=w*(sumw/sumsqw)
  ) %>% ungroup() %>% 
  select(-sqw,-sumsqw,-sumw,-nj)

cor(mydata %>% select(w,aw1,bw1))

## Level-2 (school) weights
# Total number of co-ed or all-boys/girls schools 
# Reference) https://gsis.kwdi.re.kr/statHtml/statHtml.do?orgId=338&tblId=DT_1LCB041
# w_sc <- readxl::read_excel("number_of_schools_2020.xlsx")
w_sc
scP <- w_sc %>% pivot_longer(cols=everything()) %>% 
  separate(name,c("sch_mdhg","sch_type"),sep="_") %>% 
  mutate(
    sch_mdhg=factor(sch_mdhg,labels=c("Public, High school",
                                      "Private, High school",
                                      "Middle school",
                                      "Vocational, High school",
                                      "Specialized, High school")),
    sch_mdhg=ifelse(str_detect(sch_mdhg,"Middle"),"Middle school","High school"),
    sch_type=factor(sch_type,labels=c("Co-ed school","All-boys school","All-girls school"))
  ) %>% 
  group_by(sch_mdhg,sch_type) %>% 
  summarise(n=sum(value)) %>% ungroup() %>% 
  mutate(propP=n/sum(n)) %>% 
  select(-n)
scP
scS <- mydata %>% 
  count(cluster,sch_mdhg,sch_type) %>% 
  ungroup() %>% 
  count(sch_mdhg,sch_type) %>% 
  mutate(propS=n/sum(n)) %>% 
  select(-n)
scS

sc_weight2 <- full_join(scS,scP,by=c("sch_mdhg","sch_type")) %>% 
  mutate(w2=propP/propS) %>% select(-propS,-propP)
sc_weight2

mydata <- mydata %>% 
  full_join(sc_weight2,by=c("sch_mdhg","sch_type"))

# Group-mean centering
mydata_mc <- mydata %>% 
  group_by(cluster) %>% 
  mutate(mc_ses=ses-mean(ses))

# (Unweighted) multilevel models 
M0_N <- lmer(happy~1+(1|cluster),mydata_mc)
MI_N <- lmer(happy~mc_ses+sch_mdhg+(1|cluster),mydata_mc)
MIS_N <- lmer(happy~mc_ses+sch_mdhg+(mc_ses||cluster),mydata_mc) # random intercept // slope
summary(M0_N)
summary(MI_N)

1-0.87436/0.9018  # Level-1 PRE 
1-0.01492/0.0202  # Level-2 PRE 

# Survey-weighted multilevel models
## Method A
M0_A <- mix(happy~1+(1|cluster),mydata,
            weights=c("aw1","w2"),
            cWeights=TRUE) # since they're rescaled 
MI_A <- mix(happy~ses+sch_mdhg+(1|cluster),mydata,
            weights=c("aw1","w2"),
            cWeights=TRUE,
            center_group=list("cluster"=~ses)) # group-mean centering
MIS_A <- mix(happy~ses+sch_mdhg+(ses||cluster),mydata, # random intercept // slope
             weights=c("aw1","w2"),
             cWeights=TRUE,
             center_group=list("cluster"=~ses))
summary(M0_A)
summary(MI_A)

## Method B
M0_B <- mix(happy~1+(1|cluster),mydata,
            weights=c("bw1","w2"),
            cWeights=TRUE) # since they're rescaled  
MI_B <- mix(happy~ses+sch_mdhg+(1|cluster),mydata,
            weights=c("bw1","w2"),
            cWeights=TRUE,
            center_group=list("cluster"=~ses)) # group-mean centering
MIS_B <- mix(happy~ses+sch_mdhg+(ses||cluster),mydata, #random intercept // slope
             weights=c("bw1","w2"),
             cWeights=TRUE,
             center_group=list("cluster"=~ses))

# Tidy the model outputs
## lmer(): broom.mixed::tidy()
tidy(M0_N) %>% mutate(rid=row_number()) # create row id

summ_N <- list("M0"=M0_N,"MI"=MI_N,"MIS"=MIS_N) %>% 
  map_df(~tidy(.) %>% mutate(rid=row_number()),
         .id="model")

# mix(): Create a user-defined function
tidy_mix_function <- function(object_mix){
  # Fixed effects
  FE <- summary(object_mix)$coef %>% data.frame()
  names(FE) <- c("estimate","std.error","statistic")
  # Random effects
  RE <- object_mix$vars %>% data.frame() %>% rename(estimate=".")
  # Bind them into rows
  bind_rows(FE %>% mutate(effect="fixed"),
            RE %>% mutate(effect="ran_pars")) %>% 
    as_tibble(rownames="group_term") %>% 
    select(effect,group_term,everything())
    
}
tidy_mix_function(M0_A) %>% mutate(rid=row_number())

summ_A <- list("M0"=M0_A,"MI"=MI_A,"MIS"=MIS_A) %>% 
  map_df(~tidy_mix_function(.) %>% mutate(rid=row_number()),
         .id="model")
summ_B <- list("M0"=M0_B,"MI"=MI_B,"MIS"=MIS_B) %>% 
  map_df(~tidy_mix_function(.) %>% mutate(rid=row_number()),
         .id="model")

bind_rows(summ_N %>% mutate(method="N"),
          summ_A %>% mutate(method="A"),
          summ_B %>% mutate(method="B")) %>% 
  select(-group,-term,-group_term) %>% 
  mutate(statistic=ifelse(is.na(statistic),
                          yes="",
                          no=str_c("\n(",
                                   format(round(statistic,4),nsmall=4),
                                   ")")),
         myreport=str_c(format(round(estimate,4),nsmall=4),statistic),
         name=str_c(model,"-",method)) %>% 
  left_join(summ_N %>% select(model,rid,term)) %>% 
  arrange(model) %>% 
  select(name,effect,term,myreport) %>% 
  pivot_wider(names_from=name,values_from=myreport) %>% 
  arrange(effect,term) %>% 
  mutate(across(
    .cols=everything(),
    .fns=function(x){ifelse(is.na(x)," ",x)}
  )) %>% 
  write_excel_csv("Table_Part4_Ch8_2_compare_MLM.csv")

# Fixed effect
waldTest(MIS_A,type="beta",
         coefs=c("ses","sch_mdhgMiddle school"))
waldTest(MIS_A,type="beta",
         coefs="ses",
         hypothesis=c(-0.22))

# Random effect
waldTest(MIS_A,type="Lambda")
waldTest(MIS_A,type="Lambda",
         coefs=c("cluster.ses"),
         hypothesis=c(0.06))

###############################################
# Save the results
# save(M0_A,MI_A,MIS_A,
#      M0_B,MI_B,MIS_B,
#      file="object_MLM_CS.RData")
# load("object_MLM_CS.RData")
