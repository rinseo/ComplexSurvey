####################################################################################
## Part 3. Analysis: Inferential statistics 
## Chpt 6. Linear regression models
####################################################################################

# Load libraries
library(tidyverse)  # Data manipulation 
library(survey)     # Complex survey analysis: Base-R approach
library(srvyr)      # Complex survey analysis: Tidyverse approach 
library(MASS)       # Ordinal logistic, negative binomial models
select <- dplyr::select 
library(sjstats)    # Survey-weighted negative binomial models
library(nnet)       # Multinomial logistic models
library(svyVGAM)    # Survey-weighted multinomial logistic models
library(modelr)     # Visualizing results
library(broom)      # Tidying results

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
    happy=6-pr_hd, # perceived happiness (reversed; 1[Unhappy]-5[Happy])
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

###########################################################
# For (normally distributed) continuous variables: OLS models
pred_1 <- "ses+achieve+use_smart+scale_gad+I(scale_gad^2)+covid_suffer"
formula1 <- as.formula(str_c("happy~",pred_1))
pred_2 <- "+scale_gad:covid_suffer+I(scale_gad^2):covid_suffer"
formula2 <- as.formula(str_c("happy~",pred_1,pred_2))
formula1; formula2

# Mean-centering
mydata_mc <- mydata %>% 
  mutate(across(
    .cols=c(ses,achieve,scale_gad,covid_suffer,spend_smart,scale_spaddict),
    .fns=function(x){x-mean(x,na.rm=TRUE)}
  ))
cs_design_mc <- cs_design %>%
  mutate(across(
    .cols=c(ses,achieve,scale_gad,covid_suffer,spend_smart,scale_spaddict),
    .fns=function(x){x-svymean(~x,cs_design,na.rm=TRUE)}
  ))

# Scenario 1 (NO)
ols1 <- lm(formula1,mydata_mc)
ols2 <- lm(formula2,mydata_mc)

# Scenario 2 (CS)
normal_glm1 <- svyglm(formula1, cs_design_mc)
normal_glm2 <- svyglm(formula2, cs_design_mc)

# Compare the results under two scenarios
summary(ols1)
summary(normal_glm1)

# Footnote: Weighted regression?
lm(formula1,cs_design_mc,weights=w) %>% summary()
lm(formula2,cs_design_mc,weights=w) %>% summary()

# Calculate R2 
normal_glm0 <- svyglm(happy~1, cs_design_mc)
summary(normal_glm0)$dispersion
cs_design_mc %>% summarise(var=survey_var(happy))

# Compare the model fits 
r2_ols1 <- summary(ols1)$r.squared %>% round(4)
r2_ols2 <- summary(ols2)$r.squared %>% round(4)
r2_ols1;r2_ols2;r2_ols2-r2_ols1

r2_normal_glm1 <- 1-(summary(normal_glm1)$dispersion[1]/summary(normal_glm0)$dispersion[1]) %>% round(4)
r2_normal_glm2 <- 1-(summary(normal_glm2)$dispersion[1]/summary(normal_glm0)$dispersion[1]) %>% round(4)
r2_normal_glm1;r2_normal_glm2;r2_normal_glm2-r2_normal_glm1

# User-defined function for model summary 
wrap_function <- function(object_glm){
  tidy(object_glm) %>% 
    mutate(
      sigstar=cut(p.value,c(-Inf,0.001,0.01,0.05,1),
                  labels=c("***","**","*","")),
      myreport=str_c(format(round(estimate,4),nsmall=4),
                     sigstar,"\n(",
                     format(round(std.error,4),nsmall=4),")")
    ) %>% select(-(estimate:sigstar))
}
# Tidy the model outputs
list("m1_no"=ols1,"m1_yes"=normal_glm1,
     "m2_no"=ols2,"m2_yes"=normal_glm2) %>%
  map_df(wrap_function,.id="model") %>% 
  pivot_wider(names_from=model,values_from=myreport) %>% 
  write_excel_csv("Table_Part3_Ch6_1_regression_normal.csv")

# Model comparison 
anova(ols1,ols2) 
anova(ols1,ols2)$`Pr(>F)`
anova(normal_glm1, normal_glm2, method="Wald") # Wald test (See Chapter 5)
# anova(normal_glm1, normal_glm2, method="LRT") # Rao-Scott test (See Chapter 5) 
AIC(ols1);AIC(ols2)
BIC(ols1);BIC(ols2)
AIC(normal_glm1, normal_glm2)
BIC(normal_glm1, normal_glm2, maximal=normal_glm2)

# Footnote: regTermTest() can be used
regTermTest(normal_glm2, 
            ~scale_gad:covid_suffer+I(scale_gad^2):covid_suffer) # Wald test
regTermTest(normal_glm2, 
            ~scale_gad:covid_suffer+I(scale_gad^2):covid_suffer,
            method="LRT") # Rao-Scott test

# Visualization
newdata <- mydata %>% 
  data_grid(ses=0,achieve=0,use_smart=1,
            scale_gad=(1:4)-svymean(~scale_gad,cs_design),
            covid_suffer=(1:4)-svymean(~covid_suffer,cs_design)) %>% 
  mutate(use_smart=as.factor(use_smart))
newdata

# Get predicted values
mypred_cs <- predict(normal_glm2,newdata,se.fit=TRUE) %>% as_tibble()
mypred_cs

mypred_no <- predict(ols2,newdata,se.fit=TRUE) %>% as_tibble()
mypred_no

myfig <- bind_rows(
  newdata %>% mutate(model="Complex survey analysis: Yes",
                     predy=mypred_cs$link,SE=mypred_cs$SE),
  newdata %>% mutate(model="Complex survey analysis: No",
                     predy=mypred_no$fit,SE=mypred_no$se.fit)) %>% 
  mutate(
    ll=predy-1.96*SE,ul=predy+1.96*SE,
    scale_gad=factor(scale_gad,
           labels=c("1. Very low","2. ", "3. ", "4. Very high"))
  ) 
myfig %>% 
  ggplot(aes(x=covid_suffer,y=predy,color=scale_gad,fill=scale_gad))+
  geom_line()+
  geom_ribbon(aes(ymin=ll,ymax=ul),alpha=0.2,color=NA)+
  labs(x="Level of financial hardship due to COVID-19\n(1=High, 4=Low)",
       y="Perceived happiness\n(1=Unhappy, 5=Happy)",
       color="GAD (Generalized Anxiety Disorder)",fill="GAD (Generalized Anxiety Disorder)")+
  scale_x_continuous(breaks=(1:4)-svymean(~covid_suffer,cs_design),
                     labels=c("1. Very high","2. ","3. ","4. Very low"))+
  theme_bw()+
  theme(legend.position="top",
        axis.text.x = element_text(hjust=1,angle = 90))+
  facet_wrap(~model)
ggsave("Figure_Part3_Ch6_1_interaction_covidsuffer_gad_happy.png",height=14,width=20,units='cm')

# Subpopulation(=smartphone users) analysis
# Set model formula
pred_sub0 <- "ses+achieve+(scale_gad+I(scale_gad^2))*covid_suffer"
pred_sub1 <- "+spend_smart+scale_spaddict+I(scale_spaddict^2)"
formula_sub1 <- as.formula(str_c("happy~",pred_sub0,pred_sub1))
pred_sub2 <- "+spend_smart:scale_spaddict+spend_smart:I(scale_spaddict^2)"
formula_sub2 <- as.formula(str_c("happy~",pred_sub0,pred_sub1,pred_sub2))
formula_sub1; formula_sub2

# Scenario 1 (NO)
mydata_mc_sub <- mydata_mc %>% filter(use_smart==1)
ols_sub1 <- lm(formula_sub1, mydata_mc_sub)
ols_sub2 <- lm(formula_sub2, mydata_mc_sub) 
# Scenario 2 (CS)
cs_design_mc_sub <- cs_design_mc %>% filter(use_smart==1)
normal_glm_sub1 <- svyglm(formula_sub1, cs_design_mc_sub)
normal_glm_sub2 <- svyglm(formula_sub2, cs_design_mc_sub)
# Tidy the model outputs
list("m1_no"=ols_sub1,"m1_yes"=normal_glm_sub1,
     "m2_no"=ols_sub2,"m2_yes"=normal_glm_sub2) %>%
  map_df(wrap_function,.id="model") %>% 
  pivot_wider(names_from=model,values_from=myreport) %>% 
  write_excel_csv("Table_Part3_Ch6_2_regression_normal_sub.csv")

# Compare the model fits
r2_ols1 <- summary(ols_sub1)$r.squared %>% round(4)
r2_ols2 <- summary(ols_sub2)$r.squared %>% round(4)
r2_ols1;r2_ols2;r2_ols2-r2_ols1

normal_glm_sub0 <- svyglm(happy~1, cs_design_mc %>% filter(use_smart==1))
r2_normal_glm1 <- 1-(summary(normal_glm_sub1)$dispersion[1]/summary(normal_glm_sub0)$dispersion[1]) %>% round(4)
r2_normal_glm2 <- 1-(summary(normal_glm_sub2)$dispersion[1]/summary(normal_glm_sub0)$dispersion[1]) %>% round(4)
r2_normal_glm1;r2_normal_glm2;r2_normal_glm2-r2_normal_glm1

# Model comparison 
anova(ols_sub1,ols_sub2) 
anova(ols_sub1,ols_sub2)$`Pr(>F)`
anova(normal_glm_sub1, normal_glm_sub2, method="Wald")
anova(normal_glm_sub1, normal_glm_sub2, method="LRT")
AIC(ols_sub1);AIC(ols_sub2)
BIC(ols_sub1);BIC(ols_sub2)
AIC(normal_glm_sub1,normal_glm_sub2)
BIC(normal_glm_sub1,normal_glm_sub2, maximal=normal_glm_sub2)

# Visualization
newdata <- mydata %>% 
  data_grid(
    ses=0,achieve=0,use_smart=1,scale_gad=0,covid_suffer=0,
    spend_smart=(1:5)-svymean(~spend_smart,subset(cs_design,use_smart==1)),
    scale_spaddict=(1:4)-svymean(~scale_spaddict,subset(cs_design,use_smart==1)))
# Get predicted values
mypred_cs <- predict(normal_glm_sub2,newdata,se.fit=TRUE) %>% as_tibble()
mypred_no <- predict(ols_sub2,newdata,se.fit=TRUE) %>% as_tibble()

myfig <- bind_rows(
  newdata %>% mutate(model="Complex survey analysis: Yes",
                     predy=mypred_cs$link,SE=mypred_cs$SE),
  newdata %>% mutate(model="Complex survey analysis: No",
                     predy=mypred_no$fit,SE=mypred_no$se.fit)) %>% 
  mutate(
    ll=predy-1.96*SE,ul=predy+1.96*SE,
    scale_spaddict=factor(scale_spaddict,
                     labels=c("1. Very low","2. ", "3. ", "4. Very high")),
    spend_smart=spend_smart+svymean(~spend_smart,subset(cs_design,use_smart==1))
  ) 
myfig %>% 
  ggplot(aes(x=spend_smart,y=predy,color=scale_spaddict,fill=scale_spaddict))+
  geom_line()+
  geom_ribbon(aes(ymin=ll,ymax=ul),alpha=0.2,color=NA)+
  labs(x="Hours of media use",y="Perceived happiness",
       color="Smartphone addiction",fill="Smartphone addiction")+
  scale_x_continuous(breaks=1:6)+
  theme_bw()+
  theme(legend.position="top")+
  facet_wrap(~model,nrow=1)
ggsave("Figure_Part3_Ch6_2_interaction_smartphone_happy_sub.png",height=12,width=20,units='cm')

###########################################################
# For binary variables: Binary logistic models
formula1 <- as.formula(str_c("exp_sad~",pred_1))
formula2 <- as.formula(str_c("exp_sad~",pred_1,pred_2))
formula1; formula2

# Scenario 1 (NO)
blog1 <- glm(formula1,mydata_mc,family=binomial)
blog2 <- glm(formula2,mydata_mc,family=binomial)

# Scenario 1 (NO)
# If you change into "family=binomial" an error occurs... 
binary_glm1_bin <- svyglm(formula1,cs_design_mc,family=binomial)
binary_glm1 <- svyglm(formula1,cs_design_mc,family=quasibinomial)
binary_glm2 <- svyglm(formula2,cs_design_mc,family=quasibinomial)
# ...but you can get the same results
summary(binary_glm1_bin)
summary(binary_glm1)

# Tidy the model outputs
list("m1_no"=blog1,"m1_yes"=binary_glm1,
     "m2_no"=blog2,"m2_yes"=binary_glm2) %>%
  map_df(wrap_function,.id="model") %>% 
  pivot_wider(names_from=model,values_from=myreport) %>% 
  write_excel_csv("Table_Part3_Ch6_3_regression_binary_logistic.csv")

# Model comparison
lmtest::lrtest(blog1, blog2)
anova(binary_glm1, binary_glm2, method="Wald")
anova(binary_glm1, binary_glm2, method="LRT")
c(AIC(blog1),AIC(blog2),BIC(blog1),BIC(blog2)) %>% round(0)
AIC(binary_glm1, binary_glm2) %>% round(0)
BIC(binary_glm1, binary_glm2, maximal=normal_glm2) %>% round(0)

# Visualization
newdata <- mydata %>% 
  data_grid(ses=0,achieve=0,use_smart=1,
            scale_gad=(1:4)-svymean(~scale_gad,cs_design),
            covid_suffer=0.1*(10:40)-svymean(~covid_suffer,cs_design)) %>% 
  mutate(use_smart=as.factor(use_smart))
# You can use type='response' but it's different from 95% CI
mypred_cs <- predict(binary_glm2,newdata,se.fit=TRUE) %>% as_tibble()
mypred_no <- predict(blog2,newdata,se.fit=TRUE) %>% as_tibble()

myfig <- bind_rows(
  newdata %>% mutate(model="Complex survey analysis: Yes",
                     link=mypred_cs$link,SE=mypred_cs$SE),
  newdata %>% mutate(model="Complex survey analysis: No",
                     link=mypred_no$fit,SE=mypred_no$se.fit)) %>% 
  mutate(
    predy=1/(1+exp(-link)),       # logit to probability
    ll=1/(1+exp(-link+1.96*SE)),  # logit to probability 
    ul=1/(1+exp(-link-1.96*SE)),  # logit to probability
    scale_gad=factor(scale_gad,
                     labels=c("1. Very low","2. ", "3. ", "4. Very high"))
  ) 
myfig %>% 
  ggplot(aes(x=covid_suffer,y=predy,color=scale_gad,fill=scale_gad))+
  geom_line()+
  geom_ribbon(aes(ymin=ll,ymax=ul),alpha=0.2,color=NA)+
  labs(x="Level of financial hardship due to COVID-19\n(1=High, 4=Low)",
       y="Probability of feeling depressed",
       color="GAD (Generalized Anxiety Disorder)",fill="GAD (Generalized Anxiety Disorder)")+
  scale_x_continuous(breaks=(1:4)-svymean(~covid_suffer,cs_design),
                     labels=c("1. Very high","2. ","3. ","4. Very low"))+
  theme_bw()+
  theme(legend.position="top",
        axis.text.x = element_text(hjust=1,angle = 90))+
  facet_wrap(~model)
ggsave("Figure_Part3_Ch6_3_interaction_covidsuffer_gad_sad.png",height=14,width=20,units='cm')

# Subpopulation(=smartphone users) analysis
# Set model formula
formula_sub1 <- as.formula(str_c("exp_sad~",pred_sub0,pred_sub1))
formula_sub2 <- as.formula(str_c("exp_sad~",pred_sub0,pred_sub1,pred_sub2))
formula_sub1
formula_sub2

# Scenario 1 (NO) 
blog_sub1 <- glm(formula_sub1, mydata_mc_sub, family=binomial)
blog_sub2 <- glm(formula_sub2, mydata_mc_sub, family=binomial) 
# Scenario 2 (CS)
binary_glm_sub1 <- svyglm(formula_sub1, cs_design_mc_sub, family=quasibinomial)
binary_glm_sub2 <- svyglm(formula_sub2, cs_design_mc_sub, family=quasibinomial)
# Tidy the model outputs
list("m1_no"=blog_sub1,"m1_yes"=binary_glm_sub1,
     "m2_no"=blog_sub2,"m2_yes"=binary_glm_sub2) %>%
  map_df(wrap_function,.id="model") %>% 
  pivot_wider(names_from=model,values_from=myreport) %>% 
  write_excel_csv("Table_Part3_Ch6_4_regression_binary_logistic_sub.csv")

# Model comparison
lmtest::lrtest(blog_sub1,blog_sub2)
anova(binary_glm_sub1,binary_glm_sub2, method="Wald")
anova(binary_glm_sub1,binary_glm_sub2, method="LRT")
c(AIC(blog_sub1),AIC(blog_sub2),BIC(blog_sub1),BIC(blog_sub2)) %>% round(0)
AIC(binary_glm_sub1,binary_glm_sub2) %>% round(0)
BIC(binary_glm_sub1,binary_glm_sub2, maximal=binary_glm_sub2) %>% round(0)

# Visualization
newdata <- mydata %>% 
  data_grid(ses=0,achieve=0,use_smart=1,
            scale_gad=0,covid_suffer=0,
            spend_smart=c(1,3,5),scale_spaddict=0.1*(10:40))
# Get predicted values
mypred_cs <- predict(binary_glm_sub2,newdata,se.fit=TRUE) %>% as_tibble()
mypred_no <- predict(blog_sub2,newdata,se.fit=TRUE) %>% as_tibble()
myfig <- bind_rows(
  newdata %>% mutate(model="Complex survey analysis: Yes",
                     link=mypred_cs$link,SE=mypred_cs$SE),
  newdata %>% mutate(model="Complex survey analysis: No",
                     link=mypred_no$fit,SE=mypred_no$se.fit)) %>% 
  mutate(
    predy=1/(1+exp(-link)),       # logit to probability
    ll=1/(1+exp(-link+1.96*SE)),  # logit to probability 
    ul=1/(1+exp(-link-1.96*SE)),  # logit to probability
    spend_smart=factor(spend_smart,
                       labels=c("1 hr","3 hrs","5 hrs"))
  ) 
myfig %>% 
  ggplot(aes(x=scale_spaddict,y=predy,color=spend_smart,fill=spend_smart))+
  geom_line()+
  geom_ribbon(aes(ymin=ll,ymax=ul),alpha=0.2,color=NA)+
  labs(x="Smartphone addiction",y="Probability of feeling depressed",
       color="Hours of smartphone use",fill="Hours of smartphone use")+
  scale_x_continuous(breaks=c(1,4),labels=c("1. Very low","4. Very high"))+
  theme_bw()+
  theme(legend.position="top",
        axis.text.x = element_text(angle = 90))+
  facet_wrap(~model)
ggsave("Figure_Part3_Ch6_4_interaction_smartphone_sad_sub.png",height=12,width=20,units='cm')

###########################################################
# For ordinal variables: Ordinal logistic models
pred_1 <- "ses+achieve+use_smart+scale_gad+I(scale_gad^2)+covid_suffer"
formula1 <- as.formula(str_c("happy~",pred_1))
pred_2 <- "+scale_gad:covid_suffer+I(scale_gad^2):covid_suffer"
formula2 <- as.formula(str_c("happy~",pred_1,pred_2))
formula1; formula2

# Transform the dependent variable (numeric) to factor type 
mydata_mc_ologit <- mydata_mc %>% 
  mutate(happy=as.factor(happy))
cs_design_mc_ologit <- cs_design_mc %>% 
  mutate(happy=as.factor(happy))

# Scenario 1 (NO)
olog1 <- polr(formula1,mydata_mc_ologit,Hess=TRUE)
olog2 <- polr(formula2,mydata_mc_ologit,Hess=TRUE)

# Scenario 2 (CS)
ologit_glm1 <- svyolr(formula1,cs_design_mc_ologit)
ologit_glm2 <- svyolr(formula2,cs_design_mc_ologit)

summary(olog1)
summary(ologit_glm1)

# User-defined function for model summary 
wrap_ordinal_function <- function(object_glm){
  tidy(object_glm) %>% 
    mutate(
      p.value=pt(abs(statistic),object_glm$df.residual,lower.tail=F),
      sigstar=cut(p.value,c(-Inf,0.001,0.01,0.05,1),
                  labels=c("***","**","*","")),
      myreport=str_c(format(round(estimate,4),nsmall=4),
                     sigstar,"\n(",
                     format(round(std.error,4),nsmall=4),")")
    ) %>% select(-(estimate:sigstar))
}
# Tidy the model outputs
list("m1_no"=olog1,"m1_yes"=ologit_glm1,
     "m2_no"=olog2,"m2_yes"=ologit_glm2) %>%
  map_df(wrap_ordinal_function,.id="model") %>% 
  pivot_wider(names_from=model,values_from=myreport) %>% 
  write_excel_csv("Table_Part3_Ch6_5_regression_ordinal_logistic.csv")

# Model comparison
lmtest::lrtest(olog1,olog2)
c(AIC(olog1),AIC(olog2),BIC(olog1),BIC(olog2)) %>% round(0)

# Visualization
newdata <- mydata %>% 
  data_grid(ses=0,achieve=0,use_smart=1,
            scale_gad=svyquantile(~scale_gad, cs_design, c(0.25,0.5,0.75))$scale_gad[1:3],
            covid_suffer=0.1*(10:40)-svymean(~covid_suffer,cs_design)) %>% 
  mutate(use_smart=as.factor(use_smart))
# Built-in predict() function only provides point estimates
mypred <- predict(olog2,newdata,type="probs") %>% as_tibble()
mypred

# Get predicted values with 95% CI
newdata <- newdata %>% 
  mutate(
    use_smart1=1,
    scale_gad2=scale_gad^2,
    scale_gad_covid_suffer=scale_gad*covid_suffer,
    scale_gad2_covid_suffer=scale_gad^2*covid_suffer
  ) %>% 
  select(ses,achieve,use_smart1,scale_gad,
         scale_gad2,covid_suffer,scale_gad_covid_suffer,scale_gad2_covid_suffer)
names(newdata); names(coef(olog2)) # make sure they're in the same order

# User-defined function for prediction
predict_ordinal_function <- function(newdata,object_glm){
  X <- as.matrix.data.frame(newdata)
  numY <- length(object_glm$zeta)
  predy <- LL <- UL <- matrix(0,nrow(X),numY)
  for(i in 1:numY){
    C <- matrix(0,nrow(X),ncol(X)+numY)
    C[,1:ncol(X)] <- X
    C[,ncol(X)+i] <- 1
    beta <- as.matrix(c(-object_glm$coefficients,object_glm$zeta))
    link <- as.vector(C%*%beta)
    SE <- sqrt(diag(C%*%vcov(object_glm)%*%t(C)))
    predy[,i] <- 1/(1+exp(-link)) # logit to probability
    LL[,i] <- 1/(1+exp(-link+1.96*SE)) # 95% CI lower limit
    UL[,i] <- 1/(1+exp(-link-1.96*SE)) # 95% CI upper limit
  }
  # final estimates
  list("predy"=predy,"LL"=LL,"UL"=UL) %>% 
    map(~t(apply(.,1,function(x) diff(c(0, x, 1))))) %>% # break down the cumulative prob
    map(~as.data.frame(.) %>% rownames_to_column("rid")) %>% 
    map_df(~pivot_longer(.,cols=-rid),.id="type") %>% 
    pivot_wider(names_from="type") %>% 
    # Caution: ul<->ll flipped in the last category
    rowwise() %>% mutate(ll=min(LL,UL),ul=max(LL,UL)) %>% 
    select(rid,name,predy,ll,ul) %>% ungroup()
}
mypred_cs <- predict_ordinal_function(newdata,ologit_glm2)
mypred_no <- predict_ordinal_function(newdata,olog2)
mypred_no

myfig <- newdata %>% 
  rownames_to_column("rid") %>% 
  left_join(bind_rows(
    mypred_cs %>%
      mutate(model="Complex survey analysis: Yes"),
    mypred_no %>%
      mutate(model="Complex survey analysis: No")),
    by="rid") %>% 
  mutate(name=factor(name,labels=c("Proability\n(happy = 1)",
                                   "Proability\n(happy = 2)",
                                   "Proability\n(happy = 3)",
                                   "Proability\n(happy = 4)",
                                   "Proability\n(happy = 5)")),
         scale_gad=factor(scale_gad,
                          labels=c("Low (25th quantile)",
                                   "Middle (50th quantile)",
                                   "High (75th quantile)")),
         covid_suffer=covid_suffer+svymean(~covid_suffer,cs_design)[1])
myfig %>% 
  ggplot(aes(x=covid_suffer,y=predy,color=scale_gad,fill=scale_gad))+
  geom_line()+
  geom_ribbon(aes(ymin=ll,ymax=ul),alpha=0.2,color=NA)+
  labs(x="Level of financial hardship due to COVID-19\n(1=High, 4=Low)",y="Probability",
       color="GAD (Generalized Anxiety Disorder)",fill="GAD (Generalized Anxiety Disorder)")+
  theme_bw()+
  theme(legend.position="top")+
  facet_grid(model~name)
ggsave("Figure_Part3_Ch6_5_interaction_covidsuffer_gad_order.png",height=12,width=24,units='cm')

###########################################################
# For count variables: Possion/Negative binomial models
## Poisson models
pred_1 <- "ses+achieve+use_smart+scale_gad+I(scale_gad^2)+covid_suffer"
formula1 <- as.formula(str_c("exercise~",pred_1))
pred_2 <- "+scale_gad:covid_suffer+I(scale_gad^2):covid_suffer"
formula2 <- as.formula(str_c("exercise~",pred_1,pred_2))
formula1; formula2

# Scenario 1 (NO)
pois1 <- glm(formula1, mydata_mc, family=quasipoisson)
pois2 <- glm(formula2, mydata_mc, family=quasipoisson)
# Scenario 2 (CS)
count_glm1 <- svyglm(formula1,cs_design_mc,family=quasipoisson)
count_glm2 <- svyglm(formula2,cs_design_mc,family=quasipoisson)

# Compare their test statistics 
summary(pois1)$coefficient[,3] %>% round(4)
summary(count_glm1)$coefficient[,3] %>% round(4)

# Tidy the model outputs
list("m1_no"=pois1,"m1_yes"=count_glm1,
     "m2_no"=pois2,"m2_yes"=count_glm2) %>%
  map_df(wrap_function,.id="model") %>% 
  pivot_wider(names_from=model,values_from=myreport) %>% 
  write_excel_csv("Table_Part3_Ch6_6_regression_count_poisson.csv")

# Model comparison
# lmtest::lrtest(pois1,pois2)  # not available for Quasi-poisson models 
anova(count_glm1,count_glm2,method="Wald")
anova(count_glm1,count_glm2,method="LRT")
c(AIC(pois1),AIC(pois2),BIC(pois1),BIC(pois2)) %>% round(0)
AIC(count_glm1,count_glm2) %>% round(0)
BIC(count_glm1,count_glm2, maximal=count_glm2) %>% round(0)

## Negative binomial models
# Scenario 1 (NO)
negbin1 <- glm.nb(formula1, mydata_mc)
negbin2 <- glm.nb(formula2, mydata_mc)
# Scenario 2 (CS)
count_glmnb1 <- svyglm.nb(formula1,cs_design_mc)
count_glmnb2 <- svyglm.nb(formula2,cs_design_mc)

summary(negbin1)
summary(count_glmnb1)

# User-defined function for negative binomial models
## glm.nb() 
wrap_negbin_function <- function(object_glm){
  sum_table <- tidy(object_glm)
  theta_table <- tibble(term="theta",
                        estimate=object_glm$theta,
                        std.error=object_glm$SE.theta,
                        statistic=estimate/std.error,
                        p.value=2*pnorm(abs(statistic),lower.tail=F))
  bind_rows(sum_table,theta_table) %>% 
    mutate(
      sigstar=cut(p.value,c(-Inf,0.001,0.01,0.05,1),
                  labels=c("***","**","*","")),
      myreport=str_c(format(round(estimate,4),nsmall=4),
                     sigstar,"\n(",
                     format(round(std.error,4),nsmall=4),")")
    ) %>% select(-(estimate:sigstar))
}
## svyglm.nb()
wrap_negbin_cs_function <- function(object_glm){
  sum_theta_table <- tibble(term=names(object_glm$par),
                            estimate=coef(object_glm),
                            std.error=SE(object_glm),
                            statistic=estimate/std.error,
                            p.value=2*pnorm(abs(statistic),lower.tail=F))
  sum_theta_table %>% 
    mutate(
      term=ifelse(str_starts(term,"theta"),"theta",str_remove(term,"eta\\.")),
      sigstar=cut(p.value,c(-Inf,0.001,0.01,0.05,1),
                  labels=c("***","**","*","")),
      myreport=str_c(format(round(estimate,4),nsmall=4),
                     sigstar,"\n(",
                     format(round(std.error,4),nsmall=4),")"),
    ) %>% select(-(estimate:sigstar))
}

# Tidy the model outputs
bind_rows(
  list("m1_no"=negbin1,"m2_no"=negbin2) %>% 
    map_df(wrap_negbin_function,.id="model"),
  list("m1_yes"=count_glmnb1,"m2_yes"=count_glmnb2) %>% 
    map_df(wrap_negbin_cs_function,.id="model")
) %>% 
  pivot_wider(names_from=model,values_from=myreport,names_sort=T) %>% 
  write_excel_csv("Table_Part3_Ch6_7_regression_count_negbin.csv")

# Model comparison
lmtest::lrtest(negbin1,negbin2) 
# anova(count_glmnb1,count_glmnb2,method="Wald")
c(AIC(negbin1),AIC(negbin2),BIC(negbin1),BIC(negbin2)) %>% round(0)
# AIC(count_glmnb1,count_glmnb2) %>% round(0)

tibble(
  source=rownames(summary(negbin1)$coefficients),
  CS_no=summary(negbin1)$coefficients[,3] %>% round(4),
  CS_yes=(coef(count_glmnb1)/SE(count_glmnb1))[-1] %>% round(4)
)

# Visualization
newdata <- mydata %>% 
  data_grid(ses=0,achieve=0,use_smart=1,
            scale_gad=1:4-svymean(~scale_gad,cs_design),
            covid_suffer=0.1*(10:40)-svymean(~covid_suffer,cs_design)) %>% 
  mutate(use_smart=as.factor(use_smart))
# You can use tyep='response' but it's different from 95% CI
## Poisson models
mypred_cs_P <- predict(count_glm2,newdata,se.fit=TRUE) %>% as_tibble()
mypred_no_P <- predict(pois2,newdata,se.fit=TRUE) %>% as_tibble()
## Negative binomial models
mypred_cs_N <- predict(count_glmnb2,newdata,se.fit=TRUE) %>% as_tibble()
mypred_no_N <- predict(negbin2,newdata,se.fit=TRUE) %>% as_tibble()

myfig <- bind_rows(
  newdata %>% mutate(model="Complex survey analysis: Yes,\nPossion model",
                     link=mypred_cs_P$link,SE=mypred_cs_P$SE),
  newdata %>% mutate(model="Complex survey analysis: No,\nPossion model",
                     link=mypred_no_P$fit,SE=mypred_no_P$se.fit),
  newdata %>% mutate(model="Complex survey analysis: Yes,\nNegative binomial model",
                     link=mypred_cs_N$fit,SE=mypred_cs_N$se.fit),
  newdata %>% mutate(model="Complex survey analysis: No,\nNegative binomial model",
                     link=mypred_no_N$fit,SE=mypred_no_N$se.fit)) %>% 
  separate(model,into=c("model","dist"),sep=",") %>% 
  mutate(
    predy=exp(link),  # link to probability 
    ll=exp(link-1.96*SE),
    ul=exp(link+1.96*SE),
    covid_suffer=covid_suffer+svymean(~covid_suffer,cs_design),
    scale_gad=factor(scale_gad,
                     labels=c("1. Very low","2.","3.","4. Very high"))
  ) 
myfig %>% 
  ggplot(aes(x=covid_suffer,y=predy,color=scale_gad,fill=scale_gad))+
  geom_line()+
  geom_ribbon(aes(ymin=ll,ymax=ul),alpha=0.2,color=NA)+
  labs(x="Level of financial hardship due to COVID-19\n(1=High, 4=Low)",y="Weekly workout days",
       color="GAD (Generalized Anxiety Disorder)",fill="GAD (Generalized Anxiety Disorder)")+
  theme_bw()+theme(legend.position="top")+
  facet_grid(dist~model)
ggsave("Figure_Part3_Ch6_6_interaction_covid19_gad.png",height=16,width=16,units='cm')

###########################################################
# For multi-category variables: Multinomial models
# Make a tabulation of 2[happy>3 or not] x 2[exp_sad or not] = 4 groups
mydata_mc_mlogit <- mydata_mc %>% 
  mutate(
    hpp2=ifelse(happy>3,1,0),
    type4=str_c(hpp2,"-",exp_sad),
    type4=factor(type4,labels=c("Unhappy-Nsad","Unhappy-Ysad",
                                "Happy-Nsad","Happy-Ysad"))
  )
# Set model formula 2 including interaction terms
formula2 <- type4~ses+achieve+use_smart+(scale_gad+I(scale_gad^2))*covid_suffer

# Scenario 1 (NO)
mlogit2 <- multinom(formula2,mydata_mc_mlogit,
                    Hess=TRUE, #분산추정을 위한 헤시안행렬 저장
                    trace=FALSE,model=TRUE) #모형추정과정 알림 생략
summary(mlogit2)
# Wrap-up
Table_NO_CS <- wrap_function(mlogit2) %>% 
  pivot_wider(names_from=y.level,values_from=myreport)

# Scenario 2 (CS): Parametric approach
cs_design_mc_mlogit <- cs_design_mc %>% 
  mutate(
    hpp2=ifelse(happy>3,1,0),
    type4=str_c(hpp2,"-",exp_sad),
    type4=factor(type4,labels=c("Unhappy-Nsad","Unhappy-Ysad",
                                "Happy-Nsad","Happy-Ysad"))
  )
mlogit_glm2 <- svy_vglm(formula2,cs_design_mc_mlogit,
                        # keep the reference group unchanged
                        family=multinomial(refLevel="Unhappy-Nsad"))
summary(mlogit_glm2)
# User-defined function for multinomial models
wrap_multinom_cs_function <- function(object_glm){
  summary(object_glm)$coeftable %>% 
    as_tibble(rownames="term") %>% 
    mutate(
      y.level=str_extract(term,"\\d$"),
      term=str_remove(term,":\\d$"),
      sigstar=cut(p,c(-Inf,0.001,0.01,0.05,1),
                  labels=c("***","**","*","")),
      myreport=str_c(format(round(Coef,4),nsmall=4),
                     sigstar,"\n(",
                     format(round(SE,4),nsmall=4),")"),
    ) %>%
    select(y.level,term,myreport) %>% 
    pivot_wider(names_from=y.level,values_from=myreport)
}
Table_YS_CS_YP <- wrap_multinom_cs_function(mlogit_glm2) %>% 
  rename("Unhappy-Ysad"="1","Happy-Nsad"="2","Happy-Ysad"="3")

# Scenario 2 (CS): Nonparametric approach
set.seed(1234)
cs_mlogit_boot <- cs_design_mc_mlogit %>% 
  as_survey_rep(
    type="bootstrap",
    replicates=100) # please increase to >1,000 for practical use
# Set model formula
mlogit_boot2 <- withReplicates(
  cs_mlogit_boot,
  quote(coef(multinom(type4~ses+achieve+use_smart+(scale_gad+I(scale_gad^2))*covid_suffer,
                      weights=.weights,trace=F))),
  return.replicates=TRUE # save the replicates data
)
mlogit_boot2 # you can't print out the object itself...
str(mlogit_boot2) # ...but you can print its structure!
# Resampling-based variance estimation
mlogit_vcov <- vcov(mlogit_boot2)
mlogit_vcov
# Regression coefficients
mycoef <- attr(mlogit_vcov,"means")
# Footnote: identical to the averages of replicates
apply(mlogit_boot2$replicates,2,mean) 
coef(mlogit_boot2) # this is "Theta" so don't be confused

# SE of regression coefficients
mySE <- sqrt(diag(mlogit_vcov))
mySE
# Test statistics
myZ <- mycoef/mySE
# Statistical significance tests
myP <- 2*pnorm(abs(myZ),0,1,lower.tail=FALSE)
mysigstar <- as.character(cut(myP,c(-Inf,0.001,0.01,0.05,1),labels=c("***","**","*"," ")))

# Wrap-up 
wrap_multinom_boot_function <- function(object_boot){
  # Resampling-based variance estimation
  myvcov <- vcov(object_boot)
  # Regression coefficients
  mycoef <- attr(myvcov,"means")
  # SE of regression coefficients
  mySE <- sqrt(diag(myvcov))
  # Test statistics
  myZ <- mycoef/mySE
  # Statistical significance tests
  myP <- 2*pnorm(abs(myZ),0,1,lower.tail=FALSE)
  mysigstar <- as.character(cut(myP,c(-Inf,0.001,0.01,0.05,1),labels=c("***","**","*"," ")))
  # Summary table
  sum_table <- object_boot$theta
  sum_table[] <- str_c(format(round(mycoef,4),nsmall=4),
                       mysigstar,"\n(",
                       format(round(mySE,4),nsmall=4),")")
  # into tibble data
  t(sum_table) %>% as_tibble(rownames="term")
}
Table_YS_CS_NP <- wrap_multinom_boot_function(mlogit_boot2)
  
# Tidy the model outputs
bind_rows(
  Table_NO_CS %>% mutate(model="NO"),
  Table_YS_CS_YP %>% mutate(model="YS_YP"),
  Table_YS_CS_NP %>% mutate(model="YS_NP")
) %>% 
  pivot_wider(names_from=model,values_from=c(-term,-model)) %>% 
  write_excel_csv("Table_Part3_Ch6_8_regression_multinomial_logistic.csv")

# Visualization: Estimate probability for each cell of [Covid_suffer x GAD]
newdata <- mydata %>% 
  data_grid(ses=0,achieve=0,use_smart=1,
            scale_gad=svyquantile(~scale_gad, cs_design, c(0.25,0.5,0.75))$scale_gad[1:3],
            covid_suffer=(1:4)-svymean(~covid_suffer,cs_design)) %>% 
  mutate(use_smart=as.factor(use_smart))
# Built-in predict() function only provides point estimates
mypred <- predict(mlogit2,newdata,type="probs") %>% as_tibble()
mypred

# Get predicted values with 95% CI
newdata <- newdata %>% 
  mutate(
    intercept=1,
    use_smart1=1,
    scale_gad2=scale_gad^2,
    scale_gad_covid_suffer=scale_gad*covid_suffer,
    scale_gad2_covid_suffer=scale_gad^2*covid_suffer
  ) %>% 
  select(intercept,ses,achieve,use_smart1,scale_gad,
         scale_gad2,covid_suffer,scale_gad_covid_suffer,scale_gad2_covid_suffer)
names(newdata);colnames(coef(mlogit2)) # make sure they're in the same order

# User-defined function for prediction
predict_multinom_function <- function(newdata,numY,mycoef,myvcov){
  C <- as.matrix.data.frame(newdata)
  beta <- t(matrix(mycoef,nrow=numY))
  link <- C%*%beta
  Cr <- matrix(rep(t(C),numY),nrow=nrow(C),byrow=T)
  SE <- sqrt(diag(Cr%*%myvcov%*%t(Cr)))
  predy <- cbind(1,exp(link))/(1+rowSums(exp(link))) # link to probability
  LL <- cbind(1,exp(link-1.96*SE))/(1+rowSums(exp(link-1.96*SE))) # 95% CI lower limit
  UL <- cbind(1,exp(link+1.96*SE))/(1+rowSums(exp(link+1.96*SE))) # 95% CI upper limit
  # final estimates
  list("predy"=predy,"LL"=LL,"UL"=UL) %>% 
    map(~as.data.frame(.) %>% rownames_to_column("rid")) %>% 
    map_df(~pivot_longer(.,cols=-rid),.id="type") %>% 
    pivot_wider(names_from="type") %>% 
    # Caution: ul<->ll flipped in the last category
    rowwise() %>% mutate(ll=min(LL,UL),ul=max(LL,UL)) %>% 
    select(rid,name,predy,ll,ul) %>% ungroup()
}
predict_multinom_function(newdata,numY=3,
                          coef(mlogit2),
                          vcov(mlogit2))

myfig <- newdata %>% 
  rownames_to_column("rid") %>% 
  left_join(bind_rows(
    predict_multinom_function(newdata,numY=3,
                              coef(mlogit2),
                              vcov(mlogit2)) %>% 
      mutate(model="Complex survey analysis: No"),
    predict_multinom_function(newdata,numY=3,
                              coef(mlogit_glm2),
                              vcov(mlogit_glm2)) %>% 
      mutate(model="Complex survey analysis: Yes, Parametric"),
    predict_multinom_function(newdata,numY=3,
                              attr(vcov(mlogit_boot2),"means"),
                              vcov(mlogit_boot2)) %>% 
      mutate(model="Complex survey analysis: Yes, Nonparametric")),
    by="rid") %>% 
  mutate(name=factor(name,
                     labels=c("Unhappy (happy=1,2,3)\nUndepressed",
                              "Unhappy (happy=1,2,3)\nDepressed",
                              "Happy (happy=4,5)\nUndepressed",
                              "Happy (happy=4,5)\nDepressed")),
         scale_gad=factor(scale_gad,
                          labels=c("Low (25th quantile)",
                                   "Middle (50th quantile)",
                                   "High (75th quantile)")),
         covid_suffer=covid_suffer+svymean(~covid_suffer,cs_design))
myfig %>% 
  ggplot(aes(x=covid_suffer,y=predy,color=scale_gad,fill=scale_gad))+
  geom_line()+
  geom_ribbon(aes(ymin=ll,ymax=ul),alpha=0.2,color=NA)+
  labs(x="Level of financial hardship due to COVID-19\n(1=High, 4=Low)",y="Proability",
       color="GAD (Generalized Anxiety Disorder)",fill="GAD (Generalized Anxiety Disorder)")+
  theme_bw()+
  theme(legend.position="top")+
  facet_grid(model~name)
ggsave("Figure_Part3_Ch6_7_interaction_covid19_gad_mlogit.png",height=16,width=20,units='cm')

###############################################
# Save the results
# save(mlogit_glm2,mlogit_boot2,
#      file="object_mlogit_CS.RData")
# load("object_mlogit_CS.RData")
