####################################################################################
## Part 2. Analysis: Descriptive statistics 
## Chpt 3. Complex survey data and descriptive statistics 
####################################################################################

# Load libraries
library(tidyverse)  # Data manipulation 
library(survey)     # Complex survey analysis: Base-R approach
library(srvyr)      # Complex survey analysis: Tidyverse approach 

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

# Design variables
# 1) CLUSTER: cluster
# 2) STRATA: stratification variable
# 3) W: weight variable
# 4) FPC: finite population correction variable

# Complex survey analysis: Tidyverse approach　
cs_design <- mydata %>% 
  as_survey_design(ids=cluster,     # cluster 
                   strata=strata,   # strata
                   weights=w,       # weight 
                   fpc=fpc)         # FPC 
cs_design
summary(cs_design)

# Complex survey analysis: Base-R approach
cs_design_classic <- svydesign(ids=~cluster,  # cluster 
                               strata=~strata,# strata 
                               weights=~w,    # weight 
                               fpc=~fpc,      # FPC 
                               data=mydata)
cs_design_classic

# Footnote: What does FPC variable mean?
# Cluster(=class) size
mydata %>% group_by(cluster) %>%
  summarize(n_fpc=n_distinct(fpc)) %>% count(n_fpc)

# Descriptive statistics: Categorical variables
## Tidyverse approach
# Scenario 1: Ignore the survey design (NO)
mydata %>% count(female)
mydata %>% count(female) %>% mutate(prop=n/sum(n))
# Scenario 2: Apply the survey design (CS)
cs_design %>% survey_count(female)
cs_design %>% survey_count(female,vartype="ci")
# Compare their proportions
ftable_n <- mydata %>% count(female) %>% mutate(prop=n/sum(n))
ftable_y <- cs_design %>% survey_count(female) %>% mutate(prop=n/sum(n))
ftable_n$prop
ftable_y$prop

# Footnote: Weighted means?
mydata %>% group_by(female) %>%
  summarize(n=sum(w)) %>% mutate(prop=n/sum(n))
nrow(mydata); sum(mydata$w)

# Repeat the procedure to all categorical variables
Fvars <- mydata %>% select_if(is.factor) %>% names(.)
Fvars

# Scenario 1 (NO)
Ftable_NO <- Fvars %>% 
  map_df(~mydata %>% mutate(group=!!sym(.x)) %>% 
           count(group)  %>%
           mutate(variable=.x,prop=n/sum(n)))
Ftable_NO

# Scenario 2 (CS)
Ftable_CS <- Fvars %>% 
  map_df(~cs_design %>% mutate(group=!!sym(.x)) %>% 
           survey_count(group) %>%
           mutate(variable=.x,prop=n/sum(n)))
# Compare their frequencies
Ftable <- bind_rows(Ftable_CS %>% mutate(CS="YS"),
                    Ftable_NO %>% mutate(CS="NO")) %>% 
  mutate(myreport=str_c(round(n),
                        "\n(",
                        format(round(prop,4),nsmall=4),")")) %>% 
  select(variable,group,CS,myreport) %>% 
  pivot_wider(names_from="CS",values_from="myreport")
Ftable %>% write_excel_csv("Table_Part2_Ch3_1_descriptive_freq.csv")

## Base-R approach
svytable(~sch_mdhg, cs_design_classic) # frequencies 
prop.table(svytable(~sch_mdhg, cs_design_classic)) # proportions 

# Visualization: Distribution of school types
bind_rows(
  mydata %>% count(sch_mdhg) %>%
    mutate(prop=n/sum(n),CS="No"),
  cs_design %>% survey_count(sch_mdhg) %>% select(-n_se) %>%
    mutate(prop=n/sum(n),CS="Yes")) %>% 
  ggplot(aes(x=sch_mdhg, y=prop, fill=CS))+
  geom_bar(stat='identity',position=position_dodge(width=0.8),alpha=0.7)+
  scale_fill_manual(values=c("grey70","grey10"))+
  labs(x="Middle/High school",y="Proportion",fill="Complex survey analysis")+
  coord_cartesian(ylim=c(0.4,0.55))+
  theme_bw()
ggsave("Figure_Part2_Ch3_1_descriptive_prop_school.png",width=12,height=8,units='cm')

# Descriptive statistics: Continuous variables
# Scenario 1 (NO)
mydata %>% 
  summarise(M=mean(happy))
lm(happy~1, mydata) %>% confint()

# Scenario 2 (CS)
cs_design %>% 
  summarise(happy=survey_mean(happy,vartype='ci'))
# 90% CI also available
cs_design %>% 
  summarise(happy=survey_mean(happy,vartype=c('se','ci','var','cv')),level=.90)
# Design effect
cs_design %>% 
  summarise(happy=survey_mean(happy,vartype='ci',deff=TRUE))

# Calculate their means and 95% CIs
Dvars <- mydata %>% select(female:happy) %>%
  select_if(is.double) %>% names(.)
descriptives <- cs_design %>% 
  summarise(across(
    .cols=all_of(Dvars),
    .fns=function(x){survey_mean(x,vartype='ci',
                                 deff=TRUE,na.rm=TRUE)}
  ))
descriptives

# Design effect
descriptives %>% 
  select(contains("_deff")) %>% 
  round(4)

# Mean, 95% CI lower/upper limit
descriptives %>% 
  select(-contains("_deff"),
         -contains("_low"),
         -contains("_upp")) %>% 
  round(4)
descriptives %>% 
  select(contains("_low"),
         contains("_upp")) %>% 
  round(4)

# Summarize into a table 
# Scenario 1 (NO)
Dtable_NO <- tibble(variable=Dvars) %>% 
  mutate(fit=map(variable,~lm(str_c(.x,"~1"),mydata)),
         M=map_dbl(fit,coef),
         LLUL=map(fit,confint)) %>% 
  unnest(LLUL) %>% 
  mutate(myreport=str_c(format(round(M,4),nsmall=4),
                        "\n(",
                        format(round(LLUL[,1],4),nsmall=4),
                        ", ",
                        format(round(LLUL[,2],4),nsmall=4),
                        ")")) %>% 
  select(variable,myreport)
# Scenario 2 (CS)
Dtable_CS <- descriptives %>% 
  select(!contains("_deff")) %>% 
  pivot_longer(cols=everything()) %>% 
  mutate(
    type=ifelse(str_detect(name,"_low"),"ll",
                ifelse(str_detect(name,"_upp"),"ul","mn")),
    variable=str_remove(name,"_low|_upp")
  ) %>% select(-name) %>% 
  pivot_wider(names_from="type",values_from="value") %>% 
  mutate(
    myreport=
      str_c(format(round(mn,4),nsmall=4),
            "\n(",
            format(round(ll,4),nsmall=4),
            ", ",
            format(round(ul,4),nsmall=4),
            ")")
  ) %>% select(variable,myreport) 
# Compare two scenarios
Dtable <- full_join(Dtable_NO %>% rename(no=myreport),
                    Dtable_CS %>% rename(yes=myreport))
Dtable %>% write_excel_csv("Table_Part2_Ch3_2_descriptive_M_CI.csv")

# Visualize the table
myfig <- Dtable %>%
  pivot_longer(-variable) %>% 
  mutate(value=str_remove_all(value,"\\(|\\)")) %>% 
  separate(value,into=c("mean","llul"),sep="\n") %>% 
  separate(llul,into=c("ll","ul"),sep=", ") %>% 
  arrange(name) %>% 
  mutate(across(mean:ul,as.double),
         name=ifelse(name=="no","No","Yes"),
         variable=fct_reorder(variable, row_number()))
myfig %>% 
  ggplot(aes(x=name,y=mean,color=name))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=ll,ymax=ul),width=0.2)+
  facet_wrap(~variable,scales="free")+
  coord_flip()+
  labs(x="",y="Mean with 95% CI")+
  theme_bw()+
  guides(color="none")
ggsave("Figure_Part2_Ch3_2_descriptive_M_CI.png",width=25,height=22,units='cm')

## Base-R approach
svymean(~happy, cs_design_classic, deff=TRUE) 
# 95% CI 
confint(
  svyglm(happy~1, cs_design_classic)
)

# # Visualization: Age distribution
## Base-R approach
png("Figure_Part2_Ch3_3_descriptive_agehist.png",width=960)
par(mfrow=c(1,2))
hist(mydata$agem,prob=TRUE,main="Conventional data analysis",
     col="red",xlab="Age",ylab="Proportion")
svyhist(~agem,cs_design,main="Complex survey data analysis",
        col="blue",xlab="Age",ylab="Proportion")
par(mfrow=c(1,1))
dev.off()

## Tidyverse approach
hist_no <- hist(mydata$agem,prob=TRUE) # Scenario 1 
hist_cs <- svyhist(~agem,cs_design) # Scenario 2
tibble(mids=hist_no$mids) %>% 
  ggplot(aes(x=mids))+
  geom_bar(aes(y=hist_cs$density,fill="Yes"),stat="identity",width=0.5,alpha=0.3)+
  geom_bar(aes(y=hist_no$density,fill="No"),stat="identity",width=0.5,alpha=0.3)+
  theme_bw()+
  scale_fill_manual(values=c("red","blue"))+
  labs(x="Age",y="Proportion",fill="Complex survey analysis")+
  theme(legend.position="top")
ggsave("Figure_Part2_Ch3_4_descriptive_agehist_ggplot.png",width=12,height=10,units='cm')

# Subpopulation analysis: Categorical variables
## Tidyverse approach
# Scenario 1 (NO)
prop_no <- mydata %>% 
  group_by(female,sch_type) %>% 
  count(exp_sad) %>% 
  mutate(prop=n/sum(n))
prop_no
# Scenario 2 (CS)
## Tidyverse approach
prop_yes <- cs_design %>% 
  group_by(female,sch_type) %>% 
  survey_count(exp_sad) %>% 
  mutate(prop=n/sum(n))
prop_yes

# Footnote: do NOT use filtered data
mydata %>% 
  filter(female==0) %>% # this is wrong!
  as_survey_design(ids=cluster,     # cluster 
                   strata=strata,   # strata
                   weights=w,       # weight 
                   fpc=fpc)         # FPC 

# Visualization
bind_rows(
  prop_no %>% mutate(cs="No"),
  prop_yes %>% select(-n_se) %>% mutate(cs="Yes")
) %>% 
  mutate(
    exp_sad=ifelse(exp_sad==1,"With sadness","Without sadness"),
    gender=ifelse(female==0,"Boys","Girls"),
    subpop=str_c(sch_type,",\n", gender)
  ) %>% 
  ggplot(aes(x=subpop,y=prop,fill=cs))+
  geom_bar(stat='identity',position=position_dodge(width=0.8),alpha=0.7)+
  scale_fill_manual(values=c("grey70","grey10"))+
  labs(x="Subpopulation",y="Proportion",fill="Complex survey analysis")+
  coord_cartesian(ylim=c(0.1,0.9))+
  theme_bw()+theme(legend.position="top")+
  facet_grid(~exp_sad)
ggsave("Figure_Part2_Ch3_5_descriptive_prop_subpop.png",width=20,height=12,units='cm')

# Subpopulation analysis: Continuous variables
## Tidyverse approach
# Scenario 1 (NO)
mci_no <- mydata %>% 
  group_by(female,sch_type) %>% 
  summarize(mean=mean(happy),
            ll=(lm(happy~1,cur_data()) %>% confint())[1],
            ul=(lm(happy~1,cur_data()) %>% confint())[2])
mci_no
# Scenario 2 (CS)
mci_yes <- cs_design %>% 
  group_by(female,sch_type) %>% 
  summarise(mean=survey_mean(happy,vartype='ci')) %>% 
  rename(ll=mean_low,ul=mean_upp)
mci_yes

# Visualization
bind_rows(
  mci_no %>% mutate(cs="No"),
  mci_yes %>% mutate(cs="Yes")
) %>% 
  mutate(
    gender=ifelse(female==0,"Boys","Girls"),
    subpop=str_c(sch_type,",\n", gender)
  ) %>% 
  ggplot(aes(x=subpop,y=mean,color=cs,shape=cs))+
  geom_point(size=2,position=position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=ll,ymax=ul),width=0.1,
                position=position_dodge(width=0.3))+
  labs(x="Subpopulation",y="Perceived Happiness (Mean with 95% CI)",
       shape="Complex survey analysis",color="Complex survey analysis")+
  coord_cartesian(ylim=c(2.0,2.4))+
  theme_bw()+theme(legend.position="top")
ggsave("Figure_Part2_Ch3_6_descriptive_M_CI_subpop.png",width=12,height=12,units='cm')

# Subpopulation analysis: Categorical variables
## Base-R approach
svytable(~exp_sad+sch_type+female, cs_design_classic)
prop.table(
  svytable(~exp_sad+female, subset(cs_design_classic,sch_type=="Co-ed school")),2
)
prop.table(
  svytable(~exp_sad+female, subset(cs_design_classic,sch_type!="Co-ed school")),2
)

png("Figure_Part2_Ch3_7_descriptive_subpop_fourfoldplots.png",width=720)
mytable_mix <- svytable(~exp_sad+female, subset(cs_design_classic,sch_type=="Co-ed school"))
mytable_only <- svytable(~exp_sad+female, subset(cs_design_classic,sch_type!="Co-ed school"))
colnames(mytable_mix) <- colnames(mytable_only) <- c("Boys","Girls")
rownames(mytable_mix) <- rownames(mytable_only) <- c("No","Yes")
par(mfrow=c(1,2)) # use left and right 
fourfoldplot(round(mytable_mix),main="Coeducation")
fourfoldplot(round(mytable_only),main="All-boys or All-girls")
par(mfrow=c(1,1))
dev.off()

# Subpopulation analysis: Continuous variables
## Base-R approach
svyby(~happy, by=~female+sch_type, 
      cs_design_classic, svymean) # mean
confint(
  svyby(~happy, by=~female+sch_type, 
        cs_design_classic, svymean) # # 95% CI
)
