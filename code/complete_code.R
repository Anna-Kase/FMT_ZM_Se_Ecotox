
library(tidyverse)
library(janitor)
library(brms)
library(tidybayes)
library(ggplot2)
library(dplyr)
library(tidyr)



sedata <- read.csv("data/se_data.csv")

sedata %>% 
  ggplot(aes((as.factor(Time_Point)), Se, fill=Treatment))+
  geom_boxplot()


sedata$Time_Point <- as.factor(sedata$Time_Point)
sedata$Treatment <- as.factor(sedata$Treatment)
sedata$Turtle <- as.factor(sedata$Turtle)

sedata %>% 
  group_by(Treatment, Time_Point) %>% 
  summarize(
    mean = mean(Se)
  )

sedata <- sedata %>% 
  mutate(days_exp = case_when(Time_Point == "0" ~ "0",
                              TRUE ~ "20"))

sedata$days_exp <- as.factor(sedata$days_exp)


semodel <- brm(Se ~ 1 + Treatment*days_exp +(1|Turtle) ,
               family = Gamma(link="log"),
               data = sedata,
               prior = c(prior(normal(0,10), class="Intercept"),
                         prior(normal(0,1), class="b"),
                         prior(exponential(0.1), class = "sd")),
               cores=4, chains=4, iter=4000)

saveRDS(semodel, file="models/semodel.RDS")

semodel_sensitivity <- brm(Se ~ Treatment*days_exp + (1|Turtle),
                           family = Gamma(link="log"),
                           data = sedata,
                           prior = c(prior(normal(0,20), class="Intercept"),
                                     prior(normal(0,2), class="b"),
                                     prior(exponential(0.01), class="sd")),
                           cores=4, chains=4, iter=4000) 

saveRDS(semodel_sensitivity, file="models/semodel_sensitivity.RDS")



semodel_post <- semodel$data %>% 
  distinct(Treatment, days_exp) %>% 
  add_epred_draws(semodel, re_formula = NA) 




# %>% 
#   mutate(model_type = "Model") 

#%>% 
 # mutate(days_exposure = case_when(Time_Point == "0" ~ "0",
#                                   TRUE ~ "20"))

#semodel_post$days_exposure <- as.factor(semodel_post$days_exposure)

semodel_post$.row = NULL

semodel_post %>% 
  ggplot(aes(days_exp, .epred, fill=Treatment))+
  geom_boxplot(outlier.shape = NA)+
  labs(x="Days of Zebra Mussel Selenium Exposure", 
       y="False Map Turtle Blood Selenium Content (ug/g)")+
  geom_point(data=semodel$data, aes(days_exp, Se, group=Treatment),
             position=position_dodge(width=0.75), alpha=0.25)+
  scale_fill_manual(name="Selenium Concentration", labels = c("0 ug/L", "10 ug/L"),
                    values = c("#44AA99", "#882255"))

ggsave("~/GitHub/FMT_ZM_Se_Ecotox/plots/semodel_boxplot.png")

  
  

   

semodel_post %>% 
  ungroup() %>% 
  select(-.chain, -.iteration) %>%
  distinct() %>% 
  pivot_wider(names_from = Time_Point, values_from = .epred) %>%
  mutate(diff = `4` - `0`) %>% 
  group_by(Treatment) %>% 
  summarize(prob = (sum(diff>0)/length(.draw))*100)

semodel_post %>% 
  ungroup() %>% 
  select(-.chain, -.iteration) %>% 
  distinct() %>% 
  pivot_wider(names_from = Treatment, values_from = .epred) %>% 
  mutate(diff = `0` - `10`) %>% 
  group_by(Time_Point) %>% 
  summarize(prob = (sum(diff>0)/length(.draw))*100)

semodel_post %>% 
  median_qi()




semodel_sens_post <- semodel_sensitivity$data %>% 
  distinct(Treatment, days_exp) %>% 
  add_epred_draws(semodel_sensitivity, re_formula = NA) %>% 
  mutate(model_type = "Sensitivity Analysis")

semodel_sens_post$.row = NULL

semodel_post <- semodel_post %>% 
  mutate(model_type = "Model")

both_post <- rbind(semodel_post, semodel_sens_post)

both_post %>% 
  ggplot(aes(days_exp, .epred, fill=Treatment))+
  geom_boxplot(outlier.shape=NA)+
  facet_wrap(~model_type) +
  labs(x="Days of Zebra Mussel Selenium Exposure", 
       y="False Map Turtle Blood Selenium Content (ug/g)")+
  geom_point(data=semodel$data, aes(days_exp, Se, group=Treatment),
             position=position_dodge(width=0.75), alpha=0.25)+
  scale_fill_manual(name="Selenium Concentration", labels = c("0 ug/L", "10 ug/L"),
                    values = c("#44AA99", "#882255"))
ggsave("~/GitHub/FMT_ZM_Se_Ecotox/plots/semodel_sa_boxplot.png")

  


semodel_restrictive <- brm(Se ~ Treatment*Time_Point + (1|Turtle),
                           family = Gamma(link = "log"),
                           data = sedata,
                           prior = c(prior(normal(0,1), class="Intercept"),
                                     prior(normal(0,1), class="b"),
                                     prior(exponential(0.1), class="sd")),
                           cores=4, chains=4, iter=4000)

rest_post <- semodel_restrictive$data %>% 
  distinct(Treatment, Time_Point) %>% 
  add_epred_draws(semodel_restrictive, re_formula = NA) %>% 
  mutate(model_type = "Restricted")

prior_comp <- rbind(rest_post, semodel_post)

prior_comp %>% 
  ggplot(aes(Time_Point, .epred, fill=Treatment))+
  geom_boxplot()+
  facet_wrap(~model_type)





#Testing out a super restrictive prior (normal(0,0.1)) - this produces
#pretty high estimates, unlike the normal(0,10) and normal(0,1)

semodel_super_restrictive <- brm(Se ~ Treatment*Time_Point + (1|Turtle),
                           family = Gamma(link = "log"),
                           data = sedata,
                           prior = c(prior(normal(0,0.1), class="Intercept"),
                                     prior(normal(0,1), class="b"),
                                     prior(exponential(0.1), class="sd")),
                           cores=4, chains=4, iter=4000)

super_rest_post <- semodel_super_restrictive$data %>% 
  distinct(Treatment, Time_Point) %>% 
  add_epred_draws(semodel_super_restrictive, re_formula = NA) %>% 
  mutate(model_type = "Super Restricted")

prior_three_comp <- rbind(rest_post, semodel_post, super_rest_post)

prior_three_comp %>% 
  ggplot(aes(Time_Point, .epred, fill=Treatment))+
  geom_boxplot()+
  facet_wrap(~model_type)


