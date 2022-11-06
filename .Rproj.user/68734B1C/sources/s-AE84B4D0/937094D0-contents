

library(tidyverse)
library(janitor)
library(brms)
library(tidybayes)


# average [Se] in mussels = 6 ug/g +/- 1
# average [Se] in scaup blood = ~30 ug/g  +/- 20
# https://www.researchgate.net/profile/Scott-Petrie/publication/49647821_Selenium_Concentrations_in_Greater_Scaup_and_Dreissenid_Mussels_During_Winter_on_Western_Lake_Ontario/links/00b7d5273d524319d2000000/Selenium-Concentrations-in-Greater-Scaup-and-Dreissenid-Mussels-During-Winter-on-Western-Lake-Ontario.pdf

# zm_sim <- sort(round(
#   rnorm(30, 6, 2),
#   digits=2),
#   decreasing = TRUE)
# 
# fmt_sim <- rnorm(30,15,8)
# 
# fmt_sim <- sort(round(
#   rnorm(30, 15, 8),
#   digits=2),
#   decreasing = TRUE)
# 
# sim_data <- as.data.frame(cbind(zm_sim, fmt_sim))


fmt_0_sim <- rnorm(15, 4, 2)
fmt_10_sim <- rnorm(15, 15, 8)

sim_data <- as.data.frame(cbind(fmt_0_sim,fmt_10_sim))

sim_data <- sim_data %>% 
  gather("treatment", "concentration")

sim_data %>% 
  ggplot(aes(treatment, concentration, fill=treatment))+
  geom_boxplot()+
  scale_x_discrete(labels=c("fmt_0_sim" = "0 ug/L", "fmt_10_sim" = "10 ug/L"))+
  labs(x= "Selenium Treatment Group", y="Simulated False Map Turtle Selenium Blood Concentration (ug/g)")+
  theme_bw()+
  scale_fill_manual(values=c("#009988","#EE3377"))+
  theme(legend.position = "none")


se_sim <- brm(concentration ~ treatment,
              family = Gamma(link="log"),
              data = sim_data,
              prior = c(prior(normal(4,2), class="Intercept"),
                        prior(normal(0,1), class="b")),
              cores=4, chains=1, iter=1000)

conditional_effects(se_sim)

post <- se_sim$data %>% 
  distinct(treatment) %>% 
  add_epred_draws(se_sim)

post

post$.row = NULL

post2 <- post %>% 
  select(-.chain, -.iteration) %>%
  distinct() %>% 
  pivot_wider(names_from = "treatment", values_from = ".epred") 

se_diff <- post2 %>% 
  mutate(diff=fmt_0_sim-fmt_10_sim) 

se_diff %>% 
  summarize(higher=sum(diff>0)/nrow(.))
#>99% probability that the 10 ug/L treatment group has higher blood selenium concentrations
# than the control group

post2 %>% median_qi()


tp0_0 <- rnorm(15, 3, 1)
tp0_10 <- rnorm(15, 3, 1)

tp1_0 <- rnorm(15, 4, 1)
tp1_10 <- rnorm(15, 6, 2)

tp2_0 <- rnorm(15, 4, 2)
tp2_10 <- rnorm(15, 8, 3)

tp3_0 <- rnorm(15, 5, 2)
tp3_10 <- rnorm(15, 12, 5)

tp4_0 <- rnorm(15, 5, 2)
tp4_10 <- rnorm(15, 15, 8)

sim_time <- as.data.frame(cbind(tp0_0, tp0_10, tp1_0, tp1_10, tp2_0, tp2_10, tp3_0, tp3_10, tp4_0, tp4_10))

sim_time <- sim_time %>% 
  gather("treatment", "concentration") %>% 
  mutate(trt = case_when(treatment == "tp0_0" ~ 0,
                         treatment == "tp0_10" ~ 10,
                         treatment == "tp1_0" ~ 0,
                         treatment == "tp1_10" ~ 10,
                         treatment == "tp2_0" ~ 0,
                         treatment == "tp2_10" ~ 10,
                         treatment == "tp3_0" ~ 0,
                         treatment == "tp3_10" ~ 10,
                         treatment == "tp4_0" ~ 0,
                         treatment == "tp4_10" ~ 10)) %>%
  mutate(tp = case_when(treatment == "tp0_0" ~ 0,
                         treatment == "tp0_10" ~ 0,
                         treatment == "tp1_0" ~ 1,
                         treatment == "tp1_10" ~ 1,
                         treatment == "tp2_0" ~ 2,
                         treatment == "tp2_10" ~ 2,
                         treatment == "tp3_0" ~ 3,
                         treatment == "tp3_10" ~ 3,
                         treatment == "tp4_0" ~ 4,
                         treatment == "tp4_10" ~ 4))

sim_time$tp <- as.factor(sim_time$tp)
sim_time$trt <- as.factor(sim_time$trt)

sim_time %>% 
  ggplot(aes(tp, concentration, fill=trt))+
  geom_boxplot(outlier.shape = NA)+labs(x= "Time Point", y="Simulated False Map Turtle Selenium Blood Concentration (ug/g)")+
  theme_bw()+
  scale_fill_manual(name = "Se Treatment", labels = c("0 ug/L", "10 ug/L"), 
                      values=c("#009988","#EE3377"))


time_sim <- brm(concentration ~ tp*trt,
              family = Gamma(link="log"),
              data = sim_time,
              prior = c(prior(normal(4,2), class="Intercept"),
                        prior(normal(0,1), class="b")),
              cores=4, chains=1, iter=1000)

conditional_effects(time_sim)

post <- time_sim$data %>% 
  distinct(tp, trt) %>% 
  add_epred_draws(time_sim)

post

post$.row = NULL

post2 <- post %>% 
  select(-.chain, -.iteration) %>%
  distinct() %>% 
  pivot_wider(names_from = "trt", values_from = ".epred") %>% 
  pivot_wider(id_cols=c(tp, .draw),
              names_from=tp, 
              values_from=c("0", "10"))

tp_diff <- post2 %>% 
  mutate(tp0_diff=`0_0`-`10_0`) %>%
  mutate(tp1_diff=`0_1`-`10_1`) %>%
  mutate(tp2_diff=`0_2`-`10_2`) %>%
  mutate(tp3_diff=`0_3`-`10_3`) %>%
  mutate(tp4_diff=`0_4`-`10_4`)
  

tp_diff %>% 
  summarize(higher=sum(tp0_diff>0)/nrow(.))
# 73.2% probability that the 10 ug/L treatment group has higher blood selenium concentrations
# than the control group

tp_diff %>% 
  summarize(higher=sum(tp1_diff>0)/nrow(.))
# 99.8% probability that the 10 ug/L treatment group has higher blood selenium concentrations
# than the control group

tp_diff %>% 
  summarize(higher=sum(tp2_diff>0)/nrow(.))
# 99.2% probability that the 10 ug/L treatment group has higher blood selenium concentrations
# than the control group

tp_diff %>% 
  summarize(higher=sum(tp3_diff>0)/nrow(.))
# >99.9% probability that the 10 ug/L treatment group has higher blood selenium concentrations
# than the control group

tp_diff %>% 
  summarize(higher=sum(tp4_diff>0)/nrow(.))
# >99.9% probability that the 10 ug/L treatment group has higher blood selenium concentrations
# than the control group

post2 %>% median_qi()


post3 <- post %>% 
  select(-.chain, -.iteration) %>%
  distinct() %>% 
  pivot_wider(names_from = "tp", values_from = ".epred")

tm_diff <- post3 %>% 
  mutate(tp01_diff=`0`-`1`) %>%
  mutate(tp02_diff=`0`-`2`) %>%
  mutate(tp03_diff=`0`-`3`) %>%
  mutate(tp04_diff=`0`-`4`) 


tm_diff %>% 
  summarize(higher=sum(tp01_diff>0)/nrow(.))

tm_diff %>% 
  summarize(higher=sum(tp02_diff>0)/nrow(.))

tm_diff %>% 
  summarize(higher=sum(tp03_diff>0)/nrow(.))

tm_diff %>% 
  summarize(higher=sum(tp04_diff>0)/nrow(.))
