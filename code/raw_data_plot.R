
library(tidyverse)
library(janitor)
library(brms)
library(tidybayes)

data <- read.csv("data/Raw_Data.csv")
head(data)

data %>% mutate(ZM_Num_Eaten = 20-ZM_Num_Not_Eaten) %>% 
  mutate(TP = case_when(Time_Point == "T0" ~ 0,
                        Time_Point == "T1" ~ 1,
                        Time_Point == "T2" ~ 2,
                        Time_Point == "T3" ~ 3,
                        TRUE ~ 4)) %>%
  subset(TP > 0) %>% 
  mutate(Treatment_Group = as.factor(Treatment_Group)) %>% 
  ggplot(aes(Time_Point, 
             ZM_Num_Eaten,
             group=Turtle,
             color=Treatment_Group))+
  geom_point()+
  geom_line()+
  facet_wrap(~Turtle)+
  labs(x="Experimental Time Point",
       y="Number of Zebra Mussels Consumed")+
  scale_color_manual(name = "Se Treatment", labels = c("0 ug/L", "10 ug/L"), 
                    values=c("#009988","#EE3377"))
