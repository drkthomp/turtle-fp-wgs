library(tidyverse)

stats <- read.delim("~/2021-REU/somatypus/input/bwa_flSCBLdna_CheMyd_sorted.coverage", stringsAsFactors=TRUE) %>% mutate(okay = case_when(meandepth < 80 ~ FALSE, TRUE ~ TRUE))


stats %>% summarise(mean = mean(meandepth))
ggplot(stats,aes(x=1,y=meandepth)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color=okay),width=0.2) +
  labs(title="flSCBLdna mean coverage depth by scaffold",color="Above 80x?",x="") +
  theme_minimal() +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  scale_y_continuous(breaks=c(0,20,40,60,80,100,200,300,400))

(stats %>% filter(meandepth > 80))$X.rname

primary_assembly <- (stats %>% filter(grepl("NC_05",X.rname,  fixed = TRUE)))

mean(primary_assembly$meandepth)

ggplot(primary_assembly,aes(x=1,y=meandepth)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(color=okay),width=0.2) +
  labs(title="flSCBLdna mean coverage depth by chromosome",color="Above 80x?",x="") +
  theme_minimal() +
  scale_x_discrete(labels = NULL, breaks = NULL)
