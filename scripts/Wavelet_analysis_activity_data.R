rm(list = ls())


# Load libraries ----------------------------------------------------------
library(tidyverse)
library(WaveletComp)
library(ggplot2)

# Load libraries
pacman::p_load(tidyverse, WaveletComp, ggplot2)
conflict_prefer("rename", "plyr")

# Load functions
source(file = "./functions/theme_publication.R")

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else  length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                       .fun = function(xx, col) {
                         c(N    = length2(xx[[col]], na.rm=na.rm),
                           mean = mean   (xx[[col]], na.rm=na.rm),
                           sd   = sd     (xx[[col]], na.rm=na.rm)
                         )
                       },
                       measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(as.numeric(datac$N))  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}




# Part 1 ------------------------------------------------------------------

## Format the data ------------------------------------------------------
# let's look at what the data structure looks like
dat <- read.csv("./data/TC5_Activity.csv", 
                header = T, na.strings = c(" ","NA", ""), 
                stringsAsFactors = F)
tail(dat)
glimpse(dat)

## Selecting only the columns that I am interested in & summary of data
dat <- dat %>% select(Day, Phase, ZT, FS:Total)
dat %>% head()
summary(dat) ## We have 172 data points with NAs

# First thing first, how much data do we got?
dat %>% 
  # Keep only hourly data [only keeping the ZT 24h time point]
  filter(ZT %in% c(1:24)) %>%
  # omit all the NAs
  na.omit() %>% 
  # summarize the data to see how many days of data do we have 
  group_by(Day, Phase) %>% 
  arrange(Day, Phase) %>% 
  summarise(datapoints = n())

# Let's clean up the data
dat <-
  dat %>% 
  # Keep only hourly data [only keeping the ZT 24h time point]
  filter(ZT %in% c(0:23)) %>%
  # omit all the NAs
  na.omit() %>% 
  # keep only the days post painting
  # Three phases: Entrain-II (Entrainment post disturbance), Sampling (during a disturbance), and Entrain-III (Entrainment post disturbance as well)
  filter(Phase %in% c("Entrain-I", "Painting" ,"Entrain-II", "Sampling", "Entrain-III"))
# # Two of the other phases: Blasting and Entrain-I
# filter(!Phase %in% c("Painting" ,"Entrain-II", "Sampling", "Entrain-III"))

# Summarise
dat %>% 
  group_by(Day, Phase) %>% 
  arrange(Day, Phase) %>% 
  summarise(datapoints = n())

# Correctly format the columns
dat$Day <- as.factor(dat$Day)
dat$ZT <- as.factor(dat$ZT)
dat$FA <- as.numeric(dat$FA)
dat$FS <- as.numeric(dat$FS)
dat$Total <- as.numeric(dat$Total)
dat$Phase <- factor(dat$Phase, levels=c("Entrain-I","Painting","Entrain-II","Sampling","Entrain-III"))

dat %>% head()
dat %>% glimpse()


# Entrain-I (Initial entrainment) --------------

## The colony was allowed to entrain, undisturbed for four days (Days 4-7),
## following which the colony was disturbed to change the arena setup
## Therefore, we will calculate the mean (±SE) using data from Days 4-7 only


pd <- position_dodge(0.1)

png("./results/figures/figure_1/Entrain_I_FA.png", 
    width = 1400, height = 800, res = 300)
# create summary
dat %>% 
  filter(Phase == "Entrain-I") %>% 
  filter(Day %in% c(4:7)) %>% 
    
  # make the summary table
  summarySE(.,
              # specify your measurevar (FA/FS/Total)
              measurevar= "FA", 
              groupvars=c("ZT","Phase")) %>%
    
    # make the value column
    mutate(value=FA) %>% 
    
  # Plot
  ggplot(aes(x= as.numeric(as.character(ZT)), y=value)) +
    theme_Publication() +
    scale_colour_Publication() +
    xlab("") +
    ylab("") +
    ggtitle("") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none',
          # legend.title = element_text(size = 15, colour = 'black'),
          legend.text = element_blank()) +
    ## center align the title
    theme(plot.title = element_blank()) +
    scale_x_continuous(breaks = c(0,12,24)) +
    #scale_y_continuous(limits = c(0,40)) +
    theme(panel.grid.major.x = element_line(colour = "#808080", size=0.1),
          panel.grid.major.y = element_line(colour = "#808080", size=0.2)) +
    ## if you need highlighting parts of the graph
    geom_rect(aes(xmin = 11.8, xmax = 23.8, ymin = -Inf, ymax = Inf),
              fill = "lightgrey", alpha = 0.02, color=NA) +
    geom_line(position=pd,
              # col="#F2CB05", size=2, alpha=1) + # total
              # col="#BF0404", size=2, alpha=1) + # FS
              col="#0FBF67", size=2, alpha=1) + # FA
              
    
    ## Add error bar here
    geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                  width=.4, position=pd, col="black", alpha = 0.7) +
    
    # Add the points on top of the error bars
    geom_point(position=pd, size=2.5,
               col="black", fill="black",
               show.legend = F, color="black", pch=21, alpha=0.9) +
    
    #facet_wrap(~Phase, nrow=2)
    # scale_fill_manual(values = c("black","#F20505","#F2CB05","#0FBF67")) +
    # scale_color_manual(values=c("#F20505","#F5D736","#0FBF67")) +
    theme(text = element_text(size = 25, colour = 'black'),
          legend.position = "none")
dev.off()


# Painting (Mark and recapture) --------------
pd <- position_dodge(0.1)

png("./results/figures/figure_1/Painting_Total.png", 
    width = 1400, height = 800, res = 300)
# create summary
summarySE(data=(dat %>% 
                  filter(Phase == "Painting")), 
          
          # specify your measurevar (FA/FS/Total)
          measurevar= "Total", 
          groupvars=c("ZT","Phase")) %>%
  
  # make the value column
  mutate(value=Total) %>% 
  
  # Plot
  ggplot(aes(x= as.numeric(as.character(ZT)), y=value)) +
  theme_Publication() +
  scale_colour_Publication() +
  xlab("") +
  ylab("") +
  ggtitle("") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none',
        # legend.title = element_text(size = 15, colour = 'black'),
        legend.text = element_blank()) +
  ## center align the title
  theme(plot.title = element_blank()) +
  scale_x_continuous(breaks = c(0,12,24)) +
  #scale_y_continuous(limits = c(0,40)) +
  theme(panel.grid.major.x = element_line(colour = "#808080", size=0.1),
        panel.grid.major.y = element_line(colour = "#808080", size=0.2)) +
  ## if you need highlighting parts of the graph
  geom_rect(aes(xmin = 11.8, xmax = 23.8, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.02, color=NA) +
  geom_line(position=pd,
            col="#F2CB05", size=2, alpha=1) + # total
            # col="#BF0404", size=2, alpha=1) + # FS
            # col="#0FBF67", size=2, alpha=1) + # FA
  
  
  ## Add error bar here
  geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                width=.4, position=pd, col="black", alpha = 0.7) +
  
  # Add the points on top of the error bars
  geom_point(position=pd, size=2.5,
             col="black", fill="black",
             show.legend = F, color="black", pch=21, alpha=0.9) +
  
  #facet_wrap(~Phase, nrow=2)
  # scale_fill_manual(values = c("black","#F20505","#F2CB05","#0FBF67")) +
  # scale_color_manual(values=c("#F20505","#F5D736","#0FBF67")) +
  theme(text = element_text(size = 25, colour = 'black'),
        legend.position = "none")
dev.off()


# Entrain-II (Continued entrainment/Pre-sampling entrainment) --------------
pd <- position_dodge(0.1)

png("./results/figures/figure_1/Entrain_II_FA.png", 
    width = 1400, height = 800, res = 300)
# create summary
summarySE(data=(dat %>% 
                  filter(Phase == "Entrain-II")), 
          
          # specify your measurevar (FA/FS/Total)
          measurevar= "FA", 
          groupvars=c("ZT","Phase")) %>%
  
  # make the value column
  mutate(value=FA) %>% 
  
  # Plot
  ggplot(aes(x= as.numeric(as.character(ZT)), y=value)) +
  theme_Publication() +
  scale_colour_Publication() +
  xlab("") +
  ylab("") +
  ggtitle("") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none',
        # legend.title = element_text(size = 15, colour = 'black'),
        legend.text = element_blank()) +
  ## center align the title
  theme(plot.title = element_blank()) +
  scale_x_continuous(breaks = c(0,12,24)) +
  #scale_y_continuous(limits = c(0,40)) +
  theme(panel.grid.major.x = element_line(colour = "#808080", size=0.1),
        panel.grid.major.y = element_line(colour = "#808080", size=0.2)) +
  ## if you need highlighting parts of the graph
  geom_rect(aes(xmin = 11.8, xmax = 23.8, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.02, color=NA) +
  geom_line(position=pd,
            # col="#F2CB05", size=2, alpha=1) + # total
            # col="#BF0404", size=2, alpha=1) + # FS
            col="#0FBF67", size=2, alpha=1) + # FA
  
  
  ## Add error bar here
  geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                width=.4, position=pd, col="black", alpha = 0.7) +
  
  # Add the points on top of the error bars
  geom_point(position=pd, size=2.5,
             col="black", fill="black",
             show.legend = F, color="black", pch=21, alpha=0.9) +
  
  #facet_wrap(~Phase, nrow=2)
  # scale_fill_manual(values = c("black","#F20505","#F2CB05","#0FBF67")) +
  # scale_color_manual(values=c("#F20505","#F5D736","#0FBF67")) +
  theme(text = element_text(size = 25, colour = 'black'),
        legend.position = "none")
dev.off()


# Sampling (24h of Sampling) --------------
pd <- position_dodge(0.1)

png("./results/figures/figure_1/Sampling_Total.png", 
    width = 1400, height = 800, res = 300)
# create summary
summarySE(data=(dat %>% 
                  filter(Phase == "Sampling")), 
          
          # specify your measurevar (FA/FS/Total)
          measurevar= "Total", 
          groupvars=c("ZT","Phase")) %>%
  
  # make the value column
  mutate(value=Total) %>% 
  
  # Plot
  ggplot(aes(x= as.numeric(as.character(ZT)), y=value)) +
  theme_Publication() +
  scale_colour_Publication() +
  xlab("") +
  ylab("") +
  ggtitle("") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none',
        # legend.title = element_text(size = 15, colour = 'black'),
        legend.text = element_blank()) +
  ## center align the title
  theme(plot.title = element_blank()) +
  scale_x_continuous(breaks = c(0,12,24)) +
  #scale_y_continuous(limits = c(0,40)) +
  theme(panel.grid.major.x = element_line(colour = "#808080", size=0.1),
        panel.grid.major.y = element_line(colour = "#808080", size=0.2)) +
  ## if you need highlighting parts of the graph
  geom_rect(aes(xmin = 11.8, xmax = 23.8, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.02, color=NA) +
  geom_line(position=pd,
            col="#F2CB05", size=2, alpha=1) + # total
            # col="#BF0404", size=2, alpha=1) + # FS
            # col="#0FBF67", size=2, alpha=1) + # FA
  
  
  ## Add error bar here
  geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                width=.4, position=pd, col="black", alpha = 0.7) +
  
  # Add the points on top of the error bars
  geom_point(position=pd, size=2.5,
             col="black", fill="black",
             show.legend = F, color="black", pch=21, alpha=0.9) +
  
  #facet_wrap(~Phase, nrow=2)
  # scale_fill_manual(values = c("black","#F20505","#F2CB05","#0FBF67")) +
  # scale_color_manual(values=c("#F20505","#F5D736","#0FBF67")) +
  theme(text = element_text(size = 25, colour = 'black'),
        legend.position = "none")
dev.off()


# Entrain-III (Post-sampling entrainment) --------------
pd <- position_dodge(0.1)

png("./results/figures/figure_1/Entrain_III_FA.png", 
    width = 1400, height = 800, res = 300)
# create summary
summarySE(data=(dat %>% 
                  filter(Phase == "Entrain-III")), 
          
          # specify your measurevar (FA/FS/Total)
          measurevar= "FA", 
          groupvars=c("ZT","Phase")) %>%
  
  # make the value column
  mutate(value=FA) %>% 
  
  # Plot
  ggplot(aes(x= as.numeric(as.character(ZT)), y=value)) +
  theme_Publication() +
  scale_colour_Publication() +
  xlab("") +
  ylab("") +
  ggtitle("") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position='none',
        # legend.title = element_text(size = 15, colour = 'black'),
        legend.text = element_blank()) +
  ## center align the title
  theme(plot.title = element_blank()) +
  scale_x_continuous(breaks = c(0,12,24)) +
  #scale_y_continuous(limits = c(0,40)) +
  theme(panel.grid.major.x = element_line(colour = "#808080", size=0.1),
        panel.grid.major.y = element_line(colour = "#808080", size=0.2)) +
  ## if you need highlighting parts of the graph
  geom_rect(aes(xmin = 11.8, xmax = 23.8, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.02, color=NA) +
  geom_line(position=pd,
            # col="#F2CB05", size=2, alpha=1) + # total
            # col="#BF0404", size=2, alpha=1) + # FS
            col="#0FBF67", size=2, alpha=1) + # FA
  
  
  ## Add error bar here
  geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                width=.4, position=pd, col="black", alpha = 0.7) +
  
  # Add the points on top of the error bars
  geom_point(position=pd, size=2.5,
             col="black", fill="black",
             show.legend = F, color="black", pch=21, alpha=0.9) +
  
  # set y-axis breaks
  scale_y_continuous(breaks = seq(0, 14, by = 3)) +
  
  #facet_wrap(~Phase, nrow=2)
  # scale_fill_manual(values = c("black","#F20505","#F2CB05","#0FBF67")) +
  # scale_color_manual(values=c("#F20505","#F5D736","#0FBF67")) +
  theme(text = element_text(size = 25, colour = 'black'),
        legend.position = "none") 
dev.off()



# Calculating number of observations per phase ----------------------------

# Summarise
dat %>% 
  filter(Day != c(8,9)) %>% 
  group_by(Phase) %>% 
  arrange(Phase) %>% 
  summarise(datapoints = n())


# Wavelet analysis --------------------------------------------------------

### Format the data as per Manual (p52) | not required

# tc5.data <- apply(dat[, 4:6], FUN = "diff", MAR = 2)
# tc5.data <- data.frame(Day = dat$Day[-1],
#                        Phase = dat$Phase[-1],
#                        ZT = dat$ZT[-1],
#                        tc5.data)
# head(tc5.data)
# 
# tc5.data <- tc5.data %>% filter(Phase != "Painting")
# tc5.data %>% nrow()

# Total foraging data  --------------------------------------------

set.seed(420)

# Painting - Day 11 and 12; 
### data for Day 11 starts at ZT3
dat.painting <- dat[dat$Phase=="Painting", ]
w.painting <- analyze.wavelet(dat.painting, 
                              ## Analyze which data: FS, FA and Total
                              "Total",
                              loess.span = 0,
                              method = "white.noise",
                              dt = 1, dj = 1/100,
                              lowerPeriod = 1, upperPeriod = 36,
                              make.pval = TRUE, n.sim = 100)

wt.image(w.painting, color.key = "i", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis = list(at = seq(1, 45, by = 3),
                               labels = seq(3,47, by = 3)))

## Reconstruct using all relevant periods
reconstruct(w.painting, plot.waves = F, lwd = c(1,2),
            legend.coords = "bottomright",
            spec.time.axis = list(at = seq(1,45, by=3),
                                  labels = seq(3,47, by=3)))
## Reconstruct using only 24h period
reconstruct(w.painting, sel.period = c(24), plot.waves = F, lwd = c(1,2),
            legend.coords = "bottomleft",
            spec.time.axis = list(at = seq(1,45, by=3),
                                  labels = seq(3,47, by=3)))

# Entrain - II (Undisturbed-1)
my.w <- analyze.wavelet(dat[dat$Phase == "Entrain-II",], 
                        ## Analyze which data: FS, FA and Total
                        "FA",
                        loess.span = 0,
                        method = "white.noise",
                        dt = 1, dj = 1/100,
                        lowerPeriod = 1, upperPeriod = 36,
                        make.pval = TRUE, n.sim = 100)

wt.image(my.w, color.key = "i", n.levels = 250,
         siglvl = 0.01,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis = list(at = seq(1, 96, by = 12),
                               labels = seq(0,95, by = 12)))

reconstruct(my.w, plot.waves = F, lwd = c(1,2),
            legend.coords = "bottomright", show.legend = F,
            spec.time.axis = list(at = seq(1, 96, by = 12),
                                  labels = seq(0,95, by = 12)))

reconstruct(my.w, sel.period = c(8,16,24), plot.waves = F, lwd = c(1,2),
            legend.coords = "bottomleft", show.legend = F,
            spec.time.axis = list(at = seq(1, 96, by = 12),
                                  labels = seq(0,95, by = 12)))

# Sampling (Disturbance-2)
my.w <- analyze.wavelet(dat[dat$Phase == "Sampling",], 
                        ## Analyze which data: FS, FA and Total
                        "FA",
                        loess.span = 0,
                        method = "white.noise",
                        dt = 1, dj = 1/100,
                        lowerPeriod = 1, upperPeriod = 36,
                        make.pval = TRUE, n.sim = 100)

wt.image(my.w, color.key = "i", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis = list(at = seq(1, 96, by = 12),
                               labels = seq(0,95, by = 12)))

reconstruct(my.w, plot.waves = F, lwd = c(1,2),
            legend.coords = "bottomright",
            spec.time.axis = list(at = seq(1, 96, by = 12),
                                  labels = seq(0,95, by = 12)))

reconstruct(my.w, sel.period = c(32), plot.waves = F, lwd = c(1,2),
            legend.coords = "bottomleft",
            spec.time.axis = list(at = seq(1,45, by=3),
                                  labels = seq(3,47, by=3)))


# Entrain-III (Undisturbed-2)
my.w <- analyze.wavelet(dat[dat$Phase == "Entrain-III",], 
                        ## Analyze which data: FS, FA and Total
                        "FA",
                        loess.span = 0,
                        method = "white.noise",
                        dt = 1, dj = 1/100,
                        lowerPeriod = 1, upperPeriod = 36,
                        make.pval = TRUE, n.sim = 100)

wt.image(my.w, color.key = "q", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis = list(at = seq(1, 96, by = 12),
                               labels = seq(0,95, by = 12)))

reconstruct(my.w, plot.waves = F, lwd = c(1,2),
            legend.coords = "bottomright",
            spec.time.axis = list(at = seq(1, 96, by = 12),
                                  labels = seq(0,95, by = 12)))

reconstruct(my.w, sel.period = c(24), plot.waves = F, lwd = c(1,2),
            legend.coords = "bottomleft",
            spec.time.axis = list(at = seq(1, 96, by = 12),
                                  labels = seq(0,95, by = 12)))



# In-depth analysis of undisturbed Cflo colony activity -----------------

### Data to be analyzed: Entrain-II (post-painting, pre-sampling)
tc5.data <- dat[dat$Phase == "Entrain-II",]

## compare FA and Total
my.wx <- analyze.wavelet(tc5.data, "Total", loess.span = 0,
                         method = "white.noise",
                         dt = 1, dj = 1/100,
                         lowerPeriod = 1, upperPeriod = 36,
                         make.pval = TRUE, n.sim = 100)
my.wy <- analyze.wavelet(tc5.data, "FS", loess.span = 0,
                         method = "white.noise",
                         dt = 1, dj = 1/100,
                         lowerPeriod = 1, upperPeriod = 36,
                         make.pval = TRUE, n.sim = 100)
my.wz <- analyze.wavelet(tc5.data, "FA", loess.span = 0,
                         method = "white.noise",
                         dt = 1, dj = 1/100,
                         lowerPeriod = 1, upperPeriod = 36,
                         make.pval = TRUE, n.sim = 100)

# maximum.level = 1.001*max(my.wx$Power.avg, my.wy$Power.avg, my.wz$Power.avg)

# A. Total: all extranidal activities
### Reconstruct the activity profile with period = 24h
reconstruct(my.wx, sel.period = c(24), plot.waves = F, lwd = c(1,2),
            legend.coords = "bottomleft", show.legend = F,
            spec.time.axis = list(at = seq(1, 96, by = 12),
                                  labels = seq(0,95, by = 12)),
            main = "All extranidal activities (period = 24h)")
### Dominant periods
png("./results/wavelet_analysis/Total_dominant_periods.png", width = 1000, height = 1800, res = 300)
wt.avg(my.wx, 
       # maximum.level = maximum.level, 
       spec.period.axis = list(at = c(1,2,4,8,12,24,32), las=1),
       periodtck = 1, periodtcl = NULL,
       lwd = 4, 
       sigcex = c(1,1), exponent = 1,
       show.legend = F, 
       # legend.coords = "bottom",
       main = "Total") 
dev.off()
### Wavelet power (dominant periods with accurate timeline)
png("./results/wavelet_analysis/Total_wavelet_power.png", width = 1400, height = 1000, res = 300)
wt.image(my.wx, color.key = "i", n.levels = 250,
         main = "Total")
dev.off()
# FS: feeding 
### Reconstruct the activity profile with period = 24h
png("./results/wavelet_analysis/FS_wave_recons_12_24.png", width = 1800, height = 1000, res = 300)
reconstruct(my.wy, sel.period = c(12,24), plot.waves = F, lwd = c(1,2),
            legend.coords = "bottomleft", show.legend = F,
            spec.time.axis = list(at = seq(1, 96, by = 12),
                                  labels = seq(0,95, by = 12)),
            main = "Feeding bouts (period = 24h)",
            # only.sig = F, # reconstruct using 12h waveform even for non-sig time-period
            only.coi = F)
dev.off()
### Dominant periods
png("./results/wavelet_analysis/FS_dominant_periods.png", width = 1000, height = 1800, res = 300)
wt.avg(my.wy, 
       # maximum.level = maximum.level, 
       spec.period.axis = list(at = c(1,2,4,8,12,24,32), las=1),
       periodtck = 1, periodtcl = NULL,
       # width of the line
       lwd = 4, 
       # size of the points | c(red, blue)
       sigcex = c(1,1), exponent = 1,
       show.legend = F, 
       # legend.coords = "bottomright",
       main = "FS") 
dev.off()
### Wavelet power (dominant periods with accurate timeline)
png("./results/wavelet_analysis/FS_wavelet_power.png", width = 1400, height = 1000, res = 300)
wt.image(my.wz, color.key = "i", n.levels = 250,
         main="FS")
dev.off()

# FA: non-feeding extranidal activities 
### Reconstruct the activity profile with period = 24h
reconstruct(my.wz, sel.period = c(24), plot.waves = F, lwd = c(1,2),
            legend.coords = "bottomleft", show.legend = F,
            spec.time.axis = list(at = seq(1, 96, by = 12),
                                  labels = seq(0,95, by = 12)),
            main = "Non-feeding extranidal activities (period = 24h)")
### Dominant periods
png("./results/wavelet_analysis/FA_dominant_periods.png", width = 1000, height = 1800, res = 300)
wt.avg(my.wz, 
       # maximum.level = maximum.level, 
       spec.period.axis = list(at = c(1,2,4,8,12,24,32), las=1),
       periodtck = 1, periodtcl = NULL,
       lwd = 4, sigcex = c(1,1), exponent = 1,
       show.legend = F, 
       # legend.coords = "bottomright",
       main = "FA") 
dev.off()

### Wavelet power (dominant periods with accurate timeline)
png("./results/wavelet_analysis/FA_wavelet_power.png", width = 1400, height = 1000, res = 300)
wt.image(my.wx, color.key = "i", n.levels = 250,
         main = "FA")
dev.off()

# Bivariate timeseries analysis -------------------------------------------

# FS over FA
fs.fa <- analyze.coherency(tc5.data, my.pair = c("FS", "FA"),
                           loess.span = 0,
                           dt = 1, dj = 1/100,
                           lowerPeriod = 1, upperPeriod = 36,
                           make.pval = TRUE, n.sim = 100)

# plot of the time-averaged cross-wavelet power (manual, p28)
### Not sure what this means, yet.
wc.avg(fs.fa, exponent = 1,
       siglvl = c(0.1,0.05), show.legend = F,
       sigcol = "red", sigpch = 20,
       periodlab = "period (hours), log2-scale",
       minimum.level = 0, maximum.level = 0.5,
       spec.avg.axis = list(at = seq(0, 1, by = 0.25)),
       spec.period.axis = list(at = c(1,2,4,8,24)),
       periodtck = 1, periodtcl = NULL)

## Cross-wavelet power spectrum | which.image = "wp"
wc.image(fs.fa, n.levels = 250,
         which.image = "wp", # wavelet power
         # p = 0, lvl = 0.3, # defines where to plot arrows | Not sure what p and lvl does, go figure.
         legend.params = list(lab = "cross-wavelet power levels",
                              label.digits = 2),
         siglvl.contour = 0.1, siglvl.arrow = 0.05, which.arrow.sig = "wt",
         color.key = "i",
         # maximum.level = 2.5, exponent = 0.5,
         label.time.axis = TRUE,
         show.date = F, 
         # date.format = "%F %T", date.tz = "EST5EDT",
         # spec.time.axis = list(at = at, labels = labels),
         spec.time.axis = list(at = seq(1, 96, by = 12),
                               labels = seq(0,95, by = 12)),
         # timetcl = -0.5, # outward ticks
         periodlab = "period (hours)",
         timelab = "Time elapsed (hours)",
         spec.period.axis = list(at = c(1,2,4,8,24)),
         periodtck = 1, periodtcl = NULL)

## Recompute coherency again with increased smoothing
fs.fa.2 <- analyze.coherency(tc5.data, my.pair = c("FS", "FA"),
                             loess.span = 0,
                             dt = 1, dj = 1/100,
                             window.type.t = 1, window.type.s = 1,
                             window.size.t = 5, window.size.s = 1,
                             lowerPeriod = 1, upperPeriod = 36,
                             make.pval = TRUE, n.sim = 100)

## Wavelet coherence | which.image = "wc"
wc.image(fs.fa.2, 
         which.image = "wc", 
         color.key = "i", 
         n.levels = 250,
         siglvl.contour = 0.1, siglvl.arrow = 0.05,
         legend.params = list(lab = "wavelet coherence levels",
                              label.digits=2),
         periodlab = "period (hours)",
         timelab = "Time elapsed (hours)",
         spec.time.axis = list(at = seq(1, 96, by = 12),
                               labels = seq(0,95, by = 12)),
         spec.period.axis = list(at = c(1,2,4,8,24)),
         periodtck = 1, periodtcl = NULL)


## Phase difference plot ------------------
wc.sel.phases(fs.fa.2, sel.period = 24,
              only.sig = TRUE,
              which.sig = "wc",
              siglvl = 0.05,
              phaselim = c(-pi,+pi), ## default if legend.horiz = FALSE
              legend.coords = "topright", legend.horiz = FALSE,
              main = "", sub = "", timelab = "")


at.ph <- seq(-pi, pi, by = pi/3)
labels.ph <- seq(-12, 12, by = 4)
# at.t <- seq(from = as.POSIXct("2016-10-31 00:00:00", tz = "EST5EDT"),
#             to = as.POSIXct("2016-11-06 00:00:00", tz = "EST5EDT"), by = "days")
# labels.t <- format(at.t, format = "%a")
png("./results/wavelet_analysis/FS_over_FA_phase_diff.png", width = 2400, height = 1400, res = 300)
wc.sel.phases(fs.fa.2, sel.period = 24, 
              only.sig = T,
              siglvl = 0.1,
              spec.phase.axis = list(at = at.ph, labels = labels.ph),
              timelab = "", 
              phaselab = "Phase (hours)",
              show.legend = F,
              show.date = F, 
              # date.format = "%F %T", date.tz = "EST5EDT",
              # spec.time.axis = list(at = at.t, labels = labels.t),
              spec.time.axis = list(at = seq(1, 96, by = 12),
                                    labels = seq(0,95, by = 12)))
abline(h = 0, col="black", lwd=1)
# abline(v=24, col="lightgrey", lty=2)
dev.off()


# FS vs Total -------------------------------------------------------------
fs.total <- analyze.coherency(tc5.data, my.pair = c("FS", "Total"),
                              loess.span = 0,
                              dt = 1, dj = 1/100,
                              window.type.t = 1, window.type.s = 1,
                              window.size.t = 5, window.size.s = 1,
                              lowerPeriod = 1, upperPeriod = 36,
                              make.pval = TRUE, n.sim = 100)
# time-averaged cross-wavelet power
wc.avg(fs.total, exponent = 1,
       siglvl = c(0.1,0.05), show.legend = F,
       sigcol = "red", sigpch = 20,
       periodlab = "period (hours), log2-scale",
       minimum.level = 0, 
       # maximum.level = 0.5,
       spec.avg.axis = list(at = seq(0, 1, by = 0.25)),
       spec.period.axis = list(at = c(1,2,4,8,24)),
       periodtck = 1, periodtcl = NULL)
## Cross-wavelet power spectrum | which.image = "wp"
wc.image(fs.total, n.levels = 250,
         which.image = "wp", # wavelet power
         # p = 0, lvl = 0.3, # defines where to plot arrows | Not sure what p and lvl does, go figure.
         legend.params = list(lab = "cross-wavelet power levels",
                              label.digits = 2),
         siglvl.contour = 0.1, siglvl.arrow = 0.05, which.arrow.sig = "wt",
         color.key = "i",
         # maximum.level = 2.5, exponent = 0.5,
         label.time.axis = TRUE,
         show.date = F, 
         # date.format = "%F %T", date.tz = "EST5EDT",
         # spec.time.axis = list(at = at, labels = labels),
         spec.time.axis = list(at = seq(1, 96, by = 12),
                               labels = seq(0,95, by = 12)),
         # timetcl = -0.5, # outward ticks
         periodlab = "period (hours)",
         timelab = "Time elapsed (hours)",
         spec.period.axis = list(at = c(1,2,4,8,24)),
         periodtck = 1, periodtcl = NULL)
## Wavelet coherence | which.image = "wc"
wc.image(fs.total, 
         which.image = "wc", 
         color.key = "i", 
         n.levels = 250,
         siglvl.contour = 0.1, siglvl.arrow = 0.05,
         legend.params = list(lab = "wavelet coherence levels",
                              label.digits=2),
         periodlab = "period (hours)",
         timelab = "Time elapsed (hours)",
         spec.time.axis = list(at = seq(1, 96, by = 12),
                               labels = seq(0,95, by = 12)),
         spec.period.axis = list(at = c(1,2,4,8,24)),
         periodtck = 1, periodtcl = NULL)

## Phase difference plot
at.ph <- seq(-pi, pi, by = pi/6)
labels.ph <- seq(-12, 12, by = 2)
# at.t <- seq(from = as.POSIXct("2016-10-31 00:00:00", tz = "EST5EDT"),
#             to = as.POSIXct("2016-11-06 00:00:00", tz = "EST5EDT"), by = "days")
# labels.t <- format(at.t, format = "%a")
wc.sel.phases(fs.total, sel.period = 24, 
              only.sig = T,
              siglvl = 0.1,
              spec.phase.axis = list(at = at.ph, labels = labels.ph),
              timelab = "Time elapsed (hours)", 
              phaselab = "phase (hours)",
              show.legend = F,
              show.date = F, 
              # date.format = "%F %T", date.tz = "EST5EDT",
              # spec.time.axis = list(at = at.t, labels = labels.t),
              spec.time.axis = list(at = seq(1, 96, by = 12),
                                    labels = seq(0,95, by = 12)))
abline(h = 0, col="black", lwd=1)

# FA vs Total -------------------------------------------------------------
fa.total <- analyze.coherency(tc5.data, my.pair = c("FA", "Total"),
                              loess.span = 0,
                              dt = 1, dj = 1/100,
                              window.type.t = 1, window.type.s = 1,
                              window.size.t = 5, window.size.s = 1,
                              lowerPeriod = 1, upperPeriod = 36,
                              make.pval = TRUE, n.sim = 100)
# time-averaged cross-wavelet power
wc.avg(fa.total, exponent = 1,
       siglvl = c(0.1,0.05), show.legend = F,
       sigcol = "red", sigpch = 20,
       periodlab = "period (hours), log2-scale",
       minimum.level = 0, 
       maximum.level = 0.5,
       spec.avg.axis = list(at = seq(0, 1, by = 0.25)),
       spec.period.axis = list(at = c(1,2,4,8,24)),
       periodtck = 1, periodtcl = NULL)
## Cross-wavelet power spectrum | which.image = "wp"
wc.image(fa.total, n.levels = 250,
         which.image = "wp", # wavelet power
         # p = 0, lvl = 0.3, # defines where to plot arrows | Not sure what p and lvl does, go figure.
         legend.params = list(lab = "cross-wavelet power levels",
                              label.digits = 2),
         siglvl.contour = 0.1, siglvl.arrow = 0.05, which.arrow.sig = "wt",
         color.key = "i",
         # maximum.level = 2.5, exponent = 0.5,
         label.time.axis = TRUE,
         show.date = F, 
         # date.format = "%F %T", date.tz = "EST5EDT",
         # spec.time.axis = list(at = at, labels = labels),
         spec.time.axis = list(at = seq(1, 96, by = 12),
                               labels = seq(0,95, by = 12)),
         # timetcl = -0.5, # outward ticks
         periodlab = "period (hours)",
         timelab = "Time elapsed (hours)",
         spec.period.axis = list(at = c(1,2,4,8,24)),
         periodtck = 1, periodtcl = NULL)
## Wavelet coherence | which.image = "wc"
wc.image(fa.total, 
         which.image = "wc", 
         color.key = "i", 
         n.levels = 250,
         siglvl.contour = 0.1, siglvl.arrow = 0.05,
         legend.params = list(lab = "wavelet coherence levels",
                              label.digits=2),
         periodlab = "period (hours)",
         timelab = "Time elapsed (hours)",
         spec.time.axis = list(at = seq(1, 96, by = 12),
                               labels = seq(0,95, by = 12)),
         spec.period.axis = list(at = c(1,2,4,8,24)),
         periodtck = 1, periodtcl = NULL)

## Phase difference plot
at.ph <- seq(-pi, pi, by = pi/6)
labels.ph <- seq(-12, 12, by = 2)
# at.t <- seq(from = as.POSIXct("2016-10-31 00:00:00", tz = "EST5EDT"),
#             to = as.POSIXct("2016-11-06 00:00:00", tz = "EST5EDT"), by = "days")
# labels.t <- format(at.t, format = "%a")
wc.sel.phases(fa.total, sel.period = 24, 
              only.sig = T,
              siglvl = 0.1,
              spec.phase.axis = list(at = at.ph, labels = labels.ph),
              timelab = "Time elapsed (hours)", 
              phaselab = "phase (hours)",
              show.legend = F,
              show.date = F, 
              # date.format = "%F %T", date.tz = "EST5EDT",
              # spec.time.axis = list(at = at.t, labels = labels.t),
              spec.time.axis = list(at = seq(1, 96, by = 12),
                                    labels = seq(0,95, by = 12)))
abline(h = 0, col="black", lwd=1)




# Part 2 -----------------------------------------------------------------

rm(list = ls())
# Load functions
source(file = "./functions/theme_publication.R")

# Orlando photoperiod data ------------------------------------------------
# Load data
orlando.dat <- read.csv("./data/Photoperiod_data_DateAndTime_2019.csv",
                        header = T, na.strings = c(" ","NA", ""), 
                        stringsAsFactors = F)
orlando.dat %>% glimpse()

## let's format the month column appropriately
orlando.dat$month <- factor(orlando.dat$month, levels = c("JAN","FEB","MAR","APR","MAY","JUN",
                                                          "JUL","AUG","SEP","OCT","NOV","DEC"))
## let's format the sunrise and sunset times
library(lubridate)
orlando.dat$sunrise <- as.POSIXct(strptime(orlando.dat$sunrise_time, '%I:%M %p'), format = "%H:%M:%S" )
orlando.dat$sunset <- as.POSIXct(strptime(orlando.dat$sunset_time, '%I:%M %p'), format = "%H:%M:%S" )
orlando.dat$noon <- as.POSIXct(strptime(orlando.dat$solar_noon_time, '%I:%M %p'), format = "%H:%M:%S" )

orlando.dat %>% glimpse()
orlando.dat[orlando.dat$timezone == "EST",]$sunrise <- orlando.dat[orlando.dat$timezone == "EST",]$sunrise + lubridate::hours(1)
orlando.dat[orlando.dat$timezone == "EST",]$sunset <- orlando.dat[orlando.dat$timezone == "EST",]$sunset + lubridate::hours(1)
orlando.dat[orlando.dat$timezone == "EST",]$noon <- orlando.dat[orlando.dat$timezone == "EST",]$noon + lubridate::hours(1)

# # Plot the sunrise-sunset-noon time
# ggplot(orlando.dat) +
#   geom_point(aes(x = 1:nrow(orlando.dat), y = sunrise, color = month)) +
#   geom_point(aes(x = 1:nrow(orlando.dat), y = sunset, color = month)) +
#   geom_point(aes(x = 1:nrow(orlando.dat), y = noon, color = month)) +
#   theme_Publication() +
#   scale_color_viridis_d() +
#   labs(x = "Days of 2019",
#        y = "Sunrise (time)") +
#   # label the x-axis with the last day for each month; Jan ends on 31, Feb on 59, and so on.
#   scale_x_continuous(breaks = c(1,31,59,90,120,151,181,212,
#                                 243,273,304,334,365),
#                      labels = c("Jan","Feb","Mar","Apr","May","Jun",
#                                 "Jul","Aug","Sep","Oct","Nov","Dec","Jan")) +
#   scale_y_time(minor_breaks = "1 Hour") +
#   # ylim(60,120) +
#   geom_vline(xintercept = 120, cex = 1, alpha=0.5) # Experiment starts: May (month begins on day 120)

# # Plot the sunrise-sunset angle (time of sunrise and sunset)
# ggplot(orlando.dat) +
#   geom_point(aes(x = 1:nrow(orlando.dat), y = sunrise_angle, color = month)) +
#   geom_point(aes(x = 1:nrow(orlando.dat), y = sunset_angle, color = month)) +
#   theme_Publication() +
#   scale_color_viridis_d() +
#   labs(x = "Days of 2019",
#        y = "Sunset (angle)") +
#   # label the x-axis with the last day for each month; Jan ends on 31, Feb on 59, and so on.
#   scale_x_continuous(breaks = c(1,31,59,90,120,151,181,212,
#                                 243,273,304,334,365)) +
#   # ylim(60,300) +
#   geom_vline(xintercept = 120, cex = 1, alpha=0.5) # Experiment starts: May (month begins on day 120)

# Plot the daylength (duration of sunrise to sunset); note: twilight not included
### number of hours of sunlight throughout the year
# Format the daylength as follows
orlando.dat$daylength <- strptime(orlando.dat$daylength, format = "%H:%M:%S")
# Make a column that contains the minutes data
orlando.dat$daylength_min <- as.numeric(orlando.dat$daylength$hour*60 + orlando.dat$daylength$min)


# Supp_Fig_1 --------------------------------------------------------------
png("./results/supp_fig_1/orlando_daylength.png", width = 1800, height = 800, res = 300)
ggplot(orlando.dat)+
  geom_point(aes(x = 1:nrow(orlando.dat), 
                 y = daylength_min/60), 
             color=viridis::inferno(100)[75], size = 2.5, alpha=0.3) +
  theme_Publication() +
  scale_color_viridis_d() +
  labs(x = "Days of 2019",
       y = "daylength - sunrise to sunset (hours)") +
  # label the x-axis with the last day for each month; Jan ends on 31, Feb on 59, and so on.
  scale_x_continuous(breaks = c(1,31,59,90,120,151,181,212,
                                243,273,304,334,365),
                     labels = c("Jan","Feb","Mar","Apr","May","Jun",
                                "Jul","Aug","Sep","Oct","Nov","Dec","Jan")) +
  geom_vline(xintercept = 120, cex = 1, alpha=0.5) # Experiment starts: May (month begins on day 120)
dev.off()

# Twilight ----------------------------------------------------------------

orlando.dat %>% glimpse()

# Let's convert the astro_start and astro_end to date-time format
orlando.dat$twi.start <- as.POSIXct(strptime(orlando.dat$astro_start, '%I:%M %p'), format = "%H:%M:%S" )
orlando.dat$twi.end <- as.POSIXct(strptime(orlando.dat$astro_end, '%I:%M %p'), format = "%H:%M:%S" )

orlando.dat[orlando.dat$timezone == "EST",]$twi.start <- orlando.dat[orlando.dat$timezone == "EST",]$twi.start + lubridate::hours(1)
orlando.dat[orlando.dat$timezone == "EST",]$twi.end <- orlando.dat[orlando.dat$timezone == "EST",]$twi.end + lubridate::hours(1)

## Duration between start of twilight to its end
orlando.dat$daylight_secs <- as.numeric(difftime(orlando.dat$twi.end, orlando.dat$twi.start, units = "secs"))

# # Plot the annual twilight start-end times
# ggplot(orlando.dat) +
#   geom_point(aes(x = 1:nrow(orlando.dat), y = twi.start, color = month)) +
#   geom_point(aes(x = 1:nrow(orlando.dat), y = twi.end, color = month)) +
#   geom_point(aes(x = 1:nrow(orlando.dat), y = noon, color = month)) +
#   theme_Publication() +
#   scale_color_viridis_d() +
#   labs(x = "",
#        y = "") +
#   # label the x-axis with the last day for each month; Jan ends on 31, Feb on 59, and so on.
#   scale_x_continuous(breaks = c(1,31,59,90,120,151,181,212,
#                                 243,273,304,334,365)) +
#   # scale_y_time(minor_breaks = "1 Hour") +
#   # ylim(60,120) +
#   geom_vline(xintercept = 120, cex = 1, alpha=0.5) # Experiment starts: May (month begins on day 120)

# # Plot the duration of daylight (beginning of twilight to its end);
# ggplot(orlando.dat)+
#   geom_point(aes(x = 1:nrow(orlando.dat), y = daylight_secs/3600, color = month), size = 2.5) +
#   theme_Publication() +
#   scale_color_viridis_d() +
#   labs(x = "Days of 2019",
#        y = "daylight (hours)") +
#   # label the x-axis with the last day for each month; Jan ends on 31, Feb on 59, and so on.
#   scale_x_continuous(breaks = c(1,31,59,90,120,151,181,212,
#                                 243,273,304,334,365)) +
#   geom_vline(xintercept = 120, cex = 1, alpha=0.5) # Experiment starts: May (month begins on day 120)


### Annual fluctuations in the duration of twilight -------

# # (beginning of twilight-sunrise + sunset-twilight ends);
# ggplot(orlando.dat)+
#   geom_point(aes(x = 1:nrow(orlando.dat), 
#                  y = (daylight_secs/3600) - (daylength_min/60), 
#                  color = month), size = 2.5, alpha=0.5) +
#   theme_Publication() +
#   scale_color_viridis_d() +
#   labs(x = "Days of 2019",
#        y = "twilight (hours)") +
#   # label the x-axis with the last day for each month; Jan ends on 31, Feb on 59, and so on.
#   scale_x_continuous(breaks = c(1,31,59,90,120,151,181,212,
#                                 243,273,304,334,365)) +
#   ylim(2.4,3.2) +
#   geom_vline(xintercept = 120, cex = 1, alpha=0.5) # Experiment starts: May (month begins on day 120)

## Evening twilight only
# (sunset time to end to [astro] twilight)
eve.twi.hrs <- with(orlando.dat, difftime(twi.end, sunset, units = "hours"))
png("./results/supp_fig_1/orlando_eve_twi.png", width = 1800, height = 800, res = 300)
ggplot(orlando.dat)+
  geom_point(aes(x = 1:nrow(orlando.dat), 
                 y = eve.twi.hrs*60),
             color=viridis::inferno(100)[65],
             size = 2.5, alpha=0.5) +
  theme_Publication() +
  scale_color_viridis_d() +
  labs(x = "Days of 2019",
       y = "evening twilight (minutes)") +
  # label the x-axis with the last day for each month; Jan ends on 31, Feb on 59, and so on.
  scale_x_continuous(breaks = c(1,31,59,90,120,151,181,212,
                                243,273,304,334,365),
                     labels = c("Jan","Feb","Mar","Apr","May","Jun",
                                "Jul","Aug","Sep","Oct","Nov","Dec","Jan")) +
  # ylim(1.2,1.6) +
  geom_vline(xintercept = 120, cex = 1, alpha=0.5) # Experiment starts: May (month begins on day 120)
dev.off()


## Let's split twilight into astronomical, nautical and civil twilight

# nautical
orlando.dat$nautical.start <- as.POSIXct(strptime(orlando.dat$nautical_start, '%I:%M %p'), format = "%H:%M:%S" )
orlando.dat$nautical.end <- as.POSIXct(strptime(orlando.dat$nautical_end, '%I:%M %p'), format = "%H:%M:%S" )
# civil
orlando.dat$civil.start <- as.POSIXct(strptime(orlando.dat$civil_start, '%I:%M %p'), format = "%H:%M:%S" )
orlando.dat$civil.end <- as.POSIXct(strptime(orlando.dat$civil_end, '%I:%M %p'), format = "%H:%M:%S" )

# adjust for daylight savings
orlando.dat[orlando.dat$timezone == "EST",]$nautical.start <- orlando.dat[orlando.dat$timezone == "EST",]$nautical.start + lubridate::hours(1)
orlando.dat[orlando.dat$timezone == "EST",]$nautical.end <- orlando.dat[orlando.dat$timezone == "EST",]$nautical.end + lubridate::hours(1)
orlando.dat[orlando.dat$timezone == "EST",]$civil.start <- orlando.dat[orlando.dat$timezone == "EST",]$civil.start + lubridate::hours(1)
orlando.dat[orlando.dat$timezone == "EST",]$civil.end <- orlando.dat[orlando.dat$timezone == "EST",]$civil.end + lubridate::hours(1)

# note, astronomical start and end times are the same as twi.start and twi.end
eve.twi.hrs <- as.numeric(with(orlando.dat, difftime(twi.end, sunset, units = "hours")))
# astronimal twilight = astro.end - nautical.end
eve.astro.hrs <- as.numeric(with(orlando.dat, difftime(twi.end, nautical.end, units = "hours")))
# nautical twilight = nautical.end - civil.end
eve.nautical.hrs <- as.numeric(with(orlando.dat, difftime(nautical.end, civil.end, units = "hours")))
# civil twilight = civil.end - sunset
eve.civil.hrs <- as.numeric(with(orlando.dat, difftime(civil.end, sunset, units = "hours")))

# # Let's plot all three in a plot
# ggplot(orlando.dat)+
#   # geom_point(aes(x = 1:nrow(orlando.dat), y = eve.twi.hrs*60, color = month), size = 1, alpha=0.5) +
#   # astro
#   geom_point(aes(x = 1:nrow(orlando.dat), y = eve.astro.hrs*60), color = viridis::inferno(100)[5], size = 1.5, alpha=0.9) +
#   # nautical
#   geom_point(aes(x = 1:nrow(orlando.dat), y = eve.nautical.hrs*60), color = viridis::inferno(100)[40], size = 1.5, alpha=0.9) +
#   # civil
#   geom_point(aes(x = 1:nrow(orlando.dat), y = eve.civil.hrs*60), color = viridis::inferno(100)[70], size = 1.5, alpha=0.9) +
#   theme_Publication() +
#   scale_color_viridis_d() +
#   labs(x = "Days of 2019",
#        y = "evening twilight (mins)") +
#   # label the x-axis with the last day for each month; Jan ends on 31, Feb on 59, and so on.
#   scale_x_continuous(breaks = c(1,31,59,90,120,151,181,212,
#                                 243,273,304,334,365)) +
#   ylim(22,35.5) +
#   geom_vline(xintercept = 120, cex = 1, alpha=0.5) # Experiment starts: May (month begins on day 120)

# Plot the twilight-sunrise-sunset-noon time
ggplot(orlando.dat) +
  geom_point(aes(x = 1:nrow(orlando.dat), y = sunrise), color = viridis::inferno(100)[85]) +
  geom_point(aes(x = 1:nrow(orlando.dat), y = sunset), color = viridis::inferno(100)[85]) +
  geom_point(aes(x = 1:nrow(orlando.dat), y = noon), color=viridis::inferno(100)[55]) +
  geom_point(aes(x = 1:nrow(orlando.dat), y = twi.start), color = viridis::inferno(100)[5]) +
  geom_point(aes(x = 1:nrow(orlando.dat), y = nautical.start), color = viridis::inferno(100)[40]) +
  geom_point(aes(x = 1:nrow(orlando.dat), y = civil.start), color = viridis::inferno(100)[70]) +
  geom_point(aes(x = 1:nrow(orlando.dat), y = civil.end), color = viridis::inferno(100)[70]) +
  geom_point(aes(x = 1:nrow(orlando.dat), y = nautical.end), color = viridis::inferno(100)[40]) +
  geom_point(aes(x = 1:nrow(orlando.dat), y = twi.end), color = viridis::inferno(100)[5]) +
  
  theme_Publication() +
  scale_color_viridis_d() +
  labs(x = "Days of 2019",
       y = "Time (EDT)") +
  # label the x-axis with the last day for each month; Jan ends on 31, Feb on 59, and so on.
  scale_x_continuous(breaks = c(1,31,59,90,120,151,181,212,
                                243,273,304,334,365),
                     labels = c("Jan","Feb","Mar","Apr","May","Jun",
                                "Jul","Aug","Sep","Oct","Nov","Dec","Jan")) +
  # scale_y_time(minor_breaks = "1 Hour") +
  # ylim("0","23:00") +
  geom_vline(xintercept = 120, cex = 1, alpha=0.5) # Experiment starts: May (month begins on day 120)


# Summary stats -----------------------------------------------------------

## Daylength 
mean(orlando.dat$daylength_min)/60
sd(orlando.dat$daylength_min)/60

## Daylight
mean(orlando.dat$daylight_secs)/3600
sd(orlando.dat$daylight_secs)/3600

## Dusk
mean(eve.twi.hrs)*60
sd(eve.twi.hrs)*60


# Let’s save all the files to a .RData file -------------------------------
# save(my.wx, my.wy, my.wz,
#      fs.fa.2,
#      file = "./results/TC5_wavelet_analysis.RData")

