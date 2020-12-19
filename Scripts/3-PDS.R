# First of all, clear workspace
rm(list=ls())
# Save graphical parameters
op <- par()

# Read phenotypic data
info <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Phenotypic/phenotypic.csv")

# Load 'gamlss' package for GAMM-LOESS implementation
if(require(gamlss)==0) install.packages("gamlss"); library(gamlss)

# Remove extra columns with NA's
# ('gamlss' does not allow NA's even at non-used columns)
nona <- which(!is.na(info$PDS))
info_tmp <- info[nona,c(1,4,5,9)]
info_tmp$BIDS.ID <- factor(info_tmp$BIDS.ID)
# Compute Age-Sex interaction manually
info_tmp$AgeM <- info_tmp$Age; info_tmp$AgeM[info_tmp$Sex=="F"] <- 0

# Apply GAMM-LOESS
fit <- gamlss(PDS ~ lo(~Age) + lo(~AgeM) + Sex + re(random = ~1|BIDS.ID),
              data = info_tmp)
fit_avg <- gamlss(PDS ~ lo(~Age) + lo(~AgeM) + Sex, data = info_tmp)

# Predict NA's and store fitted lines (and its standard error)
info$PDSe <- info$PDS
info$PDS_pred <- info$PDS_se <- as.numeric(NA)
info$PDS_pred[nona] <- predict(fit_avg)
info$PDS_se[nona] <- predict(getSmo(fit_avg),se=T)[[2]]
# Create newdata for prediction
info_na <- info[-nona,c(1,4,5,9)]
info_na$AgeM <- info_na$Age; info_na$AgeM[info_na$Sex=="F"] <- 0
info_na$BIDS.ID <- factor(info_na$BIDS.ID)
# First, predict the longitudinala ones
info$PDSe[-nona] <- predict(fit, newdata = info_na)
info$PDS_pred[-nona] <- predict(fit_avg, newdata = info_na)
info$PDS_se[-nona] <- predict(getSmo(fit_avg), newdata = info_na, se = T)[[2]]
# Second, predict those with one session
nona <- which(!is.na(info$PDSe))
info$PDSe[-nona] <- predict(fit_avg, newdata = info_na)[which(is.na(predict(fit, newdata = info_na)))]

# Plot PSD vs Age
if(require(ggplot2)==0) install.packages("ggplot2"); library(ggplot2)
g1 <- ggplot(data = info, mapping = aes(x = Age)) + 
  geom_line(aes(y = PDSe, group = BIDS.ID, color = Sex), size=.25) + 
  geom_point(aes(y = PDSe, color = Sex)) + 
  geom_ribbon(aes(group = Sex,
                  ymax = PDS_pred + 1.96 * PDS_se,
                  ymin = PDS_pred - 1.96 * PDS_se,
                  linetype = NA), alpha = 0.2) + 
  geom_line(aes(y = PDS_pred, color = Sex), size = 1.5) + 
  theme_light() +
  xlab("Age (years)") +
  ylab("PDS") +
  ggtitle("Pubertal Development Scale")
show(g1)
#ggsave(filename = "~/R/github/Puberty/Phenotypic/age-PDS.svg", plot = g1, device = "svg")

# Compute LRT for model terms
drop1(fit)
