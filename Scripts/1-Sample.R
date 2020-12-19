# First of all, clear workspace
rm(list=ls())

# Read phenotypic data
info <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Phenotypic/phenotypic.csv")

# Compute inter-session interval
# For session 2
ses2_idx <- which(info$Session==2)
ses2_span <- vector(mode = "numeric", length = length(ses2_idx))
for(ii in 1:length(ses2_idx)){
  # First check if the previous session has the same ID
  if(info$BIDS.ID[ses2_idx[ii]]==info$BIDS.ID[ses2_idx[ii]-1]){
    ses2_span[ii] <- info$Age[ses2_idx[ii]]-info$Age[ses2_idx[ii]-1]
  } 
}
print(paste0("Session 1-2 interval: ",mean(ses2_span),
             " (sd=",sd(ses2_span),")"))
# For session 3
ses3_idx <- which(info$Session==3)
ses3_span <- vector(mode = "numeric", length = length(ses3_idx))
for(ii in 1:length(ses3_idx)){
  # First check if the previous session has the same ID
  if(info$BIDS.ID[ses3_idx[ii]]==info$BIDS.ID[ses3_idx[ii]-1]){
    ses3_span[ii] <- info$Age[ses3_idx[ii]]-info$Age[ses3_idx[ii]-1]
  } 
}
print(paste0("Session 2-3 interval: ",mean(ses3_span),
             " (sd=",sd(ses3_span),")"))

# Compared QC sample vs. noQC sample
qc_sample <- info[which(info$QC==1),]
bq_sample <- info[which(info$QC==0),]
# Refactorize ID
qc_sample$BIDS.ID <- factor(qc_sample$BIDS.ID)
bq_sample$BIDS.ID <- factor(bq_sample$BIDS.ID)
# Removed subjects
match(bq_sample$BIDS.ID,qc_sample$BIDS.ID)
# Descriptive statistics
summary(bq_sample$Age)
table(bq_sample$Sex[which(bq_sample$Session==1)])

# Set color
if(!require(scales)) install.packages("scales"); library(scales)
huecol <- hue_pal()(2)

# Plot histogram of age by sex
# Take initial variables
edad <- info$Age
sexo <- info$Sex
id <- factor(info$BIDS.ID)
tp <- table(table(id))
# Plot histogram
seq1y <- seq(floor(min(edad)), ceiling(max(edad)), 1)
hist(edad, breaks = seq1y, col = huecol[1], axes = F,
     main = "Histogram of age by sex", xlab = "Age (years)")
hist(edad[which(sexo=="M")], breaks = seq1y, col = huecol[2], axes = F, add = T)
seq2y <- seq(floor(min(edad)),ceiling(max(edad)),2)
axis(1,seq1y,F);axis(1,seq2y)
axis(2,seq(0,25,5),las=1)
legend("topright", c("F","M"), pch = 15, col = huecol)
text(17,20,paste("N =",length(table(id)),"\n2tp =",tp[2],"\n3tp =",tp[3]),pos=4)

# Longitudinal scatter-plot
# Sort original variables
info <- info[order(info$Age),]
edad <- info$Age
id <- as.character(info$BIDS.ID)
# Extract labels and frequencies
id_lab <- unique(id)
id_num <- sapply(id_lab,function(x)sum(id==x))
n <- length(id_lab)
sexo <- info$Sex
# Plot canvas
plot(edad, 1:length(edad), type = "n", axes = F,
     ylim = c(0,n), xlim = c(floor(min(edad)), ceiling(max(edad))),
     main = "Longitudinal plot", ylab = "Subject", xlab = "Age (years)")
axis(1, seq1y, F); axis(1,seq2y)
axis(2, seq(0, 10*ceiling(n/10), 10),las=1)
# Draw points
for(ii in 1:n){
  
  # Find sessions
  pos <- which(id==id_lab[ii])
  
  # Extract sex character
  ifelse(sexo[pos[1]] == "F", sexchar <- huecol[1], sexchar <- huecol[2])
  
  # Draw points
  if(length(pos)==1) points(edad[pos], ii, pch = 21, bg = sexchar, cex = 0.8, col = "black") else{
    lines(c(edad[min(pos)],edad[max(pos)]),c(ii,ii), col = sexchar)
    for(jj in pos) points(edad[jj], ii, pch = 21, bg = sexchar, cex = 0.8, col = "black")
  }
}
# Legend
legend(min(edad), n, c("F","M"), pch = 21, pt.bg = huecol, bty = "n")
