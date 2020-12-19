# First of all, clear workspace
rm(list=ls())
# Save graphical parameters
op <- par()

# Check if whole-brain file exists
tr_file <- "https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Trajectories/connectome.csv"

# If you prefer to generate those values follow the steps within this 'if' command
# (but be aware that some computations may be take a while)
if(0){
  
  # Read phenotypic data
  info <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Phenotypic/phenotypic.csv")
  # Remove those sessions with excesive motion artifact
  info <- info[which(as.logical(info$QC)),]
  # Refactorize subjects' ID
  info$BIDS.ID <- factor(info$BIDS.ID)
  
  # Extract time-series files (those with GSR)
  ts_files <- paste0("https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Brain/GSR/sub-",
                     info$BIDS.ID,"/ses-",info$Session,"/pp150v_wGSR_NIHPD_MNI_P264_ts.txt")
  
  # Compute connectivity matrices
  ConnMx <- sapply(ts_files,
                   function(x){
                     ts <- read.table(x)
                     cmx <- cor(ts)
                     if(sum(is.na(cmx))>0) cmx[is.na(cmx)] <- 0
                     return(cmx)
                   }
  )
  # Reshape object
  nroi <- sqrt(dim(ConnMx)[1])
  ConnMx_dim <- c(nroi, nroi, nrow(info))
  ConnMx <- array(ConnMx, dim = ConnMx_dim)
  
  # Load graph theory formulas
  source("https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Scripts/Auxiliary/inet_w.R")
  # Load cost-threshold functions
  source("https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Scripts/Auxiliary/cost_th_ConnMx.R")
  # Load 'psych' package for Fisher r-to-z transformation
  if(!require(psych)) install.packages("psych"); library(psych)
  
  # Apply thresholds to Connectivity Matrix (25%)
  ConnMx_th <- apply(ConnMx,3,function(x) cost_th_ConnMx(x,0.25))
  ConnMx_th <- array(ConnMx_th,dim=ConnMx_dim)
  # Substract diagonal
  for(ii in 1:ConnMx_dim[3]) diag(ConnMx_th[,,ii]) <- 0
  # Fisher's z transformation
  ConnMx_th <- fisherz(ConnMx_th)
  
  # GRAPH THEORY
  # Weighted-Degree
  w_deg <- apply(ConnMx_th,3:2,sum)
  # Clustering Coefficient (Barrat's formula)
  w_CC <- t(sapply(1:ConnMx_dim[3], function(x) iCC_w(ConnMx_th[,,x])))
  # Global Efficiency (inverse of Dijkstra's formula)
  w_E <- t(sapply(1:ConnMx_dim[3], function(x) iEff_w(ConnMx_th[,,x])))
  # Minimum path-length distances (Dijkstra's formula)
  w_L <- t(sapply(1:ConnMx_dim[3], function(x) iL_w(ConnMx_th[,,x])))
  
  # Compute whole-brain variables
  info$D_net <- rowMeans(w_deg)
  info$C_net <- rowMeans(w_CC)
  info$E_net <- rowMeans(w_E)
  info$L_net <- rowMeans(w_L)
  
  # Extract Functional Network inferences
  atlas <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Atlas/P264/power264.csv")
  
  # Labels
  grph_mat <- c("w_deg","w_CC","w_E","w_L")
  grph_lab <- c("D","C","E","L")
  
  # Apply Functional Network averages
  for(ii in 1:length(grph_mat)){
    fn_mat <- t(aggregate(x = t(get(grph_mat[ii])), by = list(atlas$Net), FUN = mean)[,-1])
    colnames(fn_mat) <- paste0(grph_lab[ii],"_",levels(atlas$Net))
    info <- cbind(info,fn_mat[,-12]) # Remove 'Uncertain' elements
  }
  
  # Save results
  outfile <- "connectome.csv"
  #write.csv(info, outfile, quote = F, row.names = F)
}

###########################################################################
# Whole-brain
###########################################################################

# Plot result
brain_net <- read.csv(tr_file)
# Remove columns with NA's
# ('gammlss' does not allow NA's even if the column is not used)
brain_net <- brain_net[,which(apply(brain_net,2,function(x)sum(is.na(x)))==0)]
# Load 'gamlss', 'gamm4', 'ggplot2' & 'gridExtra' packages
if(!require(gamlss)) install.packages("gamlss"); library(gamlss)
if(!require(gamm4)) install.packages("gamm4"); library(gamm4)
if(!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
if(!require(gridExtra)) install.packages("gridExtra"); library(gridExtra)

# Graph labels
grph_var <- c("D_net", "C_net", "E_net", "L_net")
grph_lab <- c("Degree", "Clustering", "Efficiency", "Path-length")

# Create empty matrix to store results
DF <- matrix(data = as.numeric(NA), nrow = length(grph_var), ncol = 12)
rownames(DF) <- grph_var
colnames(DF) <- c("IPnum","IPval","slope","full_AIC",
                  "PDS_df","PDS_AIC","PDS_LRT","PDS_p",
                  "re_df","re_AIC","re_LRT","re_p")

# Generate scatterplot for each graph modality
for(ii in 1:length(grph_var)){
  
  # Select column
  grp_form <- as.formula(paste0(grph_var[ii], " ~ lo(~PDSe) + FDRMS + Coil + re(random = ~1|BIDS.ID)"))
  # Apply GAMM model
  fit <- gamlss(grp_form, data = brain_net)
  
  # Store GAMLSS results via dropping single terms
  fit_drop <- drop1(fit)
  DF[ii,4] <- fit_drop[1,2] # Full Model AIC
  DF[ii,5:8] <- unlist(fit_drop[2,1:4]) # PDS term
  DF[ii,9:12] <- unlist(fit_drop[5,1:4]) # Random effects term
  
  # Generate residuals
  res_frm <- as.formula(paste0(grph_var[ii], "~ FDRMS + Coil + (1|BIDS.ID)"))
  res_lme <- lmer(res_frm, data = brain_net, REML = F)
  brain_net$res <- scale(residuals(res_lme))
  plt_form <- as.formula("res ~ lo(~PDSe)")
  fit_plt <- getSmo(gamlss(plt_form, data = brain_net))
  brain_net$preds <- predict(fit_plt)
  brain_net$preds_se <- predict(fit_plt, se = T)[[2]]
  # Generate plot
  gSC <- ggplot(data = brain_net, mapping = aes(x = PDSe)) + 
    geom_line(aes(y = res, group = BIDS.ID, color = Sex), size=.25) + 
    geom_point(aes(y = res, color = Sex)) + 
    geom_ribbon(aes(ymax = preds + 1.96 * preds_se,
                    ymin = preds - 1.96 * preds_se,
                    linetype = NA), alpha = 0.2) + 
    geom_line(aes(y = preds), size = 1.5) +
    ylab("Residuals (z)") + xlab("PDS") + ggtitle(grph_lab[ii]) +
    theme_light()
  if(ii == 1) g1 <- gSC else if(ii == 2) g2 <- gSC else if(ii == 3) g3 <- gSC else if(ii == 4) g4 <- gSC
  
  # Store slope and point of inflection
  x <- seq(1,4,0.01)
  y <- predict(fit_plt, newdata = x)
  #plot(x,y, type = "l")
  infl <- c(FALSE, diff(diff(y)>0)!=0)
  DF[ii,1] <- sum(infl)
  DF[ii,2] <- 0
  if(DF[ii,1]>0){
    # Store the last inflection point
    DF[ii,2] <- x[infl][DF[ii,1]]
    # Compute slope after inflection point
    y <- y[x>=DF[ii,2]]
    x <- x[x>=DF[ii,2]]
  }
  # Compute slope
  DF[ii,3] <- coefficients(lm(y~x))[2]

} #for(ii in 1:length(grph_var))

# Arrange plots
g5 <- grid.arrange(g1, g2, g3, g4, nrow = 2, ncol=2)
#ggsave("brain_trajectories.pdf", g5, device = "pdf", width = 8, height = 5)

# Save results
# Compute LRT test for NA
# (drop1 does not calculate LRT-pval when the AIC is higher than full model)
DF[3,9] <- abs(DF[3,9]) # Set d.o.f. as positive value
DF[3,12] <- pchisq(DF[3,11],DF[3,9], lower.tail = F)
outfile <- "whole-brain.csv"
#write.csv(x = DF, file = outfile, quote = F)

# Generate nonlinear trajectories against age for supplementary info
# Generate scatterplot for each graph modality
for(ii in 1:length(grph_var)){
  
  # Select column
  grph_frm <- as.formula(paste0(grph_var[ii], "~ + s(Age, k=4) + FDRMS + Coil"))
  # Apply GAMM model
  fit <- gamm4::gamm4(formula = grph_frm, data = brain_net,
                      random = ~ (1|BIDS.ID), REML = F)
  
  # Generate residuals
  res_frm <- as.formula(paste0(grph_var[ii], "~ FDRMS + Coil + (1|BIDS.ID)"))
  res_lme <- lmer(res_frm, data = brain_net, REML = F)
  brain_net$res <- scale(residuals(res_lme))
  plt_form <- as.formula("res ~ + s(Age, k=4)")
  brain_net$preds <- predict(gamm4::gamm4(plt_form, data = brain_net, random = ~ (1|BIDS.ID), REML = F)$gam)
  brain_net$preds_se <- predict(gamm4::gamm4(plt_form, data = brain_net, random = ~ (1|BIDS.ID), REML = F)$gam,
                                se.fit = T)[[2]]
  # Generate plot
  gSC <- ggplot(data = brain_net, mapping = aes(x = Age)) + 
    geom_line(aes(y = res, group = BIDS.ID, color = Sex), size=.25) + 
    geom_point(aes(y = res, color = Sex)) + 
    geom_ribbon(aes(ymax = preds + 1.96 * preds_se,
                    ymin = preds - 1.96 * preds_se,
                    linetype = NA), alpha = 0.2) + 
    geom_line(aes(y = preds), size = 1.5) +
    ylab("Residuals (z)") + xlab("Age (years)") + ggtitle(grph_lab[ii]) +
    theme_light()
  if(ii == 1) g6 <- gSC else if(ii == 2) g7 <- gSC else if(ii == 3) g8 <- gSC else if(ii == 4) g9 <- gSC
  
} #for(ii in 1:length(grph_var))

# Arrange plots
g10 <- grid.arrange(g6, g7, g8, g9, nrow = 2, ncol=2)
#ggsave("brain_trajectories_age.pdf", g10, device = "pdf", width = 8, height = 5)

###########################################################################
# Functional Networks
###########################################################################

# If you prefer to generate the NIfTI files follow the steps within this 'if' command
if(0){
  # Read ROI and Network labels
  atlas <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Atlas/P264/power264.csv")
  
  # Compute GAMM models each graph theory variable and functional network
  # Find Functional Network columns
  fn_idx <- grep(paste(levels(atlas$Net), collapse = "|"), names(brain_net))
  # Create empty matrix to store results
  DF <- matrix(data = as.numeric(NA), nrow = length(fn_idx), ncol = 14)
  rownames(DF) <- names(brain_net)[fn_idx]
  colnames(DF) <- c("IPnum","IPval","slope","full_AIC",
                    "PDS_df","PDS_AIC","PDS_LRT","PDS_p","PDS_pFDR",
                    "re_df","re_AIC","re_LRT","re_p","re_pFDR")
  
  # For each functional network
  for(nn in 1:length(fn_idx)){
    
    # Select column
    fn_frm <- as.formula(paste0(names(brain_net)[fn_idx[nn]],
                         " ~ lo(~PDSe) + FDRMS + Coil + re(random = ~1|BIDS.ID)"))
    # Apply GAMM model
    fit <- gamlss(fn_frm, data = brain_net)
    
    # Store GAMLSS results via dropping single terms
    fit_drop <- drop1(fit)
    DF[nn,4] <- fit_drop[1,2] # Full Model AIC
    DF[nn,5:8] <- unlist(fit_drop[2,1:4]) # PDS term
    DF[nn,10:13] <- unlist(fit_drop[5,1:4]) # Random effects term
    
    # Extract standardized residuals
    res_frm <- as.formula(paste0(names(brain_net)[fn_idx[nn]],
                                 "~ FDRMS + Coil + (1|BIDS.ID)"))
    res_lme <- lmer(formula = res_frm, data = brain_net, REML = F)
    brain_net$res <- scale(residuals(res_lme))
    plt_form <- as.formula("res ~ lo(~PDSe)")
    fit_plt <- getSmo(gamlss(plt_form, data = brain_net))
    
    # Store slope and point of inflection
    x <- seq(1,4,0.01)
    y <- predict(fit_plt, newdata = x)
    #plot(x,y, type = "l")
    infl <- c(FALSE, diff(diff(y)>0)!=0)
    DF[nn,1] <- sum(infl)
    DF[nn,2] <- 0
    if(DF[nn,1]>0){
      # Store the last inflection point
      DF[nn,2] <- x[infl][DF[nn,1]]
      # Compute slope after inflection point
      y <- y[x>=DF[nn,2]]
      x <- x[x>=DF[nn,2]]
    }
    # Compute slope
    DF[nn,3] <- coefficients(lm(y~x))[2]
    
  }
  
  # Compute manual p-values for low LRT
  DF <- as.data.frame(DF)
  # p-values for PDS term
  for(ii in which(is.na(DF$PDS_p))){
    # Positive degrees of freedom
    DF$PDS_df[ii] <- abs(DF$PDS_df[ii])
    # Compute LRT-pvalue
    DF$PDS_p[ii] <- pchisq(q = DF$PDS_LRT[ii], df = DF$PDS_df[ii], lower.tail = F)
  }
  # p-values for Random-effects term
  for(ii in which(is.na(DF$re_p))){
    # Positive degrees of freedom
    DF$re_df[ii] <- abs(DF$re_df[ii])
    # Compute LRT-pvalue
    DF$re_p[ii] <- pchisq(q = DF$re_LRT[ii], df = DF$re_df[ii], lower.tail = F)
  }
  
  # Compute FDR correction and volume assignation at the Functional Network level
  if(!require(neurobase)) install.packages("neurobase"); library(neurobase)
  # Get the number of functional networks (excluding Uncertain ROIs)
  fn_num <- nlevels(atlas$Net)-1
  for(gg in 1:length(grph_var)){
    # Get p-value indices
    p_idx <- 1:fn_num + (gg-1)*fn_num
    
    # Apply FDR
    DF$PDS_pFDR[p_idx] <- p.adjust(p = DF$PDS_p[p_idx], method = "fdr")
    DF$re_pFDR[p_idx] <- p.adjust(p = DF$re_p[p_idx], method = "fdr")
    
    # Lastly, insert LRT values into a brain volume
    # It is necesary to download this file first
    # https://github.com/BrainMapINB/Pubertal_Connectome/blob/main/Atlas/P264/FuntionalNetworks_P264.nii.gz
    # PDS term
    nii <- readnii("FuntionalNetworks_P264.nii.gz")
    for(nn in 1:fn_num){
      # Avoid 'uncertain' label
      if(nn >= 12){
        nii[nii==nn+1] <- DF$PDS_LRT[p_idx[nn]]
      } else if(nn == 4){
        # Avoid to order DAN before DMN (fullnames Default goes before than Dorsal)
        nii[nii==nn+1] <- DF$PDS_LRT[p_idx[nn]]
      } else if(nn == 5){
        # The opposite for DAN-DMN
        nii[nii==nn-1] <- DF$PDS_LRT[p_idx[nn]]
      } else nii[nii==nn] <- DF$PDS_LRT[p_idx[nn]]
    }
    nii[nii==12] <- 0
    outfile <- paste0("RSN_",grph_lab[gg],"_PDS-LRT")
    #writenii(nim = nii, filename = outfile)
    # Random-effects term
    nii <- readnii("FuntionalNetworks_P264.nii.gz")
    for(nn in 1:fn_num){
      # Avoid 'uncertain' label
      if(nn >= 12){
        nii[nii==nn+1] <- DF$re_LRT[p_idx[nn]]
      } else if(nn == 4){
        # Avoid to order DAN before DMN (fullnames Default goes before than Dorsal)
        nii[nii==nn+1] <- DF$re_LRT[p_idx[nn]]
      } else if(nn == 5){
        # The opposite for DAN-DMN
        nii[nii==nn-1] <- DF$re_LRT[p_idx[nn]]
      } else nii[nii==nn] <- DF$re_LRT[p_idx[nn]]
    }
    nii[nii==12] <- 0
    outfile <- paste0("RSN_",grph_lab[gg],"_re-LRT")
    #writenii(nim = nii, filename = outfile)
  }
  
  # Save results
  outfile <- "functional-networks.csv"
  #write.csv(x = DF, file = outfile, quote = F)
  
  # Heat palette for BrainNet Viewer custom color palette
  aux <- t(round(col2rgb(heat.colors(64))/255,3))
  #write.table(aux, "colpalette_heat64.txt", quote = F, sep = " ", row.names = F, col.names = F)
  aux <- aux[nrow(aux):1,]
  #write.table(aux, "colpalette_heat64_rev.txt", quote = F, sep = " ", row.names = F, col.names = F)
}

# Plot Functional Networks trajectories (FDR q<0.05 corrected)
fn_df <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Trajectories/functional-networks.csv")
fn_num <- nrow(fn_df)/length(grph_lab)
fn_col <- c("white","darkorange","white","purple","red","chartreuse2","white",
            "navy","white","white","deeppink","gold","cyan3")

# Extract by graph theory measure
for(ii in 1:length(grph_lab)){
  
  # Set indices by measure
  fn_idx <- 1:fn_num+fn_num*(ii-1)
  
  # Find significant trends (FDR q<0.05)
  pFDR <- which(fn_df$PDS_pFDR[fn_idx] < 0.05)
  
  # If any FDR plot trends
  if(length(pFDR)>0){
    # Plot every trend
    new_df <- data.frame(PDSe=seq(1,4,0.1))
    gSC <- ggplot(data = new_df, mapping = aes(x = PDSe)) +
      ylab("Residuals (z)") + xlab("PDS") + ggtitle(grph_lab[ii])+
      theme_light()
    # Plot Confidence Intervals first
    for(jj in pFDR){
      # Extract residuals
      res_frm <- as.formula(paste0(fn_df$X[fn_idx[jj]], "~ FDRMS + Coil + (1|BIDS.ID)"))
      res_lme <- lmer(res_frm, data = brain_net, REML = F)
      brain_net$res <- scale(residuals(res_lme))
      plt_form <- as.formula("res ~ lo(~PDSe)")
      fit_plt <- getSmo(gamlss(plt_form, data = brain_net))
      new_df$preds <- predict(fit_plt, newdata = new_df)
      new_df$preds_se <- predict(fit_plt, newdata = new_df, se = T)[[2]]
      # Add lines and CI
      gSC <- gSC + geom_ribbon(data=new_df, aes(ymax = preds + 1.96 * preds_se,
                            ymin = preds - 1.96 * preds_se,
                            linetype = NA), alpha = 0.2, fill=fn_col[jj])
    }
    # Plot fitted line
    for(jj in pFDR){
      # Extract residuals
      res_frm <- as.formula(paste0(fn_df$X[fn_idx[jj]], "~ FDRMS + Coil + (1|BIDS.ID)"))
      res_lme <- lmer(res_frm, data = brain_net, REML = F)
      brain_net$res <- scale(residuals(res_lme))
      plt_form <- as.formula("res ~ lo(~PDSe)")
      fit_plt <- getSmo(gamlss(plt_form, data = brain_net))
      new_df$preds <- predict(fit_plt, newdata = new_df)
      new_df$preds_se <- predict(fit_plt, newdata = new_df, se = T)[[2]]
      # Add lines and CI
      gSC <- gSC + geom_line(data=new_df, aes(y = preds), size = 1.5, color=fn_col[jj])
    }
  }
  if(ii == 1) g11 <- gSC else if(ii == 3) g12 <- gSC else if(ii == 4) g13 <- gSC
}

# Arrange plots
g14 <- grid.arrange(g11, g12, g13, nrow = 1, ncol=3)
#ggsave("fnet_trajectories.pdf", g14, device = "pdf", width = 8, height = 5)
