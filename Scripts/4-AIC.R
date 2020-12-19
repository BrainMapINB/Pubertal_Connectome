# First of all, clear workspace
rm(list=ls())
# Save graphical parameters
op <- par()

# Check if AIC file exists
aic_file <- "https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/AIC/AIC_results.csv"

# If you prefer to generate those values follow the steps within this 'if' command
# (but be aware that some computations may be take a while)
if(0){
  
  # Read phenotypic data
  info <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Phenotypic/phenotypic.csv")
  # Remove those sessions with excesive motion artifact
  # And columns with NA's in order to apply 'gamlss'
  info <- info[which(as.logical(info$QC)),c(1:8,12)]
  # Refactorize subjects' ID
  info$BIDS.ID <- factor(info$BIDS.ID)
  # Compute PDS-Sex interaction manually
  info$PDSeM <- info$PDSe; info$PDSeM[info$Sex=="F"] <- 0
  
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
  source("https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Scripts/Auxiliary/cost_th.R")
  source("https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Scripts/Auxiliary/cost_th_ConnMx.R")
  # Load 'gamm4' & 'gamlss' package for GAMM analyses
  if(!require(gamm4)) install.packages("gamm4"); library(gamm4)
  if(!require(gamlss)) install.packages("gamlss"); library(gamlss)
  # Load 'psych' package for Fisher r-to-z transformation
  if(!require(psych)) install.packages("psych"); library(psych)
  
  # Graph names
  grp_var <- c("D_net","C_net","L_net","E_net")
  
  # Explore a set of thresholds
  cost <- seq(from = 0.01, to = 0.48, by = 0.01)
  
  # Store results into empty matrix
  DF <- matrix(data = 0, nrow = length(cost), ncol = 27)
  DF[,1] <- cost
  # Set column names
  colnames(DF) <- c("cost","thr_avg","thr_sd",
                    sapply(c("d","c","l","e"),
                           function(x) paste0(x, c("_LME_Age_AIC","_LME_AgeSex_AIC",
                                                  "_GAMM_sAge_AIC","_GAMM_sAgeSex_AIC",
                                                  "_GAMM_loPDS_AIC","_GAMM_loPDSSex_AIC")
                                              )
                           )
                    )
  
  for(hh in 1:length(cost)){
    
    # Print progress
    print(paste0("Cost ",cost[hh]))
    
    # Apply cost threshold
    thr <- apply(ConnMx,3,function(x) cost_th(x,cost[hh]))
    # Save threshold average
    DF[hh,2] <- mean(thr); DF[hh,3] <- sd(thr)
    # Apply thresholds to Connectivity Matrix
    ConnMx_th <- apply(ConnMx,3,function(x) cost_th_ConnMx(x,cost[hh]))
    ConnMx_th <- array(ConnMx_th,dim=ConnMx_dim)
    
    # Apply complex networks measures.
    # Substract diagonal
    for(ii in 1:ConnMx_dim[3]) diag(ConnMx_th[,,ii]) <- 0
    # Fisher's z transformation
    ConnMx_th <- fisherz(ConnMx_th)
    if(max(ConnMx_th) == Inf) ConnMx_th[ConnMx_th==Inf] <- max(ConnMx_th[ConnMx_th!=Inf])
    
    # Weighted-Degree
    w_deg <- apply(ConnMx_th,3:2,sum)
    # Clustering Coefficiente (per node) - weighted
    tic <- Sys.time()
    w_CC <- t(sapply(1:ConnMx_dim[3], function(x) iCC_w(ConnMx_th[,,x])))
    if(sum(is.na(w_CC))>0)  w_CC[is.na(w_CC)] <- 0
    toc <- Sys.time()
    print("Clustering:"); print(toc-tic)
    # Average distances (per node) - weighted
    tic <- Sys.time()
    w_L <- t(sapply(1:ConnMx_dim[3], function(x) iL_w(ConnMx_th[,,x])))
    toc <- Sys.time()
    print("Path-length:"); print(toc-tic)
    # Inverse average distances (per node) - weighted
    tic <- Sys.time()
    w_E <- t(sapply(1:ConnMx_dim[3], function(x) iEff_w(ConnMx_th[,,x])))
    toc <- Sys.time()
    print("Efficiency:"); print(toc-tic)
    
    # Whole-brain graph-theory measures along the sample
    D_net <- rowMeans(w_deg)
    C_net <- rowMeans(w_CC)
    L_net <- rowMeans(w_L)
    E_net <- rowMeans(w_E)
    
    # Apply linear models at the group level
    info$D_net <- D_net
    info$C_net <- C_net
    info$L_net <- L_net
    info$E_net <- E_net
    
    # Test different models
    for(kk in 1:length(grp_var)){

      # LME model (age slope)
      grp_form <- as.formula(paste0(grp_var[kk], " ~ Age + FDRMS + Coil + (1|BIDS.ID)"))
      fitLIN <- tryCatch(lmer(grp_form,
                              data = info,
                              REML = F),
                          error=function(e) NA)
      
      # LME model (age-sex slope)
      grp_form <- as.formula(paste0(grp_var[kk], " ~ Age*Sex + FDRMS + Coil + (1|BIDS.ID)"))
      fitLINi <- tryCatch(lmer(grp_form,
                               data = info,
                               REML = F),
                          error=function(e) NA)
      
      # GAMM model (age spline)
      grp_form <- as.formula(paste0(grp_var[kk], " ~ s(Age, k=4) + FDRMS + Coil"))
      fitGAMM <- tryCatch(gamm4(grp_form,
                                random = ~ (1|BIDS.ID),
                                data = info,
                                REML = F),
                          error=function(e) NA)
      
      # GAMM model (age*sex splines)
      grp_form <- as.formula(paste0(grp_var[kk], " ~ s(Age, k=4, by=Sex) + Sex + FDRMS + Coil"))
      fitGAMMi <- tryCatch(gamm4(grp_form,
                                 random = ~ (1|BIDS.ID),
                                 data = info,
                                 REML = F),
                           error=function(e) NA)
      
      # GAMM model (PDS spline)
      grp_form <- as.formula(paste0(grp_var[kk], " ~ lo(~PDSe) + FDRMS + Coil + re(random = ~1|BIDS.ID)"))
      fitPDS <- tryCatch(gamlss(grp_form, data = info,
                                control = gamlss.control(n.cyc = 30)),
                         error=function(e) NA)
      
      # GAMM model (PDS*sex splines)
      grp_form <- as.formula(paste0(grp_var[kk], " ~ lo(~PDSe)+lo(~PDSeM)+Sex + FDRMS + Coil + re(random = ~1|BIDS.ID)"))
      fitPDSi <- tryCatch(gamlss(grp_form, data = info,
                                 control = gamlss.control(n.cyc = 30)),
                          error=function(e) NA)
      
      # Store AIC
      DF[hh,4+(kk-1)*6] <- tryCatch(AIC(fitLIN), error=function(e) NA)
      DF[hh,5+(kk-1)*6] <- tryCatch(AIC(fitLINi), error=function(e) NA)
      DF[hh,6+(kk-1)*6] <- tryCatch(AIC(fitGAMM$mer), error=function(e) NA)
      DF[hh,7+(kk-1)*6] <- tryCatch(AIC(fitGAMMi$mer), error=function(e) NA)
      DF[hh,8+(kk-1)*6] <- tryCatch(AIC(fitPDS), error=function(e) NA)
      DF[hh,9+(kk-1)*6] <- tryCatch(AIC(fitPDSi), error=function(e) NA)
      
    }# for(kk in 1:length(grp_var)) 
    
  }# for(hh in 1:length(cost))
  
  # Save results
  outname <- "AIC_results.csv"
  write.csv(DF, outname, quote = F, row.names = F)
}

# Read AIC results
all_fits <- read.csv(aic_file)

# Plot results
# Load 'ggplot2' package
if(!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
# Model labels
mod_labs <- c("LM-Age","LM-Age.Sex","NL-Age","NL-Age.Sex","NP-PDS","NP-PDS.Sex")
# Set initial variables
grph_names <- c("d","c","l","e")
mod_names <- c("_LME_Age_AIC","_LME_AgeSex_AIC",
               "_GAMM_sAge_AIC","_GAMM_sAgeSex_AIC",
               "_GAMM_loPDS_AIC","_GAMM_loPDSSex_AIC")

# Plot best results
# Reshape data.frame
rsh_df <- as.data.frame(rep(all_fits$cost,length(mod_labs))) # cost
rsh_df <- cbind(rsh_df,rep(mod_labs, each = nrow(all_fits))) # varnames
for(gg in 1:length(grph_names)){
  col_idx <- match(paste0(grph_names[gg],mod_names),names(all_fits))
  rsh_df <- cbind(rsh_df,c(as.matrix(all_fits[,col_idx])))
}
names(rsh_df) <- c("cost","Model","Degree","Clustering","Path.Length","Efficiency")

# Plot best AIC
# Degree
aux <- rsh_df[,1:3]
g1_1 <- ggplot(data=aux[which(aux$cost<=0.1),], aes(x=cost*100, y=Degree, color=Model)) + 
    geom_line()+
    xlab("") + ylab("Degree") +
    #ggtitle("Degree") +
    scale_color_brewer(palette = "Set1", guide = F) +
    theme_light()
g1_2 <- ggplot(data=aux[which(aux$cost>0.1),], aes(x=cost*100, y=Degree, color=Model)) + 
  geom_line()+
  xlab("") + ylab("") +
  #ggtitle("Degree") +
  scale_color_brewer(palette = "Set1") +
  theme_light()
# Clustering
aux <- rsh_df[which(!is.na(rsh_df$Clustering)),c(1,2,4)]
g2_1 <- ggplot(data=aux[which(aux$cost<=0.1),], aes(x=cost*100, y=Clustering, color=Model)) + 
  geom_line()+
  xlab("") + ylab("Clustering") +
  #ggtitle("Clustering") +
  scale_color_brewer(palette = "Set1", guide = F) +
  theme_light()
g2_2 <- ggplot(data=aux[which(aux$cost>0.1),], aes(x=cost*100, y=Clustering, color=Model)) + 
  geom_line()+
  xlab("") + ylab("") +
  #ggtitle("Clustering") +
  scale_color_brewer(palette = "Set1") +
  theme_light()
# Path-Length
aux <- rsh_df[,c(1,2,5)]
g3_1 <- ggplot(data=aux[which(aux$cost<=0.1),], aes(x=cost*100, y=Path.Length, color=Model)) + 
  geom_line()+
  xlab("") + ylab("Path-Length") +
  #ggtitle("Path-Length") +
  scale_color_brewer(palette = "Set1", guide = F) +
  theme_light()
g3_2 <- ggplot(data=aux[which(aux$cost>0.1),], aes(x=cost*100, y=Path.Length, color=Model)) + 
  geom_line()+
  xlab("") + ylab("") +
  #ggtitle("Path-Length") +
  scale_color_brewer(palette = "Set1") +
  theme_light()
# Efficiency
aux <- rsh_df[which(!is.na(rsh_df$Efficiency)),c(1,2,6)]
g4_1 <- ggplot(data=aux[which(aux$cost<=0.1),], aes(x=cost*100, y=Efficiency, color=Model)) + 
  geom_line()+
  xlab("") + ylab("Efficiency") +
  #ggtitle("Efficiency") +
  scale_color_brewer(palette = "Set1", guide = F) +
  theme_light()
g4_2 <- ggplot(data=aux[which(aux$cost>0.1),], aes(x=cost*100, y=Efficiency, color=Model)) + 
  geom_line()+
  xlab("Connectivity density (%)") + ylab("") +
  #ggtitle("Efficiency") +
  scale_color_brewer(palette = "Set1") +
  theme_light()

# Combine plots
if(!require(gridExtra)) install.packages("gridExtra"); library(gridExtra)
g5_1 <- grid.arrange(g1_1, g2_1, g3_1, g4_1, ncol=1)
g5_2 <- grid.arrange(g1_2, g2_2, g3_2, g4_2, ncol=1)
ggAIC <- grid.arrange(g5_1, g5_2, ncol=2)

# Save plot
outfile <- "Connectome_AIC.pdf"
#ggsave(outfile, ggAIC, device = "pdf")
