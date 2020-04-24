# Code for Responder  in IL2 treated patient // Bosco Taddei

# General Set Up --------------------------------------------------------------------

# Set Work directory
setwd("./..")
getwd()

# Load functions
source("script_bt/iMAP_analysis.R")

#Load Packages
load_pkg()

# Load Clinical data -----------------------------------------------------------------

df_clini <- get_clini()

# Extract paper subjects
paper_subjects <- df_clini[df_clini[,"Paper"] == T,"Subject"]



# Load Cytokine data -----------------------------------------------------------------

file_cytokine <- "Lupil2_cytokines_191105/lupil2_conc_log2impnorm_20191022.csv"

df <- get_cytokine() %>%
      manage_HS()

varnum <- get_varnum(df)

# Long cytokine data frame
lg_df <- make_long(df,c("Responder","ARM"))



# Load Cytometrie data ---------------------------------------------------------------

#General Set up
all_panels <- c("01","02","03","04","05","08","09","10","12")

df <- get_cytometrie()

varnum <- get_varnum(df)

lg_df <- make_long(df,c("Responder","ARM")) %>%
         manage_panel()


# Per time point analysis -----------------------------------------------------------

#Set-up parameters
pval.max <- 0.05/length(varnum)
visit    <- "V29"
crit     <- c("Responder","ARM")
cols     <- c("#FF3333", "#00CC00")

# Get appropriate data
dat <- get_data(df, visit, crit)


dat_IL2 <- dat[dat$ARM == "ILT-101",]
dat_IL2$ARM <- NULL
dat_pcb <- dat[dat$ARM == "Placebo",]
dat_pcb$ARM <- NULL          


#t.test and boxplot if pval significant
plot_w.test(dat_IL2, 0.05, cols)
plot_w.test(dat_pcb, 0.05, cols)


#PCA
pca_IL2 <- plot_pca(dat_IL2, cols)
pca_pcb <- plot_pca(dat_pcb, cols)


# PLS-DA
plot_pls(dat_IL2,cols)
plot_pls(dat_pcb,cols)



# Evolution from baseline ----------------------------------------------------------

#Set-up parameters
visit <- "V05"
crit  <- c("Responder","ARM")
cols  <- c("#FF3333", "#00CC00")
pval.max <- 0.05/(length(varnum))

dat0 <- get_data(df, "V01", crit)
dat  <- get_data(df, "V05", crit) %>%
        dat_diff(dat0)

dat_IL2 <- dat[dat$ARM == "ILT-101",]
dat_IL2$ARM <- NULL
dat_pcb <- dat[dat$ARM == "Placebo",]
dat_pcb$ARM <- NULL 

#t.test and boxplot if pval < .05 / 62
plot_w.test(dat_IL2, 0.05, cols)
plot_w.test(dat_pcb, 0.05, cols)


#PCA
pca <- plot_pca(dat, varnum, crit, cols)





#PLS
plot_pls(dat_IL2, varnum, crit, type = "cor", 2)
