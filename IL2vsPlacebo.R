# Code for effect of IL2 // Bosco Taddei

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

# Not Normalized file
#file_cytokine <- "Lupil2_cytokines_191105/lupil2_conc_log2imp_20191022.csv"

# Normalized file
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
visit    <- "V01"
crit     <- c("ARM")
cols     <- c("#FF3333", "#00CC00")

# Get appropriate data
dat <- get_data(df, visit, crit)

#t.test and boxplot if pval significant
plot_w.test(dat, pval.max, cols)

#PCA
pca <- plot_pca(dat, cols)

# PLS-DA
plot_pls(dat,cols)



# Evolution from baseline ----------------------------------------------------------

#Set-up parameters
visit <- "V05"
crit  <- c("ARM")
cols  <- c("#FF3333", "#00CC00")
pval.max <- 0.05/(length(varnum))

dat0 <- get_data(df, "V01", crit)
dat  <- get_data(df, visit, crit) %>%
        dat_diff(dat0)


#t.test and boxplot if pval < .05 / 62
plot_w.test(dat_IL2, cols)

#PCA
pca <- plot_pca(dat, cols)

#PLS
plot_pls(dat, type = "cor", 2)