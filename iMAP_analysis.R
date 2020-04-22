# iMAP Analysis // Bosco TADDEI

# This script is analyzing cytokine and cytometric data
# collected from patient treated for lupus.


# General Set-up Functions -----------------------------------------------------------

# Packages acquisition 
load_pkg <- function(){
  library("plyr")
  library("dplyr")
  library("stringr")
  library("ggplot2")
  library("factoextra")
  library("gridExtra")
  library("ggrepel")
  library("reshape2")
  library("RGCCA")
}
load_pkg()

# Set Work directory
setwd("./..")
getwd()


# Creating the different data frames -------------------------------------------------

# Get clinical data
get_clini <- function(){
 
   df <- read.table("lupil2_populations.txt", 
                    header=TRUE, sep="\t")    %>%
         dplyr::rename(Subject = SUBJID)      %>%
         arrange(Subject)
  
  # Add quality of life columns
  df["Subj_id"]    <- factor(row.names(df))
  df["Country_id"] <- substr(df$Subject,1,3) %>%
                      as.factor()
  
  return(df)
}

# Get Cytometrie Panels data 
get_cytometrie <- function(subjects = paper_subjects, panels = all_panels){
  
  # Initialize the df
  df <- data.frame(paper_subjects) %>%
        dplyr::rename(Subject = paper_subjects)
  
  # Collect data from the panels
  for (id in panel_ids){
    file_panel <- paste0("lupilCytometrie/Panel_", id, "_norm.txt")
    panel <- read.table(file_panel, header=TRUE, sep="\t") %>%
             dplyr::rename(Subject = Sample.name)
    
    df <- join(df, panel)
  }
  
  #Reorder the df regarding Subjects
  df <- arrange(df, Subject, Visit)
  
  return(df)
}

# Get Cytokine data
get_cytokine <- function(){
  file_cytokine <- "Lupil2_cytokines_191105/lupil2_conc_log2imp_20191022.csv"
  df <- read.table(file_cytokine, header=TRUE, sep=",")
  
  # Extracting Subject id and visit id
  id <- as.vector(df[,"X"])
  id <- strsplit(id,"V")
  subject <- sapply(id, "[[", 1)
  subject <- str_pad(subject, 9, "left", "0") 
  subject <- paste0(substr(subject,1,3),"-",substr(subject,4,6),"-",substr(subject,7,9))
  visit   <- sapply(id, "[[", 2)
  visit   <- str_pad(visit, 2, "left", "0")
  visit   <- paste0("V",visit)
  df["Subject"] <- as.factor(subject)
  df["Visit"]   <- as.factor(visit)
  
  df <- df[-1] %>%
        subset(Subject %in% paper_subjects) %>%
        arrange(Subject, Visit) %>%
        dplyr::select("Subject","Visit",everything())
  
  return(df)
}


# Other information extracting functions ---------------------------------------------

# Names of the cytometrie numerical variables 
get_varnum <- function(df){
  
  varnum <- df %>%
            select_if(is.numeric) %>%
            colnames()
  
  return(varnum)
}

# Long data frame      
make_long <- function(df,crits){
  
  factors <- df %>%
             select_if(is.factor) %>%
             colnames()
  
  lg_df <- df %>%
           melt(id.vars = factors) %>%
           na.omit() %>%
           join(df_clini[,c("Subject",crits)]) %>%
           arrange(Subject,Visit)
}

# Add panel info
manage_panel <- function(lg_df){
  
  lg_df["Panel"]    <- lg_df$variable %>%
                        substr(1,8) %>%
                        as.factor()
  lg_df["variable"] <- lg_df$variable %>%
                       substring(10,1000) %>%
                       as.factor()
  
  #Re-order cols
  lg_df <- lg_df %>%
           dplyr::select("Subject","Visit","Panel",everything())
  
  return(lg_df)
}

# Remove cols that have a HS equivalent
manage_HS <- function(df){
  
  varnum  <- get_varnum(df)
  badcols <- varnum[grepl('HS',varnum)] %>%
              substring(3)
  df[,c(badcols)] <- list(NULL)
  
  return(df)
}

#Remove cols that have var == 0
drop_nul_var <- function(df){
  
  badcols <- df[,which(apply(df,2,var) == 0)] %>%
             colnames()
  print(paste("Zero variance cols removed :",badcols))
  
  df[,badcols] <- list(NULL)
  
  return(df)
}


# Analysis Functions -----------------------------------------------------------------

# Prepare data 
get_data <- function(df, varnum, visit, crit){
  
  dat <- df[df$Visit == visit,] %>%
         join(df_clini[,c("Subject", crit, "Subj_id")]) %>%
         dplyr::select(c(varnum, crit, "Subj_id"))
  row.names(dat) <- dat$Subj_id
  dat$Subj_id <- NULL
  
  return(dat)
}

#PCA
plot_pca <- function(dat, varnum, crit, cols, scale = TRUE){
  
  pca <- dat[,varnum] %>%
         na.omit() %>%
         drop_nul_var() %>%
         prcomp(scale = scale)
  
  group <- dat %>%
           na.omit() %>%
           dplyr::select(crit)
  
  eig <- fviz_eig(pca)
  
  ind <- fviz_pca_ind(pca,
                      col.ind = group[,1]) +
         scale_color_manual(values = cols)
  
  var <- fviz_pca_var(pca,
                      select.var = list(contrib = 5),
                      repel = TRUE)
  print(eig)
  print(var)
  print(ind)
  
  
  return(pca)
}

# t.test
plot_t.test <- function(dat, crit, pval.max, cols){
  
  varnum <- dat %>%
            drop_nul_var() %>%
            get_varnum()
  y <- dat[,crit]
      
  i <- 0
  vars <- c()
  for (var in varnum){
    x = dat[,var]
    
    t <- t.test(x ~ y)
    pval <- t$p.value
    
    if (pval < pval.max){
      i <- i+1
      vars[[i]] <- var
      
      bp <- ggplot(data = dat, aes(x = y, y = x, fill = y)) + 
            geom_boxplot() +
            geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1.5, color = "white", alpha = 0.4) +
            labs(title = var, y = "", 
                 x = paste("Student test pval =", pval)) +
            theme_classic() +
            scale_fill_manual(values = cols)
      print(bp)
    }
    
  }
  
  return(vars)
}

#Plot Wilcoxon test
plot_w.test <- function(dat, crit, pval.max, cols){
  
  varnum <- dat %>%
            drop_nul_var() %>%
            get_varnum()
  y <- dat[,crit]
  
  i <- 0
  vars <- c()
  for (var in varnum){
    x = dat[,var]
    
    t <- wilcox.test(x ~ y)
    pval <- t$p.value
    
    if (pval < pval.max){
      i <- i+1
      vars[[i]] <- var
      
      bp <- ggplot(data = dat, aes(x = y, y = x, fill = y)) + 
        geom_boxplot() +
        geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1.5, color = "white", alpha = 0.4) +
        labs(title = var, y = "", 
             x = paste("Wilcoxon test pval =", pval)) +
        theme_classic() +
        scale_fill_manual(values = cols)
      print(bp)
    }
    
  }
  
  return(vars)
}

# PLS DA
plot_pls <- function(dat, varnum, crit){
  
  group <- dat[,crit]; names(group) <- rownames(dat)
  y <- data.frame(model.matrix( ~  group-1, data = group))
  blocks <- list(x = dat[,varnum], y = y)
  pls <-rgcca(block = blocks, ncomp = c(3,1), type = "pls")
  
  plot(pls, i_block=1, type="both", resp = group, compx = 1, compy = 2)
}



# Load Clinical data -----------------------------------------------------------------

df_clini <- get_clini()

# Extract paper subjects
paper_subjects <- df_clini[df_clini[,"Paper"] == T,"Subject"]




# Load Cytokine data -----------------------------------------------------------------

# General Setup
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
visit    <- "V05"
crit     <- c("Responder")
cols     <- c("#FF3333", "#00CC00")

# Get appropriate data
dat <- get_data(df, varnum, visit , c("Responder","ARM"))


dat_IL2 <- dat[dat$ARM == "ILT-101",]
dat_IL2$ARM <- NULL
dat_pcb <- dat[dat$ARM == "Placebo",]
dat_pcb$ARM <- NULL          

#PCA
pca_IL2 <- plot_pca(dat_IL2,varnum, crit, cols)
pca_pcb <- plot_pca(dat_pcb,varnum, crit, cols)

#t.test and boxplot if pval significant
plot_w.test(dat_IL2, crit, 0.05, cols)
plot_w.test(dat_pcb, crit, 0.05, cols)

# Evolution from baseline ----------------------------------------------------------

#Set-up parameters
crit <- "ARM"
cols <- c("#0066CC", "#99CCFF")
pval.max <- 0.05/(length(varnum))

dat1 <- get_data(df, varnum, "V01", c("Responder","ARM"))
dat5 <- get_data(df, varnum, "V09", c("Responder","ARM"))
dat  <- dat5
dat[,varnum] <- dat5[,varnum] - dat1[,varnum]

#PCA
pca <- plot_pca(dat, varnum, crit, cols)


#t.test and boxplot if pval < .05 / 62
plot_t.test(dat, crit, pval.max, cols)


#PLS
plot_pls(dat, varnum, crit)






# Unorganized bits of script ---------------------------------------------------------


# Analyzing subject 90 "The PCA Outlier" -----------------------------------------------
subj <- dat["90",varnum]
mini <- dat[,varnum] %>%
        na.omit() %>%
        lapply("min") %>%
        as.data.frame()
maxi <- dat[,varnum] %>%
        na.omit() %>%
        lapply("max") %>%
        as.data.frame()
  
mins <- colnames(subj[, subj == mini])
maxs <- colnames(subj[, subj == maxi])

for (var in mins){
  data <- df_cytokine[, var]
  hist(data,
       main = paste("Histogram of", var, "V01"),
       xlab = var,
       col = "cadetblue3")
}

# PCA --------------------------------------------------------------------------------

# More PCA info
main_contrib <- function(df, pca, n, add_il2 = FALSE){
  contribs <- get_pca_var(pca_cytokine)$contrib[,1] %>%
    sort(decreasing = TRUE) %>%
    names()
  contribs <- contribs[1:n]
  if (add_il2){
    contribs <- append(contrib_dim1,"IL.2")
  }
  data = df[,contribs]
  if (length(contribs) == 1){
    hist(data,
         xlab = contribs)
  } else {
    plot(data)
  }
}



main_contrib(df_cytokine,pca_cytokine,0,TRUE)

contrib_dim1 <- get_pca_var(pca_cytokine)$contrib[,1] %>%
  sort(decreasing = TRUE) %>%
  names()
contrib_dim1 <- contrib_dim1[1:3]

plot(df_cytokine[,c(contrib_dim1[1:3],"IL.2")])
var_1 = contrib_dim1[1]
t.test(df_cytokine[,var_1])   
hist(df_cytokine[,var_1],
     xlab = var_1)


# Analyse monovariée
# Plot Cytokine data distribution
vars <- df_cytokine %>%
  select_if(is.numeric) %>%
  colnames() %>%
  reorder(contrib_dim1)
for (var in vars){
  data <- df_cytokine[, var]
  png(paste0("./plots/hist_",var,"_all.png"))
  hist(data,
       main = paste("Histogram of", var, "(all)"),
       xlab = var,
       col = "cadetblue3")
  dev.off()
}

# Plot Means
mean <- df_cytokine[,vars] %>%
  colMeans(na.rm = TRUE)
sem  <- df_cytokine[,vars] %>%
  apply(2,sd,na.rm = TRUE)
p <- join(df_cytokine,df_pop)[,c(vars,"Responder")] %>%
  ggplot() +
  geom_bar(aes(x=vars, y = mean), 
           stat = "identity",
           fill = "cadetblue3") +
  geom_errorbar((aes(x=vars, ymin=mean - sem, ymax = mean + sem))) +
  facet_grid( ~ Responder) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Cytokine data (Mean and SEM)")
p




# RGCCA ----------------------------------------------------------------------------------
resp <- dat_cytokine_ev5$ARM; names(resp) <- rownames(dat_cytokine_ev5)
reponse <- data.frame(model.matrix( ~  resp-1, data=resp))
blocks <- list(df_cytokine[,varnum_cytokine],df_cytometrie[,varnum_cytometrie], reponse = reponse)
resRgcca=rgcca(block = blocks,ncomp=c(3,1),type="rgcca")
plot(resRgcca,i_block=1,type="both",resp = resp,compx = 1, compy = 2)
