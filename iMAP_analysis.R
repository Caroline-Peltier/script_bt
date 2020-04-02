# iMAP Analysis // Bosco TADDEI

# This script is analyzing cytokine and cytometric data
# collected from patient treated for lupus.



# Packages acquisition --------------------------------------------------------------
library("plyr")
library("dplyr")
library("stringr")
library("ggplot2")
library("factoextra")
library("gridExtra")
library("ggrepel")
library("reshape2")
library("RGCCA")


# Loading data ----------------------------------------------------------------------
setwd("./..")
getwd()

# Get population --------------------------------------------------------------------
df_pop <- read.table("lupil2_populations.txt", 
                     header=TRUE, sep="\t") %>%
          dplyr::rename(Subject = SUBJID) %>%
          arrange(Subject)

# Add quality of life columns
df_pop["Subj_id"]    <- factor(row.names(df_pop))
df_pop["Country_id"] <- substr(df_pop$Subject,1,3) %>%
                        as.factor()

summary(df_pop)

# Extract paper subjects
paper_subjects <- df_pop[df_pop[,"Paper"] == T,"Subject"]


# Get Cytometrie Panels data --------------------------------------------------------
# Initialize the df
df_cytometrie <- data.frame(paper_subjects) %>%
                 dplyr::rename(Subject = paper_subjects)



# Collect data from the different panels
panel_ids <- c("02","03","08","09","10","12","13")
for (id in panel_ids){
  file_panel = paste0("lupilCytometrie/Panel_", id, ".txt")
  panel <- read.table(file_panel,
                      header=TRUE, sep="\t") %>%
           dplyr::rename(Subject = Sample.name)
  df_cytometrie  <- join(df_cytometrie, panel)
}

df_cytometrie <- arrange(df_cytometrie, Subject, Visit)
summary(df_cytometrie)

# Names of the cytometrie numerical variables 
varnum_cytometrie <- df_cytometrie %>%
                     select_if(is.numeric) %>%
                     colnames()
                
# Long cytometrie data frame      
lg_cytometrie <- melt(data = df_cytometrie, id.vars = c("Subject", "Visit")) %>%
                 na.omit() %>%
                 join(df_pop[,c("Subject","ARM","Responder")]) %>%
                 arrange(Subject,Visit)
lg_cytometrie["Panel"] <- lg_cytometrie$variable %>%
                          as.vector() %>%
                          strsplit(split = "_P") %>%
                          sapply("[[",1)

# Get Cytokine data ------------------------------------------------------------------
file_cytokine <- "Lupil2_cytokines_191105/lupil2_conc_log2imp_20191022.csv"
df_cytokine <- read.table(file_cytokine,
                          header=TRUE, sep=",")

# Extracting Subject id and visit id
id <- as.vector(df_cytokine[,"X"])
id <- strsplit(id,"V")
subject <- sapply(id, "[[", 1)
subject <- str_pad(subject, 9, "left", "0") 
subject <- paste0(substr(subject,1,3),"-",substr(subject,4,6),"-",substr(subject,7,9))
visit   <- sapply(id, "[[", 2)
visit   <- str_pad(visit, 2, "left", "0")
visit   <- paste0("V",visit)
df_cytokine["Subject"] <- as.factor(subject)
df_cytokine["Visit"] <- as.factor(visit)

df_cytokine <- df_cytokine[-1] %>%
               subset(Subject %in% paper_subjects) %>%
               arrange(Subject, Visit)

# Name of the cytokine numerical variables
varnum_cytokine <- df_cytokine %>%
                   select_if(is.numeric) %>%
                   colnames()

# Remove cols that have a HS equivalent
badcols <- varnum_cytokine[grepl('HS',varnum_cytokine)] %>%
           substring(3)
df_cytokine[,c(badcols)] <- list(NULL)
varnum_cytokine <- df_cytokine %>%
                   select_if(is.numeric) %>%
                   colnames()

# Remove the Outlier
df_cytokine <- df_cytokine[,"Subj_id" != 90]

summary(df_cytokine)



# Long cytokine data frame
lg_cytokine <- melt(data = df_cytokine, id.vars = c("Subject", "Visit")) %>%
               na.omit() %>%
               join(df_pop[,c("Subject","ARM","Responder")]) %>%
               arrange(Subject,Visit)



# Analyzing data ---------------------------------------------------------------------

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


# Preliminary Functions --------------------------------------------------------------

# Prepare data 
get_data <- function(df, varnum, visit, crit){
  
  dat <- df[df$Visit == visit,] %>%
         join(df_pop[,c("Subject", crit, "Subj_id")]) %>%
         dplyr::select(c(varnum, crit, "Subj_id"))
  row.names(dat) <- dat$Subj_id
  dat$Subj_id <- NULL
  
  return(dat)
}

#PCA
plot_pca <- function(dat, varnum, crit, cols){
  
  pca <- dat[,varnum] %>%
         na.omit() %>%
         prcomp( scale = TRUE)
  
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
  print(ind)
  print(var)
  
  return(pca)
}

# t.test
plot_t.test <- function(dat, varnum, crit, pval.max, cols){
  
  y = dat[,crit]
  
  i = 0
  vars = c()
  for (var in varnum){
    x = dat[,var]
    
    t <- t.test(x ~ y)
    pval = t$p.value
    
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


# PLS DA
plot_pls <- function(dat, varnum, crit){
  
  group <- dat[,crit]; names(group) <- rownames(dat)
  y <- data.frame(model.matrix( ~  group-1, data = group))
  blocks <- list(x = dat[,varnum], y = y)
  pls <-rgcca(block = blocks, ncomp = c(3,1), type = "pls")
  
  plot(pls, i_block=1, type="both", resp = group, compx = 1, compy = 2)
}


# Cytokine --------------------------------------------------------------------------

# Visit V01 : Predict Responder -----------------------------------------------------

#Set-up parameters
df       <- df_cytokine
varnum   <- varnum_cytokine
pval.max <- 0.05/length(varnum)
visit    <- "V01"
crit     <- "Responder"
cols     <- c("#FF3333", "#00CC00")

# Get appropriate data
dat <- get_data(df, varnum, visit , crit)

#PCA
pca <- plot_pca(dat, varnum, crit, cols)

dat_90 <- dat[which(rownames(dat) != "90"),]

pca_90 <- plot_pca(dat_90, varnum, crit, cols)
  
#t.test and boxplot if pval significant
plot_t.test(dat, varnum, crit, pval.max, cols)


# Visit V05 - V01 : Predict Placebo ---------------------------------------------------

#Set-up parameters
df <- df_cytokine
varnum <- varnum_cytokine
crit <- "ARM"
cols <- c("#0066CC", "#99CCFF")
pval.max <- 0.05/(length(varnum))

dat1 <- get_data(df, varnum, "V01", crit)
dat5 <- get_data(df, varnum, "V05", crit)
dat  <- dat5
dat[,varnum] <- dat5[,varnum] - dat1[,varnum]

#PCA
pca <- plot_pca(dat, varnum, crit, cols)


#t.test and boxplot if pval < .05 / 62
plot_t.test(dat, varnum, crit, pval.max, cols)


#PLS
plot_pls(dat, varnum, crit)

# ------------------------------------------------------------------------------------
resp <- dat_cytokine_ev5$ARM; names(resp) <- rownames(dat_cytokine_ev5)
reponse <- data.frame(model.matrix( ~  resp-1, data=resp))
blocks <- list(df_cytokine[,varnum_cytokine],df_cytometrie[,varnum_cytometrie], reponse = reponse)
resRgcca=rgcca(block = blocks,ncomp=c(3,1),type="rgcca")
plot(resRgcca,i_block=1,type="both",resp = resp,compx = 1, compy = 2)
