# Unorganized bits of script for iMAP Analysis // Bosco Taddei


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

for (variable in mins){
  data <- df_cytokine[, variable]
  hist(data,
       main = paste("Histogram of", variable, "V01"),
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
