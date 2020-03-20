#######################
#### Import XLS #######
#######################
rm(list = ls())
setwd("~/Desktop")
library("readxl")
library(dplyr)
library(tibble)

All_cytokine_data <- read_excel("CART_UH_combined_Data_AS_conc_table_NEW.xlsx")
All_cytokine_data <- as.data.frame(All_cytokine_data)
rownames(All_cytokine_data) <- All_cytokine_data$...1
All_cytokine_data <- All_cytokine_data[,-1]

#Subset on interest####
Tm6 <- subset(All_cytokine_data, Time == 'Tm6')
T0 <- subset(All_cytokine_data, Time == 'T0')




##############################
#### FEATURE SELECTION #######
##############################

X <- as.matrix(T0[,c(3:35)]) #(data frame of cytokines as columns and samples as rows)

#X <- scale(X, scale = T, center = T) #Complete_Response,CRS,Neurotoxicity
Y <- as.factor(T0$Complete_Response)
Y

library(glmnet)
# bionomial glmnet
fitFull <- glmnet(x = X,
                  y = Y,
                  family = "binomial",
                  alpha = 1)
set.seed(seed = 1)
fit2Full <- cv.glmnet(x = X,
                      y = Y,
                      family = "binomial",
                      lambda = fitFull$lambda,
                      type.measure = "auc",
                      nfolds = 33,
                      alpha = 1)
# print AUC
max(fit2Full$cvm)
plot(fit2Full)

# print coeffiencts
coefs <- as.matrix(coef(fit2Full, s = "lambda.min")) %>% as.data.frame()
coefs <- coefs[-1, , drop = F]
coefs <- coefs[c(abs(coefs$"1") > 0), , drop = F] %>% rownames_to_column()
colnames(coefs) <- c("cytokine", "coef")
coefs <- coefs %>% arrange(desc(coef))
coefs2 <- coefs %>% .$cytokine
coefs2




#######################################################################
#### Testing predictive model at different timepoints/cohorts #########
#######################################################################
Fit <- glm(as.factor(CRS) ~ `MCP.3`+`IL.29.IFN.1`+`IL.8`+`MCP.1`+`IL.27`+`TGF.b1`+`Fractalkine`, 
           data = x_m6, family = "binomial")
summary(Fit)




##############################
#### PCA and Heatmap #########
##############################
### PCA
library(ggfortify)
df <- All_cytokine_data[,c(4:36)]
autoplot(prcomp(df, center = T, scale. = T), data = All_cytokine_data, loadings = TRUE, loadings.label = T) +  geom_point(size = 1, aes(color = Time))

All_cytokine_data <- as.data.frame(All_cytokine_data)
rownames(All_cytokine_data) <- All_cytokine_data[,1]
All_cytokine_data <- All_cytokine_data[,-1]

annot_row <- All_cytokine_data[,c(20,19)]
mat_Clus1 <- All_cytokine_data[,c(1:17)]

hm.parameters <- list(t(mat_Clus1),
                      scale = "row",
                      cellwidth =10, cellheight=10,
                      color = colorRampPalette(c("blue1","blue1","white","red1","red1"))(100),
                      kmeans_k = NA,
                      show_rownames = T, show_colnames = T,
                      fontsize = 8, 
                      main = "",
                      clustering_method = "ward.D2",
                      cluster_rows = T, cluster_cols = F,
                      clustering_distance_rows = "canberra",
                      clustering_distance_cols = "canberra",
                      annotation_col = annot_row,
                      gaps_col = c(14, 44))
library(pheatmap)
kmean.hm <- do.call("pheatmap", hm.parameters)

#Subset on interest####
Control <- subset(All_cytokine_data, Group == 'Control')
High_Intensity <- subset(All_cytokine_data, Group == 'High Intensity')
Standard_Intensity <- subset(All_cytokine_data, Group == 'Standard Intensity')




####################################
#### Univariate Statistics #########
####################################

#Univariates for Control####
Control_V1 <- subset(Control, Visit == "V1A.0")
Control_V6 <- subset(Control, Visit == "V6.0")

datDF_Control_V1 <- Control_V1[,c(1:17)]
datDF_Control_V6 <- Control_V6[,c(1:17)]

#User vs Non.User - Wilcox test
eigen.results <- matrix(NA, ncol(datDF_Control_V1), 0)
rownames(eigen.results) <- colnames(datDF_Control_V1)
eigen.results <- as.data.frame(eigen.results)

for(i in 1:nrow(eigen.results)){
  eigen.results$Median.log2FC[i] <- log2(median(as.numeric(as.character(datDF_Control_V6[,i])), na.rm = TRUE)/median(as.numeric(as.character(datDF_Control_V1[,i])), na.rm = TRUE))
  eigen.results$Wilcox.NonPara.pval[i] <- wilcox.test(as.numeric(as.character(datDF_Control_V6[,i])), as.numeric(as.character(datDF_Control_V1[,i])), alternative = "two.sided", paired=TRUE)$p.value
  eigen.results$Para.pval[i] <- t.test(as.numeric(as.character(datDF_Control_V6[,i])), as.numeric(as.character(datDF_Control_V1[,i])), alternative = "two.sided", paired=TRUE, var.equal = TRUE)$p.value
}

write.csv(eigen.results, 'Univariate_Control_Paired.csv')




#########################################################
#### Defining Heirarichal Clusters of cytokines #########
#########################################################
# Clustering all data###
library("factoextra")
library("cluster")
library("tibble")
library("dplyr")

#Hclust
set.seed(seed = 2)
hclustFun <- function(x, k) {
  return(value = list(cluster = cutree(hclust(d = as.dist(1 - cor(t(x),
                                                                  method = "spearman"))),
                                       k)))}

fit <- clusGap(t(High_Intensity[,c(1:17)]), hclustFun, K.max = 16, B = 100)
nC <- maxSE(fit$Tab[, "gap"], fit$Tab[, "SE.sim"], method = "globalmax")

# plot gap
plot(fit, main = "Number of optimal plasma cytokine clusters")
abline(v = nC, col = "blue")

# identify clusters
clLS <- hclustFun(t(High_Intensity[,c(1:17)]), nC) %>%
  as.data.frame %>%
  rownames_to_column() %>% arrange(cluster)
clLS

#K-means
# K.max is total cyokines - 1
gap_stat <- clusGap(t(mat_Clus1), FUN = kmeans, 
                    K.max = 16, B = 500)
print(gap_stat)
fviz_gap_stat(gap_stat, maxSE = list(method = "Tibs2001SEmax"))


## Use inout from above data to optimal number of clusters
set.seed(123)
km.res <- kmeans(t(mat_Clus1), 6, nstart = 25)

write.csv(km.res[["cluster"]], 'cluster_components_meso_clusters.csv')
write.csv(km.res[["centers"]], 'cluster_centers_meso_clusters.csv')