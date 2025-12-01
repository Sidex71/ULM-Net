v#############################parameter tuning############
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
library(devtools)

devtools::install_github("Sidex71/ULM")

Idents(invitro_obj) <- invitro_obj$sorting_scheme
invit_mark <- subset(invitro_obj, idents = c("TCRb+","CD11c+"))
num_genes <- c(25, 50, 75, 100, 150, 200, 250, 300, 400, 500, 600, 800, 1000)
sig_list <- lapply(num_genes, function(x){
  ULM::GetSignature(invit_mark, ident_col = 'sorting_scheme', n= x)
})

#invitro_obj@meta.data<- invitro_obj@meta.data %>% dplyr:: select(-c("count_ulm" , "celltype_ulm", "avg_pvalue", "avg_score", 'statistic'))

objTune_list <- lapply(sig_list, function(x){
  my_obj <- invitro_obj
  my_scores <- ULM::GetCellScores(seurat_obj = my_obj, signatures = x, assay = 'RNA3', slot = 'data', layer = NULL)
  my_ass <- ULM::GetCellAssignments(score_data = my_scores)
  my_obj <- ULM::AddMetaObject(my_obj, cell_class_df = my_ass)
})  

library(caret)
confmatTune_list <- lapply(objTune_list, function(x){
  actual <- ifelse(x$sorting_scheme=='CD11c+' | 
                     x$sorting_scheme=='TCRb+', 'No', 'Yes')
  predicted <- ifelse(x$celltype_ulm=='CD11c+_TCRb+', 'Yes', 'No')
  
  my_tab <- table(actual, predicted)
  matcon <-confusionMatrix(my_tab, positive = 'Yes')
  matcon$byClass[1:4]
})

confmatTune <- do.call(rbind, confmatTune_list)
colnames(confmatTune) <- c('Precision', "Neg pred",  "Sensitivity", "Specificity")
confmatTune <- as.data.frame(confmatTune)
confmatTune <- confmatTune * 100
confmatTune$n_genes <- num_genes
confmatTune$`Neg pred` <- NULL


##############plot metrics
library(ggplot2)
library(dplyr)
library(tidyr)



df_long <- confmatTune %>% 
  pivot_longer(cols = c("Precision", "Sensitivity"),names_to = "Metric", values_to = "Value")

p1<- ggplot(df_long, aes(x = n_genes, y = Value, color = Metric)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Performance Metrics Across Increasing Gene Set Size",
    x = "Number of Genes",
    y = "Metric Value"
  ) 
png("/mnt/8TB/users/shameed/shameed/Doublet predictions/figures/parameter_tuning.png", width = 12, height = 8.5, units = 'in', res = 600)
p1
dev.off()

write_csv(confmatTune, '/mnt/8TB/users/shameed/shameed/Doublet predictions/confmatTune.csv')
saveRDS(confmatTune, '/mnt/8TB/users/shameed/shameed/Doublet predictions/confmatTune.rds')
