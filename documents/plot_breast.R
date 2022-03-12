for (i in list.files("~/Documents/main_files/AskExplain/GeneCoder/data/workflow/breast/",full.names = T)){
  load(i)
}

print(c("rotate p-value",paste(c("0vs90:    ","0vs180:    ","0vs270:    "),round(do.call('c',lapply(rotate_spot2gex,function(X){X$p.value})),5))))
print(c("displace p-value",paste(c("0vs10:    ","0vs20:    ","0vs30:    "),round(do.call('c',lapply(displace_spot2gex,function(X){X$p.value})),5))))
print(c("similar p-value",paste(c("0,0vs(-)3,(-)3:    ","0,0vs3,(-)3:    ","0,0vs(-)3,3:    ","0,0vs3,3:    "),round(do.call('c',lapply(similar_spot2gex,function(X){X$p.value})),5))))


  

sample_similarity = base_spot2gex$sample_wise
feature_similarity = base_spot2gex$gene_wise
title_name = "base_similarity"
tissue_name = "breast"

library(ggplot2)

lm <- rbind(c(1,2),
            c(1,2),
            c(1,2))

g1 <- ggplot(data.frame(Measure=sample_similarity,Metric="Predicted vs Observed \n Sample-wise Pearson Correlation"), aes(x=Metric,y=Measure)) + 
  geom_violin() + ylim(-1,1) 

g2 <- ggplot(data.frame(Measure=feature_similarity,Metric="Predicted vs Observed \n Gene-wise Pearson Correlation"), aes(x=Metric,y=Measure)) + 
  geom_violin() + ylim(-1,1) 

gg_plots <- list(g1,g2)

library(grid)
library(gridExtra)
final_plots <- arrangeGrob(
  grobs = gg_plots,
  layout_matrix = lm
)


ggsave(final_plots,filename = paste("~/Documents/main_files/AskExplain/GeneCoder/github/github_figures/jpeg_accuracy_gcode_",tissue_name,"_SPATIAL/",title_name,"_accuracy_",tissue_name,"_spatial.png",sep=""),width = 5,height=5)










first_most_variable_gene = spline_signal_adj.r.squared[1,]
second_most_variable_gene = spline_signal_adj.r.squared[2,]
third_most_variable_gene = spline_signal_adj.r.squared[3,]

title_name = "spatial_signal"
tissue_name = "breast"

library(ggplot2)

lm <- rbind(c(1,2,3),
            c(1,2,3),
            c(1,2,3))

g1 <- ggplot(data.frame(Measure=first_most_variable_gene,Metric="Adj.R.Squared of Spline fits \n 1st most variable gene"), aes(x=Metric,y=Measure)) + 
  geom_violin() + ylim(-1,1) 

g2 <- ggplot(data.frame(Measure=second_most_variable_gene,Metric="Adj.R.Squared of Spline fits \n 2nd most variable gene"), aes(x=Metric,y=Measure)) + 
  geom_violin() + ylim(-1,1) 

g3 <- ggplot(data.frame(Measure=third_most_variable_gene,Metric="Adj.R.Squared of Spline fits \n 3rd most variable gene"), aes(x=Metric,y=Measure)) + 
  geom_violin() + ylim(-1,1) 

gg_plots <- list(g1,g2,g3)

library(grid)
library(gridExtra)
final_plots <- arrangeGrob(
  grobs = gg_plots,
  layout_matrix = lm
)


ggsave(final_plots,filename = paste("~/Documents/main_files/AskExplain/GeneCoder/github/github_figures/jpeg_accuracy_gcode_",tissue_name,"_SPATIAL/",title_name,"_accuracy_",tissue_name,"_spatial.png",sep=""),width = 7,height=5)
