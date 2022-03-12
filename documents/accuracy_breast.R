# t.test for difference between observations

t.test_difference_obs <- function(true_obs,a,b){
  t.test(c(a - true_obs),c(b - true_obs))
}



library(GeneCodeR)
library(gcode)


# Input list of files

file_path_list <- GeneCodeR::extract_path_framework(F)

file_path_list$coord$path <- list.files("~/Documents/main_files/AskExplain/GeneCodeR/data/external/breast/",pattern = ".csv.gz",full.names = T)
file_path_list$gex$path <- list.files("~/Documents/main_files/AskExplain/GeneCodeR/data/external/breast/",pattern = "stdata.tsv.gz",full.names = T)
file_path_list$pixel$path <- list.files("~/Documents/main_files/AskExplain/GeneCodeR/data/external/breast/",pattern = ".jpg",full.names = T)

set.seed(1)
train_ids <- sample(c(1:length(file_path_list$coord$path)),round(length(file_path_list$coord$path)*0.6,0))
test_ids <- c(1:length(file_path_list$coord$path))[-train_ids]

train_file_path_list <- lapply(file_path_list,function(X){
  list(path=X$path[train_ids])
})
test_file_path_list <- lapply(file_path_list,function(X){
  list(path=X$path[-train_ids])
})

# Set up genecoder configuration parameters

genecoder.config <- GeneCodeR::extract_config_framework(F)

# Set up meta information

meta_info_list <- GeneCodeR::extract_meta_framework(F)
meta_info_list$coord$read_file$file_type <- "csv"
meta_info_list$gex$read_file$file_type <- "tsv"
meta_info_list$pixel$read_file$file_type <- "image"

meta_info_list$coord$read_file$meta$header <- meta_info_list$gex$read_file$meta$header <- T
meta_info_list$coord$read_file$meta$sep <- ","
meta_info_list$gex$read_file$meta$sep <- "\t"
meta_info_list$coord$read_file$meta$quote <- meta_info_list$gex$read_file$meta$quote <- ""
meta_info_list$coord$read_file$meta$row.names <- 1
meta_info_list$gex$read_file$meta$row.names <- 1

meta_info_list$coord$factor_id$labels <- 0
meta_info_list$coord$factor_id$coord_x <- 1
meta_info_list$coord$factor_id$coord_y <- 2


common_genes <- lapply(c(1:length(file_path_list$gex$path)),function(X){
  print(X)
  colnames(GeneCodeR::read_file(path = file_path_list$gex$path[X], meta_info = meta_info_list$gex$read_file)[[1]])
})

meta_info_list$gex$factor$common_genes <- Reduce("intersect",common_genes)


meta_info_list$coord$factor$labels <- lapply(c(1:length(train_file_path_list$coord$path)),function(X){
  row.names(GeneCodeR::read_file(path = train_file_path_list$coord$path[X], meta_info = meta_info_list$coord$read_file)[[1]])
})


# Extract gex data
gex_data <- GeneCodeR::prepare_gex(file_path_list = train_file_path_list,meta_info_list = meta_info_list,config = genecoder.config)
meta_info_list$gex$factor$labels <- gex_data$labels


# Extract spot data

spot_data <- GeneCodeR::prepare_spot(file_path_list = train_file_path_list,meta_info_list = meta_info_list,config = genecoder.config)

# Set up join information

join <- gcode::extract_join_framework(F)
join$complete$alpha <- c("all","all")
join$complete$data_list <- c("gex","spot")
join$complete$code <- c("all","all")
join$complete$alpha.code <- c("all","all")
join$complete$beta.code <- c("all","all")
join$complete$incode <- c("gex","spot")
join$complete$beta <- c("gex","spot")



# Set up reference information

reference <- gcode::extract_references_framework(F)
reference$data_list <- c(1,0)


# Set up data

data_list <- list(gex=gex_data$gex,spot=spot_data$spot)


# Set up gcode config

gcode.config <- gcode::extract_config(F)
gcode.config$init <- c("irlba","irlba")
gcode.config$k_dim <- 150
gcode.config$i_dim <- 300
gcode.config$j_dim <- 300


# Set up gcode model

genecoder.model <- GeneCodeR::learn_model(data_list = data_list, config = gcode.config, join = join, reference = reference)





# Evaluate on test set

meta_info_list$coord$factor$labels <- lapply(c(1:length(test_file_path_list$coord$path)),function(X){
  row.names(GeneCodeR::read_file(path = test_file_path_list$coord$path[X], meta_info = meta_info_list$coord$read_file)[[1]])
})


# Extract test gex data
test_gex_data <- GeneCodeR::prepare_gex(file_path_list = test_file_path_list,meta_info_list = meta_info_list,config = genecoder.config)
meta_info_list$gex$factor$labels <- test_gex_data$labels


# save(genecoder.model,file = "~/Documents/main_files/AskExplain/GeneCoder/data/workflow/sandbox/genecoder.RData")
# save(test_gex_data,file = "~/Documents/main_files/AskExplain/GeneCoder/data/workflow/sandbox/test_gex.RData")
# save(meta_info_list,file = "~/Documents/main_files/AskExplain/GeneCoder/data/workflow/sandbox/meta_info_list.RData")

rm(list=ls());gc()

# Set up genecoder transform information
genecoder.config <- GeneCodeR::extract_config_framework(F)
genecoder.config$transform$from <- "spot"
genecoder.config$transform$to <- "gex"

load("~/Documents/main_files/AskExplain/GeneCoder/data/workflow/sandbox/genecoder.RData")
load("~/Documents/main_files/AskExplain/GeneCoder/data/workflow/sandbox/test_gex.RData")
load("~/Documents/main_files/AskExplain/GeneCoder/data/workflow/sandbox/meta_info_list.RData")


# Base testing

# Extract test spot data

base_test_spot_data <- GeneCodeR::prepare_spot(file_path_list = test_file_path_list,meta_info_list = meta_info_list,config = genecoder.config)

base_spot2gex <- GeneCodeR::genecoder(model=genecoder.model, x = base_test_spot_data$spot, config = genecoder.config)

base_spot2gex <- list(sample_wise = do.call('c',parallel::mclapply(c(1:dim(base_spot2gex)[1]),function(X){lsa::cosine(base_spot2gex[X,],test_gex_data$gex[X,])},mc.cores = 8)),
                      gene_wise = do.call('c',parallel::mclapply(c(1:dim(base_spot2gex)[2]),function(X){lsa::cosine(base_spot2gex[,X],test_gex_data$gex[,X])},mc.cores = 8))
)

save(base_spot2gex,file = "~/Documents/main_files/AskExplain/GeneCoder/data/workflow/breast/base_spot2gex.RData")

# Spatial rotation testing

rotate_spot2gex <- list()
for (rotate_val in c(0,90,180,270)){
  genecoder.config$extract_spots$rotation <- rotate_val

  rotate_test_spot_data <- GeneCodeR::prepare_spot(file_path_list = test_file_path_list,meta_info_list = meta_info_list,config = genecoder.config)

  rotate_spot2gex[[as.character(rotate_val)]] <- GeneCodeR::genecoder(model=genecoder.model, x = rotate_test_spot_data$spot, config = genecoder.config)
}


rotate_spot2gex <- lapply(c(2:4),function(X){
  t.test_difference_obs(true_obs = test_gex_data$gex, a = rotate_spot2gex[[1]], b = rotate_spot2gex[[X]])
})

save(rotate_spot2gex,file = "~/Documents/main_files/AskExplain/GeneCoder/data/workflow/breast/rotate_spot2gex.RData")





# Spatial displacement testing

displace_spot2gex <- list()
for (displace_val in seq(0,30,10)){
  genecoder.config$extract_spots$rotation <- 0
  genecoder.config$extract_spots$displacement_x <- displace_val
  genecoder.config$extract_spots$displacement_y <- displace_val

  displace_test_spot_data <- GeneCodeR::prepare_spot(file_path_list = test_file_path_list,meta_info_list = meta_info_list,config = genecoder.config)

  displace_spot2gex[[as.character(displace_val)]] <- GeneCodeR::genecoder(model=genecoder.model, x = displace_test_spot_data$spot, config = genecoder.config)
}

displace_spot2gex <- lapply(c(2:4),function(X){
  t.test_difference_obs(true_obs = test_gex_data$gex, a = displace_spot2gex[[1]], b = displace_spot2gex[[X]])
})

save(displace_spot2gex,file = "~/Documents/main_files/AskExplain/GeneCoder/data/workflow/breast/displace_spot2gex.RData")


# Spatial similarity testing

similar_spot2gex <- list()
for (i in c(1,2,3,4)){

  if (i==1){
    similar_val.x = -3
    similar_val.y = -3
  }
  if (i==2){
    similar_val.x = 3
    similar_val.y = -3
  }
  if (i==3){
    similar_val.x = -3
    similar_val.y = 3
  }
  if (i==4){
    similar_val.x = 3
    similar_val.y = 3
  }

  genecoder.config$extract_spots$displacement_x <- similar_val.x
  genecoder.config$extract_spots$displacement_y <- similar_val.y

  similar_test_spot_data <- GeneCodeR::prepare_spot(file_path_list = test_file_path_list,meta_info_list = meta_info_list,config = genecoder.config)

  similar_spot2gex[[as.character(i)]] <- GeneCodeR::genecoder(model=genecoder.model, x = similar_test_spot_data$spot, config = genecoder.config)
}


similar_spot2gex <- lapply(c(2:4),function(X){
  t.test_difference_obs(true_obs = test_gex_data$gex, a = similar_spot2gex[[1]], b = similar_spot2gex[[X]])
})

save(similar_spot2gex,file = "~/Documents/main_files/AskExplain/GeneCoder/data/workflow/breast/similar_spot2gex.RData")





# Spatial signal testing

signal_spot2gex <- list()
for (i in c(1:10)){

  genecoder.config$extract_spots$displacement_x <- i
  genecoder.config$extract_spots$displacement_y <- 0

  signal_test_spot_data <- GeneCodeR::prepare_spot(file_path_list = test_file_path_list,meta_info_list = meta_info_list,config = genecoder.config)

  signal_spot2gex[[as.character(i)]] <- GeneCodeR::genecoder(model=genecoder.model, x = signal_test_spot_data$spot, config = genecoder.config)
}



signal_spot2gex <- lapply(c(1:4),function(X){
  t.test_difference_obs(true_obs = test_gex_data$gex, a = signal_spot2gex[[X]], b = signal_spot2gex[[X+1]])
})


plot_signal <- FALSE
if (plot_signal){

  signal_spline <- data.frame(gene_level = do.call('c',lapply(c(1:10),function(X){
    signal_spot2gex[[X]][1,1]
  })),displacement=c(1:10))

  mod.ss <- npreg::summary.ss(npreg::ss(x = signal_spline$gene_level,y = signal_spline$displacement))

  plot(signal_spline$gene_level ~ signal_spline$displacement, main = paste("Adj.R-squared:    ",mod.ss$adj.r.squared))
  lines(signal_spline$gene_level, mod.ss$y, lty = 2, col = 2, lwd = 2)

}

spline_signal_adj.r.squared <- do.call('rbind',lapply(c(order(apply(test_gex_data$gex,2,var),decreasing = T)[1:3]),function(Y){
  internal.adj.r2 <- do.call('c',pbmcapply::pbmclapply(c(1:dim(signal_spot2gex[[1]])[1]),function(Z){
    signal_spline <- data.frame(gene_level = do.call('c',lapply(c(1:10),function(X){
      signal_spot2gex[[X]][Z,Y]
    })),displacement=c(1:10))
    adj.r.squared <- try(npreg::summary.ss(npreg::ss(x = signal_spline$gene_level,y = signal_spline$displacement))$adj.r.squared,silent = F)
    if (!is.character(signal_spline)){
      return(as.numeric(adj.r.squared))
    }
  },mc.cores = 8))
}))

save(spline_signal_adj.r.squared,file = "~/Documents/main_files/AskExplain/GeneCoder/data/workflow/breast/spline_signal_adj.r.squared.RData")




