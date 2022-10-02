#' GeneCodeR - Perturbing imaging tissue by altering spatial gene codes
#'
#' @param model - the learned model via Generative Encoding
#' @param x - the dataset to be transformed
#' @param config - the GeneCodeR configuration parameters to define the modality coming from and going to
#'
#' @export
genecoder <- function(model,x,model_type,config){

  transform.info <- genecoder.transform(model,x,model_type,config$transform)

  return(transform.info)

}


#' genecoder.transform - transform between modalities
#'
#' @param model - the learned model via Generative Encoding
#' @param x - the dataset to be transformed
#' @param config.transform - the transform list within the GeneCodeR configuration parameters to define the modality coming from and going to
#'
#' @export
genecoder.transform <- function(model,x,model_type,config.transform){

  if (model_type == "gcode"){
    beta.2_signal <- (model$main.parameters$beta_signal[[config.transform$from]])
    beta.1_signal <- (model$main.parameters$beta_signal[[config.transform$to]])

    beta.2_sample <- (model$main.parameters$beta_sample[[config.transform$from]])
    beta.1_sample <- (model$main.parameters$beta_sample[[config.transform$to]])

    signal2sample_obs <- as.matrix((x) - c(model$main.parameters$intercept[[config.transform$from]]))%*%beta.2_signal%*%MASS::ginv(t(beta.2_signal)%*%(beta.2_signal))%*%t(beta.1_sample) + c(model$main.parameters$intercept[[config.transform$to]])
    sample2signal_obs <- as.matrix((x) - c(model$main.parameters$intercept[[config.transform$from]]))%*%beta.2_sample%*%MASS::ginv(t(beta.2_sample)%*%(beta.2_sample))%*%t(beta.1_signal) + c(model$main.parameters$intercept[[config.transform$to]])

    sample_obs <- as.matrix((x) - c(model$main.parameters$intercept[[config.transform$from]]))%*%beta.2_sample%*%MASS::ginv(t(beta.2_sample)%*%(beta.2_sample))%*%t(beta.1_sample) + c(model$main.parameters$intercept[[config.transform$to]])
    signal_obs <- as.matrix((x) - c(model$main.parameters$intercept[[config.transform$from]]))%*%beta.2_signal%*%MASS::ginv(t(beta.2_signal)%*%(beta.2_signal))%*%t(beta.1_signal) + c(model$main.parameters$intercept[[config.transform$to]])


    return((sample_obs + signal_obs + sample2signal_obs + signal2sample_obs)/4)
  }

}


#' learn_model - learns the Generative Encoder model
#'
#' @param data_list List of data matrices of varying dimensionality. Attempts to find similarities among all datasets with a core structure.
#' @param config Configuration parameters (required, default provided)
#' @param transfer Transferring pre-trained model parameters (not required)
#' @param recover Important information used for prediction or imputation (not required)
#' @param join Join parameters that share the same axes to be similar (not required)
#'
#'#' @export
learn_model <- function(data_list,
                        config = gcode::extract_config(F),
                        transfer = gcode::extract_transfer_framework(F),
                        recover = gcode::extract_recovery_framework(F),
                        join,
                        reference){

  gcode.model <- gcode::gcode(data_list = data_list, config = config, transfer = transfer, recover = recover, join = join, reference = reference)

  return(gcode.model)

}


#' prepare_gex - prepares gene expression data
#'
#' @param file_path_list - a list of files containing the path to all files including the gene expression path
#' @param meta_info_list - a list of meta information that helps read in files and extracts relevant information
#' @param config - the main configuration parameters to extract the image pixels for each corresponding spot
#'
#' @export
prepare_gex <- function(file_path_list,
                        meta_info_list,
                        config
){

  if (!do.call("identical",list(length(file_path_list$coord$path),
                                length(file_path_list$gex$path),
                                length(file_path_list$pixel$path)))){
    print("Error because: length of file paths not equal between coord, gex, pixel")
    break
  }

  print("Extracting gex")
  gex_data <- lapply(c(1:length(file_path_list$gex$path)),function(X){
    print(paste("Preparing gex      ",X,sep=""))

    gex_list <- read_file(file_path_list$gex$path[X],meta_info_list$gex$read_file)[[1]]

    select_ids <- unique(meta_info_list$coord$factor$labels[[X]][meta_info_list$coord$factor$labels[[X]] %in% row.names(gex_list)])
    main_data <- gex_list[select_ids,meta_info_list$gex$factor$common_genes]

    return(list(gex=main_data,labels=row.names(main_data)))
  })

  print("Done preparation!")
  return(list(gex=lapply(gex_data,function(X){X$gex}),labels=lapply(gex_data,function(X){X$labels})))

}



#' prepare_spot - prepares image information per spot
#'
#' @param file_path_list - a list of files containing the path to all files including the image path
#' @param meta_info_list - a list of meta information that helps read in files and extracts relevant information
#' @param config - the main configuration parameters to extract the image pixels for each corresponding spot
#'
#' @export
prepare_spot <- function(file_path_list,
                         meta_info_list,
                         config,
                         gex_data
){

  if (!do.call("identical",list(length(file_path_list$coord$path),
                                length(file_path_list$gex$path),
                                length(file_path_list$pixel$path)))){
    print("Error because: length of file paths not equal between coord, gex, pixel")
    break
  }

  print("Extracting spots")
  spot_data <- extract_spots(file_path_list,meta_info_list,config$extract_spots,gex_data)

  print("Done preparation!")
  return(spot_data)

}



#' @export
extract_spots <- function(file_path_list,meta_info_list,config,gex_data){
  all_meta_data <- c()
  all_coord_data <- c()
  all_spot_data <- c()
  all_gex_data <- c()
  for (i in 1:length(file_path_list$coord$path)){
    print(paste("Preparing spot      ",i,sep=""))
    meta_data <- read_file(file_path_list$meta$path[i],meta_info_list$meta$read_file)[[1]]
    coord_data <- read_file(file_path_list$coord$path[i],meta_info_list$coord$read_file)[[1]]
    image_data <- read_file(file_path_list$pixel$path[i],meta_info_list$pixel$read_file)[[1]]

    labels <- if(meta_info_list$coord$factor_id$labels==0){row.names(coord_data)}else{coord_data[,meta_info_list$coord$factor_id$labels]}
    coord_id <- c(meta_info_list$coord$factor_id$coord_x,meta_info_list$coord$factor_id$coord_y)
    spot_data <- extract_pixels(image_data,coord_data[,coord_id],displacement_x = config$displacement_x,displacement_y = config$displacement_y,rotation = config$rotation, window = config$window_size)
    row.names(spot_data) <- labels

    row.names(meta_data) <- meta_coords <- apply(do.call('rbind',lapply(strsplit(x = row.names(meta_data),split = "_"),function(X){X[2:3]})),1,function(X){paste0(X,collapse="x")})
    coord_coords <- row.names(coord_data)
    spot_coords <- row.names(spot_data)
    gex_coords <- meta_info_list$gex$factor$labels[[i]]

    all_coords <- Reduce("intersect",list(meta_coords,coord_coords,spot_coords,gex_coords))

    all_coord_data <- rbind(all_coord_data, coord_data[all_coords,])
    all_spot_data <- rbind(all_spot_data, spot_data[all_coords,])
    all_meta_data <- rbind(all_meta_data, meta_data[all_coords,])
    all_gex_data <- rbind(all_gex_data, gex_data[[i]][all_coords,])

  }

  return(list(spot = as.matrix(all_spot_data),coord = as.matrix(all_coord_data), meta = all_meta_data, gex = all_gex_data))

}

#' @export
extract_pixels <- function(image,coords,displacement_x=0,displacement_y=0,rotation=0,window=32){

  copy_coords <- coords
  copy_coords[,1] <- copy_coords[,1] + window/2 + displacement_x
  copy_coords[,2] <- copy_coords[,2] + window/2 + displacement_y

  coord_list <- do.call('rbind',pbmcapply::pbmclapply(c(1:dim(copy_coords)[1]),function(i){

    cropped_image <- magick::image_crop(image = image, geometry = paste(window,"x",window,"+",copy_coords[i,1],"+",copy_coords[i,2],sep=""))
    coord_pixels <- as.numeric(magick::image_data(cropped_image, 'rgb'))

    cropped_image <- magick::image_read(coord_pixels)
    cropped_image <- magick::image_rotate(cropped_image,degrees = rotation)
    cropped_image <- magick::image_crop(image = cropped_image, geometry = paste(window,"x",window,"+",dim(cropped_image)[1]/2,"+",dim(cropped_image)[2]/2,sep=""))

    coord_pixels <- (as.numeric(magick::image_data(cropped_image, 'rgb')))
    coord_pixels[coord_pixels==1] <- 0
    coord_pixels
  },mc.cores = 8))

  return(coord_list)
}


#' @export
read_file <- function(path,meta_info){

  if (!grepl(meta_info$file_type,pattern = "csv|tsv|txt|mtx|image")){
    print("Error because: file type exists but invalid - use one of csv,tsv,txt,mtx")
    break
  }
  if (is.null(meta_info$file_type)){
    print("Error because: file type is null")
    break
  }

  if (meta_info$file_type == "image"){

    return(lapply(path,function(X){

      magick::image_read(path)

    }))

  }

  if (meta_info$file_type == "csv"){

    return(lapply(path,function(X){

      read.csv(X,sep=",",header = meta_info$meta$header, quote = meta_info$meta$quote, row.names = meta_info$meta$row.names)

    }))

  }

  if (meta_info$file_type == "txt" | meta_info$file_type == "tsv"){

    return(lapply(path,function(X){

      read.table(X,sep="\t",header = meta_info$meta$header, quote = meta_info$meta$quote, row.names = meta_info$meta$row.names)
    }))

  }

  if (meta_info$file_type == "mtx"){

    return(lapply(path,function(X){

      Matrix::readMM(X)

    }))

  }

}



#' @export
convert_to_RGB <- function(x){
  e <- ecdf(x)
  j <- e(x)
  x <- array(j,dim(x))
  return(x)
}



#' @export
max_min_transform <- function(x){
  x <- (x+abs(min(x)))/max(x+abs(min(x)))
  return(x)
}


#' @export
plot_image <- function(vector_image){

  library(ggplot2)
  library(grid)
  library(gridExtra)


  g <- convert_to_RGB(array((vector_image),dim=c(sqrt(length(vector_image)/3),sqrt(length(vector_image)/3),3)))

  g <- magick::image_read(g)
  g <- magick::image_rotate(g,degrees = -90)
  g <- (as.numeric(magick::image_data(g, 'rgb')))

  g <- rasterGrob(g, interpolate=TRUE)
  g_plots <- ggplot2::qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
                   axis.text.y=element_blank(),axis.ticks=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),legend.position="none",
                   panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),plot.background=element_blank()) +
    theme(plot.title = element_text(size = 16, face = "bold", color="black", hjust=0.5))

  return(g_plots)

}
