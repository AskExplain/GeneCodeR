#' GeneCodeR - Perturbing imaging tissue by altering spatial gene codes
#'
#' @param model - the learned model via Generative Encoding
#' @param x - the dataset to be transformed
#' @param config - the GeneCodeR configuration parameters to define the modality coming from and going to
#'
#' @export
genecoder <- function(model,x,config){

  transform.info <- genecoder.transform(model,x,config$transform)

  return(transform.info)

}


#' genecoder.transform - transform between modalities
#'
#' @param model - the learned model via Generative Encoding
#' @param x - the dataset to be transformed
#' @param config.transform - the transform list within the GeneCodeR configuration parameters to define the modality coming from and going to
#'
#' @export
genecoder.transform <- function(model,x,config.transform){

  beta.2 <- (model$main.parameters$beta[[config.transform$from]]%*%t(model$main.code$incode[[config.transform$from]]%*%model$main.parameters$beta.code[[1]]))
  beta.1 <- (model$main.parameters$beta[[config.transform$to]]%*%t(model$main.code$incode[[config.transform$to]]%*%model$main.parameters$beta.code[[1]]))

  return(as.matrix((x) - c(model$main.parameters$intercept[[config.transform$from]]))%*%beta.2%*%MASS::ginv(t(beta.2)%*%(beta.2))%*%t(beta.1) + c(model$main.parameters$intercept[[config.transform$to]]))

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

    select_ids <- meta_info_list$coord$factor$labels[[X]][meta_info_list$coord$factor$labels[[X]] %in% row.names(gex_list)]
    main_data <- gex_list[select_ids,meta_info_list$gex$factor$common_genes]

    return(list(gex=main_data,labels=row.names(main_data)))
  })

  print("Done preparation!")
  return(list(gex=do.call('rbind',lapply(gex_data,function(X){X$gex})),labels=lapply(gex_data,function(X){X$labels})))

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
                         config
){

  if (!do.call("identical",list(length(file_path_list$coord$path),
                                length(file_path_list$gex$path),
                                length(file_path_list$pixel$path)))){
    print("Error because: length of file paths not equal between coord, gex, pixel")
    break
  }

  print("Extracting spots")
  spot_data <- extract_spots(file_path_list,meta_info_list,config$extract_spots)

  print("Done preparation!")
  return(list(spot=spot_data))

}



#' @export
extract_spots <- function(file_path_list,meta_info_list,config){

  all_spot_data <- c()
  for (i in 1:length(file_path_list$coord$path)){
    print(paste("Preparing spot      ",i,sep=""))
    coord_data <- read_file(file_path_list$coord$path[i],meta_info_list$coord$read_file)[[1]]
    image_data <- read_file(file_path_list$pixel$path[i],meta_info_list$pixel$read_file)[[1]]

    labels <- if(meta_info_list$coord$factor_id$labels==0){row.names(coord_data)}else{coord_data[,meta_info_list$coord$factor_id$labels]}
    coord_id <- c(meta_info_list$coord$factor_id$coord_x,meta_info_list$coord$factor_id$coord_y)
    spot_data <- extract_pixels(image_data,coord_data[,coord_id],displacement_x = config$displacement_x,displacement_y = config$displacement_y,rotation = config$rotation, window = config$window_size)
    row.names(spot_data) <- labels
    spot_data <- spot_data[meta_info_list$gex$factor$labels[[i]],]

    all_spot_data <- rbind(all_spot_data,spot_data)
  }

  return(all_spot_data)

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

    coord_pixels <- (as.numeric(magick::image_data(cropped_image, 'rgb')))

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

      as.matrix(read.csv(X,sep=",",header = meta_info$meta$header, quote = meta_info$meta$quote, row.names = meta_info$meta$row.names))

    }))

  }

  if (meta_info$file_type == "txt" | meta_info$file_type == "tsv"){

    return(lapply(path,function(X){

      as.matrix(read.table(X,sep="\t",header = meta_info$meta$header, quote = meta_info$meta$quote, row.names = meta_info$meta$row.names))

    }))

  }

  if (meta_info$file_type == "mtx"){

    return(lapply(path,function(X){

      as.matrix(Matrix::readMM(X))

    }))

  }

}

#' @export
plot_image <- function(vector_image){

  library(ggplot2)
  library(grid)
  library(gridExtra)

  convert_to_RGB <- function(x){
    e <- ecdf(x)
    j <- e(x)
    x <- array(j,dim(x))
    return(x)
  }

  g <- convert_to_RGB(array((vector_image),dim=c(sqrt(length(vector_image)/3),sqrt(length(vector_image)/3),3)))
  g <- rasterGrob(g, interpolate=TRUE)
  g_plots <- ggplot2::qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

  return(g_plots)

}
