#' @export
extract_config_framework <- function(verbose=F){

  config <- list(
    extract_spots = list(window_size = 64,
                         displacement_x = 0,
                         displacement_y = 0,
                         rotation = 0
    ),
    transform = list(from = NULL,
                     to = NULL
    )
  )

  if (verbose){
    print(config)
  }

  return(config)
}


#' @export
extract_meta_framework <- function(verbose=F){

  meta_info <- list(
    coord = list(
      read_file = list(
        file_type = NULL,
        meta = list(
          header = NULL,
          sep = NULL,
          quote = NULL,
          row.names = NULL
        )
      ),
      factor=list(
        labels=NULL
      ),
      factor_id=list(
        coord_x = NULL,
        coord_y = NULL,
        labels=NULL
      )
    ),
    gex = list(
      read_file = list(
        file_type = NULL,
        meta = list(
          header = NULL,
          sep = NULL,
          quote = NULL,
          row.names = NULL
        )
      ),
      factor=list(
        common_genes=NULL
      ),
      factor_id=list(
        labels=NULL
      )
    ),
    pixel = list(
      read_file = list(
        file_type = NULL,
        meta = list(
          header = NULL,
          sep = NULL,
          quote = NULL,
          row.names = NULL
        )
      ),
      factor=list(
        common_genes=NULL
      ),
      factor_id=list(
        labels=NULL
      )
    )
  )

  if (verbose){
    print(meta_info)
  }

  return(meta_info)
}



#' @export
extract_path_framework <- function(verbose=F){

  file_paths <- list(
    coord = list(path=list(NULL)),
    gex = list(path=list(NULL)),
    pixel = list(path=list(NULL))
  )

  if (verbose){
    print(file_paths)
  }

  return(file_paths)
}

