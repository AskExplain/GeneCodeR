#' extract_config_framework - the configuration framework for GeneCodeR
#'
#' @param verbose - set to TRUE to print out the config parameters, to be editted
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

#' extract_meta_framework - the meta information framework for GeneCodeR
#'
#' @param verbose - set to TRUE to print out the meta information, to be editted
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




#' extract_path_framework - the path list framework for GeneCodeR
#'
#' @param verbose - set to TRUE to print out the path list, to be editted
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

