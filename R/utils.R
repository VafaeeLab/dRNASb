check_and_create_directory <- function(dir_path){
  if(!dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE)
  }
}
