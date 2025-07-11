read_data <- function(data_path){
    logger::log_info(paste("Reading data from",data_path))
    readRDS(data_path)
}