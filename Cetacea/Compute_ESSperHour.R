# Load necessary libraries
library(readr)
library(lubridate)
library(coda)

setwd("~/Nextcloud/Recherche/1_Methods/INSANE/OBDP/Cetacea_OBD/")

# Extract timestamps and calculate duration
calculate_duration <- function(log_file) {
  # Read all lines from the log file
  log_lines <- read_lines(log_file)
  
  # Extract all lines that contain timestamps
  timestamp_lines <- grep("\\(\\d+\\.\\d+\\.\\d+\\) \\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}", log_lines, value = TRUE)
  
  # Extract the timestamps from these lines
  timestamps <- sub(".*\\((\\d+\\.\\d+\\.\\d+)\\) (\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}).*", "\\2", timestamp_lines)
  
  # Convert the timestamps to POSIXct format
  time_values <- ymd_hms(timestamps)
  
  # Check if "Job terminated" is in the log file
  if (any(grepl("Job terminated|Job exceeded its runtime limit", log_lines))) {
    # Extract the timestamp of the job termination
    end_time <- max(time_values)
  } else {
    # If the job is still running, use the current time
    end_time <- Sys.time()
  }
  
  # Calculate the duration between the first timestamp and the end time in hours
  duration <- difftime(end_time, min(time_values), units = "days")
  
  return(as.numeric(duration))
}


# Compute ESS from a log file
compute_ess <- function(trace_file) {
  # Read the trace data from the log file
  trace_data <- read.table(trace_file, header = TRUE)
  trace_data_reduced <- trace_data[!grepl("iteration|div|ne|ns", names(trace_data))]
  trace_data_reduced <- trace_data_reduced[colSums(sapply(trace_data_reduced, diff)) != 0]
  
  # Compute the ESS
  #ess_value <- effectiveSize(as.mcmc(trace_data$likelihood))
  ess_value <- mean(sapply(trace_data_reduced, effectiveSize))
  
  return(as.numeric(ess_value))
}


# Compute the ESS/day for each replicates
seeds <- 0:4

get_essperday <- function(log_path, trace_path){
  results <- sapply(seeds, function(seed) {
    log_file <- paste0(log_path, seed, ".log")
    trace_file <- paste0(trace_path, seed, ".log")
    
    ess <- compute_ess(trace_file)
    duration <- calculate_duration(log_file)
    ess_per_day <- ess / duration
    
    # Append the results to the data frame
    data.frame(Seed = seed, Duration_d = duration, ESS = ess, ESS_per_day = ess_per_day)
  })
  return (cbind(results, c("Average", apply(results[-1,], 1, function(x)(mean(as.numeric(x)))))))
}

## CBD
results_CBD <- get_essperday("Logs/CBD.14464.", "outputs/CBD_Cetaceans_2000000iter_seed")
results_CBD

## CFBD
results_CFBD <- get_essperday("Logs/CFBD.14462.", "outputs/CFBD_Cetaceans_2000000iter_seed")
results_CFBD

## COBD
results_COBD <- get_essperday("Logs/COBD.14463.", "outputs/COBD_Cetaceans_2000000iter_seed")
results_COBD

## COBD_rmFossils
results_COBD_rmFossils <- get_essperday("Logs/COBD_rmFossils.14485.", "outputs/COBD_rmFossils_Cetaceans_50000000iter_seed")
results_COBD_rmFossils

## BDD
results_BDD <- get_essperday("Logs/BDD.14472.", "outputs/BDD_Cetaceans_10000000iter_seed")
results_BDD

## FBDD
results_FBDD <- get_essperday("Logs/FBDD.14467.", "outputs/FBDD_Cetaceans_rateαprior0.1_10000000iter_seed")
results_FBDD

## OBDD
results_OBDD <- get_essperday("Logs/OBDD.14466.", "outputs/OBDD_Cetaceans_rateαprior0.1_10000000iter_seed")
results_OBDD

## OBDD_rmFossils
results_OBDD_rmFossils <- get_essperday("Logs/OBDD_rmFossils.14460.", "outputs/OBDD_rmFossils_Cetaceans_rateαprior0.1_20000000iter_seed")
results_OBDD_rmFossils

