###############################################
## 0.  Libraries & global options #############
###############################################
library(httr)          # downloads
library(fs)            # directories / paths
library(data.table)    # fast wrangling
library(haven)         # read .dta
library(progress)      # pretty download bar
library(arrow)         # parquet output
library(purrr)

DATA_DIR <- file.path(dataIn, "sipp_raw")
UNZIP_DIR <- file.path(dataIn, "sipp_unzipped")
DIR.create <- function(x) if (!dir.exists(x)) dir.create(x, recursive = TRUE)
DIR.create(DATA_DIR); DIR.create(UNZIP_DIR)

###############################################
## 1.  Years & URLs ###########################
###############################################
panels   <- 2018:2024       # 2018 panel = ref-year 2017, etc. :contentReference[oaicite:0]{index=0}
url_core <- sprintf(
  "https://www2.census.gov/programs-surveys/sipp/data/datasets/%d/pu%d_dta.zip",
  panels, panels
)
url_rw   <- sub("pu", "rw", url_core)   # replicate-weights (optional)

###############################################
## 2.  Helper: download + unzip ###############
###############################################
download_unzip <- function(url) {
  zipfile <- file.path(DATA_DIR, basename(url))
  if (!file.exists(zipfile)) {
    GET(url, write_disk(zipfile, overwrite = TRUE), timeout(600))
  }
  unzip(zipfile, exdir = UNZIP_DIR, junkpaths = TRUE)
}

pb <- progress_bar$new(total = length(c(url_core, url_rw)),
                       format = "  downloading :what [:bar] :percent")
walk(c(url_core, url_rw),
     ~{ pb$tick(tokens = list(what = basename(.x))); download_unzip(.x) })



###############################################
## 1.  Years & URLs (2008-2024 panels) ########
###############################################
panels_new  <- 2018:2024                         # modern one-file panels
panels_2008 <- 2008                              # legacy multi-wave panel
panels_2014 <- 2014                              # legacy multi-wave panel

## ---- modern panels ----------------------------------------------------------
url_core_new <- sprintf(
  "https://www2.census.gov/programs-surveys/sipp/data/datasets/%d/pu%d_dta.zip",
  panels_new, panels_new)
url_rw_new   <- sub("pu", "rw", url_core_new)

## ---- 2008 panel (waves 1-16) -------------------------------------------------
waves_08     <- 1:16
url_core_08  <- sprintf(
  "https://www2.census.gov/programs-surveys/sipp/data/datasets/2008/w%d/l08puw%d.zip",
  waves_08, waves_08)
# url_rw_08    <- sprintf(
#   "https://www2.census.gov/programs-surveys/sipp/data/datasets/2008/w%d/rw08w%d.zip",
#   waves_08, waves_08)
# 
# ## ---- 2014 panel (waves 1-4) --------------------------------------------------
# waves_14     <- 1:4                # change to 1:5, 1:6 … as more waves release
# url_core_14  <- sprintf(
#   "https://www2.census.gov/programs-surveys/sipp/data/datasets/2014/w%d/pu2014w%d_v13.zip",
#   waves_14, waves_14)
# # url_rw_14    <- sprintf(
# #   "https://www2.census.gov/programs-surveys/sipp/data/datasets/2014/w%d/rw14w%d_v13.zip",
# #   waves_14, waves_14)

## ---- full download list ------------------------------------------------------
all_urls <- c(url_core_new, url_rw_new,
              url_core_14,  url_rw_14,
              url_core_08,  url_rw_08)

all_urls <- c(url_core_14,  url_rw_14,
              url_core_08,  url_rw_08)

###############################################
## 2.  Helper: download + unzip ###############
###############################################
download_unzip <- function(url) {
  zipfile <- file.path(DATA_DIR, basename(url))
  if (!file.exists(zipfile)) {
    GET(url,
        write_disk(zipfile, overwrite = TRUE),
        timeout(600))
  }
  unzip(zipfile, exdir = UNZIP_DIR, junkpaths = TRUE)
}

pb <- progress_bar$new(
  total   = length(all_urls),
  format  = "  downloading :what [:bar] :percent  eta: :eta")

walk(all_urls, ~{
  pb$tick(tokens = list(what = basename(.x)))
  download_unzip(.x)
})



## parse older files
############################################################
##  Helper: read one legacy SIPP ASCII file  (1996-2008)  ##
############################################################
read_sipp_ascii <- function(dat_path,
                            cache_sas  = TRUE,
                            beginline  = 1,
                            buffersize = 1e5,
                            verbose    = TRUE, ...) {
  
  stopifnot(file.exists(dat_path))
  library(SAScii); library(data.table)
  
  fn <- tolower(basename(dat_path))
  
  ## Option B pattern (plain groups)
  rx <- "^l([0-9]{2})([a-z]{2})w([0-9]+)\\.dat$"
  
  pieces <- strcapture(
    rx, fn,
    data.frame(yr   = integer(),
               type = character(),
               wave = integer())
  )
  
  if (any(is.na(pieces)))
    stop("File name ‘", fn, "’ does not match the lYYttwW.dat convention.")
  
  # --- download SAS script if absent ---------------------------------------
  sas_name <- sprintf("l%02d%sw%d.sas",
                      pieces$yr, pieces$type, pieces$wave)
  sas_path <- file.path(dirname(dat_path), sas_name)
  
  if (!file.exists(sas_path)) {
    if (verbose) message("‣ SAS script not found → downloading ‘", sas_name, "’")
    base <- sprintf("https://www2.census.gov/programs-surveys/sipp/data/datasets/20%02d",
                    pieces$yr)
    url  <- file.path(base, sas_name)
    
    tf <- tempfile(fileext = ".sas")
    tryCatch(
      utils::download.file(url, tf, mode = "wb", quiet = !verbose),
      error = function(e) stop("Could not download SAS script from ", url)
    )
    if (cache_sas) file.copy(tf, sas_path)
    else           sas_path <- tf
  }
  
  # --- read the ASCII file ---------------------------------------------------
  if (verbose) message("‣ Reading ", fn)
  dat <- read.SAScii(dat_path, sas_path,
                     beginline  = beginline,
                     buffersize = buffersize) #  ...)
  
  setDT(dat)
  dat[, wave := pieces$wave]   # add wave column
  dat
}


#######################################################################
##  Batch-import: read every l08puw*.dat that is already on disk     ##
#######################################################################
read_sipp08_folder <- function(folder, file_type = "pu", ...) {
  pattern <- sprintf("^l08%sw\\d+\\.dat$", file_type)
  dat_files <- list.files(folder, pattern = pattern,
                          full.names = TRUE, ignore.case = TRUE)
  if (!length(dat_files)) stop("No matching .dat files found in ", folder)
  
  # read every wave and combine
  test <- read_sipp_ascii(dat_files[1])
  out <- lapply(dat_files, read_sipp_ascii)
  # add a wave column before binding
  out <- Map(function(dt, w) { dt[, wave := w]; dt },
             out, vapply(out, attr, 0L, "wave"))
  
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}

#######################################################################
##  Example                                                           ##
#######################################################################
folder <- file.path(dataIn, "sipp_unzipped", "older_files")

library(data.table)
sipp_2008_core <- read_sipp08_folder(folder, file_type = "pu")

# If you also downloaded replicate weights (rw) files, run:
# sipp_2008_rw <- read_sipp08_folder(folder, file_type = "rw")
