###############################################################################
#  SIPP spell-panel (with class-of-worker, RMESR status & demographics)      #
#  Modified to only process years 2018-2023 as used in the paper             #
###############################################################################
library(data.table)
library(haven)
library(tidyselect)
library(lubridate)

# ---------------------------------------------------------------------------
pick_var <- function(avail, ...) {
  cands <- tolower(c(...))
  hit   <- match(TRUE, cands %in% tolower(avail))
  if (is.na(hit))
    stop("None of {", paste(c(...), collapse = ", "), "} in ", basename(parent.frame(2)$path))
  avail[ tolower(avail) == cands[hit] ][1]
}
extract_year_from_filename_simple <- function(filename_path) {
  filename <- basename(filename_path) # Get just the filename (e.g., "pu2014w1.dta")
  
  # Check if the filename starts with "pu" and is long enough
  if (startsWith(tolower(filename), "pu") && nchar(filename) >= 6) {
    year_str <- substr(filename, 3, 6) # Extract characters from 3rd to 6th position
    year_int <- suppressWarnings(as.integer(year_str))
    
    if (!is.na(year_int)) {
      return(year_int)
    }
  }
  
  # Fallback or warning if the pattern isn't met
  warning(paste("Could not extract year from filename using simple 'puYYYY' pattern:", filename))
  return(NA_integer_)
}
# ---------------------------------------------------------------------------
###############################################################################
##  read_one_panel()  – works for
##    • “classic” panels (1996…2008)                – MONTHCODE = 1…N counter
##    • “redesign” panels 2014–2017  **with OR without SWAVE/WAVE**
##    • “annual” panels 2018+ (no wave vars, fielded once per year)
###############################################################################
read_one_panel <- function(path){
  
  avail <- names(read_dta(path, n_max = 0))
  gv    <- function(...) pick_var(avail, ...)
  
  ## ------------ (1) identifiers -------------------------------------------
  ssuid     <- gv("SSUID")                       ; pnum   <- gv("PNUM")
  monthcode <- gv("MONTHCODE","TMONTHCODE")      ; spanel <- gv("SPANEL")
  
  ## ------------ (2) main-job vars -----------------------------------------
  jobid  <- gv("EJB1_JOBID","TJB1_JOBID")
  ind    <- gv("EJB1_IND"  ,"TJB1_IND")
  occ    <- gv("EJB1_OCC"  ,"TJB1_OCC")
  
  ## ------------ (3) pay fields --------------------------------------------
  annsal    <- gv("EJB1_ANNSAL1","TJB1_ANNSAL1")
  hourly    <- gv("EJB1_HOURLY1","TJB1_HOURLY1")
  weekly    <- gv("EJB1_WKLY1"  ,"TJB1_WKLY1")
  biweekly  <- gv("EJB1_BWKLY1" ,"TJB1_BWKLY1")
  semimonth <- gv("EJB1_SMTHLY1","TJB1_SMTHLY1")
  monthly   <- gv("EJB1_MTHLY1" ,"TJB1_MTHLY1")
  
  ## ------------ (4) job-start questions (only for incumbents) -------------
  start_yr_var <- gv("EJB1_STRTYR","TJB1_STRTYR","EJB1_STRYR","TJB1_STRYR")
  start_mo_var <- gv("EJB1_STRTJAN","EJB1_STRMON","TJB1_STRMON","TJB1_STRMN")
  
  ## ------------ (5) assorted stuff ----------------------------------------
  weight    <- gv("WPFINWGT")
  union     <- gv("EJB1_UNION","TJB1_UNION")
  emp_all   <- gv("TJB1_EMPALL","EJB1_EMPALL")
  emp_size  <- gv("EJB1_EMPSIZE","TJB1_EMPSIZE")
  clwrk     <- gv("EJB1_CLWRK","TJB1_CLWRK")
  mesr      <- gv("RMESR","EMESR","MESR")
  wfh_flag  <- gv("EJB1_WSHMWRK","TJB1_WSHMWRK")
  wfh_days  <- gv("EJB1_DYSWKDH","TJB1_DYSWKDH")
  st        <- gv("TEHC_ST","TFIPSST","STATEFIP")
  age_var   <- gv("EAGE","AGE","TAGE")
  sex_var   <- gv("ESEX","SEX","TPSEX")
  edu_var   <- gv("EEDUC","EDUC","TEEDUC")
  work_state <- gv("UHHSTAT", "RJB1_COMMTYP")
  wave_var <- gv("EHFQWAVE","TPEWAVE","SWAVE","WAVE","NONE")
  
  
  keep <- c(ssuid,pnum,monthcode,spanel,
            jobid,ind,occ,
            annsal,hourly,weekly,biweekly,semimonth,monthly,
            union,emp_all,emp_size,
            clwrk,mesr,
            wfh_flag,wfh_days,
            st,age_var,sex_var,edu_var,
            weight,
            start_yr_var,start_mo_var, work_state, wave_var)
  
  dt <- as.data.table(read_dta(path, col_select = all_of(keep)))
  
  ## ------------------------------------------------------------------------
  ##  Standardise names
  ## ------------------------------------------------------------------------
  setnames(dt, tolower(names(dt)))
  setnames(dt,
           old = tolower(c(jobid,ind,occ,union,
                           annsal,hourly,weekly,biweekly,semimonth,monthly,
                           emp_all,emp_size,
                           clwrk,mesr,
                           wfh_flag,wfh_days,st,
                           age_var,sex_var,edu_var,weight,
                           start_yr_var,start_mo_var, work_state, wave_var)),
           new = c("jobid","industry","occupation","union",
                   "annsal","hourly_wage","weekly","biweekly","semimonth","monthly",
                   "empall","empsize",
                   "class_worker","emp_status",
                   "wfh_flag","wfh_days","state_fips",
                   "age","sex","educ_level","weight",
                   "job_start_year","job_start_month", "work_state", "swave"),
           skip_absent = TRUE)
  
  ## ------------------------------------------------------------------------
  ##  (6) CALENDAR YEAR & MONTH  —— three regimes
  ## ------------------------------------------------------------------------
  panel_code <- extract_year_from_filename_simple(path) # 'path' is arg to read_one_panel

  
  
  
  
  # ---- A. classic panels ≤2008 -------------------------------------------
  if (panel_code < 2014L){
    panel_start_year <- panel_code - 1L
    offset           <- dt$monthcode - 1L
    dt[, `:=`(
      year  = panel_start_year + offset %/% 12L,
      month = (offset %% 12L) + 1L
    )]
    
    # ---- B. redesign 2014–2017  --------------------------------------------
  } else if (panel_code <= 2017L){
    
    dt[, `:=`(
      year  = (panel_code - 2L) + swave,   # wave1→ −2+1 = −1 (i.e. 2013 for 2014)
      month = monthcode
    )]
    
    # ---- C. annual panels 2018+  -------------------------------------------
  } else {
    dt[, `:=`(
      year  = panel_code - 1,
      month = monthcode
    )]
  }
  
  ## ------------------------------------------------------------------------
  ##  (7) FILL IN-PANEL START-DATES  (people who took the job during panel)
  ## ------------------------------------------------------------------------
  dt[
    is.na(job_start_year) | job_start_year %in% c(-4,-9,0),
    `:=`( job_start_year  = year,
          job_start_month = month )
  ]
  
  dt[, state_fips := as.character(state_fips)]
  dt[]
}



# ---------------------------------------------------------------------------
# Only process SIPP files from 2018-2023
sipp_files <- list.files(
  file.path(dataIn, "sipp_unzipped"),
  pattern = "pu(2018|2019|2020|2021|2022|2023).*\\.dta$",   # only 2018-2023 files
  full.names = TRUE,
  ignore.case = TRUE
)

dt <- rbindlist(lapply(sipp_files, read_one_panel), use.names = TRUE, fill = TRUE)
setorder(dt, ssuid, pnum, job_start_year, month)

# annual salary
dt[, annual_salary :=
     fifelse(!is.na(annsal), annsal,
             fifelse(!is.na(hourly_wage), hourly_wage*2080,
                     fifelse(!is.na(weekly)   , weekly*52,
                             fifelse(!is.na(biweekly) , biweekly*26,
                                     fifelse(!is.na(semimonth), semimonth*24,
                                             fifelse(!is.na(monthly), monthly*12, NA_real_))))))]

# dt[, hourly_wage := fifelse(!is.na(annsal), annsal/2080,
#              fifelse(!is.na(weekly)   , weekly*52/2080,
#                      fifelse(!is.na(biweekly) , biweekly*26/2080,
#                              fifelse(!is.na(semimonth), semimonth*24/2080,
#                                      fifelse(!is.na(monthly), monthly*12/2080, NA_real_)))))]

# flags

dt[wfh_days > 5, wfh_days := 5] # some idiots report more than 5 days
dt[, employer_size := fifelse(!is.na(empall), empall, empsize)]
dt[, work_diff_state := as.numeric(work_state == 8)]

# define spells
dt[, spell_id := rleid(fifelse(is.na(jobid), -1L, jobid)), by = .(ssuid,pnum)]

spell <- dt[, .(
  start_year      = first(job_start_year),          # unchanged
  end_year      = last(job_start_year),           # last year of job
  start_month     = first(month),
  state_now       = first(state_fips),
  industry_now    = first(industry),
  occupation_now  = first(occupation),
  annual_salary   = first(annual_salary),
  annual_salary_last   = last(annual_salary),
  hourly_wage     = first(hourly_wage),
  class_worker_now= first(class_worker),
  emp_status_now  = first(emp_status),
  age_start = first(age), sex = first(sex), educ_level = first(educ_level),
  union_now       = first(union),
  weight          = first(weight),
  employer_size   = first(employer_size),
  work_diff_state = first(work_diff_state),
  wfh_0 = dplyr::coalesce(mean(wfh_days==0,na.rm=TRUE), 0),
  wfh_1 = dplyr::coalesce(mean(wfh_days==1,na.rm=TRUE), 0),
  wfh_2 = dplyr::coalesce(mean(wfh_days==2,na.rm=TRUE), 0),
  wfh_3 = dplyr::coalesce(mean(wfh_days==3,na.rm=TRUE), 0),
  wfh_4 = dplyr::coalesce(mean(wfh_days==4,na.rm=TRUE), 0),
  wfh_5 = dplyr::coalesce(mean(wfh_days==5,na.rm=TRUE), 0)
), by = .(ssuid,pnum,spell_id)]

spell[, wfh_any_now  := wfh_1 > 0 | wfh_2 > 0 | wfh_3 > 0 | wfh_4 > 0 | wfh_5 > 0]
spell[, wfh_full_now := wfh_5 == 1]

# leads
setorder(spell, ssuid,pnum,start_year,start_month)
spell[, `:=`(
  industry_next     = shift(industry_now , type="lead"),
  occupation_next   = shift(occupation_now, type="lead"),
  class_worker_next = shift(class_worker_now, type="lead"),
  emp_status_next   = shift(emp_status_now , type="lead"),
  state_next        = shift(state_now      , type="lead"),
  union_next        = shift(union_now      , type="lead"),
  wfh_any_next      = shift(wfh_any_now    , type="lead"),
  wfh_full_next     = shift(wfh_full_now   , type="lead"), 
  moved_state_next = shift(work_diff_state, type="lead"),
  work_diff_state_next = shift(work_diff_state, type="lead")
), by = .(ssuid,pnum)]

spell[, moved_state := !is.na(state_next) & state_now != state_next]

# final panel
panel_final <- spell[, .(
  ssuid, pnum, spell_id,
  start_year, end_year, start_month,
  annual_salary, annual_salary_last, hourly_wage,
  industry_now,industry_next,
  occupation_now,occupation_next,
  class_worker_now,class_worker_next,
  emp_status_now,emp_status_next,
  union_now,union_next,
  employer_size,
  wfh_any_now,wfh_any_next,
  wfh_full_now,wfh_full_next,
  state_now,state_next,
  age_start,sex,educ_level, weight, 
  work_diff_state,moved_state, work_diff_state_next, moved_state_next,
  wfh_0, wfh_1, wfh_2, wfh_3, wfh_4, wfh_5
)]

dir.create(dataOut, showWarnings = FALSE, recursive = TRUE)
fwrite(panel_final, file.path(dataOut, "sipp_job_panel.csv"))
