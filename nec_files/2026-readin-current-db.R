source("libs.R")
library(janitor)
library(lubridate)
library(readxl)
necropsy_log_master0 <- read_excel("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/necfiles0/necropsy_log-master0.xlsx")


# Function to detect and convert date columns
convert_dates <- function(x) {
  
  # Already a date
  if (inherits(x, c("Date", "POSIXct", "POSIXt"))) {
    return(as.Date(x))
  }
  
  # Numeric Excel dates
  if (is.numeric(x)) {
    
    # Typical Excel date range (~1900-2100)
    if (all(is.na(x) | (x > 20000 & x < 80000))) {
      return(as.Date(x, origin = "1899-12-30"))
    }
    
    return(x)
  }
  
  # Character dates
  if (is.character(x)) {
    
    parsed <- suppressWarnings(parse_date_time(
      x,
      orders = c(
        "ymd", "mdy", "dmy",
        "Ymd", "m/d/Y", "d/m/Y",
        "m-d-Y", "d-m-Y",
        "Y-m-d"
      )
    ))
    
    # Convert if most non-missing values parsed successfully
    non_missing <- sum(!is.na(x))
    parsed_non_missing <- sum(!is.na(parsed))
    
    if (non_missing > 0 &&
        parsed_non_missing / non_missing > 0.8) {
      return(as.Date(parsed))
    }
  }
  
  x
}
convert_mixed_dates <- function(x, threshold = 0.8) {
  
  x_trim <- trimws(as.character(x))
  
  result <- rep(as.Date(NA), length(x))
  
  # Excel serial dates
  is_num <- grepl("^\\d+(\\.\\d+)?$", x_trim)
  
  result[is_num] <- as.Date(
    as.numeric(x_trim[is_num]),
    origin = "1899-12-30"
  )
  
  # Character dates
  result[!is_num] <- suppressWarnings(
    as.Date(parse_date_time(
      x_trim[!is_num],
      orders = c(
        "ymd", "mdy", "dmy",
        "m/d/Y", "d/m/Y",
        "m-d-Y", "d-m-Y",
        "b d Y", "d b Y"
      )
    ))
  )
  
  # Only convert if enough values parsed
  non_missing <- sum(x_trim != "" & !is.na(x_trim))
  converted <- sum(!is.na(result))
  
  if (non_missing > 0 && converted / non_missing >= threshold) {
    return(result)
  }
  
  x
}
convert_excel_date_column <- function(x) {
  
  x_chr <- trimws(as.character(x))
  
  out <- suppressWarnings(as.Date(
    as.numeric(x_chr),
    origin = "1899-12-30"
  ))
  
  # Fill in failed conversions with parsed text dates
  missing <- is.na(out)
  
  out[missing] <- suppressWarnings(as.Date(
    parse_date_time(
      x_chr[missing],
      orders = c(
        "ymd", "mdy", "dmy",
        "m/d/Y", "d/m/Y",
        "m-d-Y", "d-m-Y"
      )
    )
  ))
  
  out
}

date_cols <- names(df)[
  grepl("date|dob|birth|collected|sample|visit",
        names(df),
        ignore.case = TRUE) &
    !grepl("age",
           names(df),
           ignore.case = TRUE)
]



date_cols <- grep(
  "date|dob|birth|collected|sample|visit",
  names(df),
  ignore.case = TRUE,
  value = TRUE
)





date_cols <- names(nec1)[
  grepl("date|dob|birth|collected|sample|visit|bd", names(nec1), ignore.case = TRUE) &
    !grepl("age", names(nec1), ignore.case = TRUE)
]


print(date_cols)


nec1<-necropsy_log_master0%>%
  row_to_names(1)%>%
  clean_names()
  


nec2<-nec1%>%
  mutate(across(all_of(date_cols), convert_mixed_dates))

         



saveRDS(nec2,"nec_files/2026-06-necdb.rds")
