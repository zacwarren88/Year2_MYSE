f.aggregatID_Y2 <- function(folder.IDfiles, 
                            pattern.IDfiles = "id\\.csv", 
                            multifolder = T, 
                            pattern.summaryfiles = "idsummary\\.csv") {
  
  require(parallel)
  require(readr)
  require(dplyr)
  require(stringr)
  require(lubridate)
  require(zoo)
  require(tidyr)
  
  files.IDSUMMARY <-
    list.files( #list files
      folder.IDfiles, #folder path designated above
      pattern = pattern.IDfiles, #pattern, default is idsummary
      full.names = T, 
      recursive = multifolder) #look down into folder structure
  
  spec.list <- c("EPFU","LANO","LACI","NYHU","LABO","MYSE","MYLU","PESU","MYCI","MYVO","MYTH")
  
  #remove potential .bak1 files
  files.IDSUMMARY <- files.IDSUMMARY[lapply(files.IDSUMMARY, function(x) length(grep("bak1", x, value = FALSE))) == 0]
  
  ID <- bind_rows(lapply(files.IDSUMMARY, read_csv, guess_max = 10000))
  
  names(ID) <- toupper(names(ID))
  names(ID) <- str_replace(names(ID), "-", "_")
  names(ID) <- str_replace(names(ID), " ", "_")
  names(ID) <- str_replace(names(ID), "\\*", "")
  
  ID %>% 
    separate(col = FOLDER, into = c("GRID", "SITE", "NIGHT"), sep = "\\\\") %>% 
    select(-OFFSET, 
           -OUT_FILE, 
           -(FC:QUAL), 
           -MANUAL_ID, 
           -TIME_12, 
           -HOUR, 
           -HOUR_12, 
           -ALTERNATES,
           -CHANNEL,
           -DURATION,
           -MARGIN) %>% 
    rename(NIGHT_CDT=DATE_12,
           Kpro_ID = AUTO_ID) %>% 
    #add quality value
    ##high = at least 10 pulses and 90% matching, need 2 - value =0.5
    ##medium = at least 5 pulses and 75% matching, need 3 - value = 1/3
    ##low = at least 5 pulses and 50% matching, need 4 value = 0.25
    mutate(ID_pts = case_when( Kpro_ID == "NoID" ~ 1,
                               Kpro_ID == "Noise" ~ 1,
                               MATCH_RATIO <0.5 ~ 0,                #less than 50% should be NOID
                               PULSES>=10 & MATCH_RATIO>=0.9 ~ 0.5,      #Assing High quality
                               PULSES>=5  & MATCH_RATIO>=0.75 ~ 1/3,     #Assign medium quality
                               PULSES>=5  & MATCH_RATIO>=0.5 ~ 0.25,
                               TRUE ~ -Inf )) %>%  
    group_by(SITE, NIGHT_CDT, Kpro_ID) %>% 
    mutate(DATE_SPP_score = sum(ID_pts)) %>% 
    ungroup() %>% 
    mutate(ID_CORRECTED = case_when(DATE_SPP_score>=1.0 ~ Kpro_ID,
                                    DATE_SPP_score<1.0~ "NoID")) %>% 
    ungroup() %>% 
    select(GRID, SITE, DATE, NIGHT_CDT, TIME, 
           PULSES, MATCHING, MATCH_RATIO, N, 
           ID_pts, DATE_SPP_score, Kpro_ID, ID_CORRECTED) %>% 
    mutate(Kpro_ID = case_when(Kpro_ID == "EPTFUS" ~ "EPFU",
                               Kpro_ID == "LASNOC" ~ "LANO",
                               Kpro_ID == "LASCIN" ~ "LACI",
                               Kpro_ID == "NYCHUM" ~ "NYHU",
                               Kpro_ID == "LASBOR" ~ "LABO",
                               Kpro_ID == "MYOSEP" ~ "MYSE",
                               Kpro_ID == "MYOLUC" ~ "MYLU",
                               Kpro_ID == "PERSUB" ~ "PESU",
                               Kpro_ID == "MYOCIL" ~ "MYCI",
                               Kpro_ID == "MYOVOL" ~ "MYVO",
                               Kpro_ID == "MYOTHY" ~ "MYTH",
                               TRUE ~ Kpro_ID)) %>% 
    mutate(ID_CORRECTED = case_when(ID_CORRECTED == "EPTFUS" ~ "EPFU",
                                    ID_CORRECTED == "LASNOC" ~ "LANO",
                                    ID_CORRECTED == "LASCIN" ~ "LACI",
                                    ID_CORRECTED == "NYCHUM" ~ "NYHU",
                                    ID_CORRECTED == "LASBOR" ~ "LABO",
                                    ID_CORRECTED == "MYOSEP" ~ "MYSE",
                                    ID_CORRECTED == "MYOLUC" ~ "MYLU",
                                    ID_CORRECTED == "PERSUB" ~ "PESU",
                                    ID_CORRECTED == "MYOCIL" ~ "MYCI",
                                    ID_CORRECTED == "MYOVOL" ~ "MYVO",
                                    ID_CORRECTED == "MYOTHY" ~ "MYTH",
                                    TRUE ~ ID_CORRECTED)) -> ID
  
  files.MLE <-
    list.files( #list files
      folder.IDfiles, #folder path designated above
      pattern = pattern.summaryfiles, #pattern, default is idsummary
      full.names = T, 
      recursive = multifolder #look down into folder structure
    )
  
  #remove potential .bak1 files
  files.MLE <- files.MLE[lapply(files.MLE, function(x) length(grep("bak1", x, value = FALSE))) == 0]
  
  #read file and bind
  t <- suppressWarnings(bind_rows(lapply(files.MLE, read_csv, skip = 1, col_types = cols(.default = "c"))))
  

  names(t) <- toupper(names(t))
  names(t) <- str_replace(names(t), "-", "_")
  names(t) <- str_replace(names(t), " ", "_")
  names(t) <- str_replace(names(t), "\\*", "")
  
  
  t <- t %>% 
    select(matches("X[0-9]*"), matches("^[A-Z]{6}\\_1$")) %>% #select cols with MLE and Xs
    select(which(colMeans(is.na(.))<1))
  
  
  
  
  
  
  columnstofill <- grep("X[0-9]*", names(t)) #find columns with X in header
  
  t[, columnstofill]<- sapply(t[, columnstofill],FUN = function(x) na.locf(x)) #fill NA with most recent non-NA
  t <- t[!apply(t, 1, function(r) any(r %in% c("*"))),] #remove rows with * in either site or deployment
  
  
  ##change last (in order) unnamed column to "NIGHT"
  names(t)[max( matches("X[0-9]*", vars = names(t)) )] <- "NIGHT"
  
  ###CHANGE DEPENDING ON PROJECT###
  ##change first unnamed column to "SITE"
  names(t)[min( matches("X[0-9]*", vars = names(t)) )] <- "GRID"
  
  ##change remaining unnamed column to "deployment"
  names(t)[matches("X[0-9]*", vars = names(t))] <- "SITE"
  
  names(t)[matches("EPTFUS_1", vars = names(t))] <- "EPFU"
  names(t)[matches("LASNOC_1", vars = names(t))] <- "LANO"
  names(t)[matches("LASCIN_1", vars = names(t))] <- "LACI"
  names(t)[matches("NYCHUM_1", vars = names(t))] <- "NYHU"
  names(t)[matches("LASBOR_1", vars = names(t))] <- "LABO"
  names(t)[matches("MYOSEP_1", vars = names(t))] <- "MYSE"
  names(t)[matches("MYOLUC_1", vars = names(t))] <- "MYLU"
  names(t)[matches("PERSUB_1", vars = names(t))] <- "PESU"
  names(t)[matches("MYOCIL_1", vars = names(t))] <- "MYCI"
  names(t)[matches("MYOVOL_1", vars = names(t))] <- "MYVO"
  names(t)[matches("MYOTHY_1", vars = names(t))] <- "MYTH"
  
  #convert from wide to long
  t <- t %>% 
    gather("SPECIES", "MLE_NIGHT", suppressWarnings(one_of(spec.list))) %>%   #convert from wide to long
    mutate(NIGHT = ymd(NIGHT))
  
  output <- ID %>% 
    left_join(t, by = c("GRID", "SITE", "NIGHT_CDT" = "NIGHT", "Kpro_ID" = "SPECIES")) 
  
  
  
  return(output)
}
