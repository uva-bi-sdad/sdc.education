library(tidycensus)
library(data.table)
library(dplyr)
library(tidyr)

census_api_key(Sys.getenv('census_api_key'))
source("~/git/sdc.education_dev/utils/distribution/aggregate.R")

#function for calculating ays
get_ays <- function(acs_data, tract_geoid) {
  tract_data <- acs_data[acs_data$GEOID==tract_geoid, ]

  # total population
  pop <- tract_data[tract_data$variable %in% c("B15003_001"), c("estimate")][[1]]

  # value assigned to each grade level (no schooling, nursery school, kindergarten each assigned 0)
  values <- c("005" = 1, "006" = 2, "007" = 3, "008" = 4, "009" = 5, "010" = 6,
              "011" = 7, "012" = 8, "013" = 9, "014" = 10, "015" = 11, "016" = 12,
              "017" = 12, "018" = 12, "019" = 12.5, "020" = 13, "021" = 14,
              "022" = 16, "023" = 18, "024" = 19, "025" = 20)

  ays <- 0

  # calculate years of schooling for each variable to get average years of schooling
  for (var in names(values)) {
    temp <- tract_data[tract_data$variable == paste0("B15003_", var),] %>% filter(!is.na(estimate))

    # add current grade level to average years of schooling
    if (!is.null(temp)){
      ays <- ays + ((sum(temp$estimate)/pop) * values[[var]])
    }
  }

  return (ays)

}



#function for getting different years acs data
get_acsdata <- function(geography, table, state, survey, start_year, end_year) {
  acsdata <- NULL
  for (year in start_year:end_year) {
    acsdata_temp <- get_acs(geography = geography,
                            table = table,
                            year = year,
                            state = state,
                            survey = survey) %>% mutate(year=year)

    acsdata <- rbind(acsdata, acsdata_temp)
  }

  return(acsdata)
}

#get acs data for seven years tracts
acs_yos <- get_acsdata(geography = "tract",
                       table = "B15003",
                       state = "VA",
                       survey = "acs5",
                       start_year = 2015,
                       end_year = 2021)

all_tr <- NULL
#different years calculation
for (yr in c(2021, 2020, 2019, 2018, 2017, 2016, 2015)) {
  # Get the data for the current year
  acsdata_year <- acs_yos %>% filter(year==yr)

  tract_ays <- NULL

  # Calculate AYS for all tracts
  unq_tracts <- unique(acsdata_year$GEOID)
  for (t in unq_tracts) {
    ays_t <- get_ays(acsdata_year, t)

    tract_ays <- rbind(tract_ays, data.frame(GEOID=c(t), value=c(ays_t)))
  }

  acsdata_year <- acsdata_year %>% group_by(GEOID) %>% summarise(year=first(year), moe='')
  acsdata_year <- merge(acsdata_year, tract_ays, by='GEOID')

  all_tr <- rbind(all_tr, acsdata_year)
}

all_tr <- all_tr %>% filter(!is.nan(value)) %>% rename(geoid=GEOID)

#get acs data for seven years counties
acs_yos <- get_acsdata(geography = "county",
                       table = "B15003",
                       state = "VA",
                       survey = "acs5",
                       start_year = 2015,
                       end_year = 2021)

all_ct <- NULL
#different years calculation
for (yr in c(2021, 2020, 2019, 2018, 2017, 2016, 2015)) {
  # Get the data for the current year
  acsdata_year <- acs_yos %>% filter(year==yr)

  tract_ays <- NULL

  # Calculate AYS for all tracts
  unq_tracts <- unique(acsdata_year$GEOID)
  for (t in unq_tracts) {
    ays_t <- get_ays(acsdata_year, t)

    tract_ays <- rbind(tract_ays, data.frame(GEOID=c(t), value=c(ays_t)))
  }

  acsdata_year <- acsdata_year %>% group_by(GEOID) %>% summarise(year=first(year), moe='')
  acsdata_year <- merge(acsdata_year, tract_ays, by='GEOID')

  all_ct <- rbind(all_ct, acsdata_year)
}

all_ct <- all_ct %>% filter(!is.nan(value)) %>% rename(geoid=GEOID)

# combine tract, county, health district data
tr_ct <- rbind(all_tr, all_ct)

readr::write_csv(tr_ct, xzfile('Years of Schooling/data/working/va_trct_acs5_2015_2021_years_of_schooling.csv', compression=9))

# testing <- readxl::read_excel('~/git/HOI V3_14 Variables_Raw Scores (1).xlsx') %>%
#   select('CT2', 'Education')
#
# merged <- merge(testing, all_tr %>% filter(year==2020) %>% mutate(value1=round(value, 1)), by.x='CT2', by.y='geoid')
# ggplot(merged, aes(x=Education, y=value1)) + geom_point()  + geom_abline(slope=1)






