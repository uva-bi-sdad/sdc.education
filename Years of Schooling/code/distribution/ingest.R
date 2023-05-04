base_dir <- "Years of Schooling/data"
dir.create(paste0(base_dir, "/original/reference_shapes"), FALSE, TRUE)

# get health district associations
va_id_map <- jsonlite::read_json("https://uva-bi-sdad.github.io/community/dist/shapes/VA/entity_info.json")
county_districts <- c(
  unlist(lapply(va_id_map$county, "[[", "district")),
  "11001" = "11_hd_01", "24017" = "24_hd_01", "24021" = "24_hd_01",
  "24031" = "24_hd_01", "24033" = "24_hd_01"
)
missing_districts <- county_districts[!county_districts %in% names(va_id_map$district)]
names(missing_districts) <- missing_districts

states <- c("DC", "MD", "VA")
years <- 2013:2021

# download and load maps
entity_info <- c(
  lapply(va_id_map$district, function(e) list(region_name = e$name)),
  lapply(missing_districts, function(e) list(region_name = e))
)
for (location in tolower(states)) {
  for (year in c(2020, 2010)) {
    for (level in list(c("Block%20Group", "census_block_groups"), c("Tract", "census_tracts"), c("County", "counties"))) {
      name <- paste0(location, "_geo_census_cb_", year, "_", level[[2]])
      file <- paste0(
        base_dir, "/original/reference_shapes/", location, "_",
        sub("census_", "", level[[2]], fixed = TRUE), "_", year, ".geojson"
      )
      if (!file.exists(file)) {
        tryCatch(download.file(paste0(
          "https://raw.githubusercontent.com/uva-bi-sdad/sdc.geographies/main/",
          toupper(location), "/Census%20Geographies/", level[[1]], "/", year,
          "/data/distribution/", name, ".geojson"
        ), file), error = function(e) NULL)
      }
      if (file.exists(file)) {
        d <- jsonlite::read_json(file)
        for (e in d$features) {
          if (e$properties$geoid %in% names(entity_info)) e$properties$year <- "both"
          entity_info[[e$properties$geoid]] <- e$properties
        }
      }
    }
  }
}
entity_names <- unlist(lapply(entity_info, "[[", "region_name"))
entity_names <- entity_names[!grepl(", NA", entity_names, fixed = TRUE)]
entity_year <- vapply(entity_info, function(e) if (length(e$year)) e$year else "both", "")

vars <- c(
  "1.m_0" = "B15002_003",
  "2.m_2.5" = "B15002_004",
  "3.m_5.5" = "B15002_005",
  "4.m_7.5" = "B15002_006",
  "5.m_9" = "B15002_007",
  "6.m_10" = "B15002_008",
  "7.m_11" = "B15002_009",
  "8.m_11" = "B15002_010",
  "9.m_12" = "B15002_011",
  "10.m_13" = "B15002_012",
  "11.m_14" = "B15002_013",
  "12.m_14" = "B15002_014",
  "13.m_16" = "B15002_015",
  "14.m_18" = "B15002_016",
  "15.m_19" = "B15002_017",
  "16.m_20" = "B15002_018",
  "1.f_0" = "B15002_020",
  "2.f_2.5" = "B15002_021",
  "3.f_5.5" = "B15002_022",
  "4.f_7.5" = "B15002_023",
  "5.f_9" = "B15002_024",
  "6.f_10" = "B15002_025",
  "7.f_11" = "B15002_026",
  "8.f_11" = "B15002_027",
  "9.f_12" = "B15002_028",
  "10.f_13" = "B15002_029",
  "11.f_14" = "B15002_030",
  "12.f_14" = "B15002_031",
  "13.f_16" = "B15002_032",
  "14.f_18" = "B15002_033",
  "15.f_19" = "B15002_034",
  "16.f_20" = "B15002_035"
)
base_vars <- paste0(
  rep(c("average_years_schooling", "EduGini"), each = 3), c("", "_female", "_male")
)
error_vars <- paste0(base_vars, "_error")
sv <- list(
  f = grep("f", names(vars), value = TRUE),
  m = grep("m", names(vars), value = TRUE)
)

# download ACS data, calculate variables, and aggregate
## methods from https://files.eric.ed.gov/fulltext/ED536152.pdf

y <- as.numeric(sub("^.*_", "", sv$f))
Gini <- function(m, p) {
  n <- ncol(p)
  y <- as.numeric(sub("^.*_", "", colnames(p)))
  e <- m^-1 * rowSums(p[, -1] * abs(y[-1] - y[-n]) * p[, -n])
  e[!is.finite(e)] <- 0
  e
}

data <- do.call(rbind, lapply(states, function(state) {
  do.call(rbind, lapply(years, function(year) {
    retrieve <- function(layer) {
      d <- tidycensus::get_acs(
        layer,
        variables = vars,
        year = year,
        state = state,
        output = "wide"
      )

      # overall
      d_all <- d[, paste0(sv$f, "E")] + d[, paste0(sv$m, "E")]
      total <- rowSums(d_all)
      total[total == 0] <- 1
      d$average_years_schooling <- (as.matrix(d_all) %*% y) / total
      d$EduGini <- Gini(d$average_years_schooling, d_all / total)
      d_all <- d_all + sqrt(
        (d[, paste0(sv$f, "M")] / 1.645)^2 + (d[, paste0(sv$m, "M")] / 1.645)^2
      ) * 1.645
      d$average_years_schooling_error <- (as.matrix(d_all) %*% y) / total
      d$EduGini_error <- abs(Gini(d$average_years_schooling_error, d_all / total) - d$EduGini)
      d$average_years_schooling_error <- abs(
        d$average_years_schooling_error - d$average_years_schooling
      )

      # female population
      d_female <- d[, paste0(sv$f, "E")]
      total <- rowSums(d_female)
      total[total == 0] <- 1
      d$average_years_schooling_female <- (as.matrix(d_female) %*% y) / total
      d$EduGini_female <- Gini(d$average_years_schooling_female, d_female / total)
      d_female <- d_female + d[, paste0(sv$f, "M")]
      d$average_years_schooling_female_error <- (as.matrix(d_female) %*% y) / total
      d$EduGini_female_error <- abs(Gini(
        d$average_years_schooling_female_error, d_female / total
      ) - d$EduGini_female)
      d$average_years_schooling_female_error <- abs(
        d$average_years_schooling_female_error - d$average_years_schooling_female
      )

      # male population
      d_male <- d[, paste0(sv$m, "E")]
      total <- rowSums(d_male)
      total[total == 0] <- 1
      d$average_years_schooling_male <- (as.matrix(d_male) %*% y) / total
      d$EduGini_male <- Gini(d$average_years_schooling_male, d_male / total)
      d_male <- d_male + d[, paste0(sv$m, "M")]
      d$average_years_schooling_male_error <- (as.matrix(d_male) %*% y) / total
      d$EduGini_male_error <- abs(
        Gini(d$average_years_schooling_male_error, d_male / total) - d$EduGini_male
      )
      d$average_years_schooling_male_error <- abs(
        d$average_years_schooling_male_error - d$average_years_schooling_male
      )

      d <- d[d$GEOID %in% names(entity_names), c("GEOID", base_vars, error_vars)]
      d <- d[entity_names[d$GEOID] != (if (year > 2019) 2010 else 2020), ]
      d$region_type <- layer
      d$region_name <- entity_names[d$GEOID]
      d$year <- year
      n <- nrow(d)
      list(wide = d, tall = list(cbind(
        d[, c("GEOID", "region_type", "region_name", "year")],
        measure = rep(base_vars, each = n),
        value = unlist(d[, base_vars], use.names = FALSE),
        moe = unlist(d[, error_vars], use.names = FALSE)
      )))
    }
    counties <- retrieve("county")
    do.call(rbind, c(
      retrieve("block group")$tall,
      retrieve("tract")$tall,
      counties$tall,
      lapply(
        list(county_districts[substring(counties$wide$GEOID, 1, 5)]),
        function(set) {
          counties$wide$GEOID <- set
          do.call(rbind, lapply(split(counties$wide, counties$wide$GEOID), function(e) {
            id <- e$GEOID[[1]]
            data.frame(
              GEOID = id,
              region_type = "health district",
              region_name = entity_names[[id]],
              measure = base_vars,
              year = year,
              value = colMeans(e[, base_vars]),
              moe = abs(colMeans(e[, base_vars] + e[, error_vars]) - colMeans(e[, base_vars]))
            )
          }))
        }
      )
    ))
  }))
}))
colnames(data) <- tolower(colnames(data))

vroom::vroom_write(data, paste0(base_dir, "/distribution/acs.csv.xz"), ",")
