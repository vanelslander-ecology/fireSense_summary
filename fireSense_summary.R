defineModule(sim, list(
  name = "fireSense_summary",
  description = paste("Summarizes the results of multiple fireSense simulations,",
                      "across multiple study areas, climate scenarios, and replicates."),
  keywords = "fireSense",
  authors = c(
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = c("aut", "cre")),
    person("Tati", "Micheletti", email = "tati.micheletti@gmail.com", role = "aut"),
    person("Ian MS", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = "aut")
  ),
  childModules = character(0),
  version = list(fireSense_summary = "1.0.1"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "fireSense_summary.Rmd"), ## same file
  loadOrder = list(after = c("fireSense")),
  reqdPkgs = list(
    "assertthat", "cowplot", "data.table", "fs", "ggplot2", "googledrive",
    "purrr", "qs2", "RColorBrewer", "terra", "tidyterra",
    "raster", "rasterVis", ## TODO: remove these once fireSenseUtils::plotCumulativeBurns switched to ggplot2/tidyterra
    "PredictiveEcology/fireSenseUtils@development (>= 0.1.2.9000)",
    "PredictiveEcology/SpaDES.core@development (>= 3.0.3.9003)",
    "PredictiveEcology/SpaDES.tools@development (>= 2.1.1.9000)"
  ),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("climateScenario", "character", NA, NA, NA,
                    desc = paste("name of CIMP6 climate scenarios including SSP,",
                                 "formatted as in `ClimateNA`, using underscores as separator.",
                                 "E.g., 'CanESM5_SSP370'.")),
    defineParameter("mode", "character", "single", NA, NA,
                    paste("use 'single' to run part of a simulation;",
                          "use 'multi' to run as part of postprocessing multiple runs.")),
    defineParameter("simOutputPath", "character", outputPath(sim), NA, NA,
                    desc = "Directory specifying the location of the simulation outputs."),
    defineParameter("studyAreaName", "character", NA, NA, NA,
                    desc = "name of study areas simulated."),
    defineParameter("reps", "integer", 1L:10L, 1, NA,
                    desc = paste("number of replicates/runs per study area and climate scenario.",
                                 "NOTE: `mclapply` is used internally, so you should set",
                                 "`options(mc.cores = nReps)` to run in parallel.")),
    defineParameter("years", "integer", c(2011L, 2100L), NA, NA,
                    desc = paste("Which two simulation years should be compared?",
                                 "Typically start and end years."))
  ),
  inputObjects = bindrows(
    expectsInput("burnMap", "SpatRaster",
                 desc = paste("Cumulative burn map.", "Required in single mode."),
                 sourceURL = NA),
    expectsInput("burnSummary", "data.table",
                 paste("Fire summary table from `fireSense`.", "Required in single mode."),
                 sourceURL = NA),
    expectsInput("rasterToMatch", "SpatRaster",
                 paste("template raster used by the simulations for summary reporting"),
                 sourceURL = NA)
  ),
  outputObjects = bindrows(
    createsOutput(objectName = NA, objectClass = NA, desc = NA)
  )
))

## event types
#   - type `init` is required for initialization

doEvent.fireSense_summary = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      if (P(sim)$mode == "single") {
        sim <- scheduleEvent(sim, end(sim), "fireSense_summary", "save_single", .last())
      } else if (P(sim)$mode == "multi") {
        sim <- InitMulti(sim)

        padYear <- paddedFloatToChar(time(sim), padL = ceiling(log10(end(sim) + 1)))

        f_burnSummary_plot <- fireSenseUtils::plotBurnSummary(
          climateScenario = P(sim)$climateScenario,
          studyAreaName = P(sim)$studyAreaName,
          outputDir = P(sim)$simOutputPath,
          Nreps = max(P(sim)$reps),
          years = P(sim)$years,
          pixelSize = unique(terra::res(sim$rasterToMatch))
        )
        sim <- registerOutputs(f_burnSummary_plot, sim)

        f_cumulBurn_plot <- fireSenseUtils::plotCumulativeBurns(
          climateScenario = P(sim)$climateScenario,
          studyAreaName = P(sim)$studyAreaName,
          outputDir = P(sim)$simOutputPath,
          Nreps = max(P(sim)$reps),
          years = P(sim)$years,
          rasterToMatch = sim$rasterToMatch
        )
        sim <- registerOutputs(f_cumulBurn_plot, sim)

        f_historic_plot <- fireSenseUtils::plotHistoricFires(
          climateScenario = P(sim)$climateScenario,
          studyAreaName = P(sim)$studyAreaName,
          outputDir = P(sim)$simOutputPath,
          pixelSize = unique(terra::res(sim$rasterToMatch)),
          firePolys = mod$firePolys,
          ignitionPoints = mod$ignitionFirePoints
        )
        sim <- registerOutputs(f_historic_plot, sim)
      }
    },
    save_single = {
      padYear <- paddedFloatToChar(time(sim), padL = ceiling(log10(end(sim) + 1)))

      f_burnMap <- file.path(outputPath(sim), paste0("burnMap_year", padYear, ".tif"))
      terra::writeRaster(sim$burnMap, f_burnMap, datatype = "INT2U", overwrite = TRUE)
      sim <- registerOutputs(f_burnMap, sim)

      f_burnSummary <- file.path(outputPath(sim), "fireSense_burnSummary.csv")
      data.table::fwrite(sim$burnSummary, file = f_burnSummary)
      sim <- registerOutputs(f_burnSummary, sim)
    },
    noEventWarning(sim)
  )
  return(invisible(sim))
}

InitMulti <- function(sim) {
  # # ! ----- EDIT BELOW ----- ! #

  ## check for necessary output files -----------------------------------------------
  allReps <- sprintf("rep%02d", P(sim)$reps)
  padL <- ceiling(log10(P(sim)$years[2] + 1))
  padYearStart <- paddedFloatToChar(P(sim)$years[1], padL = padL)
  padYearEnd <- paddedFloatToChar(P(sim)$years[2], padL = padL)

  checkPath(file.path(P(sim)$simOutputPath, "figures", currentModule(sim)), create = TRUE)

  bmbs <- fs::dir_ls(
    P(sim)$simOutputPath,
    regexp = "burnMap|burnSummary",
    recurse = 1,
    type = "file"
  ) |>
    grep(paste0("(", paste0(P(sim)$reps, collapse = "|"), ")"), x = _, value = TRUE) |>
    grep(paste0("_year(", paste0(P(sim)$years, collapse = "|"), ")"), x = _, value = TRUE)

  filesUserHas <- c(bmbs)

  dirsExpected <- file.path(P(sim)$simOutputPath, allReps)
  filesExpected <- as.character(sapply(dirsExpected, function(d) {
    c(
      file.path(d, sprintf("burnMap_year%04d.tif", P(sim)$years[2])),
      file.path(d, "fireSense_burnSummary.csv")
    )
  }))

  filesNeeded <- data.frame(file = filesExpected, exists = filesExpected %in% filesUserHas)

  if (!all(filesNeeded$exists)) {
    missing <- filesNeeded[filesNeeded$exists == FALSE, ]$file
    stop(
      sum(!filesNeeded$exists),
      " simulation files appear to be missing:\n",
      paste(missing, collapse = "\n")
    )
  }

  ## get historical fire points and polys -------------------------------------------

  ## TODO: use an updated/working prepInputs version (fireSenseUtils::getFirePolygons?)
  mod$firePolys <- {
    dst <- inputPath(sim)
    nfdb_url <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"
    nfdb_zip <- file.path(dst, basename(nfdb_url))

    if (!file.exists(nfdb_zip)) {
      download.file(nfdb_url, destfile = nfdb_zip)
    }

    all_nfdb_files <- fs::dir_ls(dst, regexp = "NFDB_poly_(1972to2020|2021to2024).*")

    if (length(all_nfdb_files) != 16) {
      archive::archive_extract(nfdb_zip, dst)
    }

    nfdb_shp <- fs::dir_ls(dst, regexp = "NFDB_poly_(1972to2020|2021to2024).*[.]shp$")

    purrr::map(.x = nfdb_shp, .f = function(x) {
      p <- terra::vect(x)

      ## NOTE: terra::makeValid takes so long;
      ## just drop the tiny number of invalid geometries
      p[terra::is.valid(p), ]
    }) |>
      tidyterra::bind_spat_rows() |>
      tidyterra::mutate(
        YEAR = as.integer(YEAR),
        MONTH = as.integer(MONTH),
        DAY = as.integer(DAY)
      ) |>
      terra::project(sim$rasterToMatch)
  }

  ## TODO: use an updated/working prepInputs version (fireSenseUtils::getFirePoints_NFDB_V2?)
  mod$ignitionFirePoints <- {
    dst <- inputPath(sim)

    nfdb_url <- "http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip"
    nfdb_zip <- file.path(dst, basename(nfdb_url))

    if (!file.exists(nfdb_zip)) {
      download.file(nfdb_url, destfile = nfdb_zip)
    }

    all_nfdb_files <- fs::dir_ls(dst, regexp = "NFDB_point_.*")

    if (length(all_nfdb_files) != 10) {
      archive::archive_extract(nfdb_zip, dst)
    }

    nfdb_shp <- fs::dir_ls(dst, regexp = "NFDB_point_.*[.]shp$")

    ## NOTE: using terra here because it's much faster than sf
    p <- terra::vect(nfdb_shp)

    ## NOTE: terra::makeValid takes so long;
    ## just drop the tiny number of invalid geometries
    p[terra::is.valid(p), ] |>
      tidyterra::mutate(
        YEAR = as.integer(YEAR),
        MONTH = as.integer(MONTH),
        DAY = as.integer(DAY)
      ) |>
      terra::project(sim$rasterToMatch)
  }

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #

  ## nothing here

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
