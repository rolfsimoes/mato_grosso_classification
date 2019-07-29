#' - Script citation:
#' Rolf Simoes, Adeline Maciel, Pedro Andrade, Lorena Santos, 
#' Gilberto Camara, and Michelle Picoli. Source code for: Land use and cover 
#' change maps for Mato Grosso State in Brazil. (July 2019). Available 
#' at: https://doi.org/10.5281/zenodo.3354380
#'
#' - Description:
#' This R script was used to do the classification and post processing of the
#' land use and cover change maps for Mato Grosso State, Brazil.
#'
#' - Reference paper:
#' Rolf Simoes, Michelle C. A. Picoli, Gilberto Camara, 
#' Adeline Maciel, Lorena Santos, Pedro R. Andrade,
#' Alber Sánchez, Karine Ferreira, and Alexandre Carvalho.
#' Land use and cover change maps for Mato Grosso State in Brazil: 2001-2017.
#' Originally submitted to Nature Scientific Data in Jul 2019.
#'
#' - Script Usage:
#' Before run this script, open it on any editor or R IDE (e.g. RStudio) to
#' inform the input parameters, in the section `Input Params`.
#' To install the required packages, see section `Installation` below.
#'
#' - Disclaimer:
#' This program is free software: you can redistribute it and/or modify
#' it under the terms of the GNU General Public License as published by
#' the Free Software Foundation, either version 3 of the License, or
#' (at your option) any later version.
#'
#' This program is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU General Public License for more details.
#'
#' You should have received a copy of the GNU General Public License
#' along with this program.  If not, see <https://www.gnu.org/licenses/>.


#----------------------#
##### Installation #####
#----------------------#

# install.packages('devtools')
# install.packages('rgeos')
# install.packages('sf')
# install.packages('raster')
# devtools::install_github("pedro-andrade-inpe/sits.validate")
# devtools::install_github("e-sensing/inSitu")
# devtools::install_github("e-sensing/lucCalculus")
# devtools::install_github("e-sensing/sits")


#----------------------#
##### Input Params #####
#----------------------#

# define output path
outputDir  <- "~/MT" # Output directory

# classification memory and processors to be used
mem_size   <- 2      # Max memory to be used (in GB)
processors <- 1      # For Windows: set processors = 1


#------------------------------#
##### Function definitions #####
#------------------------------#

Process <- function(outputDir) {

    library(sits)
    library(magrittr)
    library(inSitu)

    #-----------------------------------#
    ##### 1. Process classification #####
    #-----------------------------------#

    # Create directory
    if (!dir.exists(paste(outputDir, "1.Classification", sep = "/")))
        dir.create(paste(outputDir, "1.Classification", sep = "/"), recursive = TRUE)

    # Get samples from inSitu package

    data(br_mt_1_8K_9classes_6bands)

    # define bands to work with
    bands <- c("evi", "ndvi", "nir", "mir")

    # select working bands from time series
    samples.tb <- br_mt_1_8K_9classes_6bands %>%
        sits_select_bands_(bands = bands)

    # train SVM model using selected bands
    model.svm <-
        samples.tb %>%
        sits_train(ml_method = sits_svm(cost = 1, formula = sits_formula_linear()))

    # create a coverage to classify MT
    cov.tb <- sits_coverage(service = "EOCUBES",
                            name = "MOD13Q1/006",
                            bands = bands,
                            geom = sf::read_sf(system.file("extdata/MT/shape/MT.shp", package = "inSitu")))

    # classify the raster image
    rasters.tb <- sits_classify_cubes(file = paste(paste(outputDir, "1.Classification", sep = "/"), "MT", sep = "/"),
                                      coverage = cov.tb,
                                      ml_model = model.svm,
                                      memsize = mem_size,
                                      multicores = processors)

    #----------------------------------#
    ##### 2. Apply Bayesian Filter #####
    #----------------------------------#

    # Create directory
    if (!dir.exists(paste(outputDir, "2.Classification_Smooth", sep = "/")))
        dir.create(paste(outputDir, "2.Classification_Smooth", sep = "/"), recursive = TRUE)

    # smooth parameters
    noise = 10
    window = matrix(1, nrow = 3, ncol = 3, byrow = TRUE)

    # apply the Bayesian filter on the output maps
    sits_bayes_postprocess(rasters.tb,
                           window = window,
                           noise = noise,
                           file = paste(paste(outputDir, "2.Classification_Smooth", sep = "/"), "smooth", sep = "/"))

    #---------------------------#
    ##### 3. Mosaic results #####
    #---------------------------#

    # Create directory
    if (!dir.exists(paste(outputDir, "3.Classification_Mosaic", sep = "/")))
        dir.create(paste(outputDir, "3.Classification_Mosaic", sep = "/"), recursive = TRUE)

    files_input <- list.files(paste(outputDir, "2.Classification_Smooth", sep = "/"), pattern = ".*\\.tif", full.names = TRUE)
    files_years <- gsub("^.*smooth_MT_[^_]{6}_[0-9]+_[0-9]+_[0-9]+_[0-9]+_([0-9]+)_.*\\.tif", "\\1", files_input)

    for (year in unique(files_years)) {

        year_list <- files_input[files_years == year]
        res <- lapply(year_list, raster::raster)
        res$filename <- paste(paste(outputDir, "3.Classification_Mosaic", sep = "/"), sprintf("MT_%s.tif", year), sep = "/")
        do.call(raster::merge, res)
    }
}

Process2 <- function(outputDir) {

  library(sits.validate)
  library(inSitu)
  library(lucCalculus)

  setBaseDir("/")

  #----------------------------#
  ##### 4. Urban Area (11) #####
  #----------------------------#

  message("Processing 4.UrbanArea_Mask...\n")

  # Urban area new pixel value
  new_px_urban_mask <- 11

  # list input files
  inputFiles <- getTifFiles(paste(outputDir, "3.Classification_Mosaic", sep = "/"))

  # get urban area mask file path
  inputFile_mask <- system.file("extdata/MT/masks/urban_area.tif", package = "inSitu")

  # create output path
  if (!dir.exists(paste(outputDir, "4.UrbanArea_Mask", sep = "/")))
      dir.create(paste(outputDir, "4.UrbanArea_Mask", sep = "/"), recursive = TRUE)

  # copy input files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "4.UrbanArea_Mask", sep = "/")), overwrite = TRUE)

  # get output tif files
  outputFiles <- getTifFiles(paste(outputDir, "4.UrbanArea_Mask", sep = "/"))

  # for each year do...
  for (pos in 1:NROW(outputFiles)) {

    message(paste0("Processing year ", pos, "/", NROW(outputFiles), "\n"))

    mask <- raster::raster(inputFile_mask) # mask
    sits_images <- raster::raster(outputFiles[pos]) # sits images
    sits_images <- raster::crop(sits_images, raster::extent(mask)) # crop sits images

    bs <- raster::blockSize(sits_images)
    for (i in 1:bs$n) {

      message(paste0("Processing block ", i, "/", bs$n, "\n"))
      row <- bs$row[i]
      nrows <- bs$nrows[i]

      block_mask = raster::getValues(mask, row = row, nrows = nrows)
      block_sits = raster::getValues(sits_images, row = row, nrows = nrows)

      pixel_mask <- block_mask == 11  # define new pixel value

      block_sits[which(pixel_mask)]  <- new_px_urban_mask

      sits_images[row:(row + nrows - 1),] <- block_sits
    }

    message(paste0("Writing raster\n"))

    raster::writeRaster(sits_images, outputFiles[pos], overwrite = TRUE)
  }

  message("4.UrbanArea_Mask Done!\n")

  #---------------------------#
  ##### 5. Sugarcane (10) #####
  #---------------------------#

  message("Processing 5.Sugarcane_Mask...\n")

  # Sugarcane new pixel value
  new_px_sugarc_mask <- 10

  # define input files as output files of previous step
  inputFiles <- getTifFiles(paste(outputDir, "4.UrbanArea_Mask", sep = "/"))

  # get sugarcane masks path
  inputDir_mask <- system.file("extdata/MT/masks/sugarcane", package = "inSitu") # open folder mask

  # list masks files
  inputDir_mask <- getTifFiles(inputDir_mask)

  # create output path
  if (!dir.exists(paste(outputDir, "5.Sugarcane_Mask", sep = "/")))
      dir.create(paste(outputDir, "5.Sugarcane_Mask", sep = "/"), recursive = TRUE)

  # copy sits files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "5.Sugarcane_Mask", sep = "/")), overwrite = TRUE)

  # list tif files of output path
  outputFiles <- getTifFiles(paste(outputDir, "5.Sugarcane_Mask", sep = "/"))

  # for years 2003 to 2016 do...
  for (pos in 3:16) {

    message(paste0("Processing year ", pos, "/", 16, "\n"))

    mask <- raster::raster(inputDir_mask[pos - 2]) # mask
    sits_images <- raster::raster(outputFiles[pos]) # sits images

    bs <- raster::blockSize(sits_images)

    for (i in 1:bs$n) {

      message(paste0("Processing block ", i, "/", bs$n, "\n"))
      row <- bs$row[i]
      nrows <- bs$nrows[i]

      block_mask = raster::getValues(mask, row = row, nrows = nrows)
      block_sits = raster::getValues(sits_images, row = row, nrows = nrows)

      pixel_mask <- block_mask == 10  # define new pixel value

      block_sits[which(pixel_mask)]  <- new_px_sugarc_mask

      sits_images[row:(row + nrows - 1),] <- block_sits
    }

    message(paste0("Writing raster\n"))

    raster::writeRaster(sits_images, outputFiles[pos], overwrite = TRUE)
  }

  message("5.Sugarcane_Mask Done!\n")


  #-----------------------#
  ##### 6. Water (12) #####
  #-----------------------#

  message("Processing 6.Water_Mask...\n")

  # Water new pixel value
  new_px_water_mask <- 12

  # define input files as output files of previous step
  inputFiles <- getTifFiles(paste(outputDir, "5.Sugarcane_Mask", sep = "/"))

  # get water mask file path
  inputFile_mask <- system.file("extdata/MT/masks/water.tif", package = "inSitu")

  # create output path
  if (!dir.exists(paste(outputDir, "6.Water_Mask", sep = "/")))
      dir.create(paste(outputDir, "6.Water_Mask", sep = "/"))

  # copy input files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "6.Water_Mask", sep = "/")), overwrite = TRUE)

  # get output tif files
  outputFiles <- getTifFiles(paste(outputDir, "6.Water_Mask", sep = "/"))

  # for each year do...
  for (pos in 1:NROW(outputFiles)) {

    message(paste0("Processing year ", pos, "/", NROW(outputFiles), "\n"))

    mask <- raster::raster(inputFile_mask) # mask
    sits_images <- raster::raster(outputFiles[pos]) # sits images

    bs <- raster::blockSize(sits_images)
    for (i in 1:bs$n) {

      message(paste0("Processing block ", i, "/", bs$n, "\n"))
      row <- bs$row[i]
      nrows <- bs$nrows[i]

      block_mask = raster::getValues(mask, row = row, nrows = nrows)
      block_sits = raster::getValues(sits_images, row = row, nrows = nrows)

      pixel_mask <- block_mask == 12  # define new pixel value

      block_sits[which(pixel_mask)]  <- new_px_water_mask

      sits_images[row:(row + nrows - 1),] <- block_sits
    }

    message(paste0("Writing raster\n"))

    raster::writeRaster(sits_images, outputFiles[pos], overwrite = TRUE)
  }

  message("6.Water_Mask Done!\n")


  #-----------------------------#
  ##### 7. Forest Mask (13) #####
  #-----------------------------#

  message("Processing 7.Forest_Mask...\n")

  # Forest mask new pixel value
  new_px_forest_mask <- 13

  # define input files as output files of previous step
  inputFiles <- getTifFiles(paste(outputDir, "6.Water_Mask", sep = "/"))

  # get forest mask mask file path
  inputFile_mask <- system.file("extdata/MT/masks/PRODES_2001_forest.tif", package = "inSitu")

  # create output path
  if (!dir.exists(paste(outputDir, "7.Forest_Mask", sep = "/")))
      dir.create(paste(outputDir, "7.Forest_Mask", sep = "/"))

  # copy input files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "7.Forest_Mask", sep = "/")), overwrite = TRUE)

  # get output tif files
  outputFiles <- getTifFiles(paste(outputDir, "7.Forest_Mask", sep = "/"))

  # for the first year 2001 do...
  mask <- raster::raster(inputFile_mask) # mask
  sits_images <- raster::raster(outputFiles[1]) # sits images

  bs <- raster::blockSize(sits_images)
  for (i in 1:bs$n) {

    message(paste0("Processing block ", i, "/", bs$n, "\n"))
    row <- bs$row[i]
    nrows <- bs$nrows[i]

    block_mask = raster::getValues(mask, row = row, nrows = nrows)
    block_sits = raster::getValues(sits_images, row = row, nrows = nrows)

    pixel_mask <- block_mask == 13  # define new pixel value

    block_sits[which(pixel_mask)]  <- new_px_forest_mask

    sits_images[row:(row + nrows - 1),] <- block_sits
  }

  message(paste0("Writing raster\n"))

  raster::writeRaster(sits_images, outputFiles[1], overwrite = TRUE)

  message("7.Forest_Mask Done!\n")


  #-----------------------------------------------------#
  ##### 8. PRODES 2001 Mask (Cerr.=14, Sec.Veg.=15) #####
  #-----------------------------------------------------#

  # Rules:
  # (Map_PRODES == Non_Forest (pixel 4) & map_MT == Forest (pixel 3)) --> map_MT == Cerrado (pixel 14)
  # (Map_PRODES == Deforestation (pixel 7) & map_MT == Forest (pixel 3)) --> map_MT == Secondary Vegetation (pixel 15)

  message("Processing 8.PRODES_2001_Mask...\n")

  # new pixel values
  new_px_value_mask_cerrado2 <- 14
  new_px_value_mask_sec_veg <- 15

  # define input files as output files of previous step
  inputFiles <- getTifFiles(paste(outputDir, "7.Forest_Mask", sep = "/"))

  # get forest mask mask file path
  inputFile_mask <- system.file("extdata/MT/masks/PRODES_2001.tif", package = "inSitu")

  # create output path
  if (!dir.exists(paste(outputDir, "8.PRODES_2001_Mask", sep = "/")))
      dir.create(paste(outputDir, "8.PRODES_2001_Mask", sep = "/"))

  # copy input files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "8.PRODES_2001_Mask", sep = "/")), overwrite = TRUE)

  # get output tif files
  outputFiles <- getTifFiles(paste(outputDir, "8.PRODES_2001_Mask", sep = "/"))

  # for year 2001 do...
  prodes2001 <- raster::raster(inputFile_mask) # PRODES
  sits_images <- raster::raster(outputFiles[1]) # sits with forest mask

  bs <- raster::blockSize(sits_images)
  for (i in 1:bs$n) {
    message(paste0("Processing block ", i, "/", bs$n, "\n"))
    row <- bs$row[i]
    nrows <- bs$nrows[i]

    block_prodes = raster::getValues(prodes2001, row = row, nrows = nrows)
    block_sits = raster::getValues(sits_images, row = row, nrows = nrows)

    forest_mt <- block_sits == 3  # sits Forest

    # PRODES
    non_forest_prodes <- block_prodes == 4      # define pixel number - Non_Forest (PRODES)
    deforestation_prodes <- block_prodes == 7   # define pixel number - Deforestation (PRODES)

    block_sits[which(forest_mt & non_forest_prodes)]  <- new_px_value_mask_cerrado2    # forest in sits and non-forest in PRODES
    block_sits[which(forest_mt & deforestation_prodes)] <- new_px_value_mask_sec_veg   # forest in sits and deforestation in PRODES

    sits_images[row:(row + nrows - 1),] <- block_sits
  }

  raster::writeRaster(sits_images, outputFiles[1], overwrite = TRUE)

  message("8.PRODES_2001_Mask Done!\n")


  #-------------------------------------------------------------#
  ##### 9. PRODES-Cerrado 2000 Mask (Cerr.=16, Sec.Veg.=17) #####
  #-------------------------------------------------------------#

  # Rules:
  # (Map_Cerrado == Non_Anthropized_Cerrado (pixel 1) & map_MT == Forest (pixel 3)) --> map_MT == Cerrado (pixel 16)
  # (Map_Cerrado == Anthropized_Cerrado (pixel 0) & map_MT == Forest (pixel 3)) --> map_MT == Secondary Vegetation (pixel 17)

  message("Processing 9.Cerrado_2001_Mask...\n")

  # new pixel values
  new_px_value_mask_cerrado3 <- 16
  new_px_value_mask_sec_veg2 <- 17

  # define input files as output files of previous step
  inputFiles <- getTifFiles(paste(outputDir, "8.PRODES_2001_Mask", sep = "/"))

  # get forest mask mask file path
  inputFile_mask <- system.file("extdata/MT/masks/PRODES-Cerrado_2000.tif", package = "inSitu")

  # create output path
  if (!dir.exists(paste(outputDir, "9.Cerrado_2001_Mask", sep = "/")))
      dir.create(paste(outputDir, "9.Cerrado_2001_Mask", sep = "/"))

  # copy input files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "9.Cerrado_2001_Mask", sep = "/")), overwrite = TRUE)

  # get output tif files
  outputFiles <- getTifFiles(paste(outputDir, "9.Cerrado_2001_Mask", sep = "/"))

  # For year 2001 do...
  cerrado2000 <- raster::raster(inputFile_mask) # Cerrado
  sits_images <- raster::raster(outputFiles[1]) # sits with forest mask

  bs <- raster::blockSize(sits_images)
  for (i in 1:bs$n) {

    message(paste0("Processing block ", i, "/", bs$n, "\n"))
    row <- bs$row[i]
    nrows <- bs$nrows[i]

    block_cerrado <- raster::getValues(cerrado2000, row = row, nrows = nrows)
    block_sits = raster::getValues(sits_images, row = row, nrows = nrows)

    forest_mt <- block_sits == 3  # sits Forest

    # Cerrado
    non_antro_cer <- block_cerrado == 0 # define pixel number - Non_Anthropized_Cerrado (Cerrado)
    antro_cer <- block_cerrado == 14    # define pixel number - Anthropized_Cerrado (Cerrado)

    block_sits[which(forest_mt & non_antro_cer)]  <- new_px_value_mask_cerrado3 # forest in sits and Non_Anthropized_Cerrado in Cerrado
    block_sits[which(forest_mt & antro_cer)] <- new_px_value_mask_sec_veg2      # forest in sits and Anthropized_Cerrado in Cerrado

    sits_images[row:(row + nrows - 1),] <- block_sits
  }

  raster::writeRaster(sits_images, outputFiles[1], overwrite = TRUE)

  message("9.Cerrado_2001_Mask Done!\n")


  #------------------------------------------------------------#
  ##### 10. Aggregate Classes (13->3; 14,16->1; 15,17->13) #####
  #------------------------------------------------------------#

  message("Processing 10.MT_BeforeLucCalc...\n")

  # Classes' codes used in output: 1-13

  # Aggregate classes to same labels (new pixel values)
  new_px_value_forest <- 3
  new_px_value_cerrado <- 1
  new_px_value_sec_veg <- 13

  # define input files as output files of previous step
  inputFiles <- getTifFiles(paste(outputDir, "9.Cerrado_2001_Mask", sep = "/"))

  # create output path
  if (!dir.exists(paste(outputDir, "10.MT_BeforeLucCalc", sep = "/")))
      dir.create(paste(outputDir, "10.MT_BeforeLucCalc", sep = "/"))

  # copy input files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "10.MT_BeforeLucCalc", sep = "/")), overwrite = TRUE)

  # get output tif files
  outputFiles <- getTifFiles(paste(outputDir, "10.MT_BeforeLucCalc", sep = "/"))

  # For each year do...
  for (pos in 1:NROW(outputFiles)) {

    message(paste0("Processing year ", pos, "/", NROW(outputFiles), "\n"))

    sits_images <- raster::raster(outputFiles[pos]) # sits images

    bs <- raster::blockSize(sits_images)

    for (i in 1:bs$n) {

      message(paste0("Processing block ", i, "/", bs$n, "\n"))
      row <- bs$row[i]
      nrows <- bs$nrows[i]

      block_sits = raster::getValues(sits_images, row = row, nrows = nrows)

      # forest
      forest_mask <- block_sits == 13

      # Cerrado
      cerrado2_mask <- block_sits == 14
      cerrado3_mask <- block_sits == 16

      # Sec. Vegetation
      secveg1_mask <- block_sits == 15
      secveg2_mask <- block_sits == 17

      block_sits[which(forest_mask)]  <- new_px_value_forest
      block_sits[which(cerrado2_mask | cerrado3_mask)] <- new_px_value_cerrado
      block_sits[which(secveg1_mask | secveg2_mask)] <- new_px_value_sec_veg

      sits_images[row:(row + nrows - 1),] <- block_sits
    }
    message(paste0("Writing raster\n"))
    raster::writeRaster(sits_images, outputFiles[pos], overwrite = TRUE)
  }

  message("10.MT_BeforeLucCalc Done!\n")


  #-----------------------------------#
  ##### 11. LUC Calculus (CF->CC) #####
  #-----------------------------------#

  message("Processing 11.LUC_Calc_newCF...\n")

  # Rule:
  # 1. (Cerrado)(Forest) --> (Cerrado)(Cerrado)

  # define input files as output files of previous step
  inputFiles <- getTifFiles(paste(outputDir, "10.MT_BeforeLucCalc", sep = "/"))

  # create output path
  dir.create(paste(outputDir, "11.LUC_Calc_newCF", sep = "/"))

  # copy input files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "11.LUC_Calc_newCF", sep = "/")), overwrite = TRUE)

  # get output tif files
  outputFiles <- getTifFiles(paste(outputDir, "11.LUC_Calc_newCF", sep = "/"))

  # for years 2 to 17 do...
  for (pos in 2:17) {
    message(paste0("Processing year ", pos, "/", 17, "\n"))

    raster0 = raster::raster(outputFiles[pos - 1])
    raster1 = raster::raster(outputFiles[pos])      # <-- output Value

    bs <- raster::blockSize(raster0)

    for (i in 1:bs$n) {

      message(paste0("Processing block ", i, "/", bs$n, "\n"))
      row = bs$row[i]
      nrows = bs$nrows[i]

      block0 = raster::getValues(raster0, row = row, nrows = nrows)
      block1 = raster::getValues(raster1, row = row, nrows = nrows)

      block1[which(block0 == 1 & block1 == 3)] <- 1

      raster1[row:(row + nrows - 1),] <- block1
    }

    message(paste0("Writing raster\n"))

    raster::writeRaster(raster1, outputFiles[pos], overwrite = TRUE)
  }

  message("11.LUC_Calc_newCF Done!\n")


  #---------------------------------------#
  ##### 12. LUC Calculus (CCPC->CCCC) #####
  #---------------------------------------#

  message("Processing 12.LUC_Calc_newCCPC...\n")

  # Rule:
  # 1. (Cerrado)(Cerrado)(Pasture)(Cerrado) --> (Cerrado)(Cerrado)(Cerrado)(Cerrado)

  # define input files as output files of previous step
  inputFiles <- getTifFiles(paste(outputDir, "11.LUC_Calc_newCF", sep = "/"))

  # create output path
  dir.create(paste(outputDir, "12.LUC_Calc_newCCPC", sep = "/"))

  # copy input files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "12.LUC_Calc_newCCPC", sep = "/")), overwrite = TRUE)

  # get output tif files
  outputFiles <- getTifFiles(paste(outputDir, "12.LUC_Calc_newCCPC", sep = "/"))

  # for years 3 to 16 do...
  for (pos in 3:16) {

    message(paste0("Processing year ", pos, "/", 16, "\n"))

    raster0 = raster::raster(outputFiles[pos - 2])
    raster1 = raster::raster(outputFiles[pos - 1])
    raster2 = raster::raster(outputFiles[pos])     # <- output value
    raster3 = raster::raster(outputFiles[pos + 1])

    bs <- raster::blockSize(raster0)

    for (i in 1:bs$n) {

      message(paste0("Processing block ", i, "/", bs$n, "\n"))
      row = bs$row[i]
      nrows = bs$nrows[i]

      block0 = raster::getValues(raster0, row = row, nrows = nrows)
      block1 = raster::getValues(raster1, row = row, nrows = nrows)
      block2 = raster::getValues(raster2, row = row, nrows = nrows)
      block3 = raster::getValues(raster3, row = row, nrows = nrows)

      block2[which(block0 == 1 & block1 == 1 & block2 == 4 & block3 == 1)] <- 1

      raster2[row:(row + nrows - 1),] <- block2
    }

    message(paste0("Writing raster\n"))

    raster::writeRaster(raster2, outputFiles[pos], overwrite = TRUE)
  }

  message("12.LUC_Calc_newCCPC Done!\n")


  #---------------------------------------#
  ##### 13. LUC Calculus (CCSC->CCCC) #####
  #---------------------------------------#

  message("Processing 13.LUC_Calc_newCCSC...\n")

  # Rule:
  # 1. (Cerrado)(Cerrado)(Soy*)(Cerrado) --> (Cerrado)(Cerrado)(Cerrado)(Cerrado)

  # define input files as output files of previous step
  inputFiles <- getTifFiles(paste(outputDir, "12.LUC_Calc_newCCPC", sep = "/"))

  # create output path
  dir.create(paste(outputDir, "13.LUC_Calc_newCCSC", sep = "/"))

  # copy input files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "13.LUC_Calc_newCCSC", sep = "/")), overwrite = TRUE)

  # get output tif files
  outputFiles <- getTifFiles(paste(outputDir, "13.LUC_Calc_newCCSC", sep = "/"))

  # for years 3 to 16 do...
  for (pos in 3:16) {

    message(paste0("Processing year ", pos, "/", 16, "\n"))

    raster0 = raster::raster(outputFiles[pos - 2])
    raster1 = raster::raster(outputFiles[pos - 1])
    raster2 = raster::raster(outputFiles[pos])     # <- output values
    raster3 = raster::raster(outputFiles[pos + 1])

    bs <- raster::blockSize(raster0)

    for (i in 1:bs$n) {

      message(paste0("Processing block ", i, "/", bs$n, "\n"))
      row = bs$row[i]
      nrows = bs$nrows[i]

      block0 = raster::getValues(raster0, row = row, nrows = nrows)
      block1 = raster::getValues(raster1, row = row, nrows = nrows)
      block2 = raster::getValues(raster2, row = row, nrows = nrows)
      block3 = raster::getValues(raster3, row = row, nrows = nrows)

      # soy1 = 5
      block2[which(block0 == 1 & block1 == 1 & block2 == 5 & block3 == 1)] <- 1
      # soy2 = 6
      block2[which(block0 == 1 & block1 == 1 & block2 == 6 & block3 == 1)] <- 1
      # soy3 = 7
      block2[which(block0 == 1 & block1 == 1 & block2 == 7 & block3 == 1)] <- 1
      # soy4 = 8
      block2[which(block0 == 1 & block1 == 1 & block2 == 8 & block3 == 1)] <- 1
      # soy5 = 9
      block2[which(block0 == 1 & block1 == 1 & block2 == 9 & block3 == 1)] <- 1

      raster2[row:(row + nrows - 1),] <- block2
    }

    message(paste0("Writing raster\n"))

    raster::writeRaster(raster2, outputFiles[pos], overwrite = TRUE)
  }

  message("13.LUC_Calc_newCCSC Done!\n")


  #-----------------------------------------#
  ##### 14. LUC Calculus (PPCCP->PPPPP) #####
  #-----------------------------------------#

  message("Processing 14.LUC_Calc_newPPCCP...\n")

  # Rule:
  # 1. (Pasture)(Pasture)(Cerrado)(Cerrado)(Pasture) --> (Pasture)(Pasture)(Pasture)(Pasture)(Pasture)

  # define input files as output files of previous step
  inputFiles <- getTifFiles(paste(outputDir, "13.LUC_Calc_newCCSC", sep = "/"))

  # create output path
  dir.create(paste(outputDir, "14.LUC_Calc_newPPCCP", sep = "/"))

  # copy input files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "14.LUC_Calc_newPPCCP", sep = "/")), overwrite = TRUE)

  # get output tif files
  outputFiles <- getTifFiles(paste(outputDir, "14.LUC_Calc_newPPCCP", sep = "/"))

  # for years 3 to 15 do...
  for (pos in 3:15) {

    message(paste0("Processing year ", pos, "/", 16, "\n"))

    raster0 = raster::raster(outputFiles[pos - 2])
    raster1 = raster::raster(outputFiles[pos - 1])
    raster2 = raster::raster(outputFiles[pos])     # <- output values
    raster3 = raster::raster(outputFiles[pos + 1]) # <- output values
    raster4 = raster::raster(outputFiles[pos + 2])

    bs <- raster::blockSize(raster0)

    for (i in 1:bs$n) {

      message(paste0("Processing block ", i, "/", bs$n, "\n"))
      row = bs$row[i]
      nrows = bs$nrows[i]

      block0 = raster::getValues(raster0, row = row, nrows = nrows)
      block1 = raster::getValues(raster1, row = row, nrows = nrows)
      block2 = raster::getValues(raster2, row = row, nrows = nrows)
      block3 = raster::getValues(raster3, row = row, nrows = nrows)
      block4 = raster::getValues(raster4, row = row, nrows = nrows)

      # define update block
      update_mask <- block0 == 4 & block1 == 4 & block2 == 1 & block3 == 1 & block4 == 4
      # update block2
      block2[which(update_mask)] <- 4
      # update block3
      block3[which(update_mask)] <- 4

      raster2[row:(row + nrows - 1),] <- block2
      raster3[row:(row + nrows - 1),] <- block3
    }

    message(paste0("Writing raster\n"))
    raster::writeRaster(raster2, outputFiles[pos], overwrite = TRUE)
    raster::writeRaster(raster3, outputFiles[pos + 1], overwrite = TRUE)
  }

  message("14.LUC_Calc_newPPCCP Done!\n")


  #---------------------------------------#
  ##### 15. LUC Calculus (FCFF->FFFF) #####
  #---------------------------------------#

  message("Processing 15.LUC_Calc_newFCFF...\n")

  # Rule:
  # 1. (Forest)(Cerrado)(Forest)(Forest) --> (Forest)(Forest)(Forest)(Forest)

  # define input files as output files of previous step
  inputFiles <- getTifFiles(paste(outputDir, "14.LUC_Calc_newPPCCP", sep = "/"))

  # create output path
  dir.create(paste(outputDir, "15.LUC_Calc_newFCFF", sep = "/"))

  # copy input files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "15.LUC_Calc_newFCFF", sep = "/")), overwrite = TRUE)

  # get output tif files
  outputFiles <- getTifFiles(paste(outputDir, "15.LUC_Calc_newFCFF", sep = "/"))

  # for years 2 to 15 do...
  for (pos in 2:15) {

    message(paste0("Processing year ", pos, "/", 15, "\n"))

    raster0 = raster::raster(outputFiles[pos - 1])
    raster1 = raster::raster(outputFiles[pos])     # <- output values
    raster2 = raster::raster(outputFiles[pos + 1])
    raster3 = raster::raster(outputFiles[pos + 2])

    bs <- raster::blockSize(raster0)

    for (i in 1:bs$n) {

      message(paste0("Processing block ", i, "/", bs$n, "\n"))
      row = bs$row[i]
      nrows = bs$nrows[i]

      block0 = raster::getValues(raster0, row = row, nrows = nrows)
      block1 = raster::getValues(raster1, row = row, nrows = nrows)
      block2 = raster::getValues(raster2, row = row, nrows = nrows)
      block3 = raster::getValues(raster3, row = row, nrows = nrows)

      block1[which(block0 == 3 & block1 == 1 & block2 == 3 & block3 == 3)] <- 3

      raster1[row:(row + nrows - 1),] <- block1
    }

    message(paste0("Writing raster\n"))
    raster::writeRaster(raster1, outputFiles[pos], overwrite = TRUE)
  }

  message("15.LUC_Calc_newFCFF Done!\n")


  #---------------------------------------#
  ##### 16. LUC Calculus (FFCF->FFFF) #####
  #---------------------------------------#

  message("Processing 16.LUC_Calc_newFFCF...\n")

  # Rule:
  # 1. (Forest)(Forest)(Cerrado)(Forest) --> (Forest)(Forest)(Forest)(Forest)

  # define input files as output files of previous step
  inputFiles <- getTifFiles(paste(outputDir, "15.LUC_Calc_newFCFF", sep = "/"))

  # create output path
  dir.create(paste(outputDir, "16.LUC_Calc_newFFCF", sep = "/"))

  # copy input files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "16.LUC_Calc_newFFCF", sep = "/")), overwrite = TRUE)

  # get output tif files
  outputFiles <- getTifFiles(paste(outputDir, "16.LUC_Calc_newFFCF", sep = "/"))

  # for years 3 to 16 do...
  for (pos in 3:16) {

    message(paste0("Processing year ", pos, "/", 16, "\n"))

    raster0 = raster::raster(outputFiles[pos - 2])
    raster1 = raster::raster(outputFiles[pos - 1])
    raster2 = raster::raster(outputFiles[pos])     # <- output values
    raster3 = raster::raster(outputFiles[pos + 1])

    bs <- raster::blockSize(raster0)

    for (i in 1:bs$n) {

      message(paste0("Processing block ", i, "/", bs$n, "\n"))
      row = bs$row[i]
      nrows = bs$nrows[i]

      block0 = raster::getValues(raster0, row = row, nrows = nrows)
      block1 = raster::getValues(raster1, row = row, nrows = nrows)
      block2 = raster::getValues(raster2, row = row, nrows = nrows)
      block3 = raster::getValues(raster3, row = row, nrows = nrows)

      block2[which(block0 == 3 & block1 == 3 & block2 == 1 & block3 == 3)] <- 3

      raster2[row:(row + nrows - 1),] <- block2
    }

    message(paste0("Writing raster\n"))
    raster::writeRaster(raster2, outputFiles[pos], overwrite = TRUE)
  }

  message("16.LUC_Calc_newFFCF Done!\n")


  #-------------------------------------#
  ##### 17. LUC Calculus (FCF->FFF) #####
  #-------------------------------------#

  message("Processing 17.LUC_Calc_newFCF...\n")

  # Rule:
  # 1. (Forest)(Cerrado)(Forest) --> (Forest)(Forest)(Forest)

  # define input files as output files of previous step
  inputFiles <- getTifFiles(paste(outputDir, "16.LUC_Calc_newFFCF", sep = "/"))

  # create output path
  dir.create(paste(outputDir, "17.LUC_Calc_newFCF", sep = "/"))

  # copy input files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "17.LUC_Calc_newFCF", sep = "/")), overwrite = TRUE)

  # get output tif files
  outputFiles <- getTifFiles(paste(outputDir, "17.LUC_Calc_newFCF", sep = "/"))

  # for years 2 to 16 do...
  for (pos in 2:16) {

    message(paste0("Processing year ", pos, "/", 16, "\n"))

    raster0 = raster::raster(outputFiles[pos - 1])
    raster1 = raster::raster(outputFiles[pos])     # <- output value
    raster2 = raster::raster(outputFiles[pos + 1])

    bs <- raster::blockSize(raster0)

    for (i in 1:bs$n) {

      message(paste0("Processing block ", i, "/", bs$n, "\n"))
      row = bs$row[i]
      nrows = bs$nrows[i]

      block0 = raster::getValues(raster0, row = row, nrows = nrows)
      block1 = raster::getValues(raster1, row = row, nrows = nrows)
      block2 = raster::getValues(raster2, row = row, nrows = nrows)

      block1[which(block0 == 3 & block1 == 1 & block2 == 3)] <- 3

      raster1[row:(row + nrows - 1),] <- block1
    }

    message(paste0("Writing raster\n"))
    raster::writeRaster(raster1, outputFiles[pos], overwrite = TRUE)
  }

  message("17.LUC_Calc_newFCF Done!\n")


  #-----------------------------------#
  ##### 18. LUC Calculus (FC->FF) #####
  #-----------------------------------#

  message("Processing 18.LUC_Calc_newFC...\n")

  # Rule:
  # 1. (Forest)(Cerrado) --> (Forest)(Forest)

  # define input files as output files of previous step
  inputFiles <- getTifFiles(paste(outputDir, "17.LUC_Calc_newFCF", sep = "/"))

  # create output path
  dir.create(paste(outputDir, "18.LUC_Calc_newFC", sep = "/"))

  # copy input files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "18.LUC_Calc_newFC", sep = "/")), overwrite = TRUE)

  # get output tif files
  outputFiles <- getTifFiles(paste(outputDir, "18.LUC_Calc_newFC", sep = "/"))

  # for years 2 to 17 do...
  for (pos in 2:17) {

    message(paste0("Processing year ", pos, "/", 17, "\n"))

    raster0 = raster::raster(outputFiles[pos - 1])
    raster1 = raster::raster(outputFiles[pos])     # <- output values

    bs <- raster::blockSize(raster0)

    for (i in 1:bs$n) {

      message(paste0("Processing block ", i, "/", bs$n, "\n"))
      row = bs$row[i]
      nrows = bs$nrows[i]

      block0 = raster::getValues(raster0, row = row, nrows = nrows)
      block1 = raster::getValues(raster1, row = row, nrows = nrows)

      block1[which(block0 == 3 & block1 == 1)] <- 3

      raster1[row:(row + nrows - 1),] <- block1
    }

    message(paste0("Writing raster\n"))
    raster::writeRaster(raster1, outputFiles[pos], overwrite = TRUE)
  }

  message("18.LUC_Calc_newFC Done!\n")
  invisible(NULL)
}

Process3 <- function(outputDir) {

  library(sits.validate)
  library(lucCalculus)

  setBaseDir("/")

  #----------------------------------#
  ##### 19. Secondary vegetation #####
  #----------------------------------#

  message("Processing 19.NewSecVeg...\n")

  # create a RasterBrick from individual raster saved previously
  lucC_create_RasterBrick(path_open_GeoTIFFs = paste(outputDir, "18.LUC_Calc_newFC", sep = "/"), path_save_RasterBrick = paste(outputDir, "19.NewSecVeg", sep = "/"))

  # define variables to use in sits
  file <- paste0(paste(outputDir, "19.NewSecVeg", sep = "/"), "/18.LUC_Calc_newFC_brick.tif")  #  brick

  # raster object
  rb_class0 <- raster::brick(file)

  message("\nCreating blocks with RasterBrick! ")
  # blocks saved in folder - Divide and conquer
  lucC_blocks_raster_create(raster_obj = rb_class0, number_blocks_xy = 6, path_save_RasterBlocks = paste(outputDir, "19.NewSecVeg", sep = "/"))

  # always
  options(digits = 12)

  # all blocks created in previous step into folder
  all.the.files <- list.files(paste(outputDir, "19.NewSecVeg", "Blocks_RasterBrick", sep = "/"), full = TRUE, pattern = ".tif")
  # all.the.files

  gc()

  # start time
  start.time <- Sys.time()

  for (y in 1:length(all.the.files)) {
    # file
    file <- all.the.files[y]

    # define variables to plot raster
    # new variable
    my_timeline <- lubridate::as_date(c("2001-09-01", "2002-09-01", "2003-09-01", "2004-09-01",
                                        "2005-09-01", "2006-09-01", "2007-09-01", "2008-09-01",
                                        "2009-09-01", "2010-09-01", "2011-09-01", "2012-09-01",
                                        "2013-09-01", "2014-09-01", "2015-09-01", "2016-09-01",
                                        "2017-09-01"))

    # original label - see QML file, same order
    my_label <- as.character(c("Cerrado", "Fallow_Cotton", "Forest", "Pasture", "Soy_Corn", "Soy_Cotton",
                               "Soy_Fallow", "Soy_Millet", "Soy_Sunflower", "Sugarcane", "Urban Area", "Water",
                               "Secondary Vegetation"))

    # create a RasterBrick metadata file based on the information about the files
    rb_class <- raster::brick(file)

    # verify if block with raster brick has NA in values
    if (!all(is.na(raster::maxValue(rb_class)))) {

      file_name <- basename(tools::file_path_sans_ext(file))

      message(paste0("\nLoad RasterBrick! Name: ``", file_name, "'' ...\n", sep = ""))

      # 1- Discover Secondary Vegetation - LUC Calculus
      # 1.1. RECUR predicate indicates a class that appear again
      message(paste0("Start RECUR ok! RasterBrick ...! Name: ``", file_name, "'' ...\n", sep = ""))
      #system.time(
      forest_recur <- lucC_pred_recur(raster_obj = rb_class, raster_class = "Forest",
                                      time_interval1 = c("2001-09-01","2001-09-01"),
                                      time_interval2 = c("2002-09-01","2017-09-01"),
                                      label = my_label, timeline = my_timeline,
                                      remove_column = TRUE)
      #head(forest_recur)

      message(paste0("Finish RECUR! RasterBrick ...! Name: ``", file_name, "'' ...\n", sep = ""))

      # 1.2. EVOLVE to verify Forest class that occurs after a different class in 2001
      forest_evolve <- NULL

      # classes without Forest based on original label
      classes <- as.character(c("Cerrado", "Fallow_Cotton", "Pasture", "Soy_Corn", "Soy_Cotton",
                                "Soy_Fallow", "Soy_Millet", "Soy_Sunflower", "Sugarcane", "Urban Area", "Water",
                                "Secondary Vegetation"))

      message(paste0("Start EVOLVE ... in RasterBrick ...! Name: ``", file_name, "'' ...\n", sep = ""))

      #system.time(
      # percor all classes
      for (i in seq_along(classes)) {

        print(classes[i])
        temp <- lucC_pred_evolve(raster_obj = rb_class, raster_class1 = classes[i],
                                 time_interval1 = c("2001-09-01","2001-09-01"), relation_interval1 = "equals",
                                 raster_class2 = "Forest",
                                 time_interval2 = c("2002-09-01","2017-09-01"), relation_interval2 = "contains",
                                 label = my_label, timeline = my_timeline, remove_column = TRUE)

        forest_evolve <- lucC_merge(forest_evolve, temp)
      }
      #)
      message(paste0("Finish EVOLVE ... in RasterBrick ...! Name: ``", file_name, "'' ...\n", sep = ""))

      # 1.3. Merge both forest_recur and forest_evolve datas
      forest_secondary <- lucC_merge(forest_recur, forest_evolve)

      rm(forest_recur, forest_evolve)
      gc()

      # 2 - Update original raster to add new pixel value
      message(paste0("Start update pixel in RasterBrick ...! Name: ``", file_name, "'' ...\n", sep = ""))

      number_label <- 13 # secondary vegetation class
      # 2.1. update original RasterBrick with new class
      rb_class_new <- lucC_raster_update(raster_obj = rb_class,
                                         data_mtx = forest_secondary,           # without 2001
                                         timeline = my_timeline,
                                         class_to_replace = "Forest",     # only class Cerrado
                                         new_pixel_value = number_label)  # new pixel value

      # new name
      new_file_name <- paste0(dirname(file),"/", file_name, sep = "")

      message(new_file_name)
      message(paste0("Save image pixel in RasterBrick ...! Name: ``", file_name, "'' ...\n", sep = ""))

      # 2.2. save the update matrix as GeoTIFF RasterBrick
      message(paste0("Prepare image 1 -- RasterBrick ...! Name: ``", file_name, "'' ...\n", sep = ""))
      lucC_save_GeoTIFF(raster_obj = rb_class,
                        data_mtx = rb_class_new,
                        path_raster_folder = new_file_name, as_RasterBrick = TRUE ) # FALSE before

      message(paste0("Image saved pixel in RasterBrick ...! Name: ``", new_file_name, "'' ...\n", sep = ""))
      gc()

    # if have NA then is save a copy of this block
    } else {

      file_name <- basename(tools::file_path_sans_ext(file))

      # new name
      new_file_name <- paste0(dirname(file),"/", file_name, sep = "")

      message(paste0("RasterBrick ...! Name: ``", new_file_name, "'' empty! ...\n", sep = ""))
    }

    message("--------------------------------------------------\n")

    gc()
  }

  # end time
  print(Sys.time() - start.time)

  # clear environment, except these elements
  rm(list = ls()[!(ls() %in% c("outputDir"))])
  gc()

  message("19.NewSecVeg Done!\n")

}

Process4 <- function(outputDir) {

  library(sits.validate)
  library(lucCalculus)

  setBaseDir("/")

  #--------------------------------------------#
  ##### 20. Merge all blocks for each band #####
  #--------------------------------------------#

  message("Processing 20.PostProc_Mosaic...\n")

  # start time
  start.time <- Sys.time()

  # create a folder called "Allblocks" and copy results of each Raster_Block_'n' to this folder
  if (!dir.exists(paste(outputDir, "20.PostProc_Mosaic", sep = "/")))
    dir.create(paste(outputDir, "20.PostProc_Mosaic", "Allblocks", sep = "/"), recursive = TRUE)

  files <- paste(list.dirs(paste(outputDir, "19.NewSecVeg", "Blocks_RasterBrick", sep = "/"))[-1:0],
                 paste0("New_", list.dirs(paste(outputDir, "19.NewSecVeg", "Blocks_RasterBrick", sep = "/"), full.names = F)[-1:0], ".tif"), sep = "/")

  file.copy(from = files, to = paste(outputDir, "20.PostProc_Mosaic", "Allblocks", sep = "/"))

  # merge blocks into a single image
  lucC_blocks_raster_merge(path_open_GeoTIFFs = paste(outputDir, "20.PostProc_Mosaic", "Allblocks", sep = "/"), pattern_name = "New_Raster_Block_", is.rasterBrick = TRUE)

  # save each layer of brick as images
  lucC_save_rasterBrick_by_layers(path_name_GeoTIFF_Brick = paste(outputDir, "20.PostProc_Mosaic", "Allblocks", "Mosaic_New_Raster_Block_.tif", sep = "/"))

  # end time
  print(Sys.time() - start.time)

  # Rename files and save in final folder
  # copy to a new directory
  list_of_files <- list.files(path = paste(outputDir, "20.PostProc_Mosaic", "Allblocks", "Mosaic_New_Raster_Block_", sep = "/"), full.names = TRUE, recursive = TRUE)
  file.copy(from = list_of_files, to = paste(outputDir, "20.PostProc_Mosaic", sep = "/"), overwrite = TRUE)

  # list of files
  list_of_files <- list.files(paste(outputDir, "20.PostProc_Mosaic", sep = "/"), full.names = TRUE, pattern = ".tif")
  file.rename(list_of_files, sub(pattern = "New__Mosaic_New_Raster_Block_.", replacement = "Mato_Grosso_", list_of_files))

  message("20.PostProc_Mosaic Done!\n")


  #----------------------------------#
  ##### 21. Crop and mask result #####
  #----------------------------------#

  message("Processing 21.Final...\n")

  library(raster)

  # list input files
  inputFiles <- getTifFiles(paste(outputDir, "20.PostProc_Mosaic", sep = "/"))

  # get urban area mask file path
  inputFile_mask <- system.file("extdata/MT/shape/MT.shp", package = "inSitu")
  shape <- sf::read_sf(inputFile_mask)

  # create output path
  if (!dir.exists(paste(outputDir, "21.Final", sep = "/")))
      dir.create(paste(outputDir, "21.Final", sep = "/"), recursive = TRUE)

  # copy input files to output path
  file.copy(inputFiles, baseDir(paste(outputDir, "21.Final", sep = "/")), overwrite = TRUE)

  # get output tif files
  outputFiles <- getTifFiles(paste(outputDir, "21.Final", sep = "/"))

  for (file in outputFiles) {

      r <- mask(crop(raster(file), extent(shape)), shape)
      writeRaster(r, file, overwrite = TRUE)
  }

  message("21.Final Done!\n")


  invisible(NULL)
}


#-------------#
##### Run #####
#-------------#

outputDir <- path.expand(outputDir)
if (!dir.exists(outputDir))
    dir.create(outputDir, recursive = TRUE)

Process(outputDir = outputDir)

Process2(outputDir = outputDir)

Process3(outputDir = outputDir)

Process4(outputDir = outputDir)

