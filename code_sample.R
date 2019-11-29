# Author: Xavier Cama
# This is a simplified sample script used for spectral analysis  
# 27/11/19

# Packages:
library(pls)
library(fields)



# CODE:

#######################################################
# IMPORT DATA 
#######################################################


# Function to import wavelength vector from PC

import_wavelengths_PC <- function (path, filename) {
  # The name of this function is import_wavelengths_PC
  
  path = paste(path, filename,".txt", sep="")
  data = read.table(path)
  answer = data[, 1] # Wavelengths are stored in the first column
  print(length(answer)) # Check if number of wavelengths is correct
  
  answer
}


# Fetching wavelengths

wavelengths <- import_wavelengths_PC(path = "C:\\data_spectral_analysis\\"
                                    ,filename = "S1_1 Meas. (1-1)")

num_wavelengths <- length(wavelengths)



# Function to import data from PC

import_data_PC <- function (path = "C:\\", filename = "sample", num_wavelengths) {
  # The name of this function is import_data_PC
  
  # Number of measurements per sample (E.g. 25 in a 5x5 grid pattern)
  num_replicates = 1
  rows = 5
  cols = 5
  num_measurements = rows * cols
  path_filename = paste(path, filename, sep = "")
  answer = matrix(0, num_replicates * num_measurements, num_wavelengths) # Creating matrix of dim (num_replicates, num_wavelengths) before loop
  
  for (i in 1:num_replicates) {
    for (k in 1:rows) {
      for (l in 1:cols) {
        
        j = num_measurements * (i - 1) + rows * (k - 1) + l # Vector to iterate between differnt files and store them in rows of the same matrix
        replicate_filenumber = rep(c(1, 2, 3), each = num_measurements) # Vector to iterate between repeats of the same sample
        path_data = paste(path_filename, replicate_filenumber[i], " Meas. (", k, "-", l, ").txt", sep = "")
        print(path_data)
        table = as.matrix(read.table(path_data)) 
        answer[j, ] = table[, 2] # Second column contains intensity
      }
    }
  }  
  answer
}



# Fetching spectra with different ion concentrations

spectra_sample1 <- import_data_PC(path = "C:\\data_spectral_analysis\\"
                         , filename = "S1_"
                         , num_wavelengths)
spectra_sample2 <- import_data_PC(path = "C:\\data_spectral_analysis\\"
                                  , filename = "S2_"
                                  , num_wavelengths)
spectra_sample3 <- import_data_PC(path = "C:\\data_spectral_analysis\\"
                                  , filename = "S3_"
                                  , num_wavelengths)
spectra_sample4 <- import_data_PC(path = "C:\\data_spectral_analysis\\"
                                  , filename = "S4_"
                                  , num_wavelengths)
spectra_sample5 <- import_data_PC(path = "C:\\data_spectral_analysis\\"
                                  , filename = "S5_"
                                  ,num_wavelengths)

# Validation sample
spectra_sample6 <- import_data_PC(path = "C:\\data_spectral_analysis\\"
                                  , filename = "S6_"
                                  , num_wavelengths) 



#######################################################
#SPECTRAL PRE-PROCESSING 
#######################################################


# Binding spectra together by rows

all_spectra <- rbind(spectra_sample1
                     , spectra_sample2
                     , spectra_sample3
                     , spectra_sample4
                     , spectra_sample5
                     , spectra_sample6)



# Average measurements

measurements_average <- function (data = raw, num_wavelengths) {
  # This function is called measurements_average
  
  num_measurements = 25
  num_concentrations = 6 # Number of dictinc samples
  num_replicates = 1
  num_samples = num_concentrations * num_replicates
  print(dim(data)[1])
  
  answer = matrix(0, num_samples, num_wavelengths)
  
  for (i in 1:num_samples) {
    for (j in 1:num_wavelengths) {
      
      # Get average in segments (from point a to point b)
      a = 1 + ((i-1) * num_measurements)
      b = num_measurements + ((i-1) * num_measurements)
      
      answer[i,j]=mean(data[(a:b),j])
    }
  }
  answer
}



# Applying average function to spectra

av_all_spectra <- measurements_average(all_spectra, num_wavelengths)



# Function to normalise spectra using standard normal variate (SNV)

snv <- function (data = raw) {
  #The name of this function is snv
  
  samples = dim(data)[1]
  wavelengths = dim(data)[2]  
  answer = matrix(0,nrow = samples,ncol = wavelengths)
  
  for(i in 1:samples) {
    
    answer[i, ] = (data[i, ] - mean(data[i, ]))/sqrt(var(data[i, ]))
  }
  answer
}



# Applying function snv to spectra

snv_av_spectra <- snv(av_all_spectra)
snv_spectra_sample6 <- snv(spectra_sample6)


# Function to plot spectra at once

func_plot <- function (wavelengths, data, main, sub, filename) 
{
  #The name of this function is func_plot
  
  num_replicates = 1
  dim_data = dim(data)[1]
  print(dim_data)

  
  colour_palette = c("blue", "green", "cyan", "purple", "orangered", "yellow") 
  vector = seq(1, dim_data)
  legend_text = c("6 mg", "4.5 mg", "3 mg", "2 mg", "1 mg", "validation (5.5 mg)")
  
  # Save images on hard drive
  #path = ("C:\\Users\\Dairypat\\Documents\\LIBS\\LWT new chamber\\")
  #jpeg(filename = paste(path,fname,".jpeg",sep = ""), width = 1000, height = 800, quality = 100)
  
  for (i in 1:dim_data) {
    if (i == 1) {
      
      plot(x = wavelengths, y = data[1,]
           , xlim = c(300,800), ylim = c(-1,30)
           , ann = FALSE
           , type = "l", lwd = 1
           , col = colour_palette[vector][i])
    }
    else {
      
      lines(x = wavelengths, y = data[i,]
            , type="l", lwd=1
            , col=colour_palette[vector][i])
      
      title(main = main, sub = sub) #Add title
      
      legend("topleft",inset=0.02
             , legend_text
             , col = colour_palette, bty = "n"
             , lwd = 1, text.font = 1, cex = 1.5) 
    }
  }
  #graphics.off()
}

func_plot(wavelengths, snv_av_spectra, main = "Spectra", sub = "SNV")



#######################################################
#PARTIAL LEAST SQUARES REGRESSION
#######################################################


# Fetching measured ion concentrations

concentrations <- as.matrix(read.csv(file = "C:\\real_concentrations.csv", header = FALSE, sep = ","))



# Preparing data frame files

train_pls <- data.frame(ion = I(snv_av_spectra[1:5,]), ion_concentration = I(concentrations[1:5]), row.names=NULL)

test_pls <- data.frame(ion = I(snv_av_spectra[6,]), ion_concentration = I(concentrations[6]), row.names=NULL)


# Building PLSR model using "pls" package 

pls_model <- plsr(ion_concentration~ion, ncomp = 3, data = train_pls, validation = "LOO")


# PLS model plot

# Save figure on hard drive
# tiff(paste(path = "C:\\", "PLS_model.tiff", sep = "")
#      , width = 3500, height = 3000, units = "px", bg = "white", res = 1000)
# windowsFonts(A = windowsFont("Times New Roman"))
# par(ps = 9, family = "A", mar = c(3.1,3.1,1.1,1.1))

plot(x = concentrations[1:5], y = pls_model$fitted.values[,,2] # [,,2] --> ncomp = 2
     , col = "black"
     , pch = 16, cex = 0.5, lwd = 1 # Point type (16), magnification and line width 
     , asp = 1, xaxp = c(0, 7, 7), yaxp = c(0, 7, 7) # Axis ticks
     , xlim = c(0, 7), ylim = c(0, 7) # Axis limits
     , cex.axis = 0.8, ann = FALSE, mgp = c(3, 0.5, 0)) # Magnification and margins  

# Axis text
mtext(side=1,text="Measured values, mg/g",line=2)
mtext(side=2,text="Predicted values, mg/g",line=2)

# Regression line
abline(lm(pls_model$fitted.values[,,2]~concentrations[1:5], col="black"))

# Figure legend
legend("topleft", inset=0.02
       , c("Calibration","Prediction")
       , col=c("black","forestgreen")
       , bty="n", pch = c(16,17)
       , pt.cex = 0.5, cex = 0.9)

#graphics.off()


# R2, RMSE and validation of the PLS model 

summary(pls_model)
R2(pls_model, estimate = "all", newdata = test_pls) 
RMSEP(pls_model, estimate = "all", newdata = test_pls)



#######################################################
# GENERATE IMAGE FROM PLS PREDICTION
#######################################################


prediction_image <- function (data = pls_model, newdata = validation, ncomp = ncomp) {
  #This function is called prediction_image
  
  rows=5
  cols=5

  pls_prediction = array(0, dim = c(rows * cols, 1, ncomp))
  answer = matrix(0, nrow = rows, ncol=cols)
  
  for (k in 1:rows) {
    for (l in 1:cols) {
      
      j = rows * (k - 1) + l # Vector to iterate between rows and columns     
      pls_prediction = predict(data, newdata, ncomp = ncomp, type = "response") # Using function predict pls package
      #print(dim(pls_array))
      answer[k,l] = pls_prediction[j, 1, 1]
    }
  }
  answer 
}

# Predicting pixels (5x5) for validation sample using PLS model

pred_spectra_sample6 <- prediction_image(data = pls_model
                                         , newdata = snv_spectra_sample6
                                         , ncomp = 1)


# Using image.plot() from package "fields"

image.plot(pred_spectra_sample6, zlim = c(2, 7), xaxt = "n", yaxt = "n")



#######################################################