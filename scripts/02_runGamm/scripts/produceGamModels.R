# AFGR June 2017
# This script is going to be used to run the GAM models for the PS-MJ study
# The steps include:
#	1.) Running GAM's for each modality for each ROI
#		a. The model will look like:"brainData ~ s(age) + PS*MJ + sex + dataQuality"
# 	2.) Reporting significant intgeractions from PS & MJ

# Load all data and library(s)
source('/data/joy/BBL/projects/pncMJPS/scripts/01_dataPrep/scripts/dataPrep.R')
source('/data/joy/BBL/projects/pncMJPS/scripts/02_runGamm/functions/functions.R')
install_load('mgcv')

# Now run the gam's for the structural data 
# Collapse all fo the Users into MJ User bin for marcat 
strucData$marcat[strucData$marcat=='MJ Frequent User'] <- 'MJ User'

