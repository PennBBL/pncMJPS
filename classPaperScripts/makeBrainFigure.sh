#!/bin/bash 

# This script will be used to run through the production for all of the t stat images 
# For the jhu labels in the mj class paper
# It'll be a very rough script -_-

./makeJHUFigure.sh maleNvUT-KEY.csv 4 1
mv outputImage.nii.gz male/nvuT.nii.gz
fslmaths male/nvuT.nii.gz -add glassBrain male/nvuT.nii.gz
./makeJHUFigure.sh femaleNvUT-KEY.csv 4 1
mv outputImage.nii.gz female/nvuT.nii.gz
fslmaths female/nvuT.nii.gz -add glassBrain female/nvuT.nii.gz
./makeJHUFigure.sh maleUvFT-KEY.csv 4 1
mv outputImage.nii.gz male/uvfT.nii.gz
fslmaths male/uvfT.nii.gz -add glassBrain male/uvfT.nii.gz
./makeJHUFigure.sh femaleUvFT-KEY.csv 4 1
mv outputImage.nii.gz female/uvfT.nii.gz
fslmaths female/uvfT.nii.gz -add glassBrain female/uvfT.nii.gz
