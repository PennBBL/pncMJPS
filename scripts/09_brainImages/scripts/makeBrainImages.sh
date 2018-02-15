../scripts/09_brainImages/scripts/makeZScoreJLFPNCTemplateImage.sh allValsA-KEY.csv 4 1
mv outputImage.nii.gz nonVsUser/aucImg.nii.gz
../scripts/09_brainImages/scripts/makeZScoreJLFPNCTemplateImage.sh allValsT-KEY.csv 4 1
mv outputImage.nii.gz nonVsUser/tImg.nii.gz
../scripts/09_brainImages/scripts/makeZScoreJLFPNCTemplateImage.sh allValsUvFA-KEY.csv 4 1
mv outputImage.nii.gz userVsFreq/aucImg.nii.gz
../scripts/09_brainImages/scripts/makeZScoreJLFPNCTemplateImage.sh allValsUvFT-KEY.csv 4 1
mv outputImage.nii.gz userVsFreq/tImg.nii.gz
