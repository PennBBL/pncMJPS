#!/bin/bash
# f and p value maps for the null effects project
baseDir="/data/jux/BBL/projects/pncMJPS/data/"
rData="${baseDir}colorAndKey/"
modality=( volume ct gmd ) #ct reho alff tr)
contrast=( f p tno pno tfo pfo tof pof )
scriptName="/home/arosen/hiLo/scripts/05_BrainRankFigure/scripts/makeZScoreJLFPNCTemplateImage.sh"
foo() {
  baseDir="/data/jux/BBL/projects/pncMJPS/data/"
  rData="${baseDir}colorAndKey/"
  m=$1
  c=$2
  scriptName="/home/arosen/hiLo/scripts/05_BrainRankFigure/scripts/makeZScoreJLFPNCTemplateImage.sh"
  inputCSV="${rData}${m}${c}-KEY.csv"
  mkdir -p ${baseDir}/imagingFigures/${m}/${c}/ 
  cd ${baseDir}/imagingFigures/${m}/${c}/ 
  ${scriptName} ${inputCSV} 4 0 
  fslmaths outputImage.nii.gz -edge -bin tmp.nii.gz
  fslmaths outputImage.nii.gz -mul tmp.nii.gz edge.nii.gz
  fslmaths outputImage.nii.gz -sub edge.nii.gz outputImage.nii.gz
  fslmaths tmp.nii.gz -mul 5000 tmp.nii.gz
  fslmaths outputImage.nii.gz -add tmp.nii.gz outputImage.nii.gz 
  rm tmp.nii.gz edge.nii.gz
  echo "All done"
}
for s in ${contrast[*]} ; do
    for q in ${modality[*]} ; do 
      foo $q $s &  
    done
done
exit 66


modality=( volume ct gmd ) #ct reho alff tr)
contrast=( p tno pno tfo pfo tof pof )
for s in ${contrast[*]} ; do 
  for q in ${modality[*]} ; do 
    cd ${baseDir}/imagingFigures/${q}/${s}/ 
    itksnap -g /home/arosen/templateMNI/MNI152_T1_1mm_brain.nii.gz -s ${baseDir}/imagingFigures/${q}/${s}/outputImage.nii.gz -l ${rData}${q}${s}-ColorTable.txt & 
  done ; 
  #sleep 60
done
