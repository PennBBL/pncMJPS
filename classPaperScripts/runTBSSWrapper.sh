#!/bin/bash

# This script will be used to do the preprocessing, prepare the regestration, postreg, prestats,
# and finally run randomize and or Angel's voxelwise wrapper

# The first step is a for loop to put all of the images into the proper directory
# First declare the statics
inputID="/home/arosen/pncMJPS/classPaperScripts/maleIDValues.csv"
inputData="/data/joy/BBL/studies/pnc/processedData/diffusion/pncDTI_2016_04"
outputData="/data/joy/BBL/projects/pncMJPS/data/mytbss/"

## Now loop trhough each line in the id file and grab and link and FA image
## to our outpout directory
loopLength=`cat ${inputID} | wc -l`
for i in `seq 2 ${loopLength}` ; do
  bblid=`sed -n -e ${i}p ${inputID} | cut -f 1 -d ,`
  scanid=`sed -n -e ${i}p ${inputID} | cut -f 2 -d ,`

  # Now find our input data
  inputImg="${inputData}/${bblid}/*${scanid}/DTI_64/dtitk/raw/${bblid}_*x${scanid}_dti_eddy_rbvecs_dtitk_aff_diffeo_fa.nii.gz"

  # Now if the image exists link it to our output file
  if [ -f ${inputImg} ] ; then
    ln ${inputImg} ${outputData}/${bblid}_${scanid}_FA.nii.gz ;
  fi

  # Now echo something out so I know I am running
  echo "${bblid}, ${scanid} done" ;
done

# Now run the mytbss preproc script
# This step takes a lot of time for some dumb reason...
cd ${outputData}
tbss_1_preproc *nii.gz

# Now onto step 2! the registration
# This is to exciting for me to handle...
tbss_2_reg -T

## The regeistration takes about 10 min per subject
## Lets round up to 15 and wait that amount of time

