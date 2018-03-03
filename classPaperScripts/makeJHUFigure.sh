#!/bin/bash

# This script is going to be used to turn the raw pnc JLF labeled image into a z score image.
# This is going to be an experiment, but hopefully it should be a very flexible script
# The required inputs will be a csv with the 

## Declare the usage function
usage(){
echo 
echo
echo
echo "${0} <inputZScore.csv> <colOfInterest(integer)> <turnOffGlassBrain> <brainRegion>"
echo "inputZScore.csv is a csv with a column of the roi names from the jlf segmentation and the z scores"
echo "colOfInterest is the numeric value of the column that you want to replace the intensity values with"
echo
echo
exit 1
}

# First lets declare all of the static variables
jlfLookUp="/home/arosen/pncMJPS/data/outputValues.csv"
pncJlfLabelImage="/home/arosen/JHU_WM_labels_GO1_template_space.nii.gz"
if [ ! "X${4}" == "X" ] ; then 
  jlfLookUp="/data/joy/BBL/studies/pnc/template/jlf/hiLoLookup/jlf_lookup${4}.csv" ; 
fi 
echo ${jlfLookUp}
inputCSV=${1}
colOfInterest=${2}
glassBrain=${3}
workingDir=`pwd`
dateValue=`date +%y_%m_%d_%H_%M_%S`
outputImag="${workingDir}/jlfLabelImage"
tmpDir="/tmp/jlfParcel${RANDOM}/"
newValueCSV="${tmpDir}valueCSV.csv"


# Now lets make sure all of our statc variables exist 
if [ ! -f ${inputCSV}  ] ; then 
  echo
  echo 
  echo "Not all statics are present"
  echo
  usage ; 
fi

if [ "X${colOfInterest}" == "X" ] ; then 
  usage ; 
fi

# Make our tmp dir
mkdir ${tmpDir}

# Now I need to create a blank image in the output directory
fslmaths ${pncJlfLabelImage} -sub ${pncJlfLabelImage} ${tmpDir}blankImage.nii.gz

# Now I need to go through each of the values in the in the jlfLookUp and look for this roi in the
# input csv, if they are in the input csv then I will put them into the 
# ** I will need to be very careful about ensuring that each ROI returns 1 line **
loopLength=`cat ${jlfLookUp} | wc -l`
for lineValue in `seq 1 ${loopLength}` ; do 
  specLine=`sed -n -e ${lineValue}p ${jlfLookUp}` 
  roiToGrep=`echo ${specLine} | cut -f 2 -d ,` 
  echo ${roiToGrep}
  # Now lets grep our ROI in our input csv
  valueToFind=`grep "${roiToGrep}" ${inputCSV}`
  if [ ! "X${valueToFind}" == "X" ] ; then 
    # first lets ensure that only one ROI is returned
    quickCheck=`grep -c "${roiToGrep}" ${inputCSV}`
    if [ ${quickCheck} -gt 1 ] ; then 
      echo "Now fixing double match for ${roiToGrep}"
      unset arrayName
      declare -a arrayName
      for variableName in `seq 1 ${quickCheck}` ; do
        tmpString=`echo ${valueToFind} | cut -f ${variableName} -d ' '`
        newString="${tmpString//${roiToGrep}}"
        newStringLength="${#newString}"
        arrayName=("${arrayName[@]}" "${newStringLength}") ; 
      done ; 
      # Now we need to find the minimum value in our array
      minValue=`echo ${arrayName[@]} | tr ' ' '\n' | sort | head -n 1`
      fieldValue=`echo ${arrayName[@]} | tr ' ' '\n' | grep -n "${minValue}" | cut -f 1 -d:`
      if [ ${#fieldValue} -gt 1 ] ; then
        fieldValue=`echo ${fieldValue} | cut -f 1 -d ' '` ; 
      fi
      specLine=`echo ${specLine} | cut -f ${fieldValue} -d ' '`
      valueToFind=`echo ${valueToFind} | cut -f ${fieldValue} -d ' '`
    fi
    intensityValue=`echo ${specLine} | cut -f 1 -d ,` 
    newValue=`echo ${valueToFind} | cut -f ${colOfInterest} -d ","`
    echo "${intensityValue},${newValue}" >> ${newValueCSV} ; 
  elif [ ${glassBrain} -gt 0 ] ; then 
    intensityValue=`echo ${specLine} | cut -f 1 -d ,` 
    newValue="1616"
    echo "${intensityValue},${newValue}" >> ${newValueCSV} ;     
  fi 
done

# Now the slow part... we need to thrshold out all of our label intensity values and multiply them by their new value
loopLength=`cat ${newValueCSV} | wc -l`
mv ${tmpDir}blankImage.nii.gz ${tmpDir}outputImage.nii.gz
for fileLine in `seq 1 ${loopLength}` ; do 
  intensityValue=`sed -n -e ${fileLine}p ${newValueCSV} | cut -f 1 -d ,`
  newValue=`sed -n -e ${fileLine}p ${newValueCSV} | cut -f 2 -d ,`
  fslmaths ${pncJlfLabelImage} -thr ${intensityValue} -uthr ${intensityValue} -bin -mul ${newValue} ${tmpDir}tmp_${intensityValue} 
  fslmaths ${tmpDir}outputImage.nii.gz -add ${tmpDir}tmp_${intensityValue} ${tmpDir}outputImage.nii.gz
  echo -ne "${fileLine} of ${loopLength}\033[0K\r" ; 
done

# Now lets move our output image to our cwd and rm the tmp directory!
mv ${tmpDir}outputImage.nii.gz ${workingDir}/
rm -rf ${tmpDir}
exit 0
