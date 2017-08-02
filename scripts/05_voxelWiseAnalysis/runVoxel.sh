subjData="/data/joy/BBL/projects/pncMJPS/scripts/05_voxelWiseAnalysis/ctMJPSTDInteraction.rds"
outDirName="/data/joy/BBL/projects/pncMJPS/data/05_voxelWiseAnalysis"
namePath="pathVals"
maskName="/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1396_ctMask_thr9_2mm.nii.gz"
subjID="bblid"
inclusionName="inclusion"
smooth=0
covsFormula="~age+race+quality+psTD+marBin+psTD*marBin"
logFile="/data/joy/BBL/projects/pncMJPS/data/05_voxelWiseAnalysis/log/"
errFile="/data/joy/BBL/projects/pncMJPS/data/05_voxelWiseAnalysis/log/"

#Rscript /data/joy/BBL/applications/groupAnalysis/lm_voxelwise.R -c ${subjData} -o ${outDirName} -p ${namePath} -m ${maskName} -i ${inclusionName} -u ${subjID} -f ${covsFormula} -n 3 -s 0 -r FALSE -a "fdr" -d FALSE


# Now run cluster correction
subjData="/data/joy/BBL/projects/pncMJPS/scripts/05_voxelWiseAnalysis/ctMJPSTDInteraction.rds"
outDirName="/data/joy/BBL/projects/pncMJPS/data/05_voxelWiseAnalysis"
namePath="pathVals"
maskName="/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1396_ctMask_thr9_2mm.nii.gz"
subjID="bblid"
inclusionName="inclusion"
fullFormula="~age+race+quality+psTD+marBin+psTD*marBin"
redFormula="~age+race+quality+psTD+marBin"
splits=10
smooth=0
thresh=0.01
nsim=500
runCommand=FALSE


# CT w/ all MJ 
outDirName="/data/joy/BBL/projects/pncMJPS/data/05_voxelWiseAnalysis/ctAll"
mkdir -p ${outDirName}
#Rscript /data/joy/BBL/applications/groupAnalysis/gam_randomise.R -c $subjData -o $outDirName -p $namePath -m $maskName -i $inclusionName -u $subjID -f $fullFormula -r $redFormula -t 0.01 -n 1000 -e TRUE -s 0 -d FALSE

# Run CT w/ frequent MJ
inclusionName="inclusionFrequent"
outDirName="/data/joy/BBL/projects/pncMJPS/data/05_voxelWiseAnalysis/ctFrequent"
mkdir -p ${outDirName}
#Rscript /data/joy/BBL/applications/groupAnalysis/gam_randomise.R -c $subjData -o $outDirName -p $namePath -m $maskName -i $inclusionName -u $subjID -f $fullFormula -r $redFormula -t 0.01 -n 1000 -e TRUE -s 0 -d FALSE

# Now do GMD ALL
inclusionName="inclusion"
namePath="pathVals2"
outDirName="/data/joy/BBL/projects/pncMJPS/data/05_voxelWiseAnalysis/gmdAll"
mkdir -p ${outDirName}
#Rscript /data/joy/BBL/applications/groupAnalysis/gam_randomise.R -c $subjData -o $outDirName -p $namePath -m $maskName -i $inclusionName -u $subjID -f $fullFormula -r $redFormula -t 0.01 -n 1000 -e TRUE -s 0 -d FALSE

# Now do GMD w/ frequent
outDirName="/data/joy/BBL/projects/pncMJPS/data/05_voxelWiseAnalysis/gmdFrequent"
mkdir -p ${outDirName}
inclusionName="inclusionFrequent"
#Rscript /data/joy/BBL/applications/groupAnalysis/gam_randomise.R -c $subjData -o $outDirName -p $namePath -m $maskName -i $inclusionName -u $subjID -f $fullFormula -r $redFormula -t 0.01 -n 1000 -e TRUE -s 0 -d FALSE

# Now do CBF
maskName="/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/asl/n1601_PcaslCoverageMask.nii.gz"
outDirName="/data/joy/BBL/projects/pncMJPS/data/05_voxelWiseAnalysis/cbfAll"
mkdir -p ${outDirName}
subjData="/data/joy/BBL/projects/pncMJPS/scripts/05_voxelWiseAnalysis/cbfMJPSTDInteraction.rds"
inclusionName="inclusion"
namePath="pathVals3"
# All
Rscript /data/joy/BBL/applications/groupAnalysis/gam_randomise.R -c $subjData -o $outDirName -p $namePath -m $maskName -i $inclusionName -u $subjID -f $fullFormula -r $redFormula -t 0.01 -n 1000 -e TRUE -s 0 -d TRUE

# Frequent
outDirName="/data/joy/BBL/projects/pncMJPS/data/05_voxelWiseAnalysis/cbfFrequent"
mkdir -p ${outDirName}
inclusionName="inclusionFrequent"
Rscript /data/joy/BBL/applications/groupAnalysis/gam_randomise.R -c $subjData -o $outDirName -p $namePath -m $maskName -i $inclusionName -u $subjID -f $fullFormula -r $redFormula -t 0.01 -n 1000 -e TRUE -s 0 -d FALSE
