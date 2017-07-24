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

Rscript /data/joy/BBL/applications/groupAnalysis/lm_voxelwise.R -c ${subjData} -o ${outDirName} -p ${namePath} -m ${maskName} -i ${inclusionName} -u ${subjID} -f ${covsFormula} -n 3 -s 0 -r FALSE -a "fdr" -d FALSE
