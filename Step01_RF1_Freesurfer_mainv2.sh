#!/bin/bash
#$ -S /bin/bash
datadir=/mnt/i/RF1/data
atlas=/mnt/e/hipp_atlas/ashs_atlas_umcutrecht_7t_20170810_ZY
count=0
# run freesurfer on RF1 data
rm -f freesurfer_parallel.txt
for subfold in $(find $datadir -maxdepth 1 -mindepth 1 -type d ); 
do
	# echo $subfold
	# rm -r -f $subpath/freesurfer
	subid=$(basename -- $subfold)
	subpath=$subfold/7T_STUDY${subid:0:5}/nifti1
	T1nii=$(find $subpath/cstfl_* -name s202*.nii) 
	T2nii=$(find $subpath/tse_* -name s202*.nii) #$subfolder/tse.nii
	T1dir=$(dirname $T1nii)
	T2dir=$(dirname $T1nii)
	if [[ ! -f $T1dir/T1_biascorrected.nii ]] && [[ -f $T1nii ]]; then
		echo "Run bias correction for $T1nii"
		N4BiasFieldCorrection -d 3 -i $T1nii -o $T1dir/T1_biascorrected.nii
	fi
	if [[ -f $T1dir/T1_biascorrected.nii ]] && [[ ! -f $subpath/freesurfer/stats/aseg.stats ]]; then
		if  [[ -d $subpath/freesurfer ]]; then
		  rm -r $subpath/freesurfer
		fi
		if [[ -f $T2nii ]]; then
		  if [[ ! -f $T2dir/tse_biascorrected.nii ]]; then
			echo "Run bias correction for $T2nii"
			N4BiasFieldCorrection -d 3 -i $T2nii -o $T2dir/tse_biascorrected.nii
		  fi
		  echo "recon-all -sd $subpath -s freesurfer -i $T1dir/T1_biascorrected.nii  -all" >> freesurfer_parallel.txt #-T2 $subpath/tse_biascorrected.nii
		else
		  echo "recon-all -sd $subpath -s freesurfer -i $T1dir/T1_biascorrected.nii -all " >> freesurfer_parallel.txt
		fi
	fi
done
parallel < freesurfer_parallel.txt

for subfold in $(find $datadir -maxdepth 1 -mindepth 1 -type d ); 
do
	# echo $subfold
	subid=$(basename -- $subfold)
	subpath=$subfold/7T_STUDY${subid:0:5}/nifti1
	if [[ -f $subpath/freesurfer/mri/aparc+aseg.mgz ]] && [[ ! -f $subpath/freesurfer/mri/aparc+aseg.nii.gz ]] ; then
		echo "Convert .mgz to .nii.gz for $subfold"
		mri_convert $subpath/freesurfer/mri/aparc+aseg.mgz $subpath/freesurfer/mri/aparc+aseg.nii.gz
	fi
done

rm -f FSsegmentHA.txt
for subfold in $(find $datadir -maxdepth 1 -mindepth 1 -type d ); 
do
	subid=$(basename -- $subfold)
	subpath=$subfold/7T_STUDY${subid:0:5}/nifti1
	T1nii=$(find $subpath/cstfl_* -name s202*.nii) 
	T2nii=$(find $subpath/tse_* -name s202*.nii) #$subfolder/tse.nii
	T1dir=$(dirname $T1nii)
	T2dir=$(dirname $T1nii)
	if [ -f "$T1dir/T1_biascorrected.nii" ] && [ -f "$subpath/freesurfer/stats/aseg.stats" ]; then
		rm $subpath/freesurfer/scripts/IsRunningHPsubT2.T1T2.lh+rh
		if [[ -f $T2dir/tse_biascorrected.nii ]] && [[ ! -f "$subpath/freesurfer/mri/lh.amygNucVolumes-T1-T1T2.v22.txt" ]]; then
			echo "$subfold: Run Freesurfer Hippo subfield segmentation with T1 and T2."
			# segmentHA_T2.sh  bert  FILE_ADDITIONAL_SCAN   ANALYSIS_ID  USE_T1  [SUBJECTS_DIR]
			echo "segmentHA_T2.sh freesurfer $T2dir/tse_biascorrected.nii T1T2 1 $subpath " >> FSsegmentHA.txt 
		elif [[ ! -f $T2dir/tse_biascorrected.nii ]] && [[ ! -f "$subpath/freesurfer/mri/lh.amygNucVolumes-T1.v22.txt" ]]; then
			echo "$subfold: Run Freesurfer Hippo subfield segmentation with T1 only."
			echo "segmentHA_T1.sh freesurfer $subpath " >> FSsegmentHA.txt 
		fi
	fi
done
parallel < FSsegmentHA.txt

for subfold in $(find $datadir -maxdepth 1 -mindepth 1 -type d ); 
do
	# echo $subfold
	subid=$(basename -- $subfold)
	subpath=$subfold/7T_STUDY${subid:0:5}/nifti1
	if [[ -f $subpath/freesurfer/mri/lh.hippoAmygLabels-T1-T1T2.v22.HBT.mgz ]] && [[ ! -f $subpath/freesurfer/mri/lh.hippoAmygLabels-T1-T1T2.v22.HBT.nii.gz ]] ; then
		mri_convert $subpath/freesurfer/mri/rh.hippoAmygLabels-T1-T1T2.v22.HBT.mgz $subpath/freesurfer/mri/rh.hippoAmygLabels-T1-T1T2.v22.HBT.nii.gz
		mri_convert $subpath/freesurfer/mri/lh.hippoAmygLabels-T1-T1T2.v22.HBT.mgz $subpath/freesurfer/mri/lh.hippoAmygLabels-T1-T1T2.v22.HBT.nii.gz
		
		mri_convert $subpath/freesurfer/mri/lh.hippoAmygLabels-T1-T1T2.v22.CA.mgz $subpath/freesurfer/mri/lh.hippoAmygLabels-T1-T1T2.v22.CA.nii.gz
		mri_convert $subpath/freesurfer/mri/rh.hippoAmygLabels-T1-T1T2.v22.CA.mgz $subpath/freesurfer/mri/rh.hippoAmygLabels-T1-T1T2.v22.CA.nii.gz
	fi
	
	if [[ -f $subpath/freesurfer/mri/lh.hippoAmygLabels-T1.v22.HBT.mgz ]] && [[ ! -f $subpath/freesurfer/mri/lh.hippoAmygLabels-T1.v22.HBT.nii.gz ]] ; then
		mri_convert $subpath/freesurfer/mri/rh.hippoAmygLabels-T1.v22.HBT.mgz $subpath/freesurfer/mri/rh.hippoAmygLabels-T1.v22.HBT.nii.gz
		mri_convert $subpath/freesurfer/mri/lh.hippoAmygLabels-T1.v22.HBT.mgz $subpath/freesurfer/mri/lh.hippoAmygLabels-T1.v22.HBT.nii.gz
		
		mri_convert $subpath/freesurfer/mri/lh.hippoAmygLabels-T1.v22.CA.mgz $subpath/freesurfer/mri/lh.hippoAmygLabels-T1.v22.CA.nii.gz
		mri_convert $subpath/freesurfer/mri/rh.hippoAmygLabels-T1.v22.CA.mgz $subpath/freesurfer/mri/rh.hippoAmygLabels-T1.v22.CA.nii.gz
	fi
done