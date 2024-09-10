#!/bin/bash
#$ -S /bin/bash
ARCH=$(uname);
EXT_BIN=$LIB_path/ext/$ARCH/bin
PATH=/usr/local/ashs/ashs-fastashs-20170915/ext/Linux/bin:$PATH
LIB_path=/usr/local/hippodeep_pytorch-master

DataDir=/mnt/i/RF1/data
  
SUBFOLDERS=$(find $DataDir -maxdepth 1 -mindepth 1 -type d ); 
count=0
HippoSegnames=( "/ashs/final/ashs_left_lfseg_corr_nogray.nii.gz" "/ashs/final/ashs_right_lfseg_corr_nogray.nii.gz" "/freesurfer/mri/lh.hippoAmygLabels-T1-T1T2.v22.CA.FSvoxelSpace.mgz" "/freesurfer/mri/rh.hippoAmygLabels-T1-T1T2.v22.CA.FSvoxelSpace.mgz" "/freesurfer/mri/lh.hippoAmygLabels-T1-T1T2.v22.HBT.FSvoxelSpace.mgz" "/freesurfer/mri/rh.hippoAmygLabels-T1-T1T2.v22.HBT.FSvoxelSpace.mgz" "/freesurfer/mri/lh.hippoAmygLabels-T1.v22.CA.FSvoxelSpace.mgz" "/freesurfer/mri/rh.hippoAmygLabels-T1.v22.CA.FSvoxelSpace.mgz" "/freesurfer/mri/lh.hippoAmygLabels-T1.v22.HBT.FSvoxelSpace.mgz" "/freesurfer/mri/rh.hippoAmygLabels-T1.v22.HBT.FSvoxelSpace.mgz" "/HippoSeg_CCF/crop/mv_seg_L.nii.gz" "/HippoSeg_CCF/crop/mv_seg_R.nii.gz" )
HippoSegsavenames=( "CA_ashs_fMRI_lh.nii.gz" "CA_ashs_fMRI_rh.nii.gz" "CA_FS_fMRI_lh.nii.gz" "CA_FS_fMRI_rh.nii.gz" "HBT_FS_fMRI_lh.nii.gz" "HBT_FS_fMRI_rh.nii.gz" "CA_FS_fMRI_lh.nii.gz" "CA_FS_fMRI_rh.nii.gz" "HBT_FS_fMRI_lh.nii.gz" "HBT_FS_fMRI_rh.nii.gz" "CA_HippoSeg_fMRI_lh.nii.gz" "CA_HippoSeg_fMRI_rh.nii.gz" )
for subfold in $SUBFOLDERS; do
  subid=$(basename -- $subfold)
  subfolder=$subfold/7T_STUDY${subid:0:5}/nifti1
  subid=$(basename -- $subfolder)
  c3T1=$(find $subfolder/cstfl_* -name c1s20*.nii)
  c1T1=$(find $subfolder/cstfl_* -name c2s20*.nii)
  c2T1=$(find $subfolder/cstfl_* -name c3s20*.nii)
  T1nii=$(find $subfolder/cstfl_* -name STT1.nii) #$subfolder/STT1.nii
  T2nii=$(find $subfolder/tse_* -name s*.nii) #$subfolder/tse.nii
  T1Orignii=$(find $subfolder/cstfl_* -name s20*.nii) #subfolder/T1.nii
  fMRI2T1_affine=$subfolder/ants_registration/fMRI2T10GenericAffine.mat
  fMRI2T1nomask_affine=$subfolder/ants_registration/fMRI2T1nomask0GenericAffine.mat
  T22T1_affine=$subfolder/ants_registration/T22T10GenericAffine.mat
  fMRISpace_dir=$subfolder/ants_registration/fMRISpace
  fMRImean=$(find $subfolder/*FMRI* -name STmean.nii )
  # fMRI to T1 space registration
  if [ ! -f "$fMRI2T1_affine" ]; then
	if [ -f "$fMRImean" ] && [ -f "$T1nii" ] && [ -f "$c3T1" ]; then
	  echo "fMRI to T1 registration: $subfolder"
	  mkdir -p $subfolder/ants_registration
	  #rm $subfolder/ants_registration/c3_dilatedmask.nii.gz
	  masknii=$subfolder/ants_registration/fMRI2T1_mask.nii.gz
	  c3d $c1T1 $c2T1 -add -thresh 0.8 inf 1 0 -dilate 1 2x2x2vox -o $masknii
	  
	  #c3d $c3T1 -thresh 0.8 inf 1 0 -erode 1 2x2x2vox -dilate 1 4x4x4vox -o $masknii
	  ants_cmd="antsRegistration --dimensionality 3 --float 1 --output [$subfolder/ants_registration/fMRI2T1,$subfolder/ants_registration/meanfMRI_in_T1.nii.gz] --interpolation Linear --winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 --transform Affine[0.2] --masks $masknii --metric MI[$T1nii,$fMRImean,1,32,Regular,0.20] --convergence [500x250x200x100,1e-8,10] --shrink-factors 8x4x2x1 	--smoothing-sigmas 3x2x1x0vox --verbose > $subfolder/ants_registration/fMRI2T1.log"
	  eval $ants_cmd
	fi
  fi
  
  # T2 to T1 registration
  if [ ! -f $T22T1_affine ] && [ -f $T2nii ] && [ -f $T1Orignii ]; then 
	  echo "T2 to T1 registration: $subfolder"
	  mkdir -p $subfolder/ants_registration
	  	  ants_cmd="antsRegistration --dimensionality 3 --float 1 --output [$subfolder/ants_registration/T22T1,$subfolder/ants_registration/T2_in_T1.nii.gz] --interpolation Linear --winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 --transform Rigid[0.2] --metric MI[$T1Orignii,$T2nii,1,32,Regular,0.20] --convergence [500x250x200x100,1e-8,10] --shrink-factors 8x4x2x1 	--smoothing-sigmas 3x2x1x0vox --verbose > $subfolder/ants_registration/T22T1.log"
	  eval $ants_cmd
	  #echo $ants_cmd >> /mnt/c/Users/zheng/Dropbox/Neuroscience/RF1/fMRI2T1_cmd.txt
  fi
  
  # # convert Hippo subfield to fMRI space

  for ind in ${!HippoSegnames[*]}; do
	hippoSeg_path=$subfolder${HippoSegnames[$ind]}

	hippoSeg_fMRI_path=$fMRISpace_dir/${HippoSegsavenames[$ind]}
	if [[ ! -f $hippoSeg_fMRI_path ]] && [[ -f $hippoSeg_path ]]; then
	  if [[ -f $T22T1_affine ]] && [[ -f $fMRI2T1_affine ]]; then
	    # Original segmentation should be in T2 space
		echo "$hippoSeg_path: Run hippocampal subfield transform to fMRI space."
		mkdir -p $fMRISpace_dir
		transforms="-t [$fMRI2T1_affine,1] -t $T22T1_affine"
		antsApplyTransforms -d 3 -i $hippoSeg_path -r $fMRImean $transforms -u int -n MultiLabel -o $hippoSeg_fMRI_path
	  elif [[ ! -f $T22T1_affine ]] && [[ -f $fMRI2T1_affine ]]; then
	    # Original segmentation should be in T1 space
		echo "$hippoSeg_path: Run hippocampal subfield transform to fMRI space."
		mkdir -p $fMRISpace_dir
		transforms="-t [$fMRI2T1_affine,1]"
		antsApplyTransforms -d 3 -i $hippoSeg_path -r $fMRImean $transforms -u int -n MultiLabel -o $hippoSeg_fMRI_path
	  fi
	fi
  done

  # transfer Freesurfer labelling to fMRI space
  aparc_dir=$subfolder/freesurfer/mri/aparc+aseg.nii.gz
  if [[ -f $aparc_dir ]] && [[ ! -f $fMRISpace_dir/aparc+aseg_fMRI.nii.gz ]] && [[ -f $fMRI2T1_affine ]]; then
    echo "$subfolder: Transform aparc+aseg.nii.gz to fMRI space."
	mkdir -p $fMRISpace_dir
	antsApplyTransforms -d 3 -i $aparc_dir -r $fMRImean -t [$fMRI2T1_affine,1] -u int -n MultiLabel -o $fMRISpace_dir/aparc+aseg_fMRI.nii.gz
  fi
  
  if [ -f $fMRI2T1_affine ]; then
	mkdir -p $subfolder/ants_registration/fMRISpace
	transforms="-t [$fMRI2T1_affine,1]"
	if [ -f $c1T1 ] && [ ! -f $subfolder/ants_registration/fMRISpace/c1T1_fMRI.nii.gz ]; then
	  antsApplyTransforms -d 3 -i $c1T1 -r $fMRImean $transforms -o $subfolder/ants_registration/fMRISpace/c1T1_fMRI.nii.gz
	fi
	if [ -f $c2T1 ] && [ ! -f $subfolder/ants_registration/fMRISpace/c2T1_fMRI.nii.gz ]; then
	  antsApplyTransforms -d 3 -i $c2T1 -r $fMRImean $transforms -o $subfolder/ants_registration/fMRISpace/c2T1_fMRI.nii.gz
	fi
	if [ -f $c3T1 ] && [ ! -f $subfolder/ants_registration/fMRISpace/c3T1_fMRI.nii.gz ]; then
	  antsApplyTransforms -d 3 -i $c3T1 -r $fMRImean $transforms -o $subfolder/ants_registration/fMRISpace/c3T1_fMRI.nii.gz
	fi
  fi
done

R01dir=/mnt/R01_JC/data
SUBFOLDERS=$(find $R01dir/SEWARD_M_* -maxdepth 0 -mindepth 0 -type d ); 
for subfolder in $SUBFOLDERS; do
	antsdir=$subfolder/MRIdata/ants_registration
	fMRISpacedir=$antsdir/fMRISpace
	mkdir -p $fMRISpacedir
	fMRImean=$(find $subfolder/MRIdata/nifti/*bold* -name STmean.nii )
	c1T1=$(find $subfolder/MRIdata/nifti/*MPRAGE* -name c1*.nii)
	c2T1=$(find $subfolder/MRIdata/nifti/*MPRAGE* -name c2*.nii)
	c3T1=$(find $subfolder/MRIdata/nifti/*MPRAGE* -name c3*.nii)
	fMRI2T1_affine=$(find $subfolder/MRIdata/ants_registration/fMRI2T10GenericAffine.mat )
	if [ -f $c1T1 ] && [ -f $fMRI2T1_affine ] && [ -f $fMRImean ] && [ ! -f $fMRISpacedir/c1T1_fMRI.nii.gz ]; then
		echo $subfolder
		transforms="-t [$fMRI2T1_affine,1]"
		antsApplyTransforms -d 3 -i $c1T1 -r $fMRImean $transforms -o $fMRISpacedir/c1T1_fMRI.nii.gz
		antsApplyTransforms -d 3 -i $c2T1 -r $fMRImean $transforms -o $fMRISpacedir/c2T1_fMRI.nii.gz
		antsApplyTransforms -d 3 -i $c3T1 -r $fMRImean $transforms -o $fMRISpacedir/c3T1_fMRI.nii.gz
	fi
done