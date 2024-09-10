%% create CAHBT on task fMRI space
clear
clc
RF1datadir = 'I:\RF1\data';%'\\Research2new\j\RF1\data';
R01datadir = 'Y:\data';
sublist = dir([RF1datadir,'/10*']);%[dir([RF1datadir,'/10*']);dir([R01datadir,'/SEWARD_M*'])];
sides = {'lh','rh'};
%232-head; 231-body; 205-SUB, 206-CA1; 209-CA34DG
ROI_ind = [205,232;
    206,232;
    209,232;
    205,231;
    206,231;
    209,231];
ROI_name = {'Anterior SUB','Anterior CA1','Anterior CA34DG','Posterior SUB','Posterior CA1','Posterior CA34DG'};
for s = 1:numel(sublist)
    for side_ind = 1:numel(sides)
        side = sides{side_ind};
        imgdir = [sublist(s).folder,'/',sublist(s).name,'/7T_study',sublist(s).name(1:5),'/nifti1/ants_registration/fMRISpace'];
        CA_dir = [imgdir,'/CA_FS_fMRI_',side,'.nii.gz'];
        HBT_dir = [imgdir,'/HBT_FS_fMRI_',side,'.nii.gz'];
        CAHBT_dir = [imgdir,'/CAHBT_FS_fMRI_',side,'.nii.gz'];
        if exist(CA_dir,'file') && exist(HBT_dir,'file') && ~exist(CAHBT_dir,'file')
            CA_nii = load_untouch_nii(CA_dir);
            HBT_nii = load_untouch_nii(HBT_dir);
            CAHBT_nii = CA_nii;CAHBT_nii.img = zeros(size(CAHBT_nii.img));
            for r = 1:size(ROI_ind,1)
                CAHBT_nii.img(CA_nii.img==ROI_ind(r,1) & HBT_nii.img==ROI_ind(r,2)) = r;
            end
            save_untouch_nii(CAHBT_nii,CAHBT_dir);
        end
    end
end