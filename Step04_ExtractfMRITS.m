clear
clc
RF1datadir = 'I:\RF1\data';%'\\Research2new\j\RF1\data';
R01datadir = 'Y:\data';
sublist = dir([RF1datadir,'/10*']);%[dir([RF1datadir,'/10*']);dir([R01datadir,'/SEWARD_M*'])];
% T_sub = readtable('E:\RF1\Analysis/subject_list.xlsx');
fMRIatlas_set= ["CA_ashs_fMRI_lh.nii.gz","CA_ashs_fMRI_rh.nii.gz","CA_FS_fMRI_lh.nii.gz","CA_FS_fMRI_rh.nii.gz",...
    "HBT_FS_fMRI_lh.nii.gz", "HBT_FS_fMRI_rh.nii.gz",'CAHBT_FS_fMRI_lh.nii.gz','CAHBT_FS_fMRI_rh.nii.gz',...
    "CA_HippoSeg_fMRI_lh.nii.gz", "CA_HippoSeg_fMRI_rh.nii.gz",...
    "aparc+aseg_fMRI.nii.gz"];
for s = 1:numel(sublist)
    subdir = [sublist(s).folder,'/',sublist(s).name];
    niftidir = [subdir,'/7T_study',sublist(s).name(1:5),'/nifti1'];
    %% load atlas_struct file
    mkdir([subdir,'/Analysis']);
    atlas_struct_dir = [subdir,'/Analysis/fMRI_ROI.mat'];
    atlas_struct = struct;
    if exist(atlas_struct_dir,'file')
        %         load(atlas_struct_dir);
        continue;
    end

    %     if exist(atlas_struct_dir,'file')
    %         load(atlas_struct_dir,'atlas_struct');
    %         if isfield(atlas_struct,'aparc_aseg_fMRI')
    %             atlas_struct = rmfield(atlas_struct,'aparc_aseg_fMRI');
    %         end
    %     else
    %         atlas_struct = struct;
    %     end

    %% load fMRI atlas
    for a = 1:numel(fMRIatlas_set)

        atlas_dir = [niftidir,'/ants_registration/fMRIspace/',fMRIatlas_set{a}];
        atlas_key = fMRIatlas_set{a}(1:end-7);
        atlas_key = strrep(atlas_key,'+','_');
        if exist(atlas_dir,'file') && ~isfield(atlas_struct,atlas_key)
            atlas_nii = load_untouch_nii(atlas_dir);
            if strcmpi(atlas_key,'aparc_aseg_fMRI')
                aparcind = readtable('aparc_aseg_lookuptable.xlsx','ReadVariableNames',false,...
                    'Sheet','final_in_use');
                ROI_label = table2array(aparcind(:,'Var1'));
            else
                ROI_label = unique(atlas_nii.img(atlas_nii.img>0));
            end
            ROI_ind = arrayfun(@(roi) {find(atlas_nii.img(:)==roi)},ROI_label);
            atlas_struct.(atlas_key).ROI = ROI_label;
            atlas_struct.(atlas_key).ROI_ind = ROI_ind;
        end
    end
    atlas_struct_keys = fieldnames(atlas_struct);

    if isempty(atlas_struct_keys)
        continue;
    end
    %% load fMRI data
    c1T1_dir = [niftidir,'/ants_registration/fMRIspace/c1T1_fMRI.nii.gz'];
    c2T1_dir = [niftidir,'/ants_registration/fMRIspace/c2T1_fMRI.nii.gz'];
    c3T1_dir = [niftidir,'/ants_registration/fMRIspace/c3T1_fMRI.nii.gz'];
    if exist(c2T1_dir,'file') && exist(c3T1_dir,'file')
        c1T1_nii = load_untouch_nii(c1T1_dir);
        c2T1_nii = load_untouch_nii(c2T1_dir);
        c3T1_nii = load_untouch_nii(c3T1_dir);
        c2_erode = spm_erode(spm_erode(spm_erode(double(c2T1_nii.img))));
        c3_erode = spm_erode(spm_erode(spm_erode(double(c3T1_nii.img))));
    end
    fMRIsess_set = {'FMRI1','FMRI2','FMRI3','CONN1'};
    for f_id = 1:numel(fMRIsess_set)
        fMRIsession = fMRIsess_set{f_id};%['fMRI_',num2str(sess)];
        fMRI_list = dir([niftidir,'/*',fMRIsession,'*']);
        if numel(fMRI_list) == 0
            continue;
        end
        fMRIdir = [fMRI_list(1).folder,'/',fMRI_list(1).name];
        if numel(dir([fMRIdir,'/uaf*.nii']))>500
            fMRIdirectory = [fMRIdir,'/uaf*.nii'];
        elseif numel(dir([fMRIdir,'/raf*.nii']))>500
            fMRIdirectory = [fMRIdir,'/raf*.nii'];
        elseif numel(dir([fMRIdir,'/uraf*.nii']))>500
            fMRIdirectory = [fMRIdir,'/uraf*.nii'];
        else
            fMRIdirectory = [];
        end

        for a = 1:numel(atlas_struct_keys)
            if isfield(atlas_struct.(atlas_struct_keys{a}),[fMRIsession,'ROITS'])
                skip_flag = 1;
            else
                skip_flag = 0;break;
            end
        end
        if ~isfield(atlas_struct,fMRIsession) || ~isfield(atlas_struct.(fMRIsession),'mask')
            skip_flag = 0;
        end

        if isempty(fMRIdirectory) || skip_flag == 1
            continue;
        end

        fMRIdir
        [fMRIdata_orig,listings] = ZY_fmrimerge(fMRIdirectory);
        first_fmri_nii = load_untouch_nii([listings(1).folder,'/',listings(1).name]);
        res = first_fmri_nii.hdr.dime.pixdim(2:4);
        % multiple fMRI smoothing level
        sm=0;
        %         for sm = [0 2 4]
        sm_key = ['sm',num2str(sm)];
        if sm == 0
            fMRIdata = fMRIdata_orig;
        else
            [~,~,sm_filter] = steerfilter(sm,3,[res,4,4,4]);
            sm_filter = reshape(sm_filter,[1,size(sm_filter)]);
            fMRIdata = imfilter(fMRIdata_orig,sm_filter,'replicate');
        end
        for a = 1:numel(atlas_struct_keys)
            if ~isfield(atlas_struct.(atlas_struct_keys{a}),'ROI_ind')
                continue;
            end
            ROITS = cellfun(@(x) {fMRIdata(:,x)},atlas_struct.(atlas_struct_keys{a}).ROI_ind);
            % only the mean TS for the aparc+aseg is saved.
            if strcmpi(atlas_struct_keys{a},'aparc_aseg_fMRI')
                for r = 1:numel(ROITS)
                    if ~isempty(ROITS{r})
                        mean_ROITS = mean(ROITS{r});
                        std_ROITS = std(ROITS{r});

                        [~,outlier_mean] = ZY_removeoutlier(mean_ROITS');
                        [~,outlier_std] = ZY_removeoutlier(std_ROITS');
                        temp_ROITS = ROITS{r}(:,outlier_mean==0 & outlier_std==0);
                        if size(temp_ROITS,2)>=1
                            ROITS{r} = mean(temp_ROITS,2,'omitnan');
                            %                         [U,S,V] = svd(ROITS{r});
                            %                         ROITS{r} = U(:,1);
                        end
                    end
                end
            end
            atlas_struct.(atlas_struct_keys{a}).([fMRIsession,'ROITS']).(sm_key) = ROITS;
        end
        %         end
        %% extract nuisance regressors for the data.
        motionlist = dir([fMRIdir,'/rp_*af*.txt']);
        affine_parameter = readmatrix([fMRIdir,'/',motionlist(1).name]);
        atlas_struct.(fMRIsession).Affine = affine_parameter;
        if exist(c2T1_dir,'file') && exist(c3T1_dir,'file')
            TS_c2 = fMRIdata_orig(:,c2_erode>=0.5*max(c2T1_nii.img(:)));
            TS_c3 = fMRIdata_orig(:,c3_erode>=0.5*max(c2T1_nii.img(:)));
            [U,S,V] = svd(TS_c2);
            atlas_struct.(fMRIsession).WMcomp = U(:,1:5);
            [U,S,V] = svd(TS_c3);
            atlas_struct.(fMRIsession).CSFcomp = U(:,1:5);

            T1mask = c1T1_nii.img+c2T1_nii.img+c3T1_nii.img;
            T1mask = T1mask>=0.9*max(T1mask(:));
            fMRImask = squeeze(mean(fMRIdata_orig,1,'omitnan'))>=0.05*max(fMRIdata_orig(:));
            atlas_struct.(fMRIsession).GSTS = mean(fMRIdata_orig(:,fMRImask & T1mask),2);
            atlas_struct.(fMRIsession).mask = fMRImask & T1mask;
        end
    end
    save(atlas_struct_dir,'atlas_struct');
end