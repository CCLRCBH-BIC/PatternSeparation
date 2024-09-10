clear
clc
close all
datadir = 'U:\RF1\data';
sublist = dir([datadir,'/10*/Analysis/fMRI_ROI.mat']);sublist = sublist(1:end);
for s = numel(sublist):-1:1
    if numel(strfind(sublist(s).folder,'KMA'))>0
        sublist(s) = [];
    end
end
Sides = {'lh','rh'};sm_key = 'sm0';
ROI_ind = [205,232;
    206,232;
    209,232;
    205,231;
    206,231;
    209,231]; %209: CA4+DG
ROI_name = {'Anterior SUB','Anterior CA1','Anterior CA34DG','Posterior SUB','Posterior CA1','Posterior CA34DG'};
nROI = numel(ROI_name);FD_thresh = 1.5;

%% load demographical data
demo_table = readtable('S:\subject_list.xlsx','Sheet','subject_list');
perf_table = readtable('S:/subject_performance.xlsx','Sheet','contrast_v3');%_combine_lure

demo_info_name = {'age','sex','ImmediateRecall','DelayedRecall'};
demo_info = nan(numel(sublist),1);% age
perf_info = nan(numel(sublist),12);% #6: diff|lure #10: diff|foil
for s = 1:numel(sublist)
    subID = sublist(s).folder(13:17);
    %% find info in demo table
    ind = find(demo_table.subjectID==str2num(subID));
    if numel(ind)==1
        demo_info(s,1) = demo_table{ind,'age'};
        demo_info(s,2) = strcmpi(demo_table{ind,'sex'},'Female');
        demo_info(s,3) = demo_table{ind,'ImmediateRecall'};
        demo_info(s,4) = demo_table{ind,'DelayedRecall'};
    end

    %% find info in performance table
    ind = ZY_findstrINcell(perf_table.sub,sublist(s).folder(8:end-9),1);
    if numel(ind)==1
        perf_info(s,:) = table2array(perf_table(ind,2:13));
    end
end
demo_info(demo_info(:,1)<40,1) = nan;
mean_FramewiseDisp_set = nan(numel(sublist),1);
for side_ind = 1:numel(Sides)
    summary_beta = nan(nROI,2,numel(sublist));
    subROITS = cell(numel(sublist),3); % nsub x nfMRIsession
    subDesign = cell(numel(sublist),3);
    side = Sides{side_ind};
    for s = 1:numel(sublist)
        sub_dir = sublist(s).folder(1:end-9);
        subID = sub_dir(end-13:end-9);%strcmpi(subID,'10192_09292022') || 
        if strcmpi(subID,'10289_05172023') %192: signal dropout, 289: Freesurfer failed for hippocampus
            continue;
        end
        version_tag = 'v08';
        if strcmpi(subID(1:2),'10')
            cohort = 'RF1';TR=1.53;
            t_uselabel = zeros(510,1);
        else
            cohort = 'SEWRAD';TR=1.485;
            t_uselabel = zeros(524,1);
        end
        [X_transform,design_mat_version,contrast,t_useind_all,regressors_name] = RF1_GLMdesign_version(version_tag,cohort);

        t_uselabel(t_useind_all) = 1;

        fMRI_dir = [sub_dir,'/Analysis/fMRI_ROI.mat'];
        if ~exist(fMRI_dir,'file')
            continue;
        end
        load(fMRI_dir,'atlas_struct');

        CA_name = ['CA_FS_fMRI_',side];
        HBT_name = ['HBT_FS_fMRI_',side];
        if ~isfield(atlas_struct,CA_name) || ~isfield(atlas_struct,HBT_name)
            continue;
        end
        X_ensemble = [];
        Y_ensemble = [];
        for sess = 1:3
            ROITS = cell(nROI,1);
            fMRITSsession = ['FMRI',num2str(sess),'ROITS'];
            fMRISession =  ['FMRI',num2str(sess)];
            design_list = dir([sub_dir,'/eprime/*scan*/X_run',num2str(sess),'_canonical_hrf_',design_mat_version,'.mat']);
            atlas_key = CA_name;
            if numel(design_list)~=1 || ~isfield(atlas_struct,fMRISession) || ~isfield(atlas_struct,atlas_key)...
                    || ~isfield(atlas_struct.(atlas_key),fMRITSsession)
                continue;
            end
            design_dir = [design_list(1).folder,'/',design_list(1).name];
            load(design_dir,'X','ref_TR');
            X = X/50;
            subDesign{s,sess} = X;

            t_dim = size(atlas_struct.(fMRISession).WMcomp,1);


            [motion_FD,motion_rms,abs_tran,abs_rot] = ZY_motion_FDandRMS(atlas_struct.(fMRISession).Affine);
            mean_FramewiseDisp_set(s,sess) = mean(motion_FD);
            if mean(motion_FD)>=3
                continue;
            else
                exclude_ind = motion_FD>=FD_thresh;
                exclude_ind = exclude_ind + [0;exclude_ind(1:end-1)] + [exclude_ind(2:end);0];
                exclude_ind = exclude_ind(1:numel(t_uselabel));
                t_useind = find(exclude_ind==0 & t_uselabel>0);
%                 if numel(t_useind) < 100
%                     continue;
%                 end
            end
            poly_regressor = (1:t_dim)';%poly_regressor.^2,poly_regressor.^3,
            nuisance_regressors = [poly_regressor,ZY_motion_regressMotionParam(atlas_struct.(fMRISession).Affine,12),...
                atlas_struct.(fMRISession).WMcomp(:,1:3),atlas_struct.(fMRISession).CSFcomp(:,1:3)];%,atlas_struct.(fMRISession).GSTS];
            CA3_ind = find(atlas_struct.(CA_name).ROI==208); % CA3
            for r = 1:nROI
                CA_ind = find(atlas_struct.(CA_name).ROI==ROI_ind(r,1));
                HBT_ind = find(atlas_struct.(HBT_name).ROI==ROI_ind(r,2));
                [C,ia,ib] = intersect(atlas_struct.(CA_name).ROI_ind{CA_ind},atlas_struct.(HBT_name).ROI_ind{HBT_ind});

                CAROITS = atlas_struct.(atlas_key).(fMRITSsession).(sm_key){CA_ind};
                ROITS{r} = mean(CAROITS(:,ia),2);
            end

            
            ROITS_final = cellfun(@(x) {function_RegressSignalOut(zscore(nuisance_regressors(t_useind,:)),...
                x(t_useind,:))},ROITS);
            ROITS_final = cellfun(@(x) {function_detrend_LPF_motion(x,TR,0,1,0,0.1,4)},ROITS_final);
%             ROITS_final = cellfun(@(x) {zscore(x)},ROITS_final);
            ROITS_final = cellfun(@(x) {100*(x - mean(x))./mean(x)},ROITS_final);
            subROITS{s,sess} = ROITS_final;
            X_ensemble = [X_ensemble;X(t_useind,:)*X_transform];
        end

        Y_ensemble = [];
        for sess = 1:3
            if ~isempty(subROITS{s,sess})
                Y_ensemble = [Y_ensemble;cell2mat(subROITS{s,sess}')];
            end
        end
        [Cor,Beta_q,p_fstat,Const,Res,Rest] = Univaranalysisv2(X_ensemble,...
            zscore(Y_ensemble),[],contrast,0,1);
%             ind = find(p_fstat<=0.05);

%             Beta_q = Beta_q(sum(abs(Beta_q)>3,2)==0,:);
        switch version_tag
            case {'v06','v04'}
                %% for v06
                Beta_q = Beta_q(:,1:end-1)-Beta_q(:,end);
            case 'v08'
                Beta_q = Beta_q(:,1:end-2)-Beta_q(:,end-1);
        end
        Beta = Beta_q;
%             ind = find(p_fstat<1e-5);
%             if numel(ind)==0
%                 Beta = zeros(1,size(Beta_q,2));
%             elseif numel(ind)==1
%                 Beta = Beta_q(ind,:);
%             else
%                 Beta = mean(Beta_q(ind,:));
%             end
% %             Beta = mean(Beta_q);
% %             Beta = mean(Beta_q(p_fstat<0.01,:));
        summary_beta(:,1:size(Beta,2),s) = Beta;
    end

    %% post analysis
    savedir = ['U:\RF1\result\ActivationAnalysis\HippoROI\',version_tag,'/',side];

    var_struct = struct;
    var_struct.pval_sig = 0.05;
    var_struct.savedir = savedir;
    var_struct.ROInames = ROI_name;

    summary_beta_input = summary_beta;
    switch version_tag
        case 'v06'
            %% for v06
            regressors_name_input = regressors_name(1:end-1);
            var_struct.measures_anova = 'LureFAs,LureCRs';%LureFAs
        case 'v08'
            regressors_name_input = regressors_name(1:end-2);
            var_struct.measures_anova = 'LureFAs,LureCRs';
        case 'v04'
            %% for v04
            regressors_name_input = regressors_name(1:end-1);
    end
    var_struct.fig_flag = 2;
%     if side_ind==1
%         temp = squeeze(summary_beta_input(3,:,:))';
%         sublist_0919 = sublist(sum(isnan(temp),2)==0 & sum(temp==0,2)==0);
%         for j = 1:numel(sublist_0919)
%             sublist_0919(j).name = sublist_0919(j).folder(8:end-9);
%         end
%         save E:\RF1\Analysis\v8b1_sublist.mat sublist_0919
%     end
    [meanstd_beta,pval_set] = RF1TaskfMRI_PostAnalysisv2(summary_beta_input,...
        regressors_name_input,var_struct);

    for i = 1:numel(pval_set)
        temp = squeeze(summary_beta_input(i,2:3,:));
        f = figure;
        violinplot(temp(1,:)'-temp(2,:)');
        ylabel('Beta difference');
%         title('LureCR - LureFA');
        xlim([0.5 1.5]);ylim([-0.8 1.2]);
        if pval_set(i)<0.05
            text(1,0.06,['p=',sprintf('%.03f',pval_set(i))],'Units','normalized','HorizontalAlignment','right');
        end
        f.Position = [100 100 300 400];
        export_fig([var_struct.savedir,'/betadiff_',ROI_name{i},'.tif'],'-transparent','-m4');
        close
    end
    if strcmpi(version_tag,'v08')
        beta_LCR_LFA_diff = squeeze(summary_beta_input(:,2,:) - summary_beta_input(:,3,:))';
        %% association with performance accuracy
%         perf_metric = perf_info(:,6)-perf_info(:,2); % similar|lure - similar|old
        perf_metric = perf_info(:,6)-perf_info(:,2); % similar|lure - similar|new as recommended by reviewer
        [corr_perf_beta,p_perf_beta] = ZY_nancorr([perf_metric,beta_LCR_LFA_diff]);
        corr_perf_beta = corr_perf_beta(1,2:end);
        p_perf_beta = p_perf_beta(1,2:end);
        ind = find(p_perf_beta<=0.05);
        for i = 1:numel(ind)
            figure(1)
            [stats,mdl]=ZY_plot_fitwithCI(beta_LCR_LFA_diff(:,ind(i)),100*perf_metric);
            xlabel('Beta difference')
            ylabel('LDI(%)');ylim([0 100]);
            text(0.7,0.9,['R^2=',sprintf('%.02f',mdl.Rsquared.Ordinary),', p=',sprintf('%.03f',p_perf_beta(ind(i)))],'Units','normalized','HorizontalAlignment','center');
            export_fig([var_struct.savedir,'/Scatter_Acc_beta',side,ROI_name{ind(i)},'.tif'],'-transparent','-m4');
            close

            %% cubic spine regression
            figure(2)
            notnan_flag = find((~isnan(beta_LCR_LFA_diff(:,ind(i))) & ~isnan(perf_metric))>0);
            [bhat f sse1 knots]=rcspline(beta_LCR_LFA_diff(notnan_flag,ind(i)),100*perf_metric(notnan_flag),...
                [-0.3 0 0.3],1000,linspace(min(beta_LCR_LFA_diff(:,ind(i)),[],'omitnan'),max(beta_LCR_LFA_diff(:,ind(i)),[],'omitnan'),50),1);
            hold on;
            scatter(beta_LCR_LFA_diff(notnan_flag,ind(i)),100*perf_metric(notnan_flag),40,[0.5 0.5 0.5],'filled');
            hold off
            Rsq_nonlinear = 1 - sse1/sum((perf_metric(notnan_flag)-mean(perf_metric(notnan_flag))).^2)/10000;
            export_fig([var_struct.savedir,'/Scatter_Acc_beta',side,ROI_name{ind(i)},'_nonlinear',sprintf('%.02f',Rsq_nonlinear),...
                '.tif'],'-transparent','-m4');
            close
        end

       %% association with age
       demo_ind = 1;demo_name = 'age';
        [corr_perf_beta,p_perf_beta] = ZY_nancorr([demo_info(:,demo_ind),beta_LCR_LFA_diff]);
        corr_perf_beta = corr_perf_beta(1,2:end);
        p_perf_beta = p_perf_beta(1,2:end);
        ind = find(p_perf_beta<=0.05);
        for i = 1:numel(ind)
            figure(1)
            [stats,mdl]=ZY_plot_fitwithCI(beta_LCR_LFA_diff(:,ind(i)),demo_info(:,demo_ind));
            xlabel('Beta difference')
            ylabel(demo_info_name{demo_ind});%ylim([0 100]);
            text(0.7,0.9,['R^2=',sprintf('%.02f',mdl.Rsquared.Ordinary),', p=',sprintf('%.03f',p_perf_beta(ind(i)))],'Units','normalized','HorizontalAlignment','center');
            export_fig([var_struct.savedir,'/Scatter_',demo_name,'_beta',side,ROI_name{ind(i)},'.tif'],'-transparent','-m4');
            close
        end
    end
    save([var_struct.savedir,'/beta.mat'],'summary_beta');
end