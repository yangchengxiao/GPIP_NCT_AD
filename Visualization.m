%% Fig 1 -- Average/Modal controllibility brain map
clear;clc
fig_path = '/media/sf_F_DRIVE/Projects/20260515_GPIP_NCT_AD/Fig/Fig_1';

load('/media/sf_F_DRIVE/Projects/20260515_GPIP_NCT_AD/Results/Statistic_RankDiff_Controllability.mat');

% ----------- 读取 mask -----------
MASK_PATH = '/media/sf_F_DRIVE/Mask/Schaefer_2018/FreeSurfer5.3/fsaverage5/label';   % mask 文件夹
% 读取 .annot 文件
[vertices_lh, label_lh, colortable_lh] = read_annotation(fullfile(MASK_PATH,'lh.Schaefer2018_100Parcels_7Networks_order.annot'));
[vertices_rh, label_rh, colortable_rh] = read_annotation(fullfile(MASK_PATH,'rh.Schaefer2018_100Parcels_7Networks_order.annot'));

% 初始化权重向量
surf_value_lh = zeros(size(vertices_lh));
surf_value_rh = zeros(size(vertices_rh));

% 按标签赋值（跳过0背景）
values  = Results.mean_AveRank_norm;
values  = Results.mean_ModalRank_norm;

for i = 1:100
    if i <=50
        idx = find(label_lh == colortable_lh.table(i+1,5));
        surf_value_lh(idx) = values(i);

    else
        idx = find(label_rh == colortable_rh.table(i+1-50,5));
        surf_value_rh(idx) = values(i);

    end
end
% 写出 .mgh 文件
cd(fig_path);
save_mgh(surf_value_lh, ['lh.AverageControls.mgh'],eye(4));
save_mgh(surf_value_rh, ['rh.AverageControls.mgh'],eye(4));


Comm{1,1} = strcat('freeview -f $SUBJECTS_DIR/fsaverage5/surf/lh.inflated', ...
    ':annot=/media/sf_F_DRIVE/Mask/Schaefer_2018/FreeSurfer5.3/fsaverage5/label/lh.Schaefer2018_100Parcels_7Networks_order.annot:annot_outline=1',...
    ':overlay=', fig_path, '/lh.AverageControls.mgh',...
    ':overlay_threshold=0.00001,0.8:overlay_color=red', ...
    ' --viewport 3d -cam azimuth 0 elevation 0 roll 0',' ', ...
    ' --ss', {32}, fig_path, '/AverageControls_lh_Lateral.png');

Comm{2,1} = strcat('freeview -f $SUBJECTS_DIR/fsaverage5/surf/rh.inflated', ...
    ':annot=/media/sf_F_DRIVE/Mask/Schaefer_2018/FreeSurfer5.3/fsaverage5/label/rh.Schaefer2018_100Parcels_7Networks_order.annot:annot_outline=1',...
    ':overlay=', fig_path, '/rh.AverageControls.mgh',...
    ':overlay_threshold=0.00001,0.8:overlay_color=red', ...
    ' --viewport 3d -cam azimuth 0 elevation 0 roll 0',' ', ...
    ' --ss', {32}, fig_path, '/AverageControls_rh_Lateral.png');

Comm{3,1} = strcat('freeview -f $SUBJECTS_DIR/fsaverage5/surf/lh.inflated', ...
    ':annot=/media/sf_F_DRIVE/Mask/Schaefer_2018/FreeSurfer5.3/fsaverage5/label/lh.Schaefer2018_100Parcels_7Networks_order.annot:annot_outline=1',...
    ':overlay=', fig_path, '/lh.AverageControls.mgh',...
    ':overlay_threshold=0.00001,0.8:overlay_color=red', ...
    ' --viewport 3d -cam azimuth 180 elevation 0 roll 0',' ', ...
    ' --ss', {32}, fig_path, '/AverageControls_lh_Medial.png');

Comm{4,1} = strcat('freeview -f $SUBJECTS_DIR/fsaverage5/surf/rh.inflated', ...
    ':annot=/media/sf_F_DRIVE/Mask/Schaefer_2018/FreeSurfer5.3/fsaverage5/label/rh.Schaefer2018_100Parcels_7Networks_order.annot:annot_outline=1',...
    ':overlay=', fig_path, '/rh.AverageControls.mgh',...
    ':overlay_threshold=0.00001,0.8:overlay_color=red', ...
    ' --viewport 3d -cam azimuth 180 elevation 0 roll 0',' ', ...
    ' --ss', {32}, fig_path, '/AverageControls_rh_Medial.png');

%% Fig 2 -- the choosen target
clear;clc
fig_path = '/media/sf_F_DRIVE/Projects/20260515_GPIP_NCT_AD/Fig/Fig_2';
mkdir(fig_path);

load('/media/sf_F_DRIVE/Projects/20260515_GPIP_NCT_AD/Results/NCT_OPT_GPIP.mat');

% % ----------- 读取 group mask -----------
% MASK_PATH = '/media/sf_F_DRIVE/Mask/Schaefer_2018/FreeSurfer5.3/fsaverage5/label';   % mask 文件夹
% % 读取 .annot 文件
% [vertices_lh, label_lh, colortable_lh] = read_annotation(fullfile(MASK_PATH,'lh.Schaefer2018_100Parcels_7Networks_order.annot'));
% [vertices_rh, label_rh, colortable_rh] = read_annotation(fullfile(MASK_PATH,'rh.Schaefer2018_100Parcels_7Networks_order.annot'));


% indi
MASK_PATH = '/media/sf_F_DRIVE/Projects/20260515_GPIP_NCT_AD/Results/Mask/State2_AD';   % mask 文件夹
Mask_dir = dir(MASK_PATH);
Mask_dir(1:2) = [];
% 按标签赋值（跳过0背景）
k =1;
for isub = 1:28
    % load mask
    load(fullfile(MASK_PATH, Mask_dir(isub).name));
    
    % 初始化权重向量
    surf_value_lh = zeros(size(vertices_lh));
    surf_value_rh = zeros(size(vertices_rh));

    values  =Target_ROI.average_opt(isub);
    % values  =Target_ROI.modal_opt(isub);

    i = values;
    if i <=50
        %         idx = find(label_lh == colortable_lh.table(i+1,5));
        gpip_cluster_lh = gpip_cluster(1:10242);
        idx = find(gpip_cluster_lh == i);
        surf_value_lh(idx) = 1;
        
    else
        %         idx = find(label_rh == colortable_rh.table(i+1-50,5));
        gpip_cluster_rh = gpip_cluster(10243:end);
        idx = find(gpip_cluster_rh == i);
        surf_value_rh(idx) = 1;
        
    end

    % 写出 .mgh 文件
    mgh_name_lh = fullfile(fig_path, ['lh.Target_ByAverage_sub', num2str(isub),'.mgh']);
    save_mgh(surf_value_lh,mgh_name_lh ,eye(4));
    mgh_name_rh =  fullfile(fig_path, ['rh.Target_ByAverage_sub', num2str(isub),'.mgh']);
    save_mgh(surf_value_rh, mgh_name_rh,eye(4));

    fig_name_lh = fullfile(fig_path, ['lh.Target_ByAverage_sub', num2str(isub),'.png']);
    Comm{k,1} = strcat('freeview -f $SUBJECTS_DIR/fsaverage5/surf/lh.inflated', ...
        ':annot=/media/sf_F_DRIVE/Mask/Schaefer_2018/FreeSurfer5.3/fsaverage5/label/lh.Schaefer2018_100Parcels_7Networks_order.annot:annot_outline=1',...
        ':overlay=', mgh_name_lh, ...
        ':overlay_threshold=0.00001,0.8:overlay_color=red', ...
        ' --viewport 3d -cam azimuth 0 elevation 0 roll 0',' ', ...
        ' --ss', {32}, fig_name_lh);
    k=k+1;

    fig_name_rh = fullfile(fig_path, ['rh.Target_ByAverage_sub', num2str(isub),'.png']);
    Comm{k,1} = strcat('freeview -f $SUBJECTS_DIR/fsaverage5/surf/rh.inflated', ...
        ':annot=/media/sf_F_DRIVE/Mask/Schaefer_2018/FreeSurfer5.3/fsaverage5/label/rh.Schaefer2018_100Parcels_7Networks_order.annot:annot_outline=1',...
        ':overlay=', mgh_name_rh, ...
        ':overlay_threshold=0.00001,0.8:overlay_color=red', ...
        ' --viewport 3d -cam azimuth 180 elevation 0 roll 0',' ', ...
        ' --ss', {32}, fig_name_rh);
    k=k+1;

end

for isub = 1
    % 初始化权重向量
surf_value_lh = zeros(size(vertices_lh));
surf_value_rh = zeros(size(vertices_rh));

%     values  =Target_ROI.average_opt(isub);
    values  =Target_ROI.modal_opt(isub);

    i = 27;
    if i <=50
        idx = find(label_lh == colortable_lh.table(i+1,5));
        surf_value_lh(idx) = 1;

    else
        idx = find(label_rh == colortable_rh.table(i+1-50,5));
        surf_value_rh(idx) = 1;

    end

    % 写出 .mgh 文件
    mgh_name_lh = fullfile(fig_path, ['lh.Target_DLPFC_sub', num2str(isub),'.mgh']);
    save_mgh(surf_value_lh,mgh_name_lh ,eye(4));
    mgh_name_rh =  fullfile(fig_path, ['rh.Target_DLPFC_sub', num2str(isub),'.mgh']);
    save_mgh(surf_value_rh, mgh_name_rh,eye(4));

    fig_name_lh = fullfile(fig_path, ['lh.Target_DLPFC_sub', num2str(isub),'.png']);
    Comm{k,1} = strcat('freeview -f $SUBJECTS_DIR/fsaverage5/surf/lh.inflated', ...
        ':annot=/media/sf_F_DRIVE/Mask/Schaefer_2018/FreeSurfer5.3/fsaverage5/label/lh.Schaefer2018_100Parcels_7Networks_order.annot:annot_outline=1',...
        ':overlay=', mgh_name_lh, ...
        ':overlay_threshold=0.00001,0.8:overlay_color=red', ...
        ' --viewport 3d -cam azimuth 0 elevation 0 roll 0',' ', ...
        ' --ss', {32}, fig_name_lh);
    k=k+1;

    fig_name_rh = fullfile(fig_path, ['rh.Target_DLPFC_sub', num2str(isub),'.png']);
    Comm{k,1} = strcat('freeview -f $SUBJECTS_DIR/fsaverage5/surf/rh.inflated', ...
        ':annot=/media/sf_F_DRIVE/Mask/Schaefer_2018/FreeSurfer5.3/fsaverage5/label/rh.Schaefer2018_100Parcels_7Networks_order.annot:annot_outline=1',...
        ':overlay=', mgh_name_rh, ...
        ':overlay_threshold=0.00001,0.8:overlay_color=red', ...
        ' --viewport 3d -cam azimuth 180 elevation 0 roll 0',' ', ...
        ' --ss', {32}, fig_name_rh);
    k=k+1;

end

%% Fig 3
clear;clc
load('F:\Projects\20260515_GPIP_NCT_AD\Results\NCT_OPT_Schaefer.mat')

for isub = 1:28
    x_0 = zscore(State_AD{isub});
    x_f = zscore(State_CN{isub});

    x_simu = squeeze(Opt_state(isub, 101:end, Target_ROI.modal_opt(isub)));

    Origin_table1{isub,1} = [x_0, x_f, x_simu'];

    R(isub, 1) = Treatmeant_Efficacy(isub, Target_ROI.average_opt(isub));
    R(isub, 2) = Treatmeant_Efficacy(isub, Target_ROI.modal_opt(isub));
    R(isub, 3) = Treatmeant_Efficacy(isub, idx_dlpfc);


end

%%% Friedman test
% Fisher r-to-z transformation
Z = atanh(R);
[p_friedman, tbl, stats] = friedman(Z, 1, 'off');

fprintf('Friedman test: p = %.6f\n', p_friedman);

%%% Pairwise Wilcoxon signed-rank tests
pairs = [1 2; 1 3; 2 3];

p = zeros(3,1);

for i = 1:3
    p(i) = signrank(Z(:,pairs(i,1)), ...
                    Z(:,pairs(i,2)));
end

%%% FDR correction
[p_sorted,idx] = sort(p);
m = numel(p);

p_fdr_sorted = min(1, ...
    flipud(cummin(flipud(p_sorted .* m ./ (1:m)'))));

p_fdr = zeros(m,1);
p_fdr(idx) = p_fdr_sorted;

Results = table(pairs(:,1),pairs(:,2),p,p_fdr, ...
    'VariableNames', ...
    {'Condition_A','Condition_B','p','p_FDR'});

disp(Results);


