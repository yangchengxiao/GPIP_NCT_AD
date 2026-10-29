%% Keep the same SUB folders between State1_CN and  State2_AD
clear;clc
path1 = 'F:\Projects\20260515_GPIP_NCT_AD\Results\Mask\State1_CN';
path2 = 'F:\Projects\20260515_GPIP_NCT_AD\Results\Mask\State2_AD';

subjectFolders1 = dir(path1);
subjectFolders1(1:2)=[];
subjects{1} = {subjectFolders1.name};
subjectFolders2 = dir(path2);
subjectFolders2(1:2)=[];
subjects{2} = {subjectFolders2.name};
commonSubjects = intersect(subjects{1},subjects{2});

keep_common(path1, commonSubjects)
keep_common(path2, commonSubjects)

function keep_common(folder_path, common_list)
subs = dir(folder_path);
subs(1:2)=[];
for i = 1:length(subs)
    sub_name = subs(i).name;
    if ~ismember(sub_name, common_list)
        % 删除不在交集里的被试
        delete(fullfile(folder_path, sub_name));
        fprintf('⚠ 已删除 %s 下的 %s\n', folder_path, sub_name);
    end
end
end

%% == FC in Schaefer mask
% ==========================================================
% Extract ROI-wise time series -> FC
% output:
% FC = subject_num × roi_num × roi_num
% ============================================================

clear;
clc;

datapath = ...
    'F:\Projects\20260515_GPIP_NCT_AD\Data\State1_CN\FunImgRCF';

maskpath = ...
    'F:\Projects\20260515_GPIP_NCT_AD\Results\Mask\State1_CN';

fname ='F:\Projects\20260515_GPIP_NCT_AD\Results\FC_CN_Schaefer.mat';

%%% mask
addpath(genpath('F:\Projects\20250430_Indi_Segmentation_Methods\scripts\GPIP'));  % Setting Environment and Parameters
gpip_fmri = fmri_class_ycx() ;

[atlas_lh_surf.vertices, atlas_lh_surf.faces] = read_surf('F:\Mask\Schaefer_2018\FreeSurfer5.3\fsaverage5\surf\lh.inflated');
atlas_lh_surf.faces = round(atlas_lh_surf.faces + 1);
atlas_lh_surf.vertices = atlas_lh_surf.vertices - [100, 0, 0];
[atlas_rh_surf.vertices, atlas_rh_surf.faces] = read_surf('F:\Mask\Schaefer_2018\FreeSurfer5.3\fsaverage5\surf\rh.inflated');
atlas_rh_surf.faces = round(atlas_rh_surf.faces + 1);
mask_path_lh='F:\Mask\Schaefer_2018\FreeSurfer5.3\fsaverage5\label\lh.Schaefer2018_100Parcels_7Networks_order.annot';
[~, atlas_lh_surf.labels, lh_c] = read_annotation(mask_path_lh);
atlas_lh_surf.labels(atlas_lh_surf.labels == lh_c.table(1, 5)) = 0;
mask_path_rh='F:\Mask\Schaefer_2018\FreeSurfer5.3\fsaverage5\label\rh.Schaefer2018_100Parcels_7Networks_order.annot';
[~, atlas_rh_surf.labels, rh_c] = read_annotation(mask_path_rh);
atlas_rh_surf.labels(atlas_rh_surf.labels == rh_c.table(1, 5)) = 0;
atlas_rh_surf.labels = atlas_rh_surf.labels + max(atlas_lh_surf.labels) + 1;
atlas_rh_surf.labels(atlas_rh_surf.labels ==(max(atlas_lh_surf.labels) + 1)) = 0;
atlas_data = merge_surface_lh_rh(atlas_lh_surf, atlas_rh_surf);
gpip_fmri.set_cortical_surface(atlas_data);
atlas_data.original_label = atlas_data.labels ;
atlas_data.labels = gpip_fmri.pre_process_atlas_clusters_labels_nopart(atlas_data.labels, 0);


roi_num = 100;

%%% ===================== Subject list ==========================

sub_dir = dir(fullfile(maskpath, '*.mat'));

sub_num = length(sub_dir);

if sub_num == 0
    error('在 maskpath 中没有找到 .mat 文件');
end

%%% ===================== Loop subjects =========================

FC = cell(sub_num, 1);
State_amplitudes = cell(sub_num, 1);
for isub = 1:sub_num
    % Initialize signals
    % ---------------------------------------------------------

    signals = [];

    %%% --------------------------------------------------------
    % Subject name
    % ---------------------------------------------------------

    mask_name = sub_dir(isub).name;

    [~, sub_name, ~] = fileparts(mask_name);

    fprintf('[%d/%d] Processing %s\n', ...
        isub, sub_num, sub_name);


    % 转成列向量
    sub_mask = atlas_data.labels(:);


    %%% --------------------------------------------------------
    % Load lh / rh surface fMRI
    % ---------------------------------------------------------

    sub_func_dir = fullfile(datapath, sub_name);

    lh_file = fullfile(sub_func_dir, 'lh.fs5.mgh');
    rh_file = fullfile(sub_func_dir, 'rh.fs5.mgh');

    if ~exist(lh_file, 'file')
        warning('找不到: %s', lh_file);
        continue;
    end

    if ~exist(rh_file, 'file')
        warning('找不到: %s', rh_file);
        continue;
    end


    %%% --------------------------------------------------------
    % Read MGH
    %
    % MRIread is provided by FreeSurfer MATLAB toolbox
    % ---------------------------------------------------------

    lh_mri = MRIread(lh_file);
    rh_mri = MRIread(rh_file);

    lh_data = lh_mri.vol;
    rh_data = rh_mri.vol;


    %%% --------------------------------------------------------
    % Convert to vertex × time
    %
    % MRIread commonly gives:
    % 1 × vertex × 1 × time
    % ---------------------------------------------------------

    lh_data = squeeze(lh_data);
    rh_data = squeeze(rh_data);

    % 保证：
    % rows    = vertices
    % columns = time points

    if size(lh_data,1) < size(lh_data,2)
        lh_data = lh_data';
    end

    if size(rh_data,1) < size(rh_data,2)
        rh_data = rh_data';
    end


    %%% --------------------------------------------------------
    % Combine hemispheres
    % ---------------------------------------------------------

    brain_data = [
        lh_data
        rh_data
        ];

    % vertex × time


    %%% --------------------------------------------------------
    % Check dimensions
    % ---------------------------------------------------------

    if length(sub_mask) ~= size(brain_data,1)

        error(['Mask 和 fMRI 顶点数量不一致：\n' ...
            'Subject: %s\n' ...
            'Mask vertices: %d\n' ...
            'fMRI vertices: %d'], ...
            sub_name, ...
            length(sub_mask), ...
            size(brain_data,1));
    end


    %%% --------------------------------------------------------
    % Calculate mean time series within each ROI
    % ---------------------------------------------------------

    for iroi = 1:roi_num

        vertex_idx = (sub_mask == iroi);

        if any(vertex_idx)

            signals(iroi,:) = ...
                mean( ...
                brain_data(vertex_idx,:), ...
                1, ...
                'omitnan');

        else

            warning( ...
                '%s: ROI %d 没有任何 vertex', ...
                sub_name, ...
                iroi);

        end

    end
    % State_amplitudes{isub} = zscore(mean(signals,2));
    State_amplitudes{isub} = mean(signals,2);
    fc =  corr(signals');
    fc(isnan(fc)) = 0;
    fc(logical(eye(100))) = 1;

    FC{isub} =fc;

end

%%% ===================== Save ================================
save( ...
    fname, ...
    'FC', 'State_amplitudes', ...
    '-v7.3');

%% == FC in GPIP mask
% ==========================================================
% Extract ROI-wise time series -> FC
% output:
% FC = subject_num × roi_num × roi_num
% ============================================================

clear;
clc;

datapath = ...
    'F:\Projects\20260515_GPIP_NCT_AD\Data\State2_AD\FunImgRCF';

maskpath = ...
    'F:\Projects\20260515_GPIP_NCT_AD\Results\Mask\State2_AD';

fname ='F:\Projects\20260515_GPIP_NCT_AD\Results\FC_AD_GPIP.mat';

roi_num = 100;

%%% ===================== Subject list ==========================

sub_dir = dir(fullfile(maskpath, '*.mat'));

sub_num = length(sub_dir);

if sub_num == 0
    error('在 maskpath 中没有找到 .mat 文件');
end

%%% ===================== Loop subjects =========================

FC = cell(sub_num, 1);
State_amplitudes = cell(sub_num, 1);
for isub = 1:sub_num
    % Initialize signals
    % ---------------------------------------------------------

    signals = [];

    %%% --------------------------------------------------------
    % Subject name
    % ---------------------------------------------------------

    mask_name = sub_dir(isub).name;

    [~, sub_name, ~] = fileparts(mask_name);

    fprintf('[%d/%d] Processing %s\n', ...
        isub, sub_num, sub_name);


    %%% --------------------------------------------------------
    % Load individualized parcellation mask
    % ---------------------------------------------------------

    mask_file = fullfile(maskpath, mask_name);

    tmp = load(mask_file);

    % 自动获取 .mat 中的变量
    varNames = fieldnames(tmp);

    if length(varNames) ~= 1
        fprintf('Variables in %s:\n', mask_name);
        disp(varNames);

        error(['该 mat 文件中包含多个变量，' ...
            '请明确指定哪个变量是 ROI mask。']);
    end

    sub_mask = tmp.(varNames{1});

    % 转成列向量
    sub_mask = sub_mask(:);


    %%% --------------------------------------------------------
    % Load lh / rh surface fMRI
    % ---------------------------------------------------------

    sub_func_dir = fullfile(datapath, sub_name);

    lh_file = fullfile(sub_func_dir, 'lh.fs5.mgh');
    rh_file = fullfile(sub_func_dir, 'rh.fs5.mgh');

    if ~exist(lh_file, 'file')
        warning('找不到: %s', lh_file);
        continue;
    end

    if ~exist(rh_file, 'file')
        warning('找不到: %s', rh_file);
        continue;
    end


    %%% --------------------------------------------------------
    % Read MGH
    %
    % MRIread is provided by FreeSurfer MATLAB toolbox
    % ---------------------------------------------------------

    lh_mri = MRIread(lh_file);
    rh_mri = MRIread(rh_file);

    lh_data = lh_mri.vol;
    rh_data = rh_mri.vol;


    %%% --------------------------------------------------------
    % Convert to vertex × time
    %
    % MRIread commonly gives:
    % 1 × vertex × 1 × time
    % ---------------------------------------------------------

    lh_data = squeeze(lh_data);
    rh_data = squeeze(rh_data);

    % 保证：
    % rows    = vertices
    % columns = time points

    if size(lh_data,1) < size(lh_data,2)
        lh_data = lh_data';
    end

    if size(rh_data,1) < size(rh_data,2)
        rh_data = rh_data';
    end


    %%% --------------------------------------------------------
    % Combine hemispheres
    % ---------------------------------------------------------

    brain_data = [
        lh_data
        rh_data
        ];

    % vertex × time


    %%% --------------------------------------------------------
    % Check dimensions
    % ---------------------------------------------------------

    if length(sub_mask) ~= size(brain_data,1)

        error(['Mask 和 fMRI 顶点数量不一致：\n' ...
            'Subject: %s\n' ...
            'Mask vertices: %d\n' ...
            'fMRI vertices: %d'], ...
            sub_name, ...
            length(sub_mask), ...
            size(brain_data,1));
    end


    %%% --------------------------------------------------------
    % Calculate mean time series within each ROI
    % ---------------------------------------------------------

    for iroi = 1:roi_num

        vertex_idx = (sub_mask == iroi);

        if any(vertex_idx)

            signals(iroi,:) = ...
                mean( ...
                brain_data(vertex_idx,:), ...
                1, ...
                'omitnan');

        else

            warning( ...
                '%s: ROI %d 没有任何 vertex', ...
                sub_name, ...
                iroi);

        end

    end
    % State_amplitudes{isub} = zscore(mean(signals,2));
    State_amplitudes{isub} = mean(signals,2);
    fc =  corr(signals');
    fc(isnan(fc)) = 0;
    fc(logical(eye(100))) = 1;

    FC{isub} =fc;

end

%%% ===================== Save ================================
save( ...
    fname, ...
    'FC', 'State_amplitudes', ...
    '-v7.3');


%% ============================================================
%  Query Schaefer2018 atlas using an MNI coordinate
%
%  Atlas:
%  Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz
%
%  Target:
%  Left DLPFC MNI = [-38 44 26]
%
%  Requires: SPM12
% ============================================================

clear; clc;

%%% ======================= Path ================================

atlasFile = ...
    'F:\Mask\Schaefer_2018\MNI\Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz';

% Schaefer label/LUT file
labelFile = ...
    'F:\Mask\Schaefer_2018\MNI\Centroid_coordinates\Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv';

% Target MNI coordinate
mni = [-38 44 26];

% Sphere radii
sphereRadius = [3 5];   % mm


%%% ======================= Load atlas ==========================

V = spm_vol(atlasFile);
Y = spm_read_vols(V);

fprintf('\n=============================================\n');
fprintf('Target MNI coordinate: [%g %g %g]\n', ...
    mni(1), mni(2), mni(3));
fprintf('=============================================\n');


%%% ============================================================
% 1. MNI coordinate -> voxel coordinate
% =============================================================

mni_h = [mni 1]';

vox = inv(V.mat) * mni_h;

fprintf('\nVoxel coordinate (continuous):\n');
fprintf('  [%.3f %.3f %.3f]\n', ...
    vox(1), vox(2), vox(3));

% nearest voxel
voxRound = round(vox(1:3));

fprintf('Nearest voxel:\n');
fprintf('  [%d %d %d]\n', ...
    voxRound(1), voxRound(2), voxRound(3));


%%% ================= Check boundary ============================

dim = size(Y);

if any(voxRound < 1) || ...
        voxRound(1) > dim(1) || ...
        voxRound(2) > dim(2) || ...
        voxRound(3) > dim(3)

    error('Target coordinate is outside atlas volume.');

end


%%% ============================================================
% 2. Get ROI ID
% =============================================================

roiID = round(Y( ...
    voxRound(1), ...
    voxRound(2), ...
    voxRound(3)));

fprintf('\n---------------------------------------------\n');
fprintf('Exact voxel lookup\n');
fprintf('---------------------------------------------\n');

fprintf('ROI ID = %d\n', roiID);


%%% ============================================================
% Read Schaefer ROI information from centroid CSV
% ============================================================

T = readtable(labelFile, 'VariableNamingRule', 'preserve');

% ROI ID
ROI_ID = T.("ROI Label");

% ROI Name
ROI_Label = string(T.("ROI Name"));

% Centroid coordinates
ROI_R = T.R;
ROI_A = T.A;
ROI_S = T.S;


%%% ===================== Find label ============================
idx = find(ROI_ID == roiID, 1);

if roiID == 0

    fprintf('ROI ID   = 0\n');
    fprintf('Label    = Background / unlabeled voxel\n');

elseif ~isempty(idx)

    roiLabel = ROI_Label(idx);

    fprintf('ROI ID   = %d\n', roiID);
    fprintf('Label    = %s\n', roiLabel);

    fprintf('Centroid = [%d %d %d]\n', ...
        ROI_R(idx), ...
        ROI_A(idx), ...
        ROI_S(idx));

    % 提取 network
    tokens = split(roiLabel, '_');

    if numel(tokens) >= 3
        network = tokens(3);
    else
        network = "Unknown";
    end

    fprintf('Network  = %s\n', network);

else

    fprintf('ROI ID = %d\n', roiID);
    fprintf('Label not found in centroid CSV.\n');

end


%%% 对应的是ROI 4



%%      NCT
%%% ============================================================
%  NCT analysis
%  AD -> CN state transition
%
%  Based on Fang et al., NeuroImage 2022
% ============================================================

clear;
clc;

addpath('F:\Projects\20260515_GPIP_NCT_AD\scripts\NetworkControl-master');

addpath(genpath( ...
    'F:\Projects\20260515_GPIP_NCT_AD\scripts\Project_SCZ_controllability-main'));

fname = ...
    'F:\Projects\20260515_GPIP_NCT_AD\Results\NCT_OPT_Schaefer.mat';
%%% ============================================================
% Load data
% ============================================================

fc_name = ...
    'F:\Projects\20260515_GPIP_NCT_AD\Results\FC_CN_Schaefer.mat';

load(fc_name);

State_CN = State_amplitudes;
FC_CN    = FC;


fc_name = ...
    'F:\Projects\20260515_GPIP_NCT_AD\Results\FC_AD_Schaefer.mat';

load(fc_name);

State_AD = State_amplitudes;
FC_AD    = FC;


%%% ============================================================
% Basic parameters
% ============================================================

N_nodes = 100;
N_sub   = size(FC_CN,1);


% Optimal control parameters
rho = 0.3;
S   = eye(N_nodes);


% Fang et al.:
T_list = 0.05 : 0.05 : 0.75;


%%% ============================================================
% Initialize
% ============================================================

AveControl_CN   = cell(N_sub,1);
ModalControl_CN = cell(N_sub,1);

AveControl_AD   = cell(N_sub,1);
ModalControl_AD = cell(N_sub,1);


% Treatmeant_Efficacy = struct;
% Simulated_state     = struct;
% Opt_state           = struct;
% Optimal_T           = struct;
% Control_Energy      = struct;
% Target_ROI          = struct;


%%% ============================================================
% Main loop
% ============================================================

for s = 1:N_sub

    fprintf('\n========================================\n');
    fprintf('Subject %d / %d\n',s,N_sub);
    fprintf('========================================\n');


    % ========================================================
    % 1. Network controllability
    % =========================================================

    fc_CN = FC_CN{s};
    fc_AD = FC_AD{s};


    % CN
    AveControl_CN{s}   = ave_control(fc_CN);
    ModalControl_CN{s} = modal_control(fc_CN);


    % AD
    AveControl_AD{s}   = ave_control(fc_AD);
    ModalControl_AD{s} = modal_control(fc_AD);


    % ========================================================
    % Brain states
    %
    % Initial state = AD
    % Target state  = CN
    % =========================================================

    x0 = State_AD{s};
    xf = State_CN{s};


    % Make sure they are column vectors
    x0 = x0(:);
    xf = xf(:);


    % ---------------------------------------------------------
    % Construct B
    % ---------------------------------------------------------

    for iroi = 1:100

        B = zeros(N_nodes);

        diagValues = ...
            5e-4 + 5e-5 .* randn(N_nodes,1);

        diagValues(iroi) = 1;

        B(1:N_nodes+1:end) = diagValues;


        % ---------------------------------------------------------
        % Search different T
        % ---------------------------------------------------------

        corr_T   = nan(length(T_list),1);
        energy_T = nan(length(T_list),1);

        simulated_T = cell(length(T_list),1);
        opt_T       = cell(length(T_list),1);


        for iT = 1:length(T_list)

            T = T_list(iT);


            % -----------------------------------------------------
            % Optimal control
            % -----------------------------------------------------

            [X_opt,U_opt,n_err] = ...
                optim_fun( ...
                fc_AD, ...
                T, ...
                B, ...
                x0, ...
                xf, ...
                rho, ...
                S);


            % Delete t = 0
            X_opt(1,:) = [];
            U_opt(1,:) = [];


            U = zeros( ...
                N_nodes, ...
                1, ...
                size(U_opt,1));


            for it = 1:size(U_opt,1)

                U(iroi,1,it) = ...
                    U_opt(it,iroi);

            end


            % -----------------------------------------------------
            % Open-loop simulation
            % -----------------------------------------------------

            x_sim = ...
                open_loop_control( ...
                fc_AD, ...
                B, ...
                x0, ...
                U);


            x_sim = squeeze(x_sim);


            % x_sim:
            %
            % N_nodes × N_timepoints
            %
            % Compare every simulated state with TARGET state xf
            % -----------------------------------------------------

            state_corr = ...
                corr( ...
                x_sim, ...
                xf);


            % Ignore first five points
            [max_corr,max_idx] = ...
                max(state_corr(6:end));


            max_idx = max_idx + 5;


            simulated_T{iT} = ...
                x_sim(:,max_idx);


            opt_T{iT} = ...
                X_opt(max_idx,:)';


            corr_T(iT) = max_corr;


            % -----------------------------------------------------
            % Control energy
            %
            % E = integral ||B*u(t)||^2 dt
            % -----------------------------------------------------

            Bu = zeros( ...
                N_nodes, ...
                size(U_opt,1));


            for it = 1:size(U_opt,1)

                Bu(:,it) = ...
                    B * U(:,1,it);

            end


            dt = ...
                T / size(U_opt,1);


            energy_T(iT) = ...
                sum( ...
                sum(Bu.^2,1) ...
                ) * dt;


        end


        % ---------------------------------------------------------
        % Choose T giving maximum model fitting
        % ---------------------------------------------------------

        [best_corr,best_T_idx] = ...
            max(corr_T);


        Treatmeant_Efficacy(s,iroi) = ...
            best_corr;


        Optimal_T(s,iroi) = ...
            T_list(best_T_idx);


        Control_Energy(s,iroi) = ...
            energy_T(best_T_idx);


        Simulated_state(s,:,iroi) = ...
            zscore( ...
            simulated_T{best_T_idx})';


        Opt_state(s,:,iroi) = ...
            zscore( ...
            opt_T{best_T_idx})';

    end
        % ========================================================
        % Strategy 1
        %
        % Maximum Average Controllability target
        % +
        % Optimal control parameters
        % =========================================================

        [~, idx_AC] = max(AveControl_AD{s});

        Target_ROI.average_opt(s,1) = idx_AC;


        % ========================================================
        % Strategy 2
        %
        % Maximum Modal Controllability target
        % +
        % Optimal control parameters
        % =========================================================

        [~, idx_MC] = ...
            max(ModalControl_AD{s});


        Target_ROI.modal_opt(s,1) = ...
            idx_MC;



        % ========================================================
        % Strategy 3
        %
        % DLPFC target
        % +
        % Optimal control parameters
        % =========================================================

        idx_dlpfc = 27;


    fprintf( ...
        'Subject %d calculated.\n', ...
        s);

end

%%% ===================== Save ================================
save( ...
    fname, ...
    'Opt_state', 'Simulated_state', 'Control_Energy', 'Optimal_T',  'Treatmeant_Efficacy',...
    'AveControl_CN', 'ModalControl_CN', 'AveControl_AD', 'ModalControl_AD', ...
    'Target_ROI', 'idx_dlpfc', 'State_AD', 'State_CN', ...
    '-v7.3');














