%% GPIP -- Seg
clear;clc ;

runs = ["State1_CN","State2_AD"];

addpath(genpath('/media/useradmin/databoxone/YCX/matlabtoolbox/GPIP'));  % Setting Environment and Parameters

data_folder_name = fullfile('/media/useradmin/databoxone/YCX/20260515_GPIP_NCT_AD/Data', runs{1}, 'FunImgRCF');
result_folder = fullfile('/media/useradmin/databoxone/YCX/20260515_GPIP_NCT_AD/Results/Mask', runs{1});
mkdir(result_folder);


data_folder = dir(data_folder_name) ;
data_folder(1:2)=[];

if ~exist('number_iterations', 'var')
    number_iterations = 20 ;
end
if ~exist('beta', 'var')
    beta = 1e-8;
end
if ~exist('lambda', 'var')
    lambda = 1e-2;
end
if ~exist('max_batch', 'var')
    max_batch = 5;
end

tmp_folder='/media/useradmin/databoxone/YCX/20260515_GPIP_NCT_AD/TEP_foler1';
mkdir(tmp_folder)
for imask=1
    disp(strcat('Seg mask:_ ',num2str(imask))) ;

    gpip_fmri = fmri_class_ycx() ;
    %%% mask
    [atlas_lh_surf.vertices, atlas_lh_surf.faces] = read_surf('/media/useradmin/databoxone/YCX/Mask/Schaefer 2018/FreeSurfer5.3/fsaverage5/surf/lh.inflated');
    atlas_lh_surf.faces = round(atlas_lh_surf.faces + 1);
    atlas_lh_surf.vertices = atlas_lh_surf.vertices - [100, 0, 0];
    [atlas_rh_surf.vertices, atlas_rh_surf.faces] = read_surf('/media/useradmin/databoxone/YCX/Mask/Schaefer 2018/FreeSurfer5.3/fsaverage5/surf/rh.inflated');
    atlas_rh_surf.faces = round(atlas_rh_surf.faces + 1);
    mask_path_lh=fullfile('/media/useradmin/databoxone/YCX/Mask/Schaefer 2018/FreeSurfer5.3/fsaverage5/label',['lh.Schaefer2018_',num2str(imask*100),'Parcels_7Networks_order.annot']);
    [~, atlas_lh_surf.labels, lh_c] = read_annotation(mask_path_lh);
    atlas_lh_surf.labels(atlas_lh_surf.labels == lh_c.table(1, 5)) = 0;
    mask_path_rh=fullfile('/media/useradmin/databoxone/YCX/Mask/Schaefer 2018/FreeSurfer5.3/fsaverage5/label',['rh.Schaefer2018_',num2str(imask*100),'Parcels_7Networks_order.annot']);
    [~, atlas_rh_surf.labels, rh_c] = read_annotation(mask_path_rh);
    atlas_rh_surf.labels(atlas_rh_surf.labels == rh_c.table(1, 5)) = 0;
    atlas_rh_surf.labels = atlas_rh_surf.labels + max(atlas_lh_surf.labels) + 1;
    atlas_rh_surf.labels(atlas_rh_surf.labels ==(max(atlas_lh_surf.labels) + 1)) = 0;
    atlas_data = merge_surface_lh_rh(atlas_lh_surf, atlas_rh_surf);
    gpip_fmri.set_cortical_surface(atlas_data);
    atlas_data.original_label = atlas_data.labels ;
    atlas_data.labels = gpip_fmri.pre_process_atlas_clusters_labels_nopart(atlas_data.labels, 0);

    %%% load subject data
    atlas_field = 'labels' ;
    k = 1;
    for isub = 1 : length(data_folder)

        sub_name = data_folder(isub).name;


        sub_path=fullfile(data_folder_name, data_folder(isub).name);

        % 左半球功能数据
        fs_dir_lh = dir(fullfile(sub_path, 'lh*fs5*'));

        if isempty(fs_dir_lh)
            warning('被试 %s 未找到左半球 fs5 文件，跳过。', ...
                sub_name);
            continue;
        elseif length(fs_dir_lh) > 1
            warning('被试 %s 找到多个左半球文件，使用第一个：%s', ...
                sub_name, fs_dir_lh(1).name);
        end

        Subnames{k} = sub_name;
        lh_file = fullfile(fs_dir_lh(1).folder, fs_dir_lh(1).name);
        lh_data = load_mgh(lh_file);
        lh_data = reshape(lh_data, [], size(lh_data, 4));

        % 右半球功能数据
        fs_dir_rh = dir(fullfile(sub_path, 'rh*fs5*'));

        if isempty(fs_dir_rh)
            warning('被试 %s 未找到右半球 fs5 文件，跳过。', ...
                sub_name);
            continue;
        elseif length(fs_dir_rh) > 1
            warning('被试 %s 找到多个右半球文件，使用第一个：%s', ...
                sub_name, fs_dir_rh(1).name);
        end

        rh_file = fullfile(fs_dir_rh(1).folder, fs_dir_rh(1).name);
        rh_data = load_mgh(rh_file);
        rh_data = reshape(rh_data, [], size(rh_data, 4));

        % 合并双侧表面数据
        bh_data = [lh_data; rh_data];
        bh_data(isnan(bh_data)) = 0;

        gpip_fmri.set_fmri_time_series(bh_data, k) ;
        gpip_fmri.set_clusters_labels(atlas_data.(atlas_field), k) ;
        if (max_batch > 0) && (max_batch < length(data_folder))
            gpip_fmri.dump_subject_time_series(k,tmp_folder);

            k = k+1;

        end

    end
    disp(['Finish loading ' num2str(gpip_fmri.number_of_subjects) ' subjects data.']) ;

    gpip_fmri.minqi_group_clustering(number_iterations, beta, lambda, max_batch);

    for jsub = 1 : gpip_fmri.number_of_subjects
        i_sub_res_file = fullfile(result_folder,Subnames{jsub}) ;
        gpip_cluster = gpip_fmri.clusters_labels{jsub};
        save(i_sub_res_file, 'gpip_cluster');
    end
end

%% GPIP -- Seg
clear;clc ;

runs = ["State1_CN","State2_AD"];

addpath(genpath('/media/useradmin/databoxone/YCX/matlabtoolbox/GPIP'));  % Setting Environment and Parameters

data_folder_name = fullfile('/media/useradmin/databoxone/YCX/20260515_GPIP_NCT_AD/Data', runs{2}, 'FunImgRCF');
result_folder = fullfile('/media/useradmin/databoxone/YCX/20260515_GPIP_NCT_AD/Results/Mask', runs{2});
mkdir(result_folder);


data_folder = dir(data_folder_name) ;
data_folder(1:2)=[];

if ~exist('number_iterations', 'var')
    number_iterations = 20 ;
end
if ~exist('beta', 'var')
    beta = 1e-8;
end
if ~exist('lambda', 'var')
    lambda = 1e-2;
end
if ~exist('max_batch', 'var')
    max_batch = 5;
end

tmp_folder='/media/useradmin/databoxone/YCX/20260515_GPIP_NCT_AD/TEP_foler2';
mkdir(tmp_folder)
for imask=1
    disp(strcat('Seg mask:_ ',num2str(imask))) ;

    gpip_fmri = fmri_class_ycx() ;
    %%% mask
    [atlas_lh_surf.vertices, atlas_lh_surf.faces] = read_surf('/media/useradmin/databoxone/YCX/Mask/Schaefer 2018/FreeSurfer5.3/fsaverage5/surf/lh.inflated');
    atlas_lh_surf.faces = round(atlas_lh_surf.faces + 1);
    atlas_lh_surf.vertices = atlas_lh_surf.vertices - [100, 0, 0];
    [atlas_rh_surf.vertices, atlas_rh_surf.faces] = read_surf('/media/useradmin/databoxone/YCX/Mask/Schaefer 2018/FreeSurfer5.3/fsaverage5/surf/rh.inflated');
    atlas_rh_surf.faces = round(atlas_rh_surf.faces + 1);
    mask_path_lh=fullfile('/media/useradmin/databoxone/YCX/Mask/Schaefer 2018/FreeSurfer5.3/fsaverage5/label',['lh.Schaefer2018_',num2str(imask*100),'Parcels_7Networks_order.annot']);
    [~, atlas_lh_surf.labels, lh_c] = read_annotation(mask_path_lh);
    atlas_lh_surf.labels(atlas_lh_surf.labels == lh_c.table(1, 5)) = 0;
    mask_path_rh=fullfile('/media/useradmin/databoxone/YCX/Mask/Schaefer 2018/FreeSurfer5.3/fsaverage5/label',['rh.Schaefer2018_',num2str(imask*100),'Parcels_7Networks_order.annot']);
    [~, atlas_rh_surf.labels, rh_c] = read_annotation(mask_path_rh);
    atlas_rh_surf.labels(atlas_rh_surf.labels == rh_c.table(1, 5)) = 0;
    atlas_rh_surf.labels = atlas_rh_surf.labels + max(atlas_lh_surf.labels) + 1;
    atlas_rh_surf.labels(atlas_rh_surf.labels ==(max(atlas_lh_surf.labels) + 1)) = 0;
    atlas_data = merge_surface_lh_rh(atlas_lh_surf, atlas_rh_surf);
    gpip_fmri.set_cortical_surface(atlas_data);
    atlas_data.original_label = atlas_data.labels ;
    atlas_data.labels = gpip_fmri.pre_process_atlas_clusters_labels_nopart(atlas_data.labels, 0);

    %%% load subject data
    atlas_field = 'labels' ;
    k = 1;
    for isub = 1 : length(data_folder)

        sub_name = data_folder(isub).name;


        sub_path=fullfile(data_folder_name, data_folder(isub).name);

        % 左半球功能数据
        fs_dir_lh = dir(fullfile(sub_path, 'lh*fs5*'));

        if isempty(fs_dir_lh)
            warning('被试 %s 未找到左半球 fs5 文件，跳过。', ...
                sub_name);
            continue;
        elseif length(fs_dir_lh) > 1
            warning('被试 %s 找到多个左半球文件，使用第一个：%s', ...
                sub_name, fs_dir_lh(1).name);
        end

                Subnames{k} = sub_name;

        lh_file = fullfile(fs_dir_lh(1).folder, fs_dir_lh(1).name);
        lh_data = load_mgh(lh_file);
        lh_data = reshape(lh_data, [], size(lh_data, 4));

        % 右半球功能数据
        fs_dir_rh = dir(fullfile(sub_path, 'rh*fs5*'));

        if isempty(fs_dir_rh)
            warning('被试 %s 未找到右半球 fs5 文件，跳过。', ...
                sub_name);
            continue;
        elseif length(fs_dir_rh) > 1
            warning('被试 %s 找到多个右半球文件，使用第一个：%s', ...
                sub_name, fs_dir_rh(1).name);
        end

        rh_file = fullfile(fs_dir_rh(1).folder, fs_dir_rh(1).name);
        rh_data = load_mgh(rh_file);
        rh_data = reshape(rh_data, [], size(rh_data, 4));

        % 合并双侧表面数据
        bh_data = [lh_data; rh_data];
        bh_data(isnan(bh_data)) = 0;

        gpip_fmri.set_fmri_time_series(bh_data, k) ;
        gpip_fmri.set_clusters_labels(atlas_data.(atlas_field), k) ;
        if (max_batch > 0) && (max_batch < length(data_folder))
            gpip_fmri.dump_subject_time_series(k,tmp_folder);

            k = k+1;

        end

    end
    disp(['Finish loading ' num2str(gpip_fmri.number_of_subjects) ' subjects data.']) ;

    gpip_fmri.minqi_group_clustering(number_iterations, beta, lambda, max_batch);

    for jsub = 1 : gpip_fmri.number_of_subjects
        i_sub_res_file = fullfile(result_folder,Subnames{jsub}) ;
        gpip_cluster = gpip_fmri.clusters_labels{jsub};
        save(i_sub_res_file, 'gpip_cluster');
    end
end