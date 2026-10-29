%% Demographic and Clinical Characteristics
% OASIS-3 CN -> AD Converter

clear; clc;

%%% 1. 路径设置
csvFile = 'F:\Projects\20260515_GPIP_NCT_AD\Data\OASIS3_clinical_by_copied_MR_date.csv';

outFile = 'OASIS3_Demographic_Table.xlsx';

% 可选：独立人口学文件，需包含 OASISID 和 SEX 两列
% SEX: 1=Male, 2=Female（根据实际编码修改）
sexFile = '';

%%% 2. 读取数据
T = readtable(csvFile);

T.OASISID = string(T.OASISID);
T.State = string(T.State);

%%% 3. 分别提取 HC 和 AD
HC = T(T.State == "CN", :);
AD = T(T.State == "AD", :);

% 按 OASISID 进行严格配对
[commonID, ia, ib] = intersect( ...
    HC.OASISID, AD.OASISID, 'stable');

HC = HC(ia, :);
AD = AD(ib, :);

N = length(commonID);

fprintf('Paired subjects: %d\n', N);

%%% 4. 计算 MRI 扫描时年龄
% NearestCDRAge 是临床评估时年龄
% 根据日期差估计 MRI 扫描时年龄

HC.Age = HC.NearestCDRAge + ...
    (HC.MRDay - HC.NearestCDRDay) / 365.25;

AD.Age = AD.NearestCDRAge + ...
    (AD.MRDay - AD.NearestCDRDay) / 365.25;

%%% 5. 匹配性别
HC.Sex = nan(N,1);
AD.Sex = nan(N,1);

if ~isempty(sexFile) && exist(sexFile,'file')

    S = readtable(sexFile);
    S.OASISID = string(S.OASISID);

    for i = 1:N
        idx = find(S.OASISID == commonID(i),1);

        if ~isempty(idx)
            HC.Sex(i) = S.SEX(idx);
            AD.Sex(i) = S.SEX(idx);
        end
    end
end

%%% 6. 提取指标
% 指标名、HC数据、AD数据、统计方法
metrics = {
    'Age (years)', HC.Age, AD.Age, 'ttest';
    'MMSE', HC.MMSE, AD.MMSE, 'ttest';
    'CDR global', HC.CDRTOT, AD.CDRTOT, 'signrank';
    'CDR-SB', HC.CDRSUM, AD.CDRSUM, 'signrank'
};

%%% 7. 初始化人口学表
Variable = {};
HC_Value = {};
AD_Value = {};
P_Value = {};

% 添加样本量
Variable{end+1,1} = 'N';
HC_Value{end+1,1} = num2str(N);
AD_Value{end+1,1} = num2str(N);
P_Value{end+1,1} = '-';

%%% 8. 年龄
for i = 1:size(metrics,1)

    name = metrics{i,1};
    x = double(metrics{i,2}(:));
    y = double(metrics{i,3}(:));
    testMethod = metrics{i,4};

    % 均值和标准差
    HC_str = sprintf('%.2f ± %.2f', ...
        mean(x,'omitnan'), std(x,'omitnan'));

    AD_str = sprintf('%.2f ± %.2f', ...
        mean(y,'omitnan'), std(y,'omitnan'));

    % 配对检验仅保留两组均有效的被试
    valid = isfinite(x) & isfinite(y);

    xValid = x(valid);
    yValid = y(valid);

    p = NaN;

    if numel(xValid) >= 2

        switch testMethod

            case 'ttest'
                if std(xValid-yValid) > 0
                    [~,p] = ttest(xValid,yValid);
                end

            case 'signrank'
                if any(xValid ~= yValid)
                    p = signrank(xValid,yValid);
                else
                    p = 1;
                end
        end
    end

    if isnan(p)
        pStr = 'NA';
    elseif p < 0.001
        pStr = '<0.001';
    else
        pStr = sprintf('%.3f',p);
    end

    Variable{end+1,1} = name;
    HC_Value{end+1,1} = HC_str;
    AD_Value{end+1,1} = AD_str;
    P_Value{end+1,1} = pStr;

end

%%% 9. 添加性别
maleHC = sum(HC.Sex == 1);
femaleHC = sum(HC.Sex == 2);

maleAD = sum(AD.Sex == 1);
femaleAD = sum(AD.Sex == 2);

Variable{end+1,1} = 'Sex (Male/Female)';

if all(isnan(HC.Sex))
    HC_Value{end+1,1} = 'NA';
    AD_Value{end+1,1} = 'NA';
else
    HC_Value{end+1,1} = sprintf('%d/%d',maleHC,femaleHC);
    AD_Value{end+1,1} = sprintf('%d/%d',maleAD,femaleAD);
end

% 同一批被试，不需要独立组性别检验
P_Value{end+1,1} = '-';

%%% 10. 整理表格
DemographicTable = table( ...
    Variable, HC_Value, AD_Value, P_Value, ...
    'VariableNames', ...
    {'Variable','HC','AD_Converter','P_Value'});

%%% 11. 显示并保存
disp(DemographicTable);


%% Rank Difference  -- Average v.s. Model
clear;clc

result_fname = 'F:\Projects\20260515_GPIP_NCT_AD\Results\Statistic_RankDiff_Controllability.mat';


fname = ...
    'F:\Projects\20260515_GPIP_NCT_AD\Results\NCT_OPT_GPIP.mat';
load(fname);

%%% 1. Cell -> matrix
% subjects x ROIs = 28 x 100
AveControl   = cell2mat(AveControl_CN')';
ModalControl = cell2mat(ModalControl_CN')';

[nSub, nROI] = size(AveControl);


%%% ============================================================
% 2. Rank the 100 ROIs within each subject
% ============================================================

AveRank   = nan(nSub, nROI);
ModalRank = nan(nSub, nROI);

for iSub = 1:nSub

    % smallest value = rank 1
    % largest value  = rank nROI
    AveRank(iSub,:)   = tiedrank(AveControl(iSub,:));
    ModalRank(iSub,:) = tiedrank(ModalControl(iSub,:));

end


%%% ============================================================
% 3. Normalize ranks to 0-1
%
% rank = 1    -> 0
% rank = 100  -> 1
% ============================================================

AveRank_norm = (AveRank - 1) ./ (nROI - 1);

ModalRank_norm = (ModalRank - 1) ./ (nROI - 1);


%%% ============================================================
% 4. Calculate Rank Difference
%
% Positive:
% ROI has relatively higher Average controllability rank
%
% Negative:
% ROI has relatively higher Modal controllability rank
% ============================================================

RankDiff = AveRank_norm - ModalRank_norm;

%%% ============================================================
% 5. Statistical comparison for each ROI
% Wilcoxon signed-rank test
% ============================================================

p = nan(nROI,1);
z = nan(nROI,1);

mean_AveRank_norm   = nan(nROI,1);
mean_ModalRank_norm = nan(nROI,1);

meanRankDiff   = nan(nROI,1);
medianRankDiff = nan(nROI,1);

for iROI = 1:nROI

    x = AveRank_norm(:,iROI);
    y = ModalRank_norm(:,iROI);

    diff_tmp = RankDiff(:,iROI);

    % Remove NaN
    valid = ~isnan(x) & ~isnan(y);

    x = x(valid);
    y = y(valid);
    diff_tmp = diff_tmp(valid);

    % Wilcoxon signed-rank test
    [p(iROI), ~, stats] = signrank(x, y);

    % z statistic
    if isfield(stats, 'zval')
        z(iROI) = stats.zval;
    end

    % Mean normalized rank
    mean_AveRank_norm(iROI)   = mean(x);
    mean_ModalRank_norm(iROI) = mean(y);

    std_AveRank_norm(iROI)   = std(x);
    std_ModalRank_norm(iROI) = std(y);

    % Rank difference
    meanRankDiff(iROI)   = mean(diff_tmp);
    medianRankDiff(iROI) = median(diff_tmp);

end


%%% ============================================================
% 6. FDR correction across 100 ROIs
% ============================================================

p_FDR = mafdr(p, 'BHFDR', true);

h_FDR = p_FDR < 0.05;


%%% ============================================================
% 7. Create result table
% ============================================================

ROI = (1:nROI)';

Results = table( ...
    ROI, ...
    mean_AveRank_norm, ...
    mean_ModalRank_norm, ...
    meanRankDiff, ...
    medianRankDiff, ...
    z, ...
    p, ...
    p_FDR, ...
    h_FDR);


%%% ============================================================
% 8. Significant ROIs
% ============================================================

SigResults = Results(Results.h_FDR == 1, :);

fprintf('\nFDR significant ROIs: %d / %d\n', ...
    height(SigResults), nROI);


%%% ============================================================
% 9. Mean RankDiff across subjects
% Used for brain surface visualization
% ============================================================

RankDiff_mean = mean(RankDiff, 1, 'omitnan')';

% 100 x 1
%
% RankDiff_mean > 0:
% higher relative ranking in Average controllability
%
% RankDiff_mean < 0:
% higher relative ranking in Modal controllability


%%% ============================================================
% 10. Optional: only retain FDR-significant RankDiff
% ============================================================

RankDiff_sig = RankDiff_mean;

RankDiff_sig(~h_FDR) = 0;


%%% ============================================================
% 11. Plot RankDiff
% ============================================================

figure;

bar(RankDiff_mean);

hold on;

% Mark significant ROIs
sigROI = find(h_FDR);

scatter(sigROI, ...
        RankDiff_mean(sigROI), ...
        50, ...
        'filled');

yline(0, '--');

xlabel('ROI');
ylabel('Normalized Rank Difference');
title('Average vs Modal Controllability');

box off;

save(result_fname, 'Results');

%% ============================================================
% Rank-based comparison of control energy
%
% Control_Energy          : 28 x 100
% Target_ROI.average_opt  : 28 x 1
% Target_ROI.modal_opt    : 28 x 1
% DLPFC ROI               : 4
%
% Lower rank = lower required control energy
% ============================================================
clear;clc

load('F:\Projects\20260515_GPIP_NCT_AD\Results\NCT_OPT_GPIP.mat')

[nSub, nROI] = size(Control_Energy);

Target_Ave   = Target_ROI.average_opt(:);
Target_Modal = Target_ROI.modal_opt(:);

DLPFC_ROI = 27;


%%% ============================================================
% 1. Rank the control energy of 100 ROIs within each subject
%
% lowest energy  -> rank 1
% highest energy -> rank 100
% ============================================================

EnergyRank = nan(nSub, nROI);

for iSub = 1:nSub

    EnergyRank(iSub,:) = tiedrank(Control_Energy(iSub,:));

end


%%% ============================================================
% 2. Normalize rank to 0-1
%
% 0 = lowest energy
% 1 = highest energy
% ============================================================

EnergyRank_norm = ...
    (EnergyRank - 1) ./ (nROI - 1);

EnergyRank_norm = EnergyRank;


%%% ============================================================
% 3. Extract rank of the three target-selection methods
% ============================================================

subID = (1:nSub)';

% Average controllability target
idx_Ave = sub2ind( ...
    [nSub, nROI], ...
    subID, ...
    Target_Ave);

Rank_Ave = EnergyRank_norm(idx_Ave);


% Modal controllability target
idx_Modal = sub2ind( ...
    [nSub, nROI], ...
    subID, ...
    Target_Modal);

Rank_Modal = EnergyRank_norm(idx_Modal);

% Fixed DLPFC
Rank_DLPFC = EnergyRank_norm(:, DLPFC_ROI);


%%% ============================================================
% 4. Combine three methods
% ============================================================

Rank_All = [ ...
    Rank_Ave, ...
    Rank_Modal, ...
    Rank_DLPFC];

% 28 x 3


%%% ============================================================
% 5. Descriptive statistics
% ============================================================

Method = { ...
    'Average controllability'; ...
    'Modal controllability'; ...
    'DLPFC'};

MeanRank = mean(Rank_All, 1, 'omitnan')';

MedianRank = median(Rank_All, 1, 'omitnan')';

SDRank = std(Rank_All, 0, 1, 'omitnan')';

Q1 = prctile(Rank_All, 25, 1)';
Q3 = prctile(Rank_All, 75, 1)';


Descriptive = table( ...
    Method, ...
    MeanRank, ...
    SDRank, ...
    MedianRank, ...
    Q1, ...
    Q3);

disp('===== Normalized energy rank =====')
disp(Descriptive)


% ============================================================
% 6. Friedman test
% Overall difference among three methods
% ============================================================

[p_Friedman, tbl_Friedman, stats_Friedman] = ...
    friedman(Rank_All, 1, 'off');

fprintf('\nFriedman test: p = %.6g\n', p_Friedman);


%%% ============================================================
% 7. Pairwise Wilcoxon signed-rank tests
% ============================================================

PairName = { ...
    'Average vs Modal'; ...
    'Average vs DLPFC'; ...
    'Modal vs DLPFC'};

pairs = [ ...
    1 2;
    1 3;
    2 3];

p_raw = nan(3,1);
zval = nan(3,1);
effect_r = nan(3,1);

MeanRankDiff = nan(3,1);
MedianRankDiff = nan(3,1);


for i = 1:3

    x = Rank_All(:, pairs(i,1));
    y = Rank_All(:, pairs(i,2));

    valid = ~isnan(x) & ~isnan(y);

    x = x(valid);
    y = y(valid);

    % paired Wilcoxon signed-rank test
    [p_raw(i), ~, stats] = ...
        signrank(x, y, 'method', 'approximate');

    zval(i) = stats.zval;

    % effect size
    effect_r(i) = ...
        stats.zval / sqrt(length(x));

    % first method - second method
    MeanRankDiff(i) = ...
        mean(x - y);

    MedianRankDiff(i) = ...
        median(x - y);

end


%%% ============================================================
% 8. FDR correction
% ============================================================

p_FDR = mafdr(p_raw, 'BHFDR', true);

h_FDR = p_FDR < 0.05;


%%% ============================================================
% 9. Result table
% ============================================================

PairwiseResults = table( ...
    PairName, ...
    MeanRankDiff, ...
    MedianRankDiff, ...
    zval, ...
    effect_r, ...
    p_raw, ...
    p_FDR, ...
    h_FDR);

disp(' ')
disp('===== Pairwise comparison =====')
disp(PairwiseResults)


%% Correlations
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

%% Friedman test
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

%% ============================================================
% Friedman test + Pairwise Wilcoxon signed-rank tests + FDR
% R: N subjects × K conditions (currently K = 6)
%%% ============================================================

% Fisher r-to-z transformation
Z = atanh(R);

% Number of conditions
K = size(Z, 2);

%%% 1. Friedman test
[p_friedman, tbl, stats] = friedman(Z, 1, 'off');

fprintf('Friedman test: p = %.6f\n', p_friedman);

%%% 2. Pairwise Wilcoxon signed-rank tests

% Generate all pairwise comparisons
pairs = nchoosek(1:K, 2);

nPairs = size(pairs, 1);

p = zeros(nPairs, 1);
zval = zeros(nPairs, 1);


for i = 1:nPairs

     % paired Wilcoxon signed-rank test
    [p(i), ~, stats] = ...
        signrank( Z(:, pairs(i,1)), Z(:, pairs(i,2)), 'method', 'approximate');

    zval(i) = stats.zval;

end

%%% 3. FDR correction (Benjamini-Hochberg)

[p_sorted, idx] = sort(p);

m = numel(p);

p_fdr_sorted = min(1, ...
    flipud(cummin(flipud( ...
    p_sorted .* m ./ (1:m)'))));

p_fdr = zeros(m, 1);

p_fdr(idx) = p_fdr_sorted;

%%% 4. Results table

Results = table( ...
    pairs(:,1), ...
    pairs(:,2), ...
    p, ...
    p_fdr, ...
    zval, ...
    p_fdr < 0.05, ...
    'VariableNames', ...
    {'Condition_A', 'Condition_B', ...
     'p', 'p_FDR', 'Z', 'Significant'});

disp(Results);