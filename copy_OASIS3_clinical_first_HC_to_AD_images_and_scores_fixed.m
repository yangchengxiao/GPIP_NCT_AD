%% Clinical-first OASIS-3 HC-to-AD converter screening and image copying
%
% Workflow:
%   1) First screen OASIS3_UDSb4_cdr.csv for subjects converting from
%      HC/Cognitively normal to AD Dementia.
%   2) Then enter only those subjects' image folders and copy the relevant
%      MR sessions into a new folder B.
%   3) Use the dates encoded in the copied MR session folder names
%      OASxxxxx_MR_d#### to generate the clinical-scale information table.
%
% HC in the user's description corresponds to dx1 = "Cognitively normal".

clear; clc;

%% ============================ User paths ===============================
A_ROOT = 'G:\Datasets\Lifespan\OASIS\Data_OASIS\OASIS3\gaokexin-20250507_024220';
B_ROOT = 'F:\Projects\20260515_GPIP_NCT_AD\Data';

cdrFile       = 'F:\Projects\20260515_GPIP_NCT_AD\Data\OASIS3_UDSb4_cdr.csv';
legacyFile    = 'F:\Projects\20260515_GPIP_NCT_AD\Data\clinical features-纵向.xlsx';
cognitiveFile = 'F:\Projects\20260515_GPIP_NCT_AD\Data\OASIS3_UDSc1_cognitive_assessments.csv';

outputWorkbook = fullfile(B_ROOT, ...
    'OASIS3_HC_to_AD_clinical_first_image_and_clinical_table.xlsx');

%% ============================ Parameters ===============================
% "ever_HC_before_AD": keeps HC -> MCI -> AD trajectories.
% "direct_HC_to_AD": the most recent diagnosis before first AD must be HC.
converterMode = "ever_HC_before_AD";

% If no MR session is classified as CN/AD by nearest CDR, use dates around
% the estimated conversion midpoint as fallback.
allowBoundaryFallback = true;

% Classify MR sessions by nearest CDR diagnosis. Inf keeps all MR sessions
% and records the time gap. Use 365 or 730 for stricter matching.
maxDiagnosisMatchDays = Inf;

% Clinical-scale matching is based on copied MR session day.
% These do not affect subject inclusion; unmatched values become NaN.
maxLegacyScaleMatchDays = 365;
maxCognitiveScaleMatchDays = 365;

requireCNBeforeAD = true;
requireMatchingAnatFuncIndex = true;

overwriteExisting = false;
decompressNiiGz = true;
dryRun = false;        % true = no image copy; only output tables/logs

t1PositiveKeywords   = ["t1w","mprage","mpr","t1"];
t1NegativeKeywords   = ["mask","brainmask","seg","aseg","aparc","label", ...
                        "warp","fieldmap","magnitude","phase","t2","flair"];

funcPositiveKeywords = ["rest","rsfmri","rfmri","bold","func"];
funcNegativeKeywords = ["sbref","mean","mask","confound","motion", ...
                        "fieldmap","magnitude","phase"];

%% ============================ Read tables ==============================
assert(isfolder(A_ROOT), 'A_ROOT does not exist: %s', A_ROOT);
ensureFolder(B_ROOT);

assert(isfile(cdrFile), 'CDR table not found: %s', cdrFile);
assert(isfile(legacyFile), 'Legacy clinical table not found: %s', legacyFile);
assert(isfile(cognitiveFile), 'Cognitive table not found: %s', cognitiveFile);

CDR = readtable(cdrFile, 'VariableNamingRule','preserve', 'TextType','string');
LEGACY = readtable(legacyFile, 'VariableNamingRule','preserve', 'TextType','string');
COG = readtable(cognitiveFile, 'VariableNamingRule','preserve', 'TextType','string');

validateColumns(CDR, ["OASISID","OASIS_session_label", ...
    "days_to_visit","age at visit","dx1"]);
validateColumns(LEGACY, ["OASISID","OASIS_session_label","age at visit"]);
validateColumns(COG, ["OASISID","OASIS_session_label", ...
    "days_to_visit","age at visit"]);

CDR_ID   = strtrim(string(CDR.("OASISID")));
CDR_DAY  = numericVector(CDR.("days_to_visit"));
CDR_DX   = strtrim(string(CDR.("dx1")));
CDR_SESS = strtrim(string(CDR.("OASIS_session_label")));

LEGACY_ID   = strtrim(string(LEGACY.("OASISID")));
LEGACY_DAY  = parseDayFromSessionLabel(LEGACY.("OASIS_session_label"));
LEGACY_SESS = strtrim(string(LEGACY.("OASIS_session_label")));

COG_ID   = strtrim(string(COG.("OASISID")));
COG_DAY  = numericVector(COG.("days_to_visit"));
COG_SESS = strtrim(string(COG.("OASIS_session_label")));

%% ====================== Step 1: screen converters ======================
CONVERTERS = screenHCtoADConverters(CDR, converterMode);

fprintf('\n================ Clinical-first screening ================\n');
fprintf('Converter mode: %s\n', converterMode);
fprintf('HC-to-AD converters in CDR table: %d\n\n', height(CONVERTERS));

%% ====================== Step 2: scan/copy MR data ======================
allMRRows = cell(0, 25);
selectedRows = cell(0, 36);
excludedRows = cell(0, 9);
clinicalRows = cell(0, 0);  % Initialized after clinicalVariableNames is defined
copyRows = cell(0, 15);

[legacyOutputNames, legacySourceFields, legacyMin, legacyMax] = legacyScoreSpec();
[cogOutputNames, cogSourceFields, cogMin, cogMax, cogDescriptions] = cognitiveScoreSpec();

baseClinicalNames = [ ...
    "OASISID", "State", "ConverterMode", ...
    "FirstHCSession", "FirstHCDay", "LastHCSession", "LastHCDay", ...
    "FirstADSession", "FirstADDay", "EstimatedConversionDayMidpoint", ...
    "MRSession", "MRDay", "MRSelectionRule", "MRSeriesPairIndex", ...
    "T1SourceFile", "FuncSourceFile", ...
    "DiagnosisFromNearestCDR", "NearestCDRSession", "NearestCDRDay", ...
    "NearestCDRAge", "MRtoCDRGapDays", "DiagnosisMatchQuality", ...
    "MMSE", "CDRSUM", "CDRTOT", ...
    "LegacyScaleSession", "LegacyScaleDay", "LegacyScaleAge", ...
    "MRtoLegacyScaleGapDays", ...
    "CognitiveScaleSession", "CognitiveScaleDay", "CognitiveScaleAge", ...
    "MRtoCognitiveScaleGapDays", ...
    "TargetT1File", "TargetFuncFile", "CopyStatus"];

clinicalVariableNames = [baseClinicalNames, legacyOutputNames, cogOutputNames];
clinicalRows = cell(0, numel(clinicalVariableNames));

for i = 1:height(CONVERTERS)
    subID = CONVERTERS.OASISID(i);
    subjectDir = fullfile(A_ROOT, char(subID));

    if ~isfolder(subjectDir)
        excludedRows(end+1,:) = { ... %#ok<SAGROW>
            subID, CONVERTERS.FirstHCDay(i), CONVERTERS.LastHCDay(i), ...
            CONVERTERS.FirstADDay(i), 0, 0, 0, 0, ...
            "No source subject image folder"};
        continue;
    end

    MRS = listSubjectMRSessions( ...
        subjectDir, subID, requireMatchingAnatFuncIndex, ...
        t1PositiveKeywords, t1NegativeKeywords, ...
        funcPositiveKeywords, funcNegativeKeywords);

    if isempty(MRS)
        excludedRows(end+1,:) = { ... %#ok<SAGROW>
            subID, CONVERTERS.FirstHCDay(i), CONVERTERS.LastHCDay(i), ...
            CONVERTERS.FirstADDay(i), 0, 0, 0, 0, ...
            "No MR session folder found"};
        continue;
    end

    MRS = classifyMRSessionsByCDR(MRS, CDR, maxDiagnosisMatchDays);

    for r = 1:height(MRS)
        allMRRows(end+1,:) = { ... %#ok<SAGROW>
            subID, MRS.MRSession(r), MRS.MRDay(r), MRS.CompleteSession(r), ...
            MRS.SeriesPairIndex(r), MRS.MRState(r), MRS.MRSelectionBasis(r), ...
            MRS.NearestDiagnosis(r), MRS.NearestCDRSession(r), ...
            MRS.NearestCDRDay(r), MRS.NearestCDRAge(r), ...
            MRS.MRtoCDRGapDays(r), MRS.DiagnosisMatchQuality(r), ...
            MRS.AnatFolder(r), MRS.FuncFolder(r), MRS.T1File(r), ...
            MRS.FuncFile(r), CONVERTERS.FirstHCSession(i), ...
            CONVERTERS.FirstHCDay(i), CONVERTERS.LastHCSession(i), ...
            CONVERTERS.LastHCDay(i), CONVERTERS.FirstADSession(i), ...
            CONVERTERS.FirstADDay(i), ...
            CONVERTERS.EstimatedConversionDayMidpoint(i), ...
            CONVERTERS.TransitionWindowDays(i)};
    end

    [cnRow, adRow, cnRule, adRule] = selectMRPairForConverter( ...
        MRS, CONVERTERS.LastHCDay(i), CONVERTERS.FirstADDay(i), ...
        CONVERTERS.EstimatedConversionDayMidpoint(i), ...
        allowBoundaryFallback, requireCNBeforeAD);

    if isnan(cnRow) || isnan(adRow)
        excludedRows(end+1,:) = { ... %#ok<SAGROW>
            subID, CONVERTERS.FirstHCDay(i), CONVERTERS.LastHCDay(i), ...
            CONVERTERS.FirstADDay(i), height(MRS), ...
            sum(MRS.CompleteSession), sum(MRS.MRState=="CN"), ...
            sum(MRS.MRState=="AD"), ...
            "No complete chronological CN/AD MR pair after converter screening"};
        continue;
    end

    if requireCNBeforeAD && MRS.MRDay(cnRow) >= MRS.MRDay(adRow)
        excludedRows(end+1,:) = { ... %#ok<SAGROW>
            subID, CONVERTERS.FirstHCDay(i), CONVERTERS.LastHCDay(i), ...
            CONVERTERS.FirstADDay(i), height(MRS), ...
            sum(MRS.CompleteSession), sum(MRS.MRState=="CN"), ...
            sum(MRS.MRState=="AD"), ...
            "Selected CN MR day is not earlier than selected AD MR day"};
        continue;
    end

    [cnTargetT1, cnTargetFunc, cnCopyStatus, cnCopyNote] = copySelectedMR( ...
        subID, "CN", "State1_CN", MRS(cnRow,:), B_ROOT, ...
        overwriteExisting, decompressNiiGz, dryRun);

    [adTargetT1, adTargetFunc, adCopyStatus, adCopyNote] = copySelectedMR( ...
        subID, "AD", "State2_AD", MRS(adRow,:), B_ROOT, ...
        overwriteExisting, decompressNiiGz, dryRun);

    selectedRows(end+1,:) = { ... %#ok<SAGROW>
        subID, CONVERTERS.FirstHCSession(i), CONVERTERS.FirstHCDay(i), ...
        CONVERTERS.LastHCSession(i), CONVERTERS.LastHCDay(i), ...
        CONVERTERS.FirstADSession(i), CONVERTERS.FirstADDay(i), ...
        CONVERTERS.EstimatedConversionDayMidpoint(i), ...
        MRS.MRSession(cnRow), MRS.MRDay(cnRow), cnRule, ...
        MRS.MRtoCDRGapDays(cnRow), MRS.SeriesPairIndex(cnRow), ...
        MRS.T1File(cnRow), MRS.FuncFile(cnRow), ...
        cnTargetT1, cnTargetFunc, cnCopyStatus, cnCopyNote, ...
        MRS.MRSession(adRow), MRS.MRDay(adRow), adRule, ...
        MRS.MRtoCDRGapDays(adRow), MRS.SeriesPairIndex(adRow), ...
        MRS.T1File(adRow), MRS.FuncFile(adRow), ...
        adTargetT1, adTargetFunc, adCopyStatus, adCopyNote, ...
        MRS.MRDay(adRow)-MRS.MRDay(cnRow), ...
        (MRS.MRDay(adRow)-MRS.MRDay(cnRow))/365.25, ...
        height(MRS), sum(MRS.CompleteSession), ...
        sum(MRS.MRState=="CN"), sum(MRS.MRState=="AD")};

    copyRows(end+1,:) = { ... %#ok<SAGROW>
        subID, "CN", MRS.MRSession(cnRow), MRS.MRDay(cnRow), ...
        MRS.SeriesPairIndex(cnRow), MRS.T1File(cnRow), ...
        MRS.FuncFile(cnRow), cnTargetT1, cnTargetFunc, ...
        cnCopyStatus, cnCopyNote, cnRule, MRS.NearestDiagnosis(cnRow), ...
        MRS.MRtoCDRGapDays(cnRow), MRS.DiagnosisMatchQuality(cnRow)};

    copyRows(end+1,:) = { ... %#ok<SAGROW>
        subID, "AD", MRS.MRSession(adRow), MRS.MRDay(adRow), ...
        MRS.SeriesPairIndex(adRow), MRS.T1File(adRow), ...
        MRS.FuncFile(adRow), adTargetT1, adTargetFunc, ...
        adCopyStatus, adCopyNote, adRule, MRS.NearestDiagnosis(adRow), ...
        MRS.MRtoCDRGapDays(adRow), MRS.DiagnosisMatchQuality(adRow)};

    selectedStateRows = [cnRow, adRow];
    selectedStateNames = ["CN","AD"];
    selectedRules = [cnRule, adRule];
    selectedTargetsT1 = [cnTargetT1, adTargetT1];
    selectedTargetsFunc = [cnTargetFunc, adTargetFunc];
    selectedCopyStatus = [cnCopyStatus, adCopyStatus];

    for s = 1:2
        mrRow = selectedStateRows(s);
        mrDay = MRS.MRDay(mrRow);
        nearestCDRIdx = MRS.CDRRowIndex(mrRow);

        [legacyIdx, legacyGap] = nearestVisitIndex( ...
            LEGACY_ID, LEGACY_DAY, subID, mrDay, isfinite(LEGACY_DAY));
        if ~isnan(legacyIdx) && legacyGap > maxLegacyScaleMatchDays
            legacyIdx = NaN; legacyGap = NaN;
        end

        [cogIdx, cogGap] = nearestVisitIndex( ...
            COG_ID, COG_DAY, subID, mrDay, isfinite(COG_DAY));
        if ~isnan(cogIdx) && cogGap > maxCognitiveScaleMatchDays
            cogIdx = NaN; cogGap = NaN;
        end

        row = cell(1, numel(clinicalVariableNames));
        row{1}  = subID;
        row{2}  = selectedStateNames(s);
        row{3}  = converterMode;
        row{4}  = CONVERTERS.FirstHCSession(i);
        row{5}  = CONVERTERS.FirstHCDay(i);
        row{6}  = CONVERTERS.LastHCSession(i);
        row{7}  = CONVERTERS.LastHCDay(i);
        row{8}  = CONVERTERS.FirstADSession(i);
        row{9}  = CONVERTERS.FirstADDay(i);
        row{10} = CONVERTERS.EstimatedConversionDayMidpoint(i);
        row{11} = MRS.MRSession(mrRow);
        row{12} = MRS.MRDay(mrRow);
        row{13} = selectedRules(s);
        row{14} = MRS.SeriesPairIndex(mrRow);
        row{15} = MRS.T1File(mrRow);
        row{16} = MRS.FuncFile(mrRow);
        row{17} = MRS.NearestDiagnosis(mrRow);
        row{18} = MRS.NearestCDRSession(mrRow);
        row{19} = MRS.NearestCDRDay(mrRow);
        row{20} = MRS.NearestCDRAge(mrRow);
        row{21} = MRS.MRtoCDRGapDays(mrRow);
        row{22} = MRS.DiagnosisMatchQuality(mrRow);

        if isnan(nearestCDRIdx)
            row{23} = NaN; row{24} = NaN; row{25} = NaN;
        else
            row{23} = cleanRange(getNumericValue(CDR,nearestCDRIdx,"MMSE"),0,30);
            row{24} = cleanRange(getNumericValue(CDR,nearestCDRIdx,"CDRSUM"),0,18);
            row{25} = cleanRange(getNumericValue(CDR,nearestCDRIdx,"CDRTOT"),0,5);
        end

        if isnan(legacyIdx)
            row{26} = ""; row{27} = NaN; row{28} = NaN; row{29} = NaN;
        else
            row{26} = LEGACY_SESS(legacyIdx);
            row{27} = LEGACY_DAY(legacyIdx);
            row{28} = getNumericValue(LEGACY, legacyIdx, "age at visit");
            row{29} = legacyGap;
        end

        if isnan(cogIdx)
            row{30} = ""; row{31} = NaN; row{32} = NaN; row{33} = NaN;
        else
            row{30} = COG_SESS(cogIdx);
            row{31} = COG_DAY(cogIdx);
            row{32} = getNumericValue(COG, cogIdx, "age at visit");
            row{33} = cogGap;
        end

        row{34} = selectedTargetsT1(s);
        row{35} = selectedTargetsFunc(s);
        row{36} = selectedCopyStatus(s);

        legacyStart = numel(baseClinicalNames) + 1;
        for k = 1:numel(legacyOutputNames)
            if isnan(legacyIdx)
                value = NaN;
            else
                value = getNumericValue(LEGACY, legacyIdx, legacySourceFields(k));
                value = cleanRange(value, legacyMin(k), legacyMax(k));
            end
            row{legacyStart+k-1} = value;
        end

        cogStart = numel(baseClinicalNames) + numel(legacyOutputNames) + 1;
        cogValues = extractCognitiveTotals( ...
            COG, cogIdx, cogOutputNames, cogSourceFields, cogMin, cogMax);
        for k = 1:numel(cogOutputNames)
            row{cogStart+k-1} = cogValues(k);
        end

        clinicalRows(end+1,:) = row; %#ok<SAGROW>
    end
end

%% =========================== Make output tables ========================
allMRVariableNames = { ...
    'OASISID','MRSession','MRDay','CompleteSession','SeriesPairIndex', ...
    'MRState','MRSelectionBasis','NearestDiagnosis','NearestCDRSession', ...
    'NearestCDRDay','NearestCDRAge','MRtoCDRGapDays', ...
    'DiagnosisMatchQuality','AnatFolder','FuncFolder','T1File','FuncFile', ...
    'FirstHCSession','FirstHCDay','LastHCSession','LastHCDay', ...
    'FirstADSession','FirstADDay','EstimatedConversionDayMidpoint', ...
    'TransitionWindowDays'};

selectedVariableNames = { ...
    'OASISID','FirstHCSession','FirstHCDay','LastHCSession','LastHCDay', ...
    'FirstADSession','FirstADDay','EstimatedConversionDayMidpoint', ...
    'CN_MRSession','CN_MRDay','CN_MRSelectionRule','CN_MRtoCDRGapDays', ...
    'CN_SeriesPairIndex','CN_T1Source','CN_FuncSource', ...
    'CN_TargetT1','CN_TargetFunc','CN_CopyStatus','CN_CopyNote', ...
    'AD_MRSession','AD_MRDay','AD_MRSelectionRule','AD_MRtoCDRGapDays', ...
    'AD_SeriesPairIndex','AD_T1Source','AD_FuncSource', ...
    'AD_TargetT1','AD_TargetFunc','AD_CopyStatus','AD_CopyNote', ...
    'CNtoAD_MRIntervalDays','CNtoAD_MRIntervalYears', ...
    'TotalMRSessions','CompleteMRSessions','CompleteCNMRSessions', ...
    'CompleteADMRSessions'};

excludedVariableNames = { ...
    'OASISID','FirstHCDay','LastHCDay','FirstADDay', ...
    'TotalMRSessions','CompleteMRSessions','CompleteCNMRSessions', ...
    'CompleteADMRSessions','ExclusionReason'};

copyVariableNames = { ...
    'OASISID','State','MRSession','MRDay','SeriesPairIndex', ...
    'SourceT1','SourceFunc','TargetT1','TargetFunc','CopyStatus', ...
    'CopyNote','MRSelectionRule','NearestDiagnosis','MRtoCDRGapDays', ...
    'DiagnosisMatchQuality'};

ALL_MR = cell2tableOrEmpty(allMRRows, allMRVariableNames);
SELECTED = cell2tableOrEmpty(selectedRows, selectedVariableNames);
EXCLUDED = cell2tableOrEmpty(excludedRows, excludedVariableNames);
CLINICAL = cell2tableOrEmpty(clinicalRows, cellstr(clinicalVariableNames));
COPYLOG = cell2tableOrEmpty(copyRows, copyVariableNames);

DICT = buildScoreDictionary(legacyOutputNames, legacySourceFields, ...
    cogOutputNames, cogSourceFields, cogDescriptions);

README = table( ...
    ["Converter screening source"; "Converter mode"; "HC label in CDR dx1"; ...
     "AD label in CDR dx1"; "Converters found in CDR"; ...
     "Converters with selected CN/AD MR pair"; "Excluded converters"; ...
     "MR diagnosis classification"; "Clinical table date basis"; ...
     "Legacy scale match window"; "Cognitive scale match window"; "Dry run"], ...
    [string(cdrFile); converterMode; "Cognitively normal"; "AD Dementia"; ...
     string(height(CONVERTERS)); string(height(SELECTED)); string(height(EXCLUDED)); ...
     "Nearest valid CDR visit to each MR folder day"; ...
     "Copied MR folder day, e.g., OAS30001_MR_d0129 -> day 129"; ...
     string(maxLegacyScaleMatchDays) + " days"; ...
     string(maxCognitiveScaleMatchDays) + " days"; string(dryRun)], ...
    'VariableNames', {'Item','Value'});

%% ============================= Save outputs ============================
if isfile(outputWorkbook)
    delete(outputWorkbook);
end

writetable(README, outputWorkbook, 'Sheet','ReadMe');
writetable(CONVERTERS, outputWorkbook, 'Sheet','CDR_HC_to_AD_Converters');
writetable(SELECTED, outputWorkbook, 'Sheet','SelectedCopiedMR');
writetable(CLINICAL, outputWorkbook, 'Sheet','ClinicalByCopiedMRDate');
writetable(ALL_MR, outputWorkbook, 'Sheet','AllConverterMRSessions');
writetable(EXCLUDED, outputWorkbook, 'Sheet','ExcludedConverters');
writetable(COPYLOG, outputWorkbook, 'Sheet','CopyLog');
writetable(DICT, outputWorkbook, 'Sheet','ScoreDictionary');

writetable(CONVERTERS, fullfile(B_ROOT,'OASIS3_CDR_HC_to_AD_converters.csv'));
writetable(SELECTED, fullfile(B_ROOT,'OASIS3_selected_copied_MR_pairs.csv'));
writetable(CLINICAL, fullfile(B_ROOT,'OASIS3_clinical_by_copied_MR_date.csv'));
writetable(ALL_MR, fullfile(B_ROOT,'OASIS3_all_converter_MR_sessions.csv'));
writetable(EXCLUDED, fullfile(B_ROOT,'OASIS3_excluded_converters.csv'));

fprintf('\n=========================== Finished ===========================\n');
fprintf('Converters screened from CDR: %d\n', height(CONVERTERS));
fprintf('Converters copied/selected with CN and AD MR: %d\n', height(SELECTED));
fprintf('Clinical rows by copied MR date: %d\n', height(CLINICAL));
fprintf('Excluded converters: %d\n', height(EXCLUDED));
fprintf('Output workbook:\n%s\n', outputWorkbook);

if dryRun
    fprintf('\nDry run is ON: no images were copied.\n');
end

%% ============================ Local functions ==========================

function validateColumns(T, requiredNames)
    available = string(T.Properties.VariableNames);
    missing = requiredNames(~ismember(requiredNames, available));
    assert(isempty(missing), ...
        'Missing required table column(s): %s', strjoin(missing, ', '));
end

function values = numericVector(inputValues)
    if isnumeric(inputValues) || islogical(inputValues)
        values = double(inputValues);
    else
        values = str2double(strtrim(string(inputValues)));
    end
    values(~isfinite(values)) = NaN;
end

function dayValues = parseDayFromSessionLabel(labels)
    labels = string(labels);
    dayValues = nan(size(labels));
    for i = 1:numel(labels)
        token = regexp(char(labels(i)), '_d(-?\d+)$', 'tokens', 'once');
        if ~isempty(token)
            dayValues(i) = str2double(token{1});
        end
    end
end

function CONVERTERS = screenHCtoADConverters(CDR, converterMode)
    ids = strtrim(string(CDR.("OASISID")));
    days = numericVector(CDR.("days_to_visit"));
    dx = strtrim(string(CDR.("dx1")));
    sess = strtrim(string(CDR.("OASIS_session_label")));
    ages = numericVector(CDR.("age at visit"));

    uniqueIDs = unique(ids(ids~=""), 'stable');
    rows = cell(0, 16);

    for i = 1:numel(uniqueIDs)
        subID = uniqueIDs(i);
        idx = find(ids == subID & isfinite(days) & dx ~= "" & dx ~= ".");
        if isempty(idx), continue; end

        [~, order] = sort(days(idx));
        idx = idx(order);
        subDx = dx(idx);

        adPositions = find(subDx == "AD Dementia");
        if isempty(adPositions), continue; end

        chosenADPos = NaN;
        for p = 1:numel(adPositions)
            adPos = adPositions(p);
            if any(subDx(1:adPos-1) == "Cognitively normal")
                chosenADPos = adPos;
                break;
            end
        end
        if isnan(chosenADPos), continue; end

        priorValid = 1:(chosenADPos-1);
        priorHC = priorValid(subDx(priorValid) == "Cognitively normal");
        if isempty(priorHC), continue; end

        if converterMode == "direct_HC_to_AD"
            if subDx(priorValid(end)) ~= "Cognitively normal"
                continue;
            end
        elseif converterMode ~= "ever_HC_before_AD"
            error('Unknown converterMode: %s', converterMode);
        end

        firstHCIdx = idx(priorHC(1));
        lastHCIdx = idx(priorHC(end));
        firstADIdx = idx(chosenADPos);

        firstHCDay = days(firstHCIdx);
        lastHCDay = days(lastHCIdx);
        firstADDay = days(firstADIdx);

        rows(end+1,:) = { ... %#ok<AGROW>
            subID, converterMode, ...
            sess(firstHCIdx), firstHCDay, ages(firstHCIdx), ...
            sess(lastHCIdx), lastHCDay, ages(lastHCIdx), ...
            sess(firstADIdx), firstADDay, ages(firstADIdx), ...
            firstADDay - firstHCDay, ...
            (firstADDay - firstHCDay)/365.25, ...
            firstADDay - lastHCDay, ...
            (firstADDay - lastHCDay)/365.25, ...
            (lastHCDay + firstADDay)/2};
    end

    variableNames = { ...
        'OASISID','ConverterMode', ...
        'FirstHCSession','FirstHCDay','FirstHCAge', ...
        'LastHCSession','LastHCDay','LastHCAge', ...
        'FirstADSession','FirstADDay','FirstADAge', ...
        'DaysFirstHCToFirstAD','YearsFirstHCToFirstAD', ...
        'TransitionWindowDays','TransitionWindowYears', ...
        'EstimatedConversionDayMidpoint'};
    CONVERTERS = cell2tableOrEmpty(rows, variableNames);
end

function T = listSubjectMRSessions(subjectDir, subID, requireMatchingIndex, ...
        t1PositiveKeywords, t1NegativeKeywords, ...
        funcPositiveKeywords, funcNegativeKeywords)

    entries = dir(subjectDir);
    entries = entries([entries.isdir]);
    entries = entries(~ismember({entries.name},{'.','..'}));
    rows = cell(0, 10);
    pattern = ['^' regexptranslate('escape', char(subID)) '_MR_d(-?\d+)$'];

    for i = 1:numel(entries)
        token = regexp(entries(i).name, pattern, 'tokens', 'once');
        if isempty(token), continue; end

        mrDay = str2double(token{1});
        mrSession = string(entries(i).name);
        mrDir = string(fullfile(entries(i).folder, entries(i).name));

        pair = selectCorrespondingAnatFuncPair(mrDir, requireMatchingIndex, ...
            t1PositiveKeywords, t1NegativeKeywords, ...
            funcPositiveKeywords, funcNegativeKeywords);

        rows(end+1,:) = { ... %#ok<AGROW>
            string(subID), mrSession, mrDay, mrDir, ...
            pair.Complete, pair.PairIndex, pair.AnatFolder, ...
            pair.FuncFolder, pair.T1File, pair.FuncFile};
    end

    variableNames = {'OASISID','MRSession','MRDay','MRSessionPath', ...
        'CompleteSession','SeriesPairIndex','AnatFolder','FuncFolder', ...
        'T1File','FuncFile'};
    T = cell2tableOrEmpty(rows, variableNames);
    if ~isempty(T), T = sortrows(T, 'MRDay'); end
end

function pair = selectCorrespondingAnatFuncPair(mrDir, requireMatchingIndex, ...
        t1PositiveKeywords, t1NegativeKeywords, ...
        funcPositiveKeywords, funcNegativeKeywords)

    pair = struct('Complete', false, 'PairIndex', NaN, ...
        'AnatFolder', "", 'FuncFolder', "", ...
        'T1File', "", 'FuncFile', "", 'PairScore', -Inf);

    entries = dir(char(mrDir));
    entries = entries([entries.isdir]);
    entries = entries(~ismember({entries.name},{'.','..'}));

    anatIdx = [];
    funcIdx = [];
    anatFolders = strings(0,1);
    funcFolders = strings(0,1);

    for i = 1:numel(entries)
        anatToken = regexp(entries(i).name, '^anat(\d+)$', ...
            'tokens','once','ignorecase');
        funcToken = regexp(entries(i).name, '^func(\d+)$', ...
            'tokens','once','ignorecase');

        if ~isempty(anatToken)
            anatIdx(end+1,1) = str2double(anatToken{1}); %#ok<AGROW>
            anatFolders(end+1,1) = string(fullfile(entries(i).folder, entries(i).name)); %#ok<AGROW>
        elseif ~isempty(funcToken)
            funcIdx(end+1,1) = str2double(funcToken{1}); %#ok<AGROW>
            funcFolders(end+1,1) = string(fullfile(entries(i).folder, entries(i).name)); %#ok<AGROW>
        end
    end

    if isempty(anatIdx) || isempty(funcIdx), return; end

    if requireMatchingIndex
        candidateIdx = intersect(anatIdx, funcIdx);
    else
        candidateIdx = unique([anatIdx; funcIdx]);
    end

    for c = 1:numel(candidateIdx)
        current = candidateIdx(c);
        ai = find(anatIdx == current, 1, 'first');
        fi = find(funcIdx == current, 1, 'first');
        if isempty(ai) || isempty(fi), continue; end

        [t1File, t1Score] = chooseBestNiftiInFolder(anatFolders(ai), ...
            "T1", t1PositiveKeywords, t1NegativeKeywords);
        [funcFile, funcScore] = chooseBestNiftiInFolder(funcFolders(fi), ...
            "FUNC", funcPositiveKeywords, funcNegativeKeywords);
        if t1File == "" || funcFile == "", continue; end

        pairScore = t1Score + funcScore;
        if pairScore > pair.PairScore || ...
                (pairScore == pair.PairScore && current < pair.PairIndex)
            pair.Complete = true;
            pair.PairIndex = current;
            pair.AnatFolder = anatFolders(ai);
            pair.FuncFolder = funcFolders(fi);
            pair.T1File = t1File;
            pair.FuncFile = funcFile;
            pair.PairScore = pairScore;
        end
    end
end

function [bestFile, bestScore] = chooseBestNiftiInFolder(folderPath, modality, ...
        positiveKeywords, negativeKeywords)

    bestFile = "";
    bestScore = -Inf;
    niiFiles = dir(fullfile(char(folderPath), '**', '*.nii'));
    gzFiles = dir(fullfile(char(folderPath), '**', '*.nii.gz'));
    files = [niiFiles; gzFiles];
    if isempty(files), return; end

    for i = 1:numel(files)
        filePath = string(fullfile(files(i).folder, files(i).name));
        lowerName = lower(string(files(i).name));
        score = log10(max(double(files(i).bytes),1)) / 100;

        for k = 1:numel(positiveKeywords)
            if contains(lowerName, lower(positiveKeywords(k)))
                score = score + 100 - 5*(k-1);
            end
        end
        for k = 1:numel(negativeKeywords)
            if contains(lowerName, lower(negativeKeywords(k)))
                score = score - 500;
            end
        end

        try
            info = niftiinfo(char(filePath));
            imageSize = double(info.ImageSize);
            if modality == "T1"
                if numel(imageSize)==3 || ...
                        (numel(imageSize)>=4 && imageSize(4)==1)
                    score = score + 500;
                else
                    score = score - 500;
                end
            else
                if numel(imageSize)>=4 && imageSize(4)>1
                    score = score + 1000 + min(imageSize(4),1000)/100;
                else
                    score = score - 1000;
                end
            end
        catch
        end

        if score > bestScore
            bestScore = score;
            bestFile = filePath;
        end
    end
end

function MRS = classifyMRSessionsByCDR(MRS, CDR, maxDiagnosisMatchDays)
    ids = strtrim(string(CDR.("OASISID")));
    days = numericVector(CDR.("days_to_visit"));
    dx = strtrim(string(CDR.("dx1")));
    sess = strtrim(string(CDR.("OASIS_session_label")));
    ages = numericVector(CDR.("age at visit"));

    n = height(MRS);
    MRS.MRState = repmat("Unclassified", n, 1);
    MRS.MRSelectionBasis = repmat("", n, 1);
    MRS.NearestDiagnosis = repmat("", n, 1);
    MRS.NearestCDRSession = repmat("", n, 1);
    MRS.NearestCDRDay = nan(n,1);
    MRS.NearestCDRAge = nan(n,1);
    MRS.MRtoCDRGapDays = nan(n,1);
    MRS.DiagnosisMatchQuality = repmat("", n,1);
    MRS.CDRRowIndex = nan(n,1);

    for i = 1:n
        [idx, gap] = nearestVisitIndex(ids, days, MRS.OASISID(i), ...
            MRS.MRDay(i), isfinite(days) & dx ~= "" & dx ~= ".");
        if isnan(idx)
            MRS.DiagnosisMatchQuality(i) = "No diagnostic visit";
            continue;
        end

        MRS.NearestDiagnosis(i) = dx(idx);
        MRS.NearestCDRSession(i) = sess(idx);
        MRS.NearestCDRDay(i) = days(idx);
        MRS.NearestCDRAge(i) = ages(idx);
        MRS.MRtoCDRGapDays(i) = gap;
        MRS.CDRRowIndex(i) = idx;

        if gap > maxDiagnosisMatchDays
            MRS.MRState(i) = "Unclassified";
            MRS.MRSelectionBasis(i) = "Nearest CDR beyond threshold";
        else
            if dx(idx) == "Cognitively normal"
                MRS.MRState(i) = "CN";
            elseif dx(idx) == "AD Dementia"
                MRS.MRState(i) = "AD";
            else
                MRS.MRState(i) = "Other";
            end
            MRS.MRSelectionBasis(i) = "Nearest CDR diagnosis";
        end
        MRS.DiagnosisMatchQuality(i) = gapQuality(gap);
    end
end

function quality = gapQuality(gap)
    if ~isfinite(gap)
        quality = "No match";
    elseif gap <= 180
        quality = "<=180 days";
    elseif gap <= 365
        quality = "181-365 days";
    elseif gap <= 730
        quality = "366-730 days";
    else
        quality = ">730 days";
    end
end

function [cnRow, adRow, cnRule, adRule] = selectMRPairForConverter( ...
        MRS, lastHCDay, firstADDay, midpointDay, ...
        allowBoundaryFallback, requireCNBeforeAD)

    cnRow = NaN; adRow = NaN;
    cnRule = ""; adRule = "";
    complete = find(MRS.CompleteSession);
    if isempty(complete), return; end

    cnPrimary = complete(MRS.MRState(complete)=="CN" & ...
        MRS.MRDay(complete) < firstADDay);
    adPrimary = complete(MRS.MRState(complete)=="AD" & ...
        MRS.MRDay(complete) > lastHCDay);

    if ~isempty(cnPrimary)
        cnRow = chooseNearestMR(cnPrimary, MRS.MRDay, lastHCDay, "latestTie");
        cnRule = "Nearest-to-LastHC among MR sessions classified CN";
    elseif allowBoundaryFallback
        cnFallback = complete(MRS.MRDay(complete) <= midpointDay);
        if ~isempty(cnFallback)
            cnRow = chooseNearestMR(cnFallback, MRS.MRDay, lastHCDay, "latestTie");
            cnRule = "Fallback: nearest-to-LastHC before conversion midpoint";
        end
    end

    if ~isempty(adPrimary)
        adRow = chooseNearestMR(adPrimary, MRS.MRDay, firstADDay, "earliestTie");
        adRule = "Nearest-to-FirstAD among MR sessions classified AD";
    elseif allowBoundaryFallback
        adFallback = complete(MRS.MRDay(complete) >= midpointDay);
        if ~isempty(adFallback)
            adRow = chooseNearestMR(adFallback, MRS.MRDay, firstADDay, "earliestTie");
            adRule = "Fallback: nearest-to-FirstAD after conversion midpoint";
        end
    end

    if isnan(cnRow) || isnan(adRow), return; end
    if cnRow == adRow || (requireCNBeforeAD && MRS.MRDay(cnRow) >= MRS.MRDay(adRow))
        cnRow = NaN; adRow = NaN; cnRule = ""; adRule = "";
    end
end

function chosen = chooseNearestMR(candidateRows, mrDays, targetDay, tieRule)
    gaps = abs(mrDays(candidateRows) - targetDay);
    minGap = min(gaps);
    tied = candidateRows(gaps == minGap);
    if numel(tied) == 1
        chosen = tied;
    elseif tieRule == "latestTie"
        [~, k] = max(mrDays(tied)); chosen = tied(k);
    else
        [~, k] = min(mrDays(tied)); chosen = tied(k);
    end
end

function [targetT1, targetFunc, copyStatus, copyNote] = copySelectedMR( ...
        subID, shortState, stateFolder, mrRow, B_ROOT, ...
        overwriteExisting, decompressNiiGz, dryRun)

    t1Dir = fullfile(B_ROOT, char(stateFolder), 'T1Img', char(subID));
    funcDir = fullfile(B_ROOT, char(stateFolder), 'FunImg', char(subID));
    ensureFolder(t1Dir); ensureFolder(funcDir);

    session = string(mrRow.MRSession);
    pairIndex = string(mrRow.SeriesPairIndex);
    targetT1 = string(fullfile(t1Dir, subID + "_" + shortState + "_" + ...
        session + "_anat" + pairIndex + "_T1.nii"));
    targetFunc = string(fullfile(funcDir, subID + "_" + shortState + "_" + ...
        session + "_func" + pairIndex + "_BOLD.nii"));

    try
        [targetT1, statusT1, noteT1] = copyNiftiFile(mrRow.T1File, ...
            targetT1, overwriteExisting, decompressNiiGz, dryRun);
        [targetFunc, statusFunc, noteFunc] = copyNiftiFile(mrRow.FuncFile, ...
            targetFunc, overwriteExisting, decompressNiiGz, dryRun);

        if statusT1 == "ERROR" || statusFunc == "ERROR"
            copyStatus = "ERROR";
        elseif dryRun
            copyStatus = "DRY_RUN";
        elseif statusT1 == "EXISTING" && statusFunc == "EXISTING"
            copyStatus = "EXISTING";
        else
            copyStatus = "COPIED";
        end
        copyNote = noteT1 + "; " + noteFunc;
    catch ME
        copyStatus = "ERROR";
        copyNote = string(ME.message);
    end
end

function [actualTarget, status, note] = copyNiftiFile(sourceFile, requestedNiiTarget, ...
        overwriteExisting, decompressNiiGz, dryRun)

    sourceFile = string(sourceFile);
    requestedNiiTarget = string(requestedNiiTarget);
    actualTarget = requestedNiiTarget;
    status = "ERROR"; note = "";

    if ~isfile(sourceFile)
        note = "Source file does not exist";
        return;
    end

    sourceIsGz = endsWith(lower(sourceFile), ".nii.gz");
    if sourceIsGz && ~decompressNiiGz
        actualTarget = requestedNiiTarget + ".gz";
    end

    existingCandidates = [requestedNiiTarget, requestedNiiTarget + ".gz"];
    existingMask = arrayfun(@(x) isfile(x), existingCandidates);
    if any(existingMask) && ~overwriteExisting
        actualTarget = existingCandidates(find(existingMask,1,'first'));
        status = "EXISTING";
        note = "Destination file already exists";
        return;
    end

    if dryRun
        status = "DRY_RUN";
        note = "Validated source and destination path";
        return;
    end

    if overwriteExisting
        for i = 1:numel(existingCandidates), deleteIfExists(existingCandidates(i)); end
    end

    ensureFolder(fileparts(char(actualTarget)));

    if sourceIsGz
        if decompressNiiGz
            temporaryGz = requestedNiiTarget + ".gz";
            copyfile(char(sourceFile), char(temporaryGz));
            gunzip(char(temporaryGz), fileparts(char(requestedNiiTarget)));
            deleteIfExists(temporaryGz);
            generatedNii = erase(temporaryGz, ".gz");
            if generatedNii ~= requestedNiiTarget && isfile(generatedNii)
                movefile(char(generatedNii), char(requestedNiiTarget));
            end
            actualTarget = requestedNiiTarget;
        else
            copyfile(char(sourceFile), char(actualTarget));
        end
    else
        copyfile(char(sourceFile), char(actualTarget));
    end

    if isfile(actualTarget)
        status = "COPIED"; note = "Copied";
    else
        status = "ERROR"; note = "Destination file was not created";
    end
end

function [idx, gap] = nearestVisitIndex(allIDs, allDays, targetID, targetDay, mask)
    idx = NaN; gap = NaN;
    candidates = find(allIDs == targetID & isfinite(allDays) & mask);
    if isempty(candidates), return; end

    gaps = abs(allDays(candidates) - targetDay);
    minGap = min(gaps);
    tied = candidates(gaps == minGap);
    if numel(tied) == 1
        idx = tied;
    else
        previous = tied(allDays(tied) <= targetDay);
        if ~isempty(previous)
            [~, k] = max(allDays(previous)); idx = previous(k);
        else
            [~, k] = min(allDays(tied)); idx = tied(k);
        end
    end
    gap = abs(allDays(idx) - targetDay);
end

function value = getNumericValue(T, rowIndex, variableName)
    value = NaN;
    if isnan(rowIndex) || ~ismember(variableName, string(T.Properties.VariableNames))
        return;
    end
    raw = T{rowIndex, char(variableName)};
    if iscell(raw), raw = raw{1}; end
    if isnumeric(raw) || islogical(raw)
        value = double(raw(1));
    else
        value = str2double(strtrim(string(raw)));
        if numel(value) > 1, value = value(1); end
    end
    if ~isfinite(value), value = NaN; end
end

function value = cleanRange(value, minimumValue, maximumValue)
    if ~isfinite(value) || value < minimumValue || value > maximumValue
        value = NaN;
    end
end

function [names, fields, minimums, maximums] = legacyScoreSpec()
    names = ["Legacy_LOGIMEM","Legacy_MEMUNITS","Legacy_BOSTON", ...
        "Legacy_WAIS","Legacy_MEMTIME_Minutes"];
    fields = ["LOGIMEM","MEMUNITS","BOSTON","WAIS","MEMTIME"];
    minimums = [0,0,0,0,0];
    maximums = [25,25,60,93,180];
end

function [names, fields, minimums, maximums, descriptions] = cognitiveScoreSpec()
    names = [ ...
        "WAIS_Info_Early", "WAIS_BlockDesign_Early", "WAIS_Similarities", ...
        "AnimalFluency", "VegetableFluency", "MentalControl", ...
        "AssociateMemory", "SRT_TotalRecall", "DigitSpan_Early_Total", ...
        "DigitSymbol", "TrailMakingA_Early", "TrailMakingB_Early", ...
        "BostonNaming", "LetterNumberSequencing", "Simon_OverallAccuracy", ...
        "Switch_OverallAccuracy", "WAISIII_Information", ...
        "WMSIII_LogicalMemory_Immediate", "WMSIII_VerbalPairedAssociates", ...
        "WAISIII_BlockDesign", "WMSIII_LogicalMemory_Delayed", ...
        "LineOrientation", "SlossonReading", "UDS2_DigitSpan_Total", ...
        "TrailMakingA_UDS", "Benson_Copy", "Benson_DelayedRecall", ...
        "Benson_Recognition", "UDS_VerbalFluency_Total", ...
        "CraftStory_Immediate", "CraftStory_Delayed", "MINT_Total", ...
        "UDS3_NumberSpan_Total", "NumberSymbol", "MoCA_Total", ...
        "AuditoryConsonantTrigrams", "SpatialRelations", "TrailMakingB_UDS"];

    fields = [ ...
        "PSY019","PSY021","SIM","ANIMALS","VEG","mentcont", ...
        "asscmem","srttotal","","digsym","tma","tmb","bnt", ...
        "lettnum","simon","switch","inform","logmem","pairs", ...
        "block","lmdelay","line","slosson","","TRAILA", ...
        "udsbentc","udsbentd","udsbenrs","udsvertn","crafturs", ...
        "craftdre","minttots","","numsym","mocatots","trigram", ...
        "spatial","trailb"];

    minimums = zeros(1, numel(names));
    maximums = [29,48,33,100,100,9,21,48,20,93,180,180,60,21, ...
        100,100,28,50,32,68,50,30,200,24,150,17,17,1, ...
        200,25,25,32,28,200,30,60,81,300];

    descriptions = [ ...
        "Early WAIS Information total", "Early WAIS Block Design total", ...
        "WAIS Similarities total", "Animal semantic fluency total", ...
        "Vegetable semantic fluency total", "Mental Control total", ...
        "Associate Memory total", "Selective Reminding Test total recall", ...
        "Early digit span forward plus backward total", "Digit Symbol total", ...
        "Early Trail Making A time", "Early Trail Making B time", ...
        "Boston Naming Test total", "Letter-Number Sequencing total", ...
        "Simon task overall accuracy", "Switch task overall accuracy", ...
        "WAIS-III Information total", "WMS-III Logical Memory immediate total", ...
        "WMS-III Verbal Paired Associates total", "WAIS-III Block Design total", ...
        "WMS-III Logical Memory delayed total", "Judgment of Line Orientation total", ...
        "Slosson reading total", "UDS2 digit span forward plus backward total", ...
        "UDS Trail Making A time", "Benson figure copy total", ...
        "Benson figure delayed recall total", "Benson figure recognition", ...
        "UDS letter fluency correct total", "Craft Story immediate total", ...
        "Craft Story delayed total", "Multilingual Naming Test total", ...
        "UDS3 number span forward plus backward correct total", "Number-Symbol total", ...
        "MoCA total", "Auditory Consonant Trigrams total", ...
        "Spatial Relations total", "UDS Trail Making B time"];
end

function values = extractCognitiveTotals(COG, rowIndex, names, fields, minimums, maximums)
    values = nan(1, numel(names));
    if isnan(rowIndex), return; end

    for i = 1:numel(names)
        if fields(i) ~= ""
            value = getNumericValue(COG, rowIndex, fields(i));
        else
            switch names(i)
                case "DigitSpan_Early_Total"
                    value = strictSum(cleanRange(getNumericValue(COG,rowIndex,"digfor"),0,12), ...
                                      cleanRange(getNumericValue(COG,rowIndex,"digback"),0,12));
                case "UDS2_DigitSpan_Total"
                    value = strictSum(cleanRange(getNumericValue(COG,rowIndex,"DIGIF"),0,12), ...
                                      cleanRange(getNumericValue(COG,rowIndex,"DIGIB"),0,12));
                case "UDS3_NumberSpan_Total"
                    value = strictSum(cleanRange(getNumericValue(COG,rowIndex,"digforct"),0,14), ...
                                      cleanRange(getNumericValue(COG,rowIndex,"digbacct"),0,14));
                otherwise
                    value = NaN;
            end
        end
        values(i) = cleanRange(value, minimums(i), maximums(i));
    end
end

function out = strictSum(a, b)
    if isfinite(a) && isfinite(b), out = a + b; else, out = NaN; end
end

function DICT = buildScoreDictionary(legacyNames, legacyFields, cogNames, cogFields, cogDescriptions)
    rows = cell(0, 5);
    rows(end+1,:) = {"MMSE","MMSE","OASIS3_UDSb4_cdr.csv", ...
        "Mini-Mental State Examination total","Higher is better"}; %#ok<AGROW>
    rows(end+1,:) = {"CDRSUM","CDRSUM","OASIS3_UDSb4_cdr.csv", ...
        "Clinical Dementia Rating Sum of Boxes","Higher is worse"}; %#ok<AGROW>
    rows(end+1,:) = {"CDRTOT","CDRTOT","OASIS3_UDSb4_cdr.csv", ...
        "Clinical Dementia Rating global score","Higher is worse"}; %#ok<AGROW>

    for i = 1:numel(legacyNames)
        rows(end+1,:) = {legacyNames(i), legacyFields(i), ...
            "clinical features-纵向.xlsx", "Legacy clinical/psychometric field", ...
            "Mostly higher is better; MEMTIME is timing metadata"}; %#ok<AGROW>
    end

    for i = 1:numel(cogNames)
        field = cogFields(i);
        if field == ""
            if cogNames(i) == "DigitSpan_Early_Total"
                field = "digfor + digback";
            elseif cogNames(i) == "UDS2_DigitSpan_Total"
                field = "DIGIF + DIGIB";
            elseif cogNames(i) == "UDS3_NumberSpan_Total"
                field = "digforct + digbacct";
            end
        end
        if contains(cogNames(i), "TrailMaking")
            direction = "Lower is better";
        else
            direction = "Higher is better";
        end
        rows(end+1,:) = {cogNames(i), field, ...
            "OASIS3_UDSc1_cognitive_assessments.csv", ...
            cogDescriptions(i), direction}; %#ok<AGROW>
    end

    DICT = cell2table(rows, 'VariableNames', ...
        {'OutputVariable','SourceField','SourceTable','Description','Direction'});
end

function T = cell2tableOrEmpty(rows, variableNames)
    if isempty(rows)
        T = cell2table(cell(0,numel(variableNames)), 'VariableNames', variableNames);
    else
        T = cell2table(rows, 'VariableNames', variableNames);
    end
end

function ensureFolder(folderPath)
    if ~isfolder(folderPath), mkdir(folderPath); end
end

function deleteIfExists(filePath)
    filePath = string(filePath);
    if isfile(filePath), delete(filePath); end
end
