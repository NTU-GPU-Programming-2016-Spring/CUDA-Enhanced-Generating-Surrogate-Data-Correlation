function [] = convert2csv(strPath, mappingExcel, rowBased, removeCounts, selectCounts)
    % Output path.
    outputPath = strcat(strPath, '_csv/');
    % All .stc files.
    files = getAllFiles(strPath);
    % Load the mapping info.
    [xlsNum, xlsTxt, xlsRaw] = xlsread(mappingExcel, 'A2:E19');
    for i = 1:length(files)
        % Generate the shorter name of stc file.
        string = files{i};
        pattern = 'lol_stc\\([^\\]+).*\\(\d+)\\.*-(\wh)\.stc$';
        [tokens, matches] = regexp(string, pattern, 'tokens', 'match');
        % Infomation in pathname.
        subjectNo = tokens{1}(1);
        testTime  = tokens{1}(2);
        halfBrain = tokens{1}(3);
        % Mapping the data to find the view of game replay.
        [findRows, findCols]   = find(ismember(xlsTxt, subjectNo));
        [findRows2, findCols2] = find(ismember(xlsTxt(findRows, :), testTime));
        replayView = xlsTxt(1, findCols2);
        % Generate file name.
        newFileName = strcat(outputPath, subjectNo, '-', testTime, '-', replayView, '-', halfBrain, '.csv');
        % Create the output folder, then write the csv.
        if ~exist(outputPath, 'dir')
            mkdir(outputPath);
        end
        stc = inverse_read_stc(files{i});
        %%% Remove some columns from start.
        stc(:, 1:removeCounts) = [];
        %%% Select some columns from start.
        stc = stc(:, 1:selectCounts);
        if rowBased
            stc = stc.';
        end
        disp(strcat('Now writing: ', newFileName));
        csvwrite(newFileName{1}, stc);
    end
end

function fileList = getAllFiles(dirName)
    dirData = dir(dirName);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    if ~isempty(fileList)
        fileList = cellfun(@(x) fullfile(dirName,x), ...
                           fileList, 'UniformOutput', false);
    end
    subDirs = {dirData(dirIndex).name};
    validIndex = ~ismember(subDirs, {'.', '..'});
    
    for iDir = find(validIndex)
        nextDir = fullfile(dirName, subDirs{iDir});
        fileList = [fileList; getAllFiles(nextDir)];
    end
end