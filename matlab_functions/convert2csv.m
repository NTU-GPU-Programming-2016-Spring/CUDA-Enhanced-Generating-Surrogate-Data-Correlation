function [] = convert2csv(strPath, rowBased)
    files = getAllFiles(strPath);
    outputPath = strcat(strPath, '_csv/');
    for i = 1:length(files)
        % Generate the shorter name of stc file.
        string = files{i};
        pattern = '^..\\lol_stc\\([^_]+)_.*\\(\d+)\\.*-(\wh)\.stc$';
        [tokens, matches] = regexp(string, pattern, 'tokens', 'match');
        newFileName = strcat(outputPath, tokens{1}(1), '-', tokens{1}(2), '-', tokens{1}(3), '.csv');
        % Create the output folder, then write the csv.
        if ~exist(outputPath, 'dir')
            mkdir(outputPath);
        end
        stc = inverse_read_stc(files{i});
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