function [] = render_all_brains(strPath, mappingExcel)
    % Output path.
    outputPathRoot = strcat(strPath, '_brain_image_render/');
    % All .stc files.
    files = getAllFiles(strPath);
    % Load the mapping info.
    [xlsNum, xlsTxt, xlsRaw] = xlsread(mappingExcel, 'A2:E19');

    for i = 1:length(files)
        % Generate the shorter name of stc file.
        string = files{i};
        pattern = '^..\\lol_stc\\([^\\]+).*\\(\d+)\\.*-(\wh)\.stc$';
        [tokens, matches] = regexp(string, pattern, 'tokens', 'match');
        % Infomation in pathname.
        subjectNo = tokens{1}(1);
        testTime  = tokens{1}(2);
        halfBrain = tokens{1}(3);
        % Mapping the data to find the view of game replay.
        [findRows, findCols]   = find(ismember(xlsTxt, subjectNo));
        [findRows2, findCols2] = find(ismember(xlsTxt(findRows, :), testTime));
        replayView = xlsTxt(1, findCols2);
        % Create the output folder.
        outputPath = strcat(outputPathRoot, subjectNo, '/', testTime, '_', replayView, '_', halfBrain, '/');
        outputPath = outputPath{1};
        if ~exist(outputPath, 'dir')
            mkdir(outputPath);
        end
        % Export the fMRI image.
        [stc, vec] = inverse_read_stc(files{i});
        for j = 1:size(stc, 2)
            figure;
            outputImagePath = strcat(outputPath, int2str(j));
            set(gcf, 'Visible', 'off');
            etc_render_fsbrain('overlay_value', stc(:,j), 'overlay_vertex', vec, 'overlay_threshold', [250 550], 'view_angle', [-90 0], 'surf', 'pial');
            hgexport(gcf, outputImagePath, hgexport('factorystyle'), 'Format', 'png');
        end
        % Log.
        disp(strcat('Now processing --- subjectNo: ', subjectNo, ', testTime: ', testTime, ', halfBrain: ', halfBrain, ' replayView: ', replayView));
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