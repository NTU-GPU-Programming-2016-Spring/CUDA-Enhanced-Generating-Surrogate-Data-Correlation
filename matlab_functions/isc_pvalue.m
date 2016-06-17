function [stc, surr_data, mean_value, p_value] = isc_pvalue(input_dir, substr)

REPEAT_TIMES = 200;
stc = [];

% Reading data from files with substr in filename under input_dir
all_files = getAllFiles(input_dir);
pattern = strcat('^.*', substr, '.*$');
match_files = regexp(all_files, pattern, 'match');
match_files = [match_files{:}];

for i=1:length(match_files)
  fprintf('Reading %s ...\r', char(match_files(i)));
  stc(:, :, i) = dlmread(char(match_files(i)));
end
fprintf('Reading %s ... done\n', char(match_files(i)));

% Surrogate testing
tic;
surr_data = [];
p_value = [];
mean_value = [];
for v_idx=1:size(stc, 1)
  fprintf('Processing position ... %5d\r', v_idx);
  data = squeeze(stc(v_idx, :, :));

  for time=1:REPEAT_TIMES
    surr = reorderingData(data);
    corr_mean = calCorrCoefMean(surr);
    % surr_data(v_idx, time, :, :) = surr;
    mean_value(v_idx, time) = corr_mean;
  end

  % Calculate P-value
  corr_mean = calCorrCoefMean(data);
  p = length(find(mean_value(v_idx, :) > corr_mean)) / REPEAT_TIMES;
  p_value(v_idx) = p;
end
fprintf('Processing position ... %5s\n', 'done');
toc

% dlmwrite('./AD1-pos.csv', mean_value, 'precision', 10);

end

function corr_mean = calCorrCoefMean(matrix)
tmp = tril(corrcoef(matrix), -1);
rr=tmp(find(tmp));
z=0.5.*log((1+rr)./(1-rr))./(1/sqrt(size(matrix,1)/2.34-3));
zm=mean(rr);
corr_mean= (exp(2.*zm)-1)./(exp(2.*zm)+1);
end

function surr = reorderingData(data)
surr = [];
for column=1:size(data, 2)
  surr(:, column) = surr_ft_algorithm3(data(:, column));
end
end

function fileList = getAllFiles(dirName)
dirData = dir(dirName);
dirIndex = [dirData.isdir];
fileList = {dirData(~dirIndex).name}';
if ~isempty(fileList)
  fileList = cellfun(@(x) fullfile(dirName,x), fileList, 'UniformOutput', false);
end
subDirs = {dirData(dirIndex).name};
validIndex = ~ismember(subDirs, {'.', '..'});

for iDir = find(validIndex)
  nextDir = fullfile(dirName, subDirs{iDir});
  fileList = [fileList; getAllFiles(nextDir)];
end
end
