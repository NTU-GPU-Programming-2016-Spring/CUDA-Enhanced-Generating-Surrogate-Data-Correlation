function [] = preprocessing(input_dir, output_dir)
confound_polynomial_order=2;
confound_sinusoidal_order=3;

brains = {'lh'; }; % 'rh';};
views = {'AD1'; }; % 'AD2'; 'SU1'; 'SU2';};

all_files = getAllFiles(input_dir);

for brain_idx=1:length(brains)
  for view_idx=1:length(views)
    fprintf('Processing %s-%s ...\n', char(views(view_idx)), char(brains(brain_idx)));
    pattern = strcat('^.*', views(view_idx), '-', brains(brain_idx), '.*$');
    files = regexp(all_files, pattern, 'match');
    subjects = [files{:}];

    %reading data
    stc=[];
    for subj_idx=1:length(subjects)
      fprintf('Reading %s ...\r', char(subjects(subj_idx)));
      stc(:, :, subj_idx) = csvread(char(subjects(subj_idx)));
    end
    fprintf('Reading %s ... done\n', char(subjects(subj_idx)));

    %remove global mean
    for s_idx=1:size(stc,3)
        tmp=squeeze(stc(:,:,s_idx));
        ga=mean(tmp,1);
        tmp=tmp-tmp*ga'*inv(ga*ga')*ga;
        stc(:,:,s_idx)=tmp;
    end;

    %prepare nuisance effects to be removed by regression
    timeVec=[1:size(stc,2)]';
    D_poly=[];
    D_sinu=[];
    D_poly=ones(length(timeVec),1);
    for i=1:confound_polynomial_order
        tmp=timeVec.^(i);
        D_poly(:,i+1)=fmri_scale(tmp(:),1,0);
    end;
    for i=1:confound_sinusoidal_order
        D_sinu(:,i*2-1)=sin(timeVec.*i./timeVec(end).*pi);
        D_sinu(:,i*2)=cos(timeVec.*i./timeVec(end).*pi);
    end;
    D=cat(2,D_poly,D_sinu);
    if(~isempty(D))
        D_prep=D*inv(D'*D)*D';
    else
        D_prep=[];
    end;

    for v_idx=1:size(stc,1)
        if(mod(v_idx,100)==0)
            fprintf('Preprocessing [%1.1f%%]...\r',v_idx/size(stc,1)*100);
        end;

        %remove nuisance effects
        if(~isempty(D_prep))
            stc(v_idx, :, :) = squeeze(stc(v_idx,:,:))-D_prep*squeeze(stc(v_idx,:,:));
        else
            stc(v_idx, :, :) = squeeze(stc(v_idx,:,:));
        end;
    end;
    fprintf('Preprocessing [100.0%%]... done\n');

    for subj_idx=1:length(subjects)
      if(isunix)
        filename = strsplit(char(subjects(subj_idx)), '/');
      else
        filename = strsplit(char(subjects(subj_idx)), '\');
      end
      filename = filename(2);
      path = sprintf('%s/%s', output_dir, char(filename));

      fprintf('Writing %s ...\r', path);
      dlmwrite(path, stc(:, :, subj_idx), 'precision', 15);
    end
    fprintf('Writing %s ... done\n\n', path);
    end
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
