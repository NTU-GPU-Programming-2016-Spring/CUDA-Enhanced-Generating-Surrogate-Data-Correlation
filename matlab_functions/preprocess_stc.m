function []=preprocess_stc(input_folder, output_folder)
subject_all={
    's01_120415';
    's02_120415';
    's03_012116';
    's04_012116';
    's05_012116';
    's06_012116';
    's07_012116';
    's08_012816';
    's09_012816';
    's10_012816';
    's11_032516';
    's12_032516';
    's13_032816';
    's14_032816';
    's15_042416';
    's16_042416';
    's17_042416';
    
    };

AD1=[
    1 0 %s01
    1 0 %s02
    0 1 %s03
    0 0 %s04
    0 0 %s05
    0 0 %s06
    0 0 %s07
    0 0 %s08
    0 1 %s09
    0 0 %s10
    0 1 %s11
    0 0 %s12
    0 0 %s13
    1 0 %s14
    0 0 %s15
    1 0 %s16
    1 0 %s17
    ];


AD2=[
    0 1 %s01
    0 1 %s02
    1 0 %s03
    1 0 %s04
    1 0 %s05
    1 0 %s06
    0 1 %s07
    0 1 %s08
    1 0 %s09
    0 1 %s10
    0 0 %s11
    0 0 %s12
    0 0 %s13
    0 0 %s14
    0 0 %s15
    0 0 %s16
    0 0 %s17
    ];

SU1=[
    0 0 %s01
    0 0 %s02
    0 0 %s03
    0 0 %s04
    0 0 %s05
    0 0 %s06
    0 0 %s07
    0 0 %s08
    0 0 %s09
    0 0 %s10
    1 0 %s11
    1 0 %s12
    1 0 %s13
    0 1 %s14
    0 1 %s15
    0 1 %s16
    0 1 %s17
    ];


SU2=[
    0 0 %s01
    0 0 %s02
    0 0 %s03
    0 1 %s04
    0 1 %s05
    0 1 %s06
    1 0 %s07
    1 0 %s08
    0 0 %s09
    1 0 %s10
    0 0 %s11
    0 1 %s12
    0 1 %s13
    0 0 %s14
    1 0 %s15
    0 0 %s16
    0 0 %s17
    ];

FILE_NAME = {
    {'002' '0' '006' '0'};
    {'002' '0' '006' '0'};
    {'004' '0' '003' '0'};
    {'0' '0' '003' '004'};
    {'0' '0' '003' '004'};
    {'0' '0' '003' '004'};
    {'0' '0' '004' '003'};
    {'0' '0' '004' '003'};
    {'004' '0' '003' '0'};
    {'0' '0' '004' '003'};
    {'003' '002' '0' '0'};
    {'0' '004' '0' '005'};
    {'0' '002' '0' '003'};
    {'002' '003' '0' '0'};
    {'0' '009' '0' '007'};
    {'002' '004' '0' '0'};
    {'002' '005' '0' '0'};
};

stc_folder_all={
    'AD1';
    'SU1';
    'AD2';
    'SU2';
    };

hemi={
    'lh';
    'rh';
    };


output_stem='isc_051816_z';

confound_polynomial_order=2;
confound_sinusoidal_order=3;

%%%%%% AD1

%%%%%% AD1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for stc_idx=1:length(stc_folder_all)
    stc_type = stc_folder_all{stc_idx};
    idx = find(sum(eval(stc_type),2)>eps); %<----change this one among 'AD1', 'AD2', 'SU1', or 'SU2'
    subject = {subject_all{idx}};
    stc_folder = {sprintf('%s', stc_type)}; %<----change this one among 'AD1' (1), 'AD2' (2), 'SU1' (3), or 'SU2' (4)
    disp(stc_folder);

    oo=sprintf('%s_%s',output_stem,stc_folder_all{stc_idx});
    for hemi_idx=1:2
        %reading data
        %clear tmp variable
        clearvars stc;
        for subj_idx=1:length(subject)
            stc_folder = FILE_NAME{idx(subj_idx)}{stc_idx};
            d=dir(sprintf('%s/%s/unpack/bold/%s/*_2_fsaverage*-%s.stc', input_folder, subject{subj_idx}, stc_folder, hemi{hemi_idx}));
            [stc(:,:,subj_idx),a,b,c]=inverse_read_stc(sprintf('%s/%s/unpack/bold/%s/%s', input_folder, subject{subj_idx}, stc_folder, d(1).name));
            fprintf('[%s]::<%s>::<%s>\t%05d time ponits\n',stc_folder, subject{subj_idx}, hemi{hemi_idx}, size(stc,2));
        end;
        
        fprintf('\n');
        
        %remove the first 10 time points
        stc=stc(:,6:445,:);
        
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
                fprintf('[%1.1f%%]...\r',v_idx/size(stc,1)*100);
            end;
            %remove nuisance effects
            if(~isempty(D_prep))
                stc(v_idx,:,:)=squeeze(stc(v_idx,:,:))-D_prep*squeeze(stc(v_idx,:,:));
            else
                stc(v_idx,:,:)=squeeze(stc(v_idx,:,:));
            end;
        end;
        for subj_idx=1:length(subject)
            filename = sprintf('%s-%s-%s.csv', stc_type, subject{subj_idx}, hemi{hemi_idx});
            path_name = sprintf('%s/%s', output_folder, char(filename));
            fprintf('Writing %s ...\r', path_name);

            dlmwrite(path_name, stc(:, :, subj_idx), 'precision', 15);
            fprintf('Writing %s ... done\n\n', path);
        end;
    end;
end;
fprintf('Preprocessing [100.0%%]... done\n');
