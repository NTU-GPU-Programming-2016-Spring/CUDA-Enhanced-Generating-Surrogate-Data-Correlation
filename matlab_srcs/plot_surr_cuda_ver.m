% surrogate data drawing
% draw both median of correlation and mean of correlation

close all; clear all;
%------ set save path and read path ------
pic_folder='result_pics/cuda';
cor_folder_name = 'stcs';

%------ end ------

roles={'AD1','AD2','SU1','SU2'};


% load stc files

% brain graph threshold
threshold = [0.95 0.99];
tic;
for role_idx=1:1
	mean_cor_file_name = sprintf('isc_051816_z_%s_gmr_cor_mean',roles{role_idx});
	median_cor_file_name = sprintf('isc_051816_z_%s_gmr_cor_median',roles{role_idx});

	[stc_lh,v_lh]=inverse_read_stc(sprintf('./%s/%s-lh.stc', cor_folder_name, mean_cor_file_name));
	[stc_rh,v_rh]=inverse_read_stc(sprintf('./%s/%s-rh.stc', cor_folder_name, mean_cor_file_name));

	[md_stc_lh,v_lh]=inverse_read_stc(sprintf('./%s/%s-lh.stc', cor_folder_name, median_cor_file_name));
	[md_stc_rh,v_rh]=inverse_read_stc(sprintf('./%s/%s-rh.stc', cor_folder_name, median_cor_file_name));

	lh_cor = stc_lh(:,1);
	rh_cor = stc_rh(:,1);
	m_cor = cat(1, lh_cor, rh_cor);

	md_lh_cor = md_stc_lh(:,1);
	md_rh_cor = md_stc_rh(:,1);
	md_cor = cat(1, lh_cor, rh_cor);

	folder_name = sprintf('%s-pos', roles{role_idx});
	lh_p_val_series = {};
	rh_p_val_series = {};
	md_lh_p_val_series = {};
	md_rh_p_val_series = {};

	for brain_idx=1:10242
		file_name = sprintf('%s-pos-%s.csv', roles{role_idx}, int2str(brain_idx));

		file_dest = sprintf('./%s/%s', folder_name, file_name);
		dest = csvread(file_dest);
		surr_size = size(dest, 2);
		p_val = sum(find(dest >= m_cor(brain_idx)))/surr_size;
		md_p_val = sum(find(dest >= md_cor(brain_idx)))/surr_size;

		if brain_idx <= 10242
			lh_p_val_series = [lh_p_val_series;p_val];
			md_lh_p_val_series = [md_lh_p_val_series;md_p_val];
		else
			rh_p_val_series = [rh_p_val_series;p_val];
			md_rh_p_val_series = [md_rh_p_val_series;md_p_val];
		end	
		fprintf('No of brain point : %s \n', int2str(brain_idx));
	end;

	lh_p_val_series = cell2mat(lh_p_val_series);
	rh_p_val_series = cell2mat(rh_p_val_series);	

	lh_inv_p_val_series = 1 - lh_p_val_series.*2.*20484;
	rh_inv_p_val_series = 1 - rh_p_val_series.*2.*20484;
	
	% convert negative values to 0
	lh_inv_p_val_series(find(lh_inv_p_val_series < 0)) = 0.001;
	rh_inv_p_val_series(find(rh_inv_p_val_series < 0)) = 0.001;

	md_lh_p_val_series = cell2mat(md_lh_p_val_series);
	md_rh_p_val_series = cell2mat(md_rh_p_val_series);	

	md_lh_inv_p_val_series = 1 - md_lh_p_val_series.*2.*20484;
	md_rh_inv_p_val_series = 1 - md_rh_p_val_series.*2.*20484;

	% convert negative values to 0
	md_lh_inv_p_val_series(find(md_lh_inv_p_val_series < 0)) = 0.001;
	md_rh_inv_p_val_series(find(md_rh_inv_p_val_series < 0)) = 0.001;

	% left brain figures
	figure;
    etc_render_fsbrain('overlay_value',lh_inv_p_val_series,'overlay_vertex',v_lh,'overlay_threshold', threshold,'view_angle',[-90 0]);
    hgexport(gcf, sprintf('./%s/%s/%s_mean_%s',pic_folder, roles{role_idx}, roles{role_idx}, 'lh'), hgexport('factorystyle'),'Format','png');
    figure;
    etc_render_fsbrain('overlay_value',md_lh_inv_p_val_series,'overlay_vertex',v_lh,'overlay_threshold', threshold,'view_angle',[-90 0]);
    hgexport(gcf, sprintf('./%s/%s/%s_median_%s',pic_folder, roles{role_idx}, roles{role_idx}, 'lh'), hgexport('factorystyle'),'Format','png');

    % right brain figures
    % figure;
    % etc_render_fsbrain('overlay_value',rh_inv_p_val_series,'overlay_vertex',v_rh,'overlay_threshold', threshold,'view_angle',[90 0]);
    % hgexport(gcf, sprintf('./%s/%s/%s_mean_%s',pic_folder, roles{role_idx}, roles{role_idx}, 'rh'), hgexport('factorystyle'),'Format','png');
    % etc_render_fsbrain('overlay_value',md_rh_inv_p_val_series,'overlay_vertex',v_rh,'overlay_threshold', threshold,'view_angle',[90 0]);
    % hgexport(gcf, sprintf('./%s/%s/%s_median_%s',pic_folder, roles{role_idx}, roles{role_idx}, 'rh'), hgexport('factorystyle'),'Format','png');
end;

toc