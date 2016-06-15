% draw original matlab methods approach
close all; clear all;
% stc folder name
stc_folder = 'stcs';
pic_folder = 'pics';
hemi = {'left', 'right'};
time_points = 439;

fstem={
    'isc_051816_z_AD1_gmr_median';
    'isc_051816_z_AD1_gmr_mean';
    'isc_051816_z_AD2_gmr_median';
    'isc_051816_z_AD2_gmr_mean';
    'isc_051816_z_SU1_gmr_median';
    'isc_051816_z_SU1_gmr_mean';
    'isc_051816_z_SU2_gmr_median';
    'isc_051816_z_SU2_gmr_mean';
    };

threshold=[0.95 0.99];


% fisher z correlation deviation
dev = 1/sqrt(time_points-3);

for f_idx=1:length(fstem)
    
    [stc_lh,v_lh]=inverse_read_stc(sprintf('./%s/%s-lh.stc', stc_folder, fstem{f_idx}));
    [stc_rh,v_rh]=inverse_read_stc(sprintf('./%s/%s-rh.stc', stc_folder, fstem{f_idx}));
        
    stc=cat(1,stc_lh(:,1),stc_rh(:,1));
    
    % calculate Z value & p value
    stc_lh = 1 - 2*(1-normcdf(stc_lh/dev));
    stc_rh = 1 - 2*(1-normcdf(stc_rh/dev));   

    Zid_l=2;
    Zid_h=3;

    figure;
    etc_render_fsbrain('overlay_value',stc_lh(:,1),'overlay_vertex',v_lh,'overlay_threshold', threshold,'view_angle',[-90 0]);
    hgexport(gcf, sprintf('./%s/%s/%s',pic_folder, hemi{1}, fstem{f_idx}), hgexport('factorystyle'),'Format','png');

    figure;
    etc_render_fsbrain('overlay_value',stc_rh(:,1),'overlay_vertex',v_rh,'overlay_threshold', threshold,'view_angle',[90 0]);
    hgexport(gcf, sprintf('./%s/%s/%s', pic_folder, hemi{2}, fstem{f_idx}), hgexport('factorystyle'),'Format','png');
    
end;

return;