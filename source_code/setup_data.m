function M2 = setup_data(M2, above_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load injection rate, wellhead pressure data. 
% Input:
% M2 = structure containing setup for pressure diffusion problem

% Output:
% M2 = same as input, but now with injection history and measured well head
% pressure

% Taiyi Wang 09/03/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
H4_r_struct = load(fullfile(above_dir, 'input/Cooper_Basin_HAB_4_Injection_Rate.mat'), 'd');
H4_p_struct = load(fullfile(above_dir, 'input/Cooper_Basin_HAB_4_Wellhead_Pressure.mat'), 'd');

% Injection rate converted to m^3/s
H4_rate = H4_r_struct.d.Injection_rate./60; 

% Wellhead pressure Pa
H4_wh_p = H4_p_struct.d.Wellhead_pressure.*1e6;

% Injection time converted to seconds
H4_cal_time = H4_r_struct.d.Date;    % caldenar time
H4_time = (H4_cal_time-H4_cal_time(1)) * 86400;

% make sure the time is unique
[H4_time, unique_idcs] = unique(H4_time);
H4_rate = H4_rate(unique_idcs);
H4_wh_p = H4_wh_p(unique_idcs);
H4_cal_time = H4_cal_time(unique_idcs);
clear idcs H4_r_struct H4_p_struct unique_idcs
    
M2.H4_time = H4_time;
M2.H4_rate = H4_rate;
M2.H4_wh_p = H4_wh_p;
M2.H4_cal_time = H4_cal_time;



end
