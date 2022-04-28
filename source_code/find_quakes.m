function quakes = find_quakes(t_ss_s, Vx_ss_s, Dx_ss_s, M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the characteristics of earthquakes from off-fault seismic patches

% Input:
% t_ss_s = Nevent x 1 cell of Ntime x 1 time from spring slider simulations
% Vx_ss_s = Nevent x 1 cell of Ntime x 1 velocity from spring slider simulations
% Dx_ss_s = Nevent x 1 cell of Ntime x 1 displacement from spring slider simulations
% M       = M1 structure of parameters

% Output:
% quakes = Nevent with eartuquakes x 4 cell: 
% col1 = index of seismic patch
% col2 = vector of timing for earthquakes
% col3 = vector of maximum velocities 
% col4 = vector of moment magnitudes for each event
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unpack variables
ss_dx = M.ss_dx; ss_dz = M.ss_dz;
N_ss = M.N_ss; mu = M.mu;

quakes = cell(1, 3);      % store 1. timing of seismic events 2. maximum velocity during seismic events 2. moment magnitude of seismic events
A = ss_dx * ss_dz;

count = 0;
Nseismic_events = 0;
% determine the time when the spring sliders 'pop' for the first time based on velocities
for i = 1:N_ss 
    [pks, locs] = findpeaks(Vx_ss_s{i}, t_ss_s{i}, 'MinPeakProminence', 0.1);

    times = locs./86400;    % times at which this spring slider popped
    max_vels = pks;         % maximum velocities of these events
    Nevents = length(times);
    Nseismic_events = Nseismic_events+Nevents;
    
    if Nevents > 0
        count = count + 1;
        quakes{count, 1} = i;                % index of the seismic patch
        quakes{count, 2} = times;            % timing in days
        quakes{count, 3} = max_vels;         % maximum velocity      
        
        % Compute displacement of each seismic event
        slip_history = Dx_ss_s{i}; % slip history of this patch
        s_list = zeros(Nevents, 1);
        t_s = t_ss_s{i};
        
        for j = 1:Nevents
            [~, idx] = min(abs(t_s./86400 - (times(j)+0.1))); % add 0.1 days to ensure that displacement is calculated 
            s_list(j) = slip_history(idx);
        end
        
        if length(s_list) > 1
            if length(s_list) > 2
                ds_list = [s_list(1); diff(s_list)];
            else
                ds_list = [s_list(1); s_list(2) - s_list(1)];
            end
        else
            ds_list = s_list;
        end
        
        % mark seismic events in space-time plot
        for k = 1:Nevents
        % Compute the moment magnitude of each event
            
            S = ds_list(k);
            M0 = mu * A * S;
            Mw = 2/3*(log10(M0) - 9.1);
            
            quakes{count, 4} = Mw;
        end
    end
end

end