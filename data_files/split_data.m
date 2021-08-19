% Jonathan Garcia
% 2021-08-17

% All files referenced in paper_figures.m:
% load('data_files/flux_array_Ca.mat');
% load('data_files/flux_array_Ca_no_Cb.mat');
% load('data_files/flux_array_hists.mat');
% load('data_files/flux_array_runs.mat');
% load('data_files/flux_array_Sv_Av.mat');
% load('data_files/impulse_response.mat');
% load('data_files/mixed_Ca.mat');
% load('data_files/mixed_runs.mat');
% load('data_files/mixed_sim_VDCC.mat');
% load('data_files/RNP_Ca.mat');
% load('data_files/RNP_fits.mat');
% load('data_files/RNP_runs.mat');
% load('data_files/RNP_Sv_Av.mat');
% load('data_files/SNARE_0.mat');
% load('data_files/SNARE_states.mat');
% load('data_files/true_facil.mat');

% All files that need splitting to < 25 MB for GitHub:
% load('data_files/flux_array_runs.mat');
% load('data_files/RNP_Ca.mat');
% load('data_files/RNP_fits.mat');
% load('data_files/RNP_runs.mat');

%%

clear all; %#ok
close all;
load('data_files/flux_array_runs.mat');

new_dir = 'flux_array_runs';
mkdir(sprintf('data_files/%s', new_dir));

% Save all data in smaller data files.
TODO = [...
    sprintf(['\nsave(''data_files/%s/presyn.mat'', ' ...
        '''time'', ''locs'');'], new_dir), ...
    sprintf(['\nsave(''data_files/%s/synchr.mat'', ' ...
        '''S_avg'', ''S_rel'', ''sr0'');'], new_dir), ...
    sprintf(['\nsave(''data_files/%s/asynch.mat'', ' ...
        '''A_avg'', ''A_rel'', ''ar0'');'], new_dir)];
eval(TODO);


%%

clear all; %#ok
close all;
load('data_files/RNP_Ca.mat');

new_dir = 'RNP_Ca';
mkdir(sprintf('data_files/%s', new_dir));

% Create cell arrays to contain subcomponents of larger data structures.
TODO = '';
for n = 1:8
    TODO = [TODO, ...
    	sprintf('\nCa_avg_%d = cell(size(Ca_avg));', n)]; %#ok
end
eval(TODO);

% Store data from large structures into smaller structures.
TODO = '';
for r = 1:5
for c = 1:4
    N = length(Ca_avg{r,c});
    for n = 1:N
        TODO = [TODO, ...
            sprintf('\nCa_avg_%d{%d,%d} = Ca_avg{%d,%d}{%d};', ...
                n, r, c, r, c, n)]; %#ok
    end
end
end
eval(TODO);

% Save all data in smaller data files.
TODO = sprintf(['\nsave(''data_files/%s/RNP.mat'', ' ...
    '''rISI'', ''nrAP'', ''prob'');'], new_dir);
for n = 1:8
    TODO = [TODO, ...
        sprintf(['\nsave(''data_files/%s/Ca_avg_%d.mat'', ' ...
            '''Ca_avg_%d'');'], new_dir, n, n)]; %#ok
end
eval(TODO);


%%

clear all; %#ok
close all;
load('data_files/RNP_fits.mat');

new_dir = 'RNP_fits';
mkdir(sprintf('data_files/%s', new_dir));

% % Create cell arrays to contain subcomponents of larger data structures.
% TODO = '';
% for n = 1:8
%     TODO = [TODO, ...
%     	sprintf('\nS_F_%d = cell(size(S_F));', n), ...
%     	sprintf('\nA_F_%d = cell(size(A_F));', n)]; %#ok
% end
% eval(TODO);
% 
% % Store data from large structures into smaller structures.
% TODO = '';
% for r = 1:5
% for c = 1:4
%     N = length(S_F{r,c});
%     for n = 1:N
%         TODO = [TODO, ...
%             sprintf('\nS_F_%d{%d,%d} = S_F{%d,%d}{%d};', ...
%                 n, r, c, r, c, n), ...
%             sprintf('\nA_F_%d{%d,%d} = A_F{%d,%d}{%d};', ...
%                 n, r, c, r, c, n)]; %#ok
%     end
% end
% end
% eval(TODO);

% Save all data in smaller data files.
TODO = [...
    sprintf(['\nsave(''data_files/%s/data.mat'', ''dat'', ...' ...
        '\n    ''S400'', ''S_Pv'', ''S_Pvs'', ' ...
        '''S_Ts'', ''S_Ns'', ''S_Ls'', ''S_FVU'', ...' ...
        '\n    ''A400'', ''A_Pv'', ''A_Pvs'', ' ...
        '''A_Ts'', ''A_Ns'', ''A_Ls'', ''A_FVU'');'], new_dir), ...
    sprintf('\nsave(''data_files/%s/S_F.mat'', ''S_F'');', new_dir), ...
    sprintf('\nsave(''data_files/%s/A_F.mat'', ''A_F'');', new_dir)];
% for n = 1:8
%     TODO = [TODO, ...
%         sprintf(['\nsave(''data_files/%s/S_F_%d.mat'', ' ...
%             '''S_F_%d'');'], new_dir, n, n), ...
%         sprintf(['\nsave(''data_files/%s/A_F_%d.mat'', ' ...
%             '''A_F_%d'');'], new_dir, n, n)]; %#ok
% end
eval(TODO);


%%

clear all; %#ok
close all;
load('data_files/RNP_runs.mat');

new_dir = 'RNP_runs';
mkdir(sprintf('data_files/%s', new_dir));

str = {'A_avg', 'S_rel', 'A_rel', 'C_rel', 'P_rel'};

% Create cell arrays to contain subcomponents of larger data structures.
TODO = '';
for n = 1:8
    TODO = [TODO, ...
        sprintf('\nS_avg_%da = cell(size(S_avg));', n), ...
        sprintf('\nS_avg_%db = cell(size(S_avg));', n)]; %#ok
    for s = 1:length(str)
    	TODO = [TODO, ...
            sprintf('\n%s_%d  = cell(size(%s));', str{s}, n, str{s})]; %#ok
    end
end
eval(TODO);

% Store data from large structures into smaller structures.
TODO = '';
for r = 1:5
for c = 1:4
    N = length(S_avg{r,c});
    for n = 1:N
        TODO = [TODO, ...
            sprintf('\nS_avg_%da{%d,%d} = S_avg{%d,%d}{%d}(1:3,:);', ...
                n, r, c, r, c, n), ...
            sprintf('\nS_avg_%db{%d,%d} = S_avg{%d,%d}{%d}(4:6,:);', ...
                n, r, c, r, c, n)]; %#ok
        for s = 1:length(str)
            TODO = [TODO, ...
                sprintf('\n%s_%d{%d,%d}  = %s{%d,%d}{%d};', ...
                    str{s}, n, r, c, str{s}, r, c, n)]; %#ok
        end
    end
end
end
eval(TODO);

% Save all data in smaller data files.
TODO = sprintf('\nsave(''data_files/%s/time.mat'', ''time'');', new_dir);
for n = 1:8
    TODO = [TODO, ...
        sprintf(['\nsave(''data_files/%s/S_avg_%da.mat'', ' ...
            '''S_avg_%da'');'], new_dir, n, n), ...
        sprintf(['\nsave(''data_files/%s/S_avg_%db.mat'', ' ...
            '''S_avg_%db'');'], new_dir, n, n)]; %#ok
    for s = 1:length(str)
    	TODO = [TODO, ...
            sprintf(['\nsave(''data_files/%s/%s_%d.mat'',  ' ...
            '''%s_%d'');'], new_dir, str{s}, n, str{s}, n)]; %#ok
    end
end
eval(TODO);
