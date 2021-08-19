% Jonathan Garcia
% 2021-08-18

% All directories that need merging for running paper_figures.m:
% 'data_files/flux_array_runs'
% 'data_files/RNP_Ca'
% 'data_files/RNP_fits'
% 'data_files/RNP_runs'


%%

clear all; %#ok
close all;

% Load all of the smaller data files.
src_dir = 'flux_array_runs';
mats = dir(sprintf('data_files/%s', src_dir));
mats = {mats.name};
for m = 1:length(mats)
    if length(mats{m}) > 4 && all(mats{m}(end-3:end) == '.mat')
        load(sprintf('data_files/%s/%s', src_dir, mats{m}));
    end
end

% Save all data in a larger data file.
TODO = sprintf(['\nsave(''data_files/test_%s.mat'', ' ...
        '''time'', ''locs'', ...' ...
        '\n    ''S_avg'', ''S_rel'', ''sr0'', ...' ...
        '\n    ''A_avg'', ''A_rel'', ''ar0'');'], src_dir);
eval(TODO);


%%

clear all; %#ok
close all;

% Load all of the smaller data files.
src_dir = 'RNP_Ca';
mats = dir(sprintf('data_files/%s', src_dir));
mats = {mats.name};
for m = 1:length(mats)
    if length(mats{m}) > 4 && all(mats{m}(end-3:end) == '.mat')
        load(sprintf('data_files/%s/%s', src_dir, mats{m}));
    end
end

% Create cell arrays to contain subcomponents in larger data structures.
Ca_avg = cell(size(Ca_avg_1));

% Store data from smaller structures into large structures.
TODO = '';
for r = 1:5
for c = 1:4
    if isempty(Ca_avg_1{r,c})
        continue
    end
    TODO = [TODO, ...
        sprintf('\nCa_avg{%d,%d} = cell(8, 1)', r, c)]; %#ok
    for n = 1:8
        TODO = [TODO, ...
            sprintf('\nCa_avg{%d,%d}{%d} = Ca_avg_%d{%d,%d};', ...
                r, c, n, n, r, c)]; %#ok
    end
end
end
eval(TODO);

% Save all data in a larger data file.
TODO = sprintf(['\nsave(''data_files/test_%s.mat'', ' ...
    '''rISI'', ''nrAP'', ''prob'', ''Ca_avg'');'], src_dir);
eval(TODO);


%%

clear all; %#ok
close all;

% Load all of the smaller data files.
src_dir = 'RNP_fits';
mats = dir(sprintf('data_files/%s', src_dir));
mats = {mats.name};
for m = 1:length(mats)
    if length(mats{m}) > 4 && all(mats{m}(end-3:end) == '.mat')
        load(sprintf('data_files/%s/%s', src_dir, mats{m}));
    end
end

% Save all data in a larger data file.
TODO = sprintf(['\nsave(''data_files/test_%s.mat'', ''dat'', ...' ...
        '\n    ''S400'', ''S_Pv'', ''S_Pvs'', ' ...
        '''S_Ts'', ''S_Ns'', ''S_Ls'', ''S_FVU'', ...' ...
        '\n    ''A400'', ''A_Pv'', ''A_Pvs'', ' ...
        '''A_Ts'', ''A_Ns'', ''A_Ls'', ''A_FVU'', ...' ...
        '\n    ''S_F'', ''A_F'');'], src_dir);
eval(TODO);


%%

clear all; %#ok
close all;

% Load all of the smaller data files.
src_dir = 'RNP_runs';
mats = dir(sprintf('data_files/%s', src_dir));
mats = {mats.name};
for m = 1:length(mats)
    if length(mats{m}) > 4 && all(mats{m}(end-3:end) == '.mat')
        load(sprintf('data_files/%s/%s', src_dir, mats{m}));
    end
end

str = {'A_avg', 'S_rel', 'A_rel', 'C_rel', 'P_rel'};

% Create cell arrays to contain subcomponents in larger data structures.
TODO = sprintf('\nS_avg = cell(size(S_avg_1a));');
for s = 1:length(str)
    TODO = [TODO, ...
        sprintf('\n%s = cell(size(%s_1));', str{s}, str{s})]; %#ok
end
eval(TODO);

% Store data from large structures into smaller structures.
TODO = '';
for r = 1:5
for c = 1:4
    if isempty(S_avg_1a{r,c})
        continue
    end
    TODO = [TODO, ...
        sprintf('\nS_avg{%d,%d} = cell(8, 1);', r, c)]; %#ok
    for s = 1:length(str)
        TODO = [TODO, ...
            sprintf('\n%s{%d,%d} = cell(8, 1);', str{s}, r, c)]; %#ok
    end
    for n = 1:8
        TODO = [TODO, ...
            sprintf('\nS_avg{%d,%d}{%d} = [S_avg_%da{%d,%d}; S_avg_%db{%d,%d}];', ...
                r, c, n, n, r, c, n, r, c)]; %#ok
        for s = 1:length(str)
            TODO = [TODO, ...
                sprintf('\n%s{%d,%d}{%d} = %s_%d{%d,%d};', ...
                    str{s}, r, c, n, str{s}, n, r, c)]; %#ok
        end
    end
end
end
eval(TODO);

% Save all data in smaller data files.
TODO = sprintf(['\nsave(''data_files/test_%s.mat'', ''time'', ...' ...
    '\n    ''S_avg'''], src_dir);
for s = 1:length(str)
    TODO = [TODO, sprintf([', ''%s'''], str{s})]; %#ok
end
TODO = [TODO, ');'];
eval(TODO);
