% Jonathan Garcia
% 2015-11-25

function [f_prev, spk_dt] = facil_inp(facil_F, dat)
% Facilitation of release parameter depends only on the facilitation
% of the previous spike and on the time interval since the last spike.
% Determine these independent parameters: f_next(f_prev, spk_dt)

% Make room for all the variables.
f_prev = zeros(length(dat.prob), length(dat.nrAP), length(dat.rISI));
spk_dt = zeros(length(dat.prob), length(dat.nrAP), length(dat.rISI));


% Determine each point's previous facilitation value.

% Probes after single spike facilitate from equilibrium (PPF).
f_prev(:,1,:) = 1;

for R = 1:length(dat.rISI)
    
    % Spikes on ramp facilitate from the previous ramp spikes.
    f_prev(1,2:end,R) = facil_F(1,1:end-1,R);
    
    % Probe spikes facilitate from last respective ramp spike.
    for P = 2:length(dat.prob)
        f_prev(P,:,R) = facil_F(1,:,R);
    end
    
end


% Label every point with its respective time delay from last spike.

% Base case starts at equilibrium, having settled over infinite time.
spk_dt(1,1,:) = Inf;

% All ramping spikes use their respective ISIs.
for R = 1:length(dat.rISI)
    spk_dt(1,2:end,R) = dat.rISI(R);
end

% All probe spikes use their respective ISIs.
for P = 2:length(dat.prob)
    spk_dt(P,:,:) = dat.prob(P);
end

end
