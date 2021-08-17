% Jonathan Garcia
% 2016-04-14

function [vars_F, param] = facil_fit(vars0F, sgn, param, dat, vlims)
% Fit facilitation function parameters for a particular fitting parameter.

if nargin < 4
    vlims = [];
%     vlims = ones(5, 1) * [0 Inf];
%     vlims(5,1) = -Inf;
%     vlims = ones(4, 1) * [0 Inf];
%     vlims(4,1) = -Inf;
%     vlims = ones(6, 1) * [0 Inf];
%     vlims([3,6],1) = -Inf;
end

facil_M = facil_map(param, dat);

% Find best-fit facilitation function.
[vars_F, err_F] = fminsearch(@(vars) ...
    facil_err(vars, sgn, dat, facil_M, vlims), [1; vars0F(:)], ...
    optimset('Display', 'off', 'MaxFunEvals', 1e6, ...
             'MaxIter',  1e6,  'TolFun',      1e-6));

% Don't let fitting algorithms push any variable outside its limits.
vars_F(:) = vlim_bounce(vars_F(:), [[0 Inf]; vlims]);

scaleF = vars_F(1);
vars_F = reshape(vars_F(2:end), 3, length(vars0F(:))/3);

% Due to the nature of the N-variable (number of steps to saturation),
% error surfaces tend to have a long minimum trough in the N-direction.
% Find the smallest N for each component that maintains the error minimum.
for j = 1:size(vars_F,2)
    N_bound = [0, vars_F(1,j)]; % "best" N-value is really upper bound
    
    test_VF = vars_F;           % all variables the same ...
    while diff(N_bound) > 1e-6
        test_VF(1,j) = mean(N_bound);   % ... except one
        
        % Find error incurred by new set of parameters.
        new_err = facil_err([1; test_VF(:)], sgn, dat, facil_M, vlims);
        
        % Update upper and lower bounds.
        if new_err < err_F + 1e-12      % lower N adds < 1 ppt to error
            N_bound(2) = mean(N_bound); % --> decrease upper bound
        else                            % lower N increases error too much
            N_bound(1) = mean(N_bound); % --> increase lower bound
        end
    end
    
    % Replace original N-value with the new lower bound.
    vars_F(1,j) = N_bound(2);
end

% Replace actual facilitated parameter values with those
% approximated by the best-fit facilitation function.
facil = scaleF * facil_gen(vars_F, sgn, dat);
facil(:,1,2:end) = NaN; facil = facil(~isnan(facil(:)));

param(:) = facil(:) * param(1);

end
