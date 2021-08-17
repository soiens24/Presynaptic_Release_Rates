% Jonathan Garcia
% 2017-04-09

function [Pv, Ls, Ns, Ts, Pvs, FVU, F] = fit_conv_facil(T, Y, r0, tL, ...
    V, dat, f_num, RNPs, fit_Ps, dispTF, linVlog)

% exponential rate of release convolved with
% ex-Gaussian distribution of start times
XXG = @(t, T, P, K, mu, sg) P.*K./(K.*T-1) .* ...
    (exp(-(t-mu-sg.^2./(2*T))./T).*normcdf(t, mu+sg.^2./T, sg) - ...
     exp(-(t-mu-sg.^2.* K/2 ).*K).*normcdf(t, mu+sg.^2.*K, sg));

% probability that the latest spike has begun
% causing release for a given release component
CXG = @(t, K, mu, sg) normcdf(t, mu, sg) - ...
    exp(-(t-mu-sg^2* K/2 )*K).*normcdf(t, mu+sg^2*K, sg);

Nc = length(V.Ts);

if nargin < 11
    linVlog = [1 1];
if nargin < 10
    dispTF = true;
if nargin < 9
    fit_Ps = true;
if nargin < 8
    RNPs = [];
end
end
end
end
if isempty(RNPs)
    % Make room to store the fitted histogram areas (P).
    Pv = cell(length(dat.nrAP), length(dat.rISI));
    
    % ... and the facilitation metaparameters.
    Ls = cell(Nc, 1); Ns = cell(Nc, 1); Ts = cell(Nc, 1);
else
    % Use previously fitted areas and facilitation metaparameters
    % as starting points for new fits.
    Pv = RNPs.Pv;
    Ls = RNPs.Ls;     Ns = RNPs.Ns;     Ts = RNPs.Ts;
end

if dispTF
    figure;
end
j = 0;
FVU = zeros(2, ...
    length(dat.prob) * (1 + length(dat.rISI)*(length(dat.nrAP)-1)));
F = Y;
for R = 1:length(dat.rISI)
for N = 1:length(dat.nrAP)
    
    % Skip redundant cases.
    if R > 1 && N == 1
        continue
    end
    
    if isempty(Pv{N,R})
        % Make room to store the fitted histogram areas (P).
        Pv{N,R} = cell(length(dat.prob), 1);
    end
    
for P = 1:length(dat.prob)
    
    % Find case that includes all but the latest spike.
    if P == 1 && N == 1
        last_Y = r0 + 0*T{1,1}{1};
        last_P = V.Ps;
    elseif P == 1 && N == 2
        last_Y = Y{1,1}{1};
        last_P = Pv{1,1}{1};
    elseif P == 1
        last_Y = Y{N-1,R}{1};
        last_P = Pv{N-1,R}{1};
    else
        last_Y = Y{N,R}{1};
        last_P = Pv{N,R}{1};
    end
    
    % Get timing of the latest spike to know where to begin measuring.
    spk_T = dat.rISI(R) * (dat.nrAP(N) - 1) + dat.prob(P);
    
    % Find the time since the last spike.
    if P > 1
        ISI = dat.prob(P);
    elseif N > 1
        ISI = dat.rISI(R);
    else
        ISI = Inf;
    end
    
    % Find the new best-fit release "probabilities" for the current spike.
    ix = find((T{N,R}{P} >= spk_T) & (T{N,R}{P} <= spk_T + tL));
    t = T{N,R}{P}(ix) - spk_T;
    y = Y{N,R}{P}(ix) - r0;
    y0 = last_Y(ix) - r0;
    if isempty(Pv{N,R}{P})
        Pv0 = last_P;
    else
        Pv0 = Pv{N,R}{P};
    end
    if fit_Ps
        Pv{N,R}{P} = fminsearch(@(v) ...
            conv_err(t, y, v, y0, ISI, last_P), Pv0);
        
        % Set the base-case parameters as the lower-bound.
        Pv{N,R}{P} = vlim_bounce(Pv{N,R}{P}(:), V.Ps(:) * [1 Inf]);
    end
    
    % Calculate errors of linear and logarithmic fits.
    [~, f] = conv_err(t, y, Pv{N,R}{P}, y0, ISI, last_P);
    j = j + 1;
    FVU(1,j) = NMSE(    y(:)+r0 ,     f(:)+r0 );
    FVU(2,j) = NMSE(log(y(:)+r0), log(f(:)+r0));
    F{N,R}{P} = f;
    
    if dispTF
        plot(t, y+r0, 'k');
        hold on;
        plot(t, y0+r0, 'r');
        plot(t, f+r0, 'b');
        hold off;
        axis([0 tL, 1e-9 10]);
        set(gca, 'XScale', 'lin', 'YScale', 'log');
        drawnow;
    end
    
end
end
end

% Uncover facilitation metaparameters for each release component.
F_maps = cell(Nc, 1);
Pvs = [Pv{:}]; Pvs = [Pvs{:}];
if dispTF
    figure; hold all; plot(Pvs');
    set(gca, 'YScale', 'log'); drawnow;
end

% Facilitation is current parameter value relative to initial.
facil = Pvs ./ Pvs(:,1);

for C = 1:Nc
    % Reshape facilitation values into a map for easy access.
    F_maps{C} = ...
        zeros(length(dat.prob), length(dat.nrAP), length(dat.rISI));
    j = 1;
    for R = 1:length(dat.rISI)
    for N = 1:length(dat.nrAP)
        % Skip the cases that are redundant from (R, N) == (1, 1).
        if R > 1 && N == 1 && ...
                size(facil,2) < ...
                length(dat.rISI) * length(dat.nrAP) * length(dat.prob)
            F_maps{C}(:,N,R) = F_maps{C}(:,1,1);
            continue
        end
        
        % Place facilitation value in its appropriate space.
        for P = 1:length(dat.prob)
            F_maps{C}(P,N,R) = facil(C,j);
            j = j + 1;
        end
    end
    end
    
    % number of facilitation components for this release component
    if ~isempty(f_num)
        Nf = f_num(C);
    else
        Nf = 1;
    end
    
    % Find best-fit facilitation function for this component.
    if Nf > 0
        if isempty(Ls{C})
            L0 = 1 * ones(Nf, 1);
            N0 = 2 * ones(Nf, 1);
            T0 = 3 * ones(Nf, 1);
        else
            L0 = Ls{C};
            N0 = Ns{C};
            T0 = Ts{C};
        end
        [vars_F, err_F] = fminsearch(@(v) facil_err(1, ...
            v(0*Nf+(1:Nf)), v(1*Nf+(1:Nf)), v(2*Nf+(1:Nf)), ...
            F_maps{C}), [L0, N0, T0], ...
            optimset('MaxFunEvals', 1e6, 'MaxIter', 1e6));
        Ls{C} = vars_F(0*Nf+(1:Nf));
        Ns{C} = vars_F(1*Nf+(1:Nf));
        Ts{C} = vars_F(2*Nf+(1:Nf));

        % Keep metaparameters within limits.
        Ls{C} = vlim_bounce(Ls{C}(:), ones(Nf,1) * [0 7]);
        Ns{C} = vlim_bounce(Ns{C}(:), ones(Nf,1) * [0 7]);
        Ts{C} = vlim_bounce(Ts{C}(:), ones(Nf,1) * [0 7]);
    else
        Ls{C} = 0;
        Ns{C} = 0;
        Ts{C} = -10;
    end
    
    % Due to the nature of the Ns (number of steps to saturation),
    % error surfaces tend to have a long minimum trough in the
    % N-direction. Find the smallest N for each component that
    % maintains the error minimum.
    XsC = Ls{C} ./ Ns{C};   % Saturation levels (L) are proportional to N.
    for j = 1:length(Ts{C})
        N_bound = [0, Ns{C}(j)];    % "best" N-value is really upper bound
        
        test_Ns = Ns{C};            % all variables the same ...
        while diff(N_bound) > 1e-6
            test_Ns(j) = mean(N_bound); % ... except one
            
            % Find error incurred by new set of parameters.
            new_err = ...
                facil_err(1, XsC.*test_Ns, test_Ns, Ts{C}, F_maps{C});
            
            % Update upper and lower bounds.
            if new_err < err_F + 1e-12      % lower N adds < 1 ppt to error
                N_bound(2) = mean(N_bound); % --> decrease upper bound
            else                            % lower N increases error more
                N_bound(1) = mean(N_bound); % --> increase lower bound
            end
        end
        
        % Replace original N-value with the new lower bound.
        Ns{C}(j) = N_bound(2);
    end
    Ls{C} = XsC .* Ns{C};   % Restore true saturation levels.
    Ls{C}(isnan(Ls{C})) = 0;
    
    % Replace actual facilitated parameter values with those
    % approximated by the best-fit facilitation function.
    ff = facil_gen(Ls{C}, Ns{C}, Ts{C});
    ff(:,1,2:end) = NaN;
    facil(C,:) = ff(~isnan(ff(:)));
end
Pvs = Pvs(:,1) .* facil;
if dispTF
    plot(Pvs'); drawnow;
end

% Determine error in the shape of a release profile based on the current
% spike and on the shape of the previous profile.
    function [err, f] = conv_err(t, y, Ps, y0, ISI, last_P)
        % Set the base-case parameters as the lower-bound.
        Ps = vlim_bounce(Ps(:), V.Ps(:) * [1 Inf]);
        
        % Assume what the shape of the last profile would have looked like
        % without the most recent spike interfering (after time 0).
        last_F = conv_fit(t+ISI, V.Ts, last_P, V.Ks, V.Ms, V.Ss);
        
        % Find out how much this profile differs from the actual to
        % estimate the latent component from any previous activity.
        prev_D = y0(:) - sum(last_F, 2);
        
        % Introduce the latest spike's interference component-by-component.
        for c = 1:Nc
            last_F(:,c) = last_F(:,c) .* ...
                (1 - CXG(t(:), V.Ks(c), V.Ms(c), V.Ss(c)));
        end
        
        % Now get the shape of the current profile.
        curr_F = conv_fit(t, V.Ts, Ps, V.Ks, V.Ms, V.Ss);
        
        % Combine profiles and find the error of the fit.
        f = sum(curr_F, 2) + sum(last_F, 2) + prev_D;
        err = linVlog(1) * NMSE(    y(:)+r0 ,     f(:)+r0 ) ...
            + linVlog(2) * NMSE(log(y(:)+r0), log(f(:)+r0));
    end

% Create a profile of multiple components of exponentials (Ts) convolved
% with the Ca influx Gaussian, further convolved with another exponential
% (Ks) to account for diffusion. Relative contributions of each component
% (Ps) are normalized to ensure a correct total area under the curve (yp).
    function f = conv_fit(t, Ts, Ps, Ks, Ms, Ss)
        f = t(:) * (1:length(Ts));
        
        for c = 1:Nc
            % Ensure that there are no numerical problems.
            if Ks(c) > 1/Ts(c) - 1e-9
                Ks(c) = Ks(c) + 2e-9;
            end
            
            % Construct this component.
            f(:,c) = XXG(t(:), Ts(c), Ps(c), Ks(c), Ms(c), Ss(c));
            f(isnan(f(:,c)),c) = 0;
            
            % Restore Ks.
            if Ks(c) > 1/Ts(c)
                Ks(c) = Ks(c) - 2e-9;
            end
        end
    end

% Given a matrix of facilitation values for increasing numbers of spikes
% for a given set of ISIs, find the error when attempting to replicate the
% facilitation from scratch.
    function err = facil_err(f0, Ls, Ns, Ts, facil_M)
        % Keep metaparameters within limits.
        f0 = vlim_bounce(f0, [0 Inf]);
        Ls = vlim_bounce(Ls(:), ones(Nf,1) * [0 7]);
        Ns = vlim_bounce(Ns(:), ones(Nf,1) * [0 7]);
        Ts = vlim_bounce(Ts(:), ones(Nf,1) * [0 7]);
        
        % Build expected facilitation values.
        facil_F = f0 * facil_gen(Ls, Ns, Ts);
        
        % Determine error from expected.
        err = linVlog(1) * NMSE(    facil_M(:) ,     facil_F(:) ) ...
            + linVlog(2) * NMSE(log(facil_M(:)), log(facil_F(:)));
    end

% Generate the expected facilitation values of each
% spike for a given set of ISIs and number of spikes.
    function facil_F = facil_gen(Ls, Ns, Ts)
        
        Xs = Ls ./ Ns;   % facilitation exponent
        Xs(isnan(Xs)) = 1;
        
        Ns = exp(Ns);    % exploring steps in log space
        Ts = exp(Ts);    % exploring time constants in log space
        
        fs = zeros(length(dat.prob), length(dat.nrAP), ...
            length(dat.rISI), Nf);
        
        for c = 1:Nf
        for r = 1:length(dat.rISI)
            % First faciltate along the ramps, with no probes.
            f0 = 0;
            for n = 1:length(dat.nrAP)
                d_T = dat.rISI(r);
                fv = f0 * exp(-d_T / Ts(c));    % left-over facilitation
                d_f = 1 - (fv / Ns(c)).^Ns(c);  % step size
                fs(1,n,r,c) = (fv + d_f);
                f0 = fs(1,n,r,c);
            end
            
            % Now, facilitate the probes after the ramps.
            f0 = fs(1,:,r,c);
            for p = 2:length(dat.prob)
                d_T = dat.prob(p);
                fv = f0 * exp(-d_T / Ts(c));    % left-over facilitation
                d_f = 1 - (fv / Ns(c)).^Ns(c);  % step size
                fs(p,:,r,c) = (fv + d_f);
            end
        end
            fs(:,:,:,c) = fs(:,:,:,c).^Xs(c);
        end
        
        facil_F = prod(fs, 4);
    end

end
