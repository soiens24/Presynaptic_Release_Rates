% Jonathan Garcia
% 2017-02-16

function [V, Fs] = fit_conv_profiles(t, Y, r0, locs, V, C, fit_TF)
% V: {Ts, Ps, Ks, Ms, Ss}

% exponential rate of release convolved with
% ex-Gaussian distribution of start times
XXG = @(t, T, P, K, mu, sg) P*K/(K*T-1) * ...
    (exp(-(t-mu-sg^2/(2*T))/T).*normcdf(t, mu+sg^2/T, sg) - ...
     exp(-(t-mu-sg^2* K/2 )*K).*normcdf(t, mu+sg^2*K, sg));

% maximum K that will yield results within numerical range for doubles
max_K = @(sg, N) (sqrt(2*log(realmax) + N.^2) - N) ./ sg;

% number of components that constitute a profile
Nc = length(V.Ts);
% N1 = length(V.T1);
% N0 = length(V.T0);

if nargin < 7
    fit_TF = true;
if nargin < 6
    C = false(Nc, 4);   % no constant parameters by default
end
end
if isempty(C)
    C = false(Nc, 4);   % no constant parameters by default
end

% time and time steps
t = t(:); % dt = t(2) - t(1);

% Display initial approximation to the true functions.
if fit_TF
    figure; hold all;
    set(gcf, 'Units', 'inches', 'Position', [1 1, 8, 6]);
    plot(t, Y, 'b'); set(gca, 'YScale', 'log'); xlim([0 100]);
    plot(t, conv_fit_all(V.Ts, V.Ps, V.Ks, V.Ms, V.Ss), 'k');
    ax = gca;
    
    figure;
    set(gcf, 'Units', 'inches', 'Position', [10 1, 8, 6]);
    sp = subplot(221); hold all;
    plot(locs, V.Ps, 'k'); set(gca, 'YScale', 'log');
    sk = subplot(222); hold all;
    plot(locs, V.Ks, 'k'); set(gca, 'YScale', 'log');
    sm = subplot(223); hold all;
    plot(locs, V.Ms, 'k'); set(gca, 'YScale', 'log');
    ss = subplot(224); hold all;
    plot(locs, V.Ss, 'k'); set(gca, 'YScale', 'log');
    
    drawnow;
end

% Initialize and run fits.
Ps = V.Ps; Ks = V.Ks; Ms = V.Ms; Ss = V.Ss;
if fit_TF
    
    P0 = median(Ps(C(:,1),:), 2); K0 = median(Ks(C(:,2),:), 2);
    M0 = median(Ms(C(:,3),:), 2); S0 = median(Ss(C(:,4),:), 2);
    
    % Run fits for one distance at a time.
    i1 = (1:Nc-length(P0)) + 0;
    i2 = (1:Nc-length(K0)) + length(i1);
    i3 = (1:Nc-length(M0)) + length([i1,i2]);
    i4 = (1:Nc-length(S0)) + length([i1,i2,i3]);
    for D = 1:length(locs)
%         vars = fminsearch(@(v) conv_err(V.Ts, ...
%             v(i1), P0, v(i2), K0, ...
%             Ms(~C(:,3),D), M0, Ss(~C(:,4),D), S0, Y(:,D)), ...
%             [Ps(~C(:,1),D), Ks(~C(:,2),D)], ...
%             optimset('MaxFunEvals', 1e6, 'MaxIter', 1e6));
%         Ps(~C(:,1),D) = vars(i1);
%         Ks(~C(:,2),D) = vars(i2);
        vars = fminsearch(@(v) conv_err(V.Ts, ...
            v(i1), P0, v(i2), K0, ...
            v(i3), M0, v(i4), S0, Y(:,D)), ...
            [Ps(~C(:,1),D); Ks(~C(:,2),D); ...
             Ms(~C(:,3),D); Ss(~C(:,4),D)], ...
            optimset('MaxFunEvals', 1e6, 'MaxIter', 1e6));
        Ps(~C(:,1),D) = vars(i1);
        Ks(~C(:,2),D) = vars(i2);
        Ms(~C(:,3),D) = vars(i3);
        Ss(~C(:,4),D) = vars(i4);
        
        % Make sure all parameters continue to fall within bounds.
        [Ps(:,D), Ks(:,D), Ms(:,D), Ss(:,D)] = ...
            bound_params(Ps(:,D), Ks(:,D), Ms(:,D), Ss(:,D));
        
        % Display fitted approximation to the true functions.
        if runs == 1
            plot(ax, t, conv_fit(V.Ts, ...
                Ps(:,D), Ks(:,D), Ms(:,D), Ss(:,D)), 'r');
        else
            ax.Children(length(locs)-D+1).YData = ...
                conv_fit(V.Ts, Ps(:,D), Ks(:,D), Ms(:,D), Ss(:,D));
        end
        drawnow;
    end
    
    % Find the values of the distance-independent parameters.
    i1 = (1:length(P0)) + 0;
    i2 = (1:length(K0)) + length(i1);
    i3 = (1:length(M0)) + length([i1,i2]);
    i4 = (1:length(S0)) + length([i1,i2,i3]);
    if ~isempty([i1,i2,i3,i4])
        vars = fminsearch(@(v) conv_all_err(V.Ts, ...
            Ps(~C(:,1),:), v(i1), Ks(~C(:,2),:), v(i2), ...
            Ms(~C(:,3),:), v(i3), Ss(~C(:,4),:), v(i4), Y), ...
            [P0; K0; M0; S0], ...
            optimset('MaxFunEvals', 1e6, 'MaxIter', 1e6));
        P0 = vars(i1); K0 = vars(i2);
        M0 = vars(i3); S0 = vars(i4);
        for D = 1:length(locs)
            Ps(C(:,1),D) = P0; Ks(C(:,2),D) = K0;
            Ms(C(:,3),D) = M0; Ss(C(:,4),D) = S0;
            [Ps(:,D), Ks(:,D), Ms(:,D), Ss(:,D)] = ...
                bound_params(Ps(:,D), Ks(:,D), Ms(:,D), Ss(:,D));
        end
    end
    
    % Display parameters across distances.
    if runs == 1
        plot(sp, locs, Ps); plot(sk, locs, Ks);
        plot(sm, locs, Ms); plot(ss, locs, Ss);
    else
        for cc = 1:Nc
            sp.Children(Nc-cc+1).YData = Ps(cc,:);
            sk.Children(Nc-cc+1).YData = Ks(cc,:);
            sm.Children(Nc-cc+1).YData = Ms(cc,:);
            ss.Children(Nc-cc+1).YData = Ss(cc,:);
        end
    end
    
%     % Keep iterating until all parameters come out sorted.
%     if is_sorted(Ps, 2, 'descend') && is_sorted(Ks, 2, 'descend') && ...
%             is_sorted(Ms, 2, 'ascend') && is_sorted(Ss, 2, 'ascend')
%         break
%     else
%         Ps = sort(Ps, 2, 'descend');    % shrinking with distance
%         Ks = sort(Ks, 2, 'descend');    % slowing with distance
%         Ms = sort(Ms, 2, 'ascend');     % later start with distance
%         Ss = sort(Ss, 2, 'ascend');     % more spread with distance
%     end
end

% Store fitted parameters in the output structure.
V.Ps = Ps; V.Ks = Ks;
V.Ms = Ms; V.Ss = Ss;

Fs = conv_fit_all(V.Ts, V.Ps, V.Ks, V.Ms, V.Ss);

% Print out the parameters for easy copy/paste.
if fit_TF
    names = {'Ps', 'Ks', 'Ms', 'Ss'};
    for n = 1:length(names)
        fprintf('V.%s = [ ...\n\t', names{n});
        for cc = 1:Nc
            for D = 1:length(locs)
                fprintf('%.2e', V.(names{n})(cc,D));
                if D < length(locs)
                    fprintf(', ');
                    if mod(D, 5) == 0
                        fprintf('...\n\t');
                    end
                end
            end
            if cc < Nc
                fprintf('; ...\n\t');
            end
        end
        fprintf('];\n');
    end
end

% Compute error across all distance cases.
    function err = conv_all_err(Ts, Pd, P0, Kd, K0, Md, M0, Sd, S0, Y)
        % Combine distance-dependent with distance-independent parameters.
        P1s = zeros(Nc, length(locs));
        P1s(~C(:,1),:) = Pd;
        P1s(C(:,1),:) = P0(:) * ones(1, length(locs));
        
        K1s = zeros(Nc, length(locs));
        K1s(~C(:,2),:) = Kd;
        K1s(C(:,2),:) = K0(:) * ones(1, length(locs));
        
        M1s = zeros(Nc, length(locs));
        M1s(~C(:,3),:) = Md;
        M1s(C(:,3),:) = M0(:) * ones(1, length(locs));
        
        S1s = zeros(Nc, length(locs));
        S1s(~C(:,4),:) = Sd;
        S1s(C(:,4),:) = S0(:) * ones(1, length(locs));
        
        % Make sure all parameters continue to fall within bounds.
        for d = 1:length(locs)
            [P1s(:,d), K1s(:,d), M1s(:,d), S1s(:,d)] = ...
                bound_params(P1s(:,d), K1s(:,d), M1s(:,d), S1s(:,d));
        end
        
        F = conv_fit_all(Ts, P1s, K1s, M1s, S1s);
        err = NMSE(Y(:), F(:)) + NMSE(log(Y(:)), log(F(:)));
    end

% Figure out a measure of error between the fit with
% the current set of parameters and the actual profile.
    function err = conv_err(Ts, Pd, P0, Kd, K0, Md, M0, Sd, S0, y)
        % Combine distance-dependent with distance-independent parameters.
        P1s = zeros(Nc, 1); P1s(~C(:,1)) = Pd; P1s(C(:,1)) = P0;
        K1s = zeros(Nc, 1); K1s(~C(:,2)) = Kd; K1s(C(:,2)) = K0;
        M1s = zeros(Nc, 1); M1s(~C(:,3)) = Md; M1s(C(:,3)) = M0;
        S1s = zeros(Nc, 1); S1s(~C(:,4)) = Sd; S1s(C(:,4)) = S0;
        
        % Make sure all parameters continue to fall within bounds.
        [P1s, K1s, M1s, S1s] = bound_params(P1s, K1s, M1s, S1s);
        
        f = conv_fit(Ts, P1s, K1s, M1s, S1s);
        err = NMSE(y, f) + NMSE(log(y), log(f));
    end

% Fit profiles at all distances simultaneously.
    function F = conv_fit_all(Ts, Ps, Ks, Ms, Ss)
        F = 0*Y;
        
        % Generate fits at all positions.
        for d = 1:length(locs)
            F(:,d) = conv_fit(Ts, Ps(:,d), Ks(:,d), Ms(:,d), Ss(:,d));
        end
    end

% Create a profile of multiple components of exponentials (Ts) convolved
% with the Ca influx Gaussian, further convolved with another exponential
% (Ks) to account for diffusion. Relative contributions of each component
% (Ps) are normalized to ensure a correct total area under the curve (yp).
    function f = conv_fit(Ts, Ps, Ks, Ms, Ss)
        f = r0 + 0*t;
        
        Nc = length(Ts);
        for c = 1:Nc
            % Ensure that there are no numerical problems.
            if Ks(c) > 1/Ts(c) - 1e-9
                Ks(c) = Ks(c) + 2e-9;
            end
            
            % Construct this component.
            comp = XXG(t, Ts(c), Ps(c), Ks(c), Ms(c), Ss(c));
            comp(isnan(comp)) = 0;
            f = f + comp;
            
            % Restore Ks.
            if Ks(c) > 1/Ts(c)
                Ks(c) = Ks(c) - 2e-9;
            end
        end
    end

% Keep parameters within reasonable bounds.
    function [Ps, Ks, Ms, Ss] = bound_params(Ps, Ks, Ms, Ss)
        Ps = vlim_bounce(Ps(:), V.PL);
        Ms = vlim_bounce(Ms(:), V.ML);
        Ss = vlim_bounce(Ss(:), V.SL); %[V.SL(:,1), min(V.SL(:,2),Ms)]);
        Ks = vlim_bounce(Ks(:), [V.KL(:,1), min(V.KL(:,2),max_K(Ss,5))]);
    end

% Determine whether a matrix is sorted.
    function TF = is_sorted(A, dim, dir)
        if strcmp(dir, 'ascend')
            B = (diff(A, 1, dim) >= 0);
        elseif strcmp(dir, 'descend')
            B = (diff(A, 1, dim) <= 0);
        end
        TF = all(B(:));
    end



end
