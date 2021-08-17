% Jonathan Garcia
% 2017-04-14

function [f_t, f_n, idx] = facil_response(t, spks, Ls, Ns, Ts)
% Plot response of facilitation function to an arbitrary spike train.

% number of facilitation components
Nf = length(Ls);
if length(Ns) ~= Nf || length(Ts) ~= Nf
    error('Lengths of Ls, Ns, and Ts must be equal.');
end

Xs = Ls ./ Ns;          % facilitation exponent
Xs(isnan(Xs)) = 1;

Ns = exp(Ns);           % exploring steps in log space
Ts = exp(Ts);           % exploring time constants in log space

% Make vectors horizontal.
t = t(:)'; spks = spks(:)';

f_t = zeros(Nf, 1) * t;     % internal facilitation parameter over time
f_n = zeros(Nf, 1) * spks;  % facilitation observed at each spike

for c = 1:Nf
    % First, deterministically calculate the exact facilitation sequence.
    f0 = 0; dT = Inf;
    for n = 1:length(spks)
        if n > 1
            f0 = f_n(c,n-1);
            dT = spks(n) - spks(n-1);
        end
        fv = f0 * exp(-dT / Ts(c));     % left-over facilitation
        dF = 1 - (fv / Ns(c))^Ns(c);    % step size
        f_n(c,n) = (fv + dF);           % facilitation on current spike
    end
    
    % Now use the facilitation sequence to plot what
    % happens to the internal state between spikes.
    for n = 1:length(spks)
        f_t(c,t>=spks(n)) = ...
            f_n(c,n) * exp(-(t(t>=spks(n)) - spks(n)) / Ts(c));
    end
end

% Apply nonlinearities and combine components.
f_t = prod(f_t.^Xs(:), 1);
f_n = prod(f_n.^Xs(:), 1);

% Know where to look to plot stems.
idx = [0; 0] * spks;
for n = 1:length(spks)
    b = find(t >= spks(n), 1, 'first');
    if b == 1
        idx(:,n) = 1;
    elseif isempty(b)
        idx(:,n) = length(t);
    else
        idx(1,n) = b - 1;
        idx(2,n) = b;
    end
end

end
