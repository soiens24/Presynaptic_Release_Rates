% Jonathan Garcia
% 2015-10-29

function err = facil_err(vars, sgn, dat, facil_V, vlims)
% Given a matrix of facilitation values for increasing numbers of spikes
% for a given set of ISIs, find the error when attempting to replicate the
% facilitation from scratch.

if nargin < 4
    vlims = [];
end

% Don't let fitting algorithms push any variable outside its limits.
vars(:) = vlim_bounce(vars(:), [[0 Inf]; vlims]);

% Build expected facilitation values.
facil_F = vars(1) * facil_gen(vars(2:end), sgn, dat);

% % Get ranges in log space for fit and original data,
% % and use them to normalize the mean squared error.
% logYR = (max(log(facil_V(:))) - min(log(facil_V(:))));
% logFR = (max(log(facil_F(:))) - min(log(facil_F(:))));

% Determine error from expected.
err = NMSE(log(facil_V), log(facil_F));
% err = mean((log(facil_V(:)) - log(facil_F(:))).^2 / (logYR * logFR));

end
