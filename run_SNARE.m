% Jonathan Garcia
% 2015-09-29

function [S_avg, A_avg, S_rel, A_rel] = run_SNARE(time, Ca_wave, S_0, A_0)

if length(time) ~= length(Ca_wave)
    error('Time must have the same length as Ca trace.');
end

% SNARE Dynamics (Nadkarni, Bartol, Sejnowski, Levine - 2010)
ks_p =   0.612e-1;  % 1/(uM*ms) - synchronous association rate
ks_n =   2.32;      %  msec^-1  - synchronous dissociation rate
ka_p =   3.822e-3;  % 1/(uM*ms) - asynchronous association rate
ka_n =  13e-3;      %  msec^-1  - asynchronous dissociation rate
b    =   0.25;      % Ca-unbinding nonlinearity factor
relS =   2*3;       %  msec^-1  - synchronous-activated release rate
relA =   2*0.025;   %  msec^-1  - asynchronous-activated release rate

dt = time(2) - time(1);     % assuming constant steps

% Track synchronous, asynchronous, and release states over time.
S_avg = zeros(6, length(time));     % 6 sync states [0-5 Ca bound]
S_rel = zeros(1, length(time));     % release from sync mechanism

A_avg = zeros(3, length(time));     % 3 async states [0-2 Ca bound]
A_rel = zeros(1, length(time));     % release from async mechanism

% Make state transitions easier with time-dependent transition matrices.
S_n = @(n) n * b^(n-1) * ks_n;
S_p = @(n, Ca) Ca * (5-n) * ks_p;
S_trans = @(Ca) ...
    [               -S_p(0,Ca), +S_n(1), 0, 0, 0, 0; ...
    +S_p(0,Ca), -S_n(1)-S_p(1,Ca), +S_n(2), 0, 0, 0; ...
    0, +S_p(1,Ca), -S_n(2)-S_p(2,Ca), +S_n(3), 0, 0; ...
    0, 0, +S_p(2,Ca), -S_n(3)-S_p(3,Ca), +S_n(4), 0; ...
    0, 0, 0, +S_p(3,Ca), -S_n(4)-S_p(4,Ca), +S_n(5); ...
    0, 0, 0, 0, +S_p(4,Ca), -S_n(5)];

A_n = @(n) n * b^(n-1) * ka_n;
A_p = @(n, Ca) Ca * (2-n) * ka_p;
A_trans = @(Ca) ...
    [              -A_p(0,Ca), +A_n(1), 0; ...
    +A_p(0,Ca), -A_n(1)-A_p(1,Ca), +A_n(2); ...
    0, +A_p(1,Ca), -A_n(2)];

% Track probabilities of being in states/releasing over time.
S_avg(:,1) = S_0; S_rel(1) = S_0(6) * relS;
A_avg(:,1) = A_0; A_rel(1) = A_0(3) * relA;

% For numerical stability, do calculations at finer time resolution.
N = ceil(dt / 0.005);     % Split into time steps of 0.005 msec. % 50;

for t = 2:length(time)
    
    % Initialize new states as previous states.
    S_avg(:,t) = S_avg(:,t-1);
    A_avg(:,t) = A_avg(:,t-1);
    
for j = 1:N     % stability
    
    % -------- Determine release rates. -------- %
    
    % Synchronous SNARE
    S_last = S_avg(:,t);
    S_rel(t) = S_last(6) * relS;   % Hz
    
    % Asynchronous SNARE
    A_last = A_avg(:,t);
    A_rel(t) = A_last(3) * relA;   % Hz
    
    
    % -------- Remove released vesicles. -------- %
    
    % Probabilities of different types of release
    p_S = S_rel(t) * dt/N;
    p_A = A_rel(t) * dt/N;
    
    % Remove synchronously released vesicles from releasible state.
    S_last(6) = S_last(6) - p_S;
    
    % Remove asynchronously released vesicles from releasible state.
    A_last(3) = A_last(3) - p_A;
    
    
    % -------- Apply state transitions. -------- %
    
    % Synchronous SNARE
    S_avg(:,t) = S_last + S_trans(Ca_wave(t)) * S_last * dt/N;
    
    % Asynchronous SNARE
    A_avg(:,t) = A_last + A_trans(Ca_wave(t)) * A_last * dt/N;
    
        
    % -------- Renormalize states by docked probability. -------- %
    
    % Synchronous SNARE
    S_avg(:,t) = S_avg(:,t) / sum(S_avg(:,t));
    
    % Asynchronous SNARE
    A_avg(:,t) = A_avg(:,t) / sum(A_avg(:,t));
    
end

end

end
