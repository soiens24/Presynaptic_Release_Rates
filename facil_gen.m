% Jonathan Garcia
% 2015-10-29

function facil_F = facil_gen(vars, sgn, dat)
% Generate the expected facilitation values of each
% spike for a given set of ISIs and number of spikes.

vars = vars(:);

if ~any(abs(sgn) == 1 | sgn == 0)
    error('Signs must be +/-1 or 0.');
end

% Sort out variables.
Ls = vars(1:3:end-2);   % saturation level
Ns = vars(2:3:end-1);   % number of steps to saturation
Ts = vars(3:3:end-0);   % (ms) decay time constant

Xs = sgn .* Ls ./ Ns;   % facilitation exponent
Xs(isnan(Xs)) = 1;

% Ns = vars(1:3:end-2);   % saturation points % number of steps to saturation
% Ts = vars(2:3:end-1);   % (ms) decay time constant
% Xs = vars(3:3:end-0);   % magnitude and direction of facilitation
% Xs = Xs .* sgn;         % attaching sign to exponent

Ns = exp(Ns);           % exploring steps in log space
Ts = exp(Ts);           % exploring time constants in log space

Nf = length(vars) / 3;

fs = zeros(length(dat.prob), length(dat.nrAP), length(dat.rISI), Nf);

for c = 1:Nf
for R = 1:length(dat.rISI)
    % First faciltate along the ramps, with no probes.
    f0 = 0;
    for N = 1:length(dat.nrAP)
        d_T = dat.rISI(R);
        fv = f0 * exp(-d_T / Ts(c));    % left-over facilitation
        d_f = 1 - (fv / Ns(c)).^Ns(c);  % step size
        fs(1,N,R,c) = (fv + d_f);
        f0 = fs(1,N,R,c);
    end
    
    % Now, facilitate the probes after the ramps.
    f0 = fs(1,:,R,c);
    for P = 2:length(dat.prob)
        d_T = dat.prob(P);
        fv = f0 * exp(-d_T / Ts(c));    % left-over facilitation
        d_f = 1 - (fv / Ns(c)).^Ns(c);  % step size
        fs(P,:,R,c) = (fv + d_f);
    end
end
    fs(:,:,:,c) = fs(:,:,:,c).^Xs(c);
end

facil_F = prod(fs, 4);

end

% % Sort out variables.
% if length(vars) ~= 6
%     error(sprintf(['Facilitation fitting parameters:\n' ...
%         'ha, Ta, hb, Tr, Tf, Fx'])); %#ok
% end
% % h = vars(1);
% % Tf = vars(2);
% % Tr = vars(3);
% % Fx = vars(4);
% 
% % ha = vars(1); Ta = vars(2); %ha(Ta>1e+4) = 0;
% % hb = vars(3); Tb = vars(4); %hb(Tb<1e-4) = 0;
% % Fx = vars(5);
% 
% ha = vars(1); Ta = vars(2);
% hb = vars(3); Xb = vars(4); Tr = vars(5);
% Fx = vars(6);
% 
% % ha = vars(1); Ta = vars(2); xa = vars(3);
% % hb = vars(1); Tb = vars(2); xb = vars(3);
% 
% facil_F = zeros(length(dat.prob), length(dat.nrAP), length(dat.rISI), 2);
% 
% for R = 1:length(dat.rISI)
%     % First faciltate along the ramps, with no probes.
%     for N = 2:length(dat.nrAP)
%         delT = dat.rISI(R);
%         rise = exp(-delT / Tr);
%         fall = delT^(-Xb); %exp(-delT / Tf);
% %         facil_F(1,N,R) = (facil_F(1,N-1,R) + h) * fall * (1 - rise) ...
% %                         + facil_F(1,N-1,R) * rise;
%         
%         facil_F(1,N,R,1) = (facil_F(1,N-1,R,1) + ha) * exp(-delT / Ta);
% %         facil_F(1,N,R,2) = (facil_F(1,N-1,R,2) + hb) * delT^(-Tb);
%         facil_F(1,N,R,2) = (facil_F(1,N-1,R,2) + hb) * fall * (1-rise) ...
%                           + facil_F(1,N-1,R,2) * rise;
%         
% %         facil_F(1,N,R,2) = (facil_F(1,N-1,R,2) + hb) * exp(-delT / Tb);
%     end
%     
%     % Now, facilitate the probes after the ramps.
%     for P = 2:length(dat.prob)
%         delT = dat.prob(P);
%         rise = exp(-delT / Tr);
%         fall = delT^(-Xb); %exp(-delT / Tf);
% %         facil_F(P,:,R) = (facil_F(1,:,R) + h) * fall * (1 - rise) ...
% %                         + facil_F(1,:,R) * rise;
%         
%         facil_F(P,:,R,1) = (facil_F(1,:,R,1) + ha) * exp(-delT / Ta);
% %         facil_F(P,:,R,2) = (facil_F(1,:,R,2) + hb) * delT^(-Tb);
%         facil_F(P,:,R,2) = (facil_F(1,:,R,2) + hb) * fall * (1-rise) ...
%                           + facil_F(1,:,R,2) * rise;
%         
% %         facil_F(P,:,R,2) = (facil_F(1,:,R,2) + hb) * exp(-delT / Tb);
%     end
% end
% 
% % Combine the facilitation components.
% facil_F = (sum(facil_F, 4) + 1).^Fx;
% % facil_F = (facil_F + 1).^Fx;
% % facil_F(:,:,:,1) = (facil_F(:,:,:,1) + 1).^xa;
% % facil_F(:,:,:,2) = (facil_F(:,:,:,2) + 1).^xb;
% % facil_F = prod(facil_F, 4);
