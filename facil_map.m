% Jonathan Garcia
% 2015-11-24

function facil_F = facil_map(param, dat)
% Reshape facilitation values into a map for easy access.

% Facilitation is current parameter value relative to initial.
facil = param / param(1);

facil_F = zeros(length(dat.prob), length(dat.nrAP), length(dat.rISI));
c = 1;
for R = 1:length(dat.rISI)
    for N = 1:length(dat.nrAP)
        
        % Skip the cases that are redundant from (R, N) == (1, 1).
        if R > 1 && N == 1
            facil_F(:,N,R) = facil_F(:,1,1);
            continue
        end
        
        % Place facilitation value in its appropriate space.
        for P = 1:length(dat.prob)
            facil_F(P,N,R) = facil(c);
            c = c + 1;
        end
        
    end
end

end



% for R = 1:length(dat.rISI)
% figure(3); %subplot(4, 1, R);
% hold all; plot(dat.prob(2:end), facil_V1(2:end,:,R) - 1);
% axis([1 1000, 0.1 1000]);
% set(gca, 'XScale', 'log', 'YScale', 'log');
% end
% 
% for R = 1:length(dat.rISI)
% figure(2); %subplot(4, 1, R);
% hold all; plot(facil_V1(1,:,R), facil_V1(2:end,:,R)-1);
% axis([1 1000, 0.1 1000]);
% set(gca, 'XScale', 'log', 'YScale', 'log');
% end
