% Jonathan Garcia
% 2019-04-21
% 2021-08-16
% Multiscale modeling of presynaptic dynamics from molecular to mesoscale
% https://doi.org/10.1101/2021.04.16.440229

Ca_color = [0.0 0.0 0.0];
Sr_color = [0.0 0.0 1.0]; S0_color = [0.5 0.5 0.8];
Ar_color = [1.0 0.0 0.0]; A0_color = [0.8 0.5 0.5];
Cr_color = [0.5 0.8 0.5];
color_grad = @(c1, c2, N) sqrt([linspace(c1(1)^2, c2(1)^2, N); ...
                                linspace(c1(2)^2, c2(2)^2, N); ...
                                linspace(c1(3)^2, c2(3)^2, N)]');


%% Fig 1. SNARE Complex Structure and Dynamics.
% (A) SNAREpins prior to vesicle fusion.
% (B) Binding of Ca2+ to synaptotagmin (Syt-1 here, Syt-7 attaches to
% target membrane [12, 13]) triggers full zippering of SNARE complex and,
% in turn, vesicle fusion [14, 15].


%% Fig 2. Spatial Modeling of Spike-Evoked Ca2+ Transients.
% (A) Ca2+ sensors (dark yellow through dark blue filled circles) at vesicle
% cluster centers, displaced linearly from cluster of Ca2+ channels (blue
% half-disk on the left); distance in ?m, dn=0.160+0.105n for n?{0,…,16}.
% (B) [Ca2+]i measured over time in MCell (dark yellow through dark blue)
% and in the deterministic well-mixed model (maroon). MCell traces averaged
% from 2000 trials of MCell simulations with ?t=0.1 ms. Color transitions
% from yellow for vesicles proximal to the VDCC Ca2+ source to blue for
% vesicles far away, as in A. Proximally (distally) measured [Ca2+]i
% displays more (fewer) components of decay than are evident in the
% deterministic model.
% (C) Logarithmic plots of peak [Ca2+]i (blue) and peak time (red) as a
% function of distance from Ca2+ source; peak [Ca2+]i drops off
% exponentially with distance from VDCC cluster; amplitude of latent Ca2+
% dominates over the initial action-potential-evoked influx after 1.4 ?m.

load('data_files/mixed_Ca.mat');
load('data_files/flux_array_Ca.mat');

new_time = logspace(log10(0.1), log10(450), 1000);
new_Ca = interp1(time{1}, Ca_avg{1}, new_time);
kern = exp(-0.5 * (-25:25).^2 / 5^2);
kern = kern' / sum(kern);
new_Ca = conv2(new_Ca - Ca_0, kern, 'same') + Ca_0;

% (A) Ca2+ sensors (dark yellow through dark blue filled circles) at
% vesicle cluster centers, displaced linearly from cluster of Ca2+ channels
% (blue half-disk on the left); distance in ?m, d_n=0.160+0.105n for
% n?{0,…,16}.

% Plot calcium transients in response to a single action potential,
% coming from different distances from the VDCC cluster.

% Find VDCC locations.
rC =  0.050;    % 50 nm - VDCC cluster radius
nC = 50;        % number of VDCC channels in the cluster
RC = rC * sqrt(rand(1, nC));
TC = pi * (rand(1, nC) - 0.5);
xC = RC .* cos(TC);
yC = RC .* sin(TC);

% Get vesicle shape.
rV = 0.0175;    % 35-nm diameters
tV = linspace(0, 2*pi, 96);
xV = rV * cos(tV);
yV = rV * sin(tV);

% relative arrangement of vesicles within a cluster
xVC = rV * [0, -2, 2,       -1,      -1,        1,       1];
yVC = rV * [0,  0, 0, -sqrt(3), sqrt(3), -sqrt(3), sqrt(3)];

figure; fig_sz = [5 4];
set(gcf, 'Units', 'inches', 'PaperUnits', 'inches');
set(gcf, 'Position', [1 1, fig_sz], 'PaperPosition', [1 1, fig_sz]);

hv = axes('Units', 'inches', 'Position', [0.1+(5-4)/2 2.2+(2-1)/2, 4 1]);
hold all;

% Draw rectangle of the axon and the VDCC cluster on one side.
fill([0 2 2 0 0], [-0.25 -0.25 0.25 0.25 -0.25], [0.9 0.8 0.7]);
plot(xC, yC, '.', 'Color', [0 0.3 0.7]);

% Now draw all the vesicles.
cmap = parula(2*length(locs));
for j = 1:length(locs)
    % Plot all 7 vesicles in a cluster.
    fill(xV + xVC(1) + locs(j)+1, yV + yVC(1), ...
        0.7*cmap(2*(length(locs)-j)+1,:), 'EdgeColor', ...
        0.3*cmap(2*(length(locs)-j)+1,:));
    for v = 2:7
        fill(xV + xVC(v) + locs(j)+1, yV + yVC(v), ...
            [0.7 0.7 0.7], 'EdgeColor', [0.3 0.3 0.3]);
    end
end
axis([0 2, -0.25 0.25]);
xlabel('Vesicle Distance from VDCC Cluster (\mum)');
ylabel('\mum');

set(hv, 'XTick', [0 locs(1:2:end)+1 2], 'YTick', 0.2*[-1 0 1], ...
    'FontSize', 10);

% (B) Profiles of [Ca2+]i averaged from 2000 trials of MCell simulations
% with ?t=0.1 ms; color corresponds to distance from VDCC cluster (yellow
% proximal to blue distal).

% Plot calcium traces for all distances over 10 msec.
hc = subplot(223); hold all;
cmap = parula(2*length(locs));
samp = [1,5,9,13,17];
plts = zeros(length(samp), 1);
lbls = cell(length(samp), 1);
for j = 1:length(samp)
%     plts(j) = plot(new_time, new_Ca(:,samp(j)), ...
%         'LineWidth', 1.5, 'Color', 0.7*cmap(2*(length(locs)-samp(j))+1,:));
    plts(j) = plot(time{1}, Ca_avg{1}(:,samp(j)), ...
        'LineWidth', 1.5, 'Color', 0.7*cmap(2*(length(locs)-samp(j))+1,:));
    lbls{j} = sprintf(' %.3f \\mum', 1+locs(samp(j)));
end
wm = plot(t, Ca_Cb45, 'LineWidth', 2.5, 'Color', [0.5 0 0]);
axis([0.7 450, 3e-2 1e2]);
xlabel('Time (ms)'); ylabel('[Ca^{2+}]_i (\muM)');
set(hc, 'XScale', 'log', 'YScale', 'log', 'FontSize', 10, ...
    'XTick', [0.1 0.3 1 3 10 30 100 300], 'YTick', 10.^(-2:2), ...
    'Position', [0.13 0.12, 0.37 0.34]);
legend([plts; wm], [lbls; ' well-mixed'], ...
    'FontSize', 8, 'Box', 'off', 'Position', [0.3 0.25, 0.2 0.2]);

% (C) Logarithmic plots of peak [Ca2+]i (blue) and peak time (red) as a
% function of distance from Ca2+ source; latent Ca2+ component dominates
% over the initial action-potential-evoked influx after 1.4 ?m.

% Plot peak [Ca2+]i as a function of distance from the VDCC cluster.
ax_pos = [0.67 0.13, 0.22 0.34];
[peakC, it] = max(new_Ca, [], 1);  % peak [Ca]
[peakCM, itm] = max(Ca_Cb45);
peakT = new_time(it);              % peak time
peakTM = t(itm);

% Plot peak [Ca] as a function of distance.
hd = axes('Position', ax_pos); hold all;
plot(locs+1, peakC, '.-', 'Color', [0 0 0.7], ...
    'LineWidth', 1, 'MarkerSize', 10);
plot([0 2], Ca_0 * [1 1], '--', 'Color', [0 0 0.7], 'LineWidth', 1);
plot([0 2], peakCM * [1 1], ':', 'Color', [0 0 0.7], 'LineWidth', 2);
axis([0 2, 0.5e-1 1e2]);
xlabel('Distance (\mum)'); ylabel('Peak [Ca^{2+}]_i (\muM)');

set(hd, 'XTick', 1+locs([1,9,17]), 'YTick', 10.^(-1:2), ...
    'YScale', 'log', 'FontSize', 10, 'YColor', [0 0 0.7]);

% Display length constant.
LC = log(peakC(1:10)) / [ones(1,10); 1+locs(1:10)];
LC = -1 / LC(2);
disp(LC);

% Plot peak time as a function of distance.
ht = axes('Position', ax_pos); hold all;
plot(locs+1, peakT, '.-', 'Color', [0.7 0 0], ...
    'LineWidth', 1, 'MarkerSize', 10);
plot([0 2], peakTM * [1 1], ':', 'Color', [0.7 0 0], 'LineWidth', 2);
axis([0 2, 0.5 1e2]);
ylabel('Peak Time (ms)');

set(ht, 'XTick', [], 'YTick', 10.^(0:2), 'YScale', 'log', 'FontSize', 10, ...
    'Color', 'none', 'Box', 'off', ...
    'YAxisLocation', 'right', 'YColor', [0.7 0 0]);


%% Fig 3. Model of Ca2+-Evoked, Synaptotagmin-Mediated Neurotransmitter Release.
% (A, B) Model adapted from Sun et al. [53]. ?S and ?A represent rates of
% vesicle fusion from the releasable states of the synchronous and
% asynchronous mechanisms, respectively.
% (A) Ca2+-bound states for Syt-1 (synchronous release); Sn indicates n
% Ca2+ ions bound to the synchronous release mechanism.
% (B) Ca2+-bound states for Syt-7 (asynchronous release); An indicates n
% Ca2+ ions bound to the asynchronous release mechanism.
% (C, D) Action-potential-like stimulus delivered to model axon starting at
% 0 ms. Diffusion is assumed to be instantaneous, and molecular state
% probabilities are tracked deterministically over time.
% (C) Free [Ca2+]i in response to single action potential.
% (D) Instantaneous vesicle release rate in response to buffered Ca2+ from
% both synaptotagmin-mediated release mechanisms.

% (A, B) Model adapted from Sun et al. [53]. ?S and ?A represent rates of
% vesicle fusion from the releasable states of the synchronous and
% asynchronous mechanisms, respectively.

% Load SNARE binding/unbinding rates and state probabilities,
% and convert them into normalized values.
load('data_files/SNARE_states.mat');
to_range = @(data, mxD, p, mx) round(10*(data / mxD).^p * mx) / 10;

% normalizing linewidths for state transitions
mnDT = min([S_n, S0p, S1p, A_n, A0p, A1p]);
mxDT = max([S_n, S0p, S1p, A_n, A0p, A1p]);
mxT  =  8;  pT = log(0.3/ mxT) / log(mnDT/mxDT);
mxAT = 16; pAT = log(3.0/mxAT) / log(mnDT/mxDT);

% normalizing linewidths for release transitions
mnDR = min([S0rel, S1rel, A0rel, A1rel]);
mxDR = max([S0rel, S1rel, A0rel, A1rel]);
mxR  =  4;  pR = log(0.3/ mxR) / log(mnDR/mxDR);
mxAR =  8; pAR = log(3.0/mxAR) / log(mnDR/mxDR);

font_max = 12;
font_med = 10;
font_min =  8;

% data for drawing state circles
radius = 0.17;                  % inches
rel_xr = 0.35; rel_yr = 0.15;   % inches
Drel   = 0.35;                  % inches

tt = linspace(0, 2*pi, 1+4*20);
xs = radius*cos(tt); ys = radius*sin(tt);   % state circle
xr = rel_xr*cos(tt); yr = rel_yr*sin(tt);   % release ellipse

% Determine where things go within the figure.
AB_sz = [5.9, 1.25];                % inches
fig_sz = [AB_sz(1), 2*AB_sz(2)];    % inches

L_buf = 0.2; R_buf = 0.2;
xcA = L_buf + radius;                               % start
xcB = AB_sz(1) - 2*rel_xr - Drel - radius - R_buf;  % end
xcR = AB_sz(1) - rel_xr - R_buf; yc = AB_sz(2) / 2; % positions

% Draw a figure that show both state diagrams together, with labels
% on the arrows to represent how transition rates are calculated.
figure;
set(gcf, 'Units', 'inches', 'PaperUnits', 'inches');
set(gcf, 'Position', [1 1, fig_sz], 'PaperPosition', [1 1, fig_sz]);
for diag = 1:2

% Drawing parameters depend on which diagram.
switch diag
    % (A) Ca2+-bound states for Syt-1 (synchronous release); Sn indicates n
    % Ca2+ ions bound to the synchronous release mechanism.
    case 1
        xc = linspace(xcA, xcB, 6); D_xc = xc(2) - xc(1) - 2*radius;
        color = color_grad(S0_color, Cr_color, 6);
        type = 'S'; NC = 5;
        neg_L = to_range(S_n, mxDT, pT, mxT);   % unbinding linewidths
        neg_A = to_range(S_n, mxDT, pAT, mxAT); % unbinding arrowheads
        
        pos_L = to_range(S0p, mxDT, pT, mxT);   % binding linewidths 0
        pos_A = to_range(S0p, mxDT, pAT, mxAT); % binding arrowheads 0
        
        rel_L = to_range(S0rel, mxDR, pR, mxR);     % release linewidth 0
        rel_A = to_range(S0rel, mxDR, pAR, mxAR);   % release arrowhead 0
        
    % (B) Ca2+-bound states for Syt-7 (asynchronous release); An indicates n
    % Ca2+ ions bound to the asynchronous release mechanism.
    case 2
        xc = linspace(xcA, xcB, 3); D_xc = xc(2) - xc(1) - 2*radius;
        color = color_grad(A0_color, Cr_color, 3);
        type = 'A'; NC = 2;
        neg_L = to_range(A_n, mxDT, pT, mxT);   % unbinding linewidths
        neg_A = to_range(A_n, mxDT, pAT, mxAT); % unbinding arrowheads
        
        pos_L = to_range(A0p, mxDT, pT, mxT);   % binding linewidths 0
        pos_A = to_range(A0p, mxDT, pAT, mxAT); % binding arrowheads 0
        
        rel_L = to_range(A0rel, mxDR, pR, mxR);     % release linewidth 0
        rel_A = to_range(A0rel, mxDR, pAR, mxAR);   % release arrowhead 0
end

% Set up subplot axis for this state transition diagram.
base_y = (2-diag)*AB_sz(2);
ax_pos = [0 base_y, AB_sz] ./ [fig_sz, fig_sz];
ax = axes('Position', ax_pos); hold all;

% Draw the circles along the state diagram.
for j = 1:length(xc)
    fill(xc(j) + xs, yc + ys, color(j,:));
    ht = text(xc(j), yc, sprintf('%s_%d', type, j-1));
    set(ht, 'FontSize', font_max, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
end

% Draw the arrows between the states.
for j = 1:length(pos_L)
    % Make positive and negative arrows.
    arP = annotation('arrow'); arN = annotation('arrow');
    set([arP arN], 'HeadStyle', 'vback1', 'Units', 'inches');
    
    % Ca-binding (positive) arrow size
    arP.LineWidth  = pos_L(j);
    arP.HeadWidth  = pos_A(j);
    arP.HeadLength = pos_A(j);
    
    % Ca-unbinding (negative) arrow size
    arN.LineWidth  = neg_L(j);
    arN.HeadWidth  = neg_A(j);
    arN.HeadLength = neg_A(j);
    
    xL = xc( j ) + radius;
    xR = xc(j+1) - radius;
    
    arP.Position = [xL yc+0.4*radius+base_y,  D_xc 0];
    arN.Position = [xR yc-0.4*radius+base_y, -D_xc 0];
    
    % Add labels.
    p_lbl = text(xL+D_xc/2, yc+radius, ...
        sprintf('%d\\cdotk_{%s+}\\cdot[Ca^{2+}]_i', NC-j+1, type), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
        'FontSize', font_min);
    n_lbl = text(xR-D_xc/2, yc-radius, ...
        sprintf('%d\\cdotb^%d\\cdotk_{%s-}', j, j-1, type), ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
        'FontSize', font_min);
end

% Draw release state.
LR = 1;
xrc = xc(end) + (1+1.5)*radius + LR;
fill(xcR + xr, yc + yr, [0.9 0.9 0.9]);
ht = text(xcR, yc, 'release');
set(ht, 'FontSize', font_med, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

% Draw release arrow.
arR = annotation('arrow');
set(arR, 'HeadStyle', 'vback1', 'Units', 'inches');
arR.LineWidth  = rel_L;
arR.HeadWidth  = rel_A;
arR.HeadLength = rel_A;
arR.Position = [xc(end)+radius yc+base_y, Drel 0];

% Add label.
r_lbl = text(xc(end)+radius+Drel/2, yc+0.6*radius, ...
    sprintf('\\gamma_%s', type), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
    'FontSize', font_min, 'FontWeight', 'bold');

% After plots, size the axes and make them invisible.
axis(ax, [0 AB_sz(1), 0 AB_sz(2)]); axis off;

end

% (C, D) Action-potential-like stimulus delivered to model axon starting at
% 0 ms. Diffusion is assumed to be instantaneous, and molecular state
% probabilities are tracked deterministically over time.

load('data_files/mixed_Ca.mat');
load('data_files/mixed_runs.mat');
% dt = 0.1;       % ms
% Ca_0 = 0.100;   % uM

figure; set(gcf, 'Units', 'inches', 'Position', [1 1, 6 2]);

% (C) Free [Ca2+]i in response to single action potential.

sp3 = subplot(121); %axes('Position', [0.19 0.32, 0.75 0.16]);
plot(sp3, t, Ca_Cb45, 'Color', Ca_color, 'LineWidth', 2);
xlabel(sp3, 'Time (ms)');
ylabel(sp3, '[Ca^{2+}]_i (\muM)');
axis(sp3, [0 120, 0 5]);
set(sp3, 'XTick', 0:20:120, 'YTick', 0:5, ...
    'Box', 'off', 'FontSize', 10);

% (D) Instantaneous vesicle release rate in response to buffered Ca2+ from
% both synaptotagmin-mediated release mechanisms.

sp4 = subplot(122); %axes('Position', [0.19 0.08, 0.75 0.16]);
plot(sp4, t, Sr_Cb45+Ar_Cb45, 'Color', Cr_color, 'LineWidth', 2);
xlabel(sp4, 'Time (ms)');
ylabel(sp4, sprintf('Release\nRate (ms^{-1})'));
axis(sp4, [0 120, 0 3e-4]);
set(sp4, 'XTick', 0:20:120, 'YTick', 0:1e-4:3e-4, ...
    'Box', 'off', 'FontSize', 10);


%% Fig 4. Synchronous and Asynchronous Release in MCell.
% Color indicates distance from VDCC source, with yellow representing a
% nearby Ca2+ sensor and dark blue a distant one (as in Fig 2 A-B).
% Action-potential-like stimulus delivered at 0 ms (left), followed by
% another at 20 ms (center) and 100 ms (right).
% (A) Spike-evoked Ca2+ traces that drive release.
% (B) Synchronous release raster.
% (C) Asynchronous release raster.
% (D) Synchronous (tall, thin bars) and asynchronous (short, wide bars)
% release stacked histogram. Most synchronous releases happen close to the
% Ca2+ source; asynchronous releases distributed across all distances.

% Plot release histograms resulting from calcium traces at different PPFs.

load('data_files/flux_array_Ca.mat');
load('data_files/flux_array_hists.mat');

figure; fig_sz = [6 6];
set(gcf, 'Units', 'inches', 'PaperUnits', 'inches');
set(gcf, 'Position', [1 1, fig_sz], 'PaperPosition', [1 1, fig_sz]);

cmap = parula(2*length(locs));
pp = [1 5 7]; axs = zeros(4, length(pp));
Lbuf = 0.12;  Rbuf = 0.02; Xbuf = 0.01;
Xsz = (1 - Lbuf - Rbuf - (length(pp)-1)*Xbuf) / length(pp);
Tbuf = 0.05; Bbuf = 0.1; Ybuf = 0.05;
Ysz = (1 - Tbuf - Bbuf - 3*Ybuf) / 4;
for ip = 1:length(pp)
    
    % (A) Spike-evoked Ca2+ traces that drive release.
    
    % Plot [Ca2](t).
    axs(1,ip) = axes('Position', ...
        [Lbuf+(ip-1)*(Xsz+Xbuf) 1-Tbuf-1*Ysz-0*Ybuf, Xsz Ysz]);
    hold all;
    for j = 1:length(locs)
        plot(time{pp(ip)}, Ca_avg{pp(ip)}(:,j), ...
            'LineWidth', 1, 'Color', 0.7*cmap(2*(length(locs)-j)+1,:));
    end
    axis([0 250, 0.5e-1 1e2]);
    set(gca, 'YScale', 'log', 'YTick', [], 'XTick', 0:100:250);
    if ip == 1
        ylabel('[Ca^{2+}]_i (\muM)');
        set(gca, 'YTick', 10.^(-1:2));
    end
    
    % (B) Synchronous release raster.
    
    % Plot synchronous release raster.
    axs(2,ip) = axes('Position', ...
        [Lbuf+(ip-1)*(Xsz+Xbuf) 1-Tbuf-2*Ysz-1*Ybuf, Xsz Ysz]);
    hold all;
    for j = 1:length(locs)
        si = (Sh{pp(ip)}(:,2) == j) & (Sh{pp(ip)}(:,1) <= 250);
        plot(Sh{pp(ip)}(si,1), Sh{pp(ip)}(si,3), '.', ...
            'Color', 0.7*cmap(2*(length(locs)-j)+1,:), 'MarkerSize', 8);
    end
    axis([0 250, 0 2001]);
    set(gca, 'YTick', [], 'XTick', 0:100:250);
    if ip == 1
        ylabel('Trial #');
        set(gca, 'YTick', [1 1000 2000], 'YTickLabel', {'1', '1k', '2k'});
    end
    
    % (C) Asynchronous release raster.
    
    % Plot asynchronous release raster.
    axs(3,ip) = axes('Position', ...
        [Lbuf+(ip-1)*(Xsz+Xbuf) 1-Tbuf-3*Ysz-2*Ybuf, Xsz Ysz]);
    hold all;
    for j = 1:length(locs)
        ai = (Ah{pp(ip)}(:,2) == j) & (Ah{pp(ip)}(:,1) <= 250);
        plot(Ah{pp(ip)}(ai,1), Ah{pp(ip)}(ai,3), '.', ...
            'Color', 0.7*cmap(2*(length(locs)-j)+1,:), 'MarkerSize', 8);
    end
    axis([0 250, 0 2001]);
    set(gca, 'YTick', [], 'XTick', 0:100:250);
    if ip == 1
        ylabel('Trial #');
        set(gca, 'YTick', [1 1000 2000], 'YTickLabel', {'1', '1k', '2k'});
    end
    
    % (D) Synchronous (tall, thin bars) and asynchronous (short, wide bars)
    % release stacked histogram. Most synchronous releases happen close to the
    % Ca2+ source; asynchronous releases distributed across all distances.
    
    % Plot combined release histogram.
    axs(4,ip) = axes('Position', ...
        [Lbuf+(ip-1)*(Xsz+Xbuf) 1-Tbuf-4*Ysz-3*Ybuf, Xsz Ysz]);
    hold all;
    Sbar = bar(Sht, 1000*S_hist{pp(ip)}, ...
        1, 'stacked', 'EdgeColor', 'none');
    Abar = bar(Aht, 1000*A_hist{pp(ip)}, ...
        1, 'stacked', 'EdgeColor', 'none');
    axis([0 250, 0.03 300]);
    xlabel('Time (ms)');
    set(gca, 'YTick', [], 'XTick', 0:100:250, 'YScale', 'log');
    if ip == 1
        ylabel('Releases / sec');
        set(gca, 'YTick', [0.1, 1, 10, 100]);
    end
    
    % Fix histogram colors.
    for j = 1:length(locs)
        set([Sbar(j) Abar(j)], ...
            'FaceColor', 0.7*cmap(2*(length(locs)-j)+1,:));
    end
end


%% Fig 5. Multi-Exponential Shape of Ca2+-Driven Vesicle Release Rate.
% (A, B) Plots given as semi-log to highlight exponential decay components
% (straight line segments of profiles).
% (A) A single, spike-evoked [Ca2+]i transient, which drives
% (B) the synchronous and asynchronous release rates.
% (C) Instantaneous time constants for Ca2+, synchronous, and asynchronous
% curves, calculated from the well-mixed model (see Eq (2)). Long release
% rate time constants (around 80 ms and 1000 ms; dashed lines) follow Ca2+
% curve due to slow un-buffering of latent Ca2+. Asynchronous starts high
% because fast and slow components have comparable magnitude and become
% conflated; it goes up to infinity where additive effects cause the curve
% to flatten.

load('data_files/flux_array_Ca.mat');
load('data_files/flux_array_runs.mat');
load('data_files/mixed_Ca.mat');
load('data_files/mixed_runs.mat');

figure; set(gcf, 'Units', 'inches', 'Position', [1 1, 6 4]);

% (A) A single, spike-evoked [Ca2+]i transient, which drives

subplot(221); hold all;
plot(time{1}, mean(Ca_avg{1}(:,3:4), 2), ...
    'LineWidth', 3, 'Color', Ca_color);
axis([0 100, 5e-2 2e1]);
xlabel('Time (ms)'); ylabel('[Ca^{2+}]_i (\muM)');
set(gca, 'YScale', 'log', 'YTick', 10.^(-2:2), 'XTick', 0:20:100, ...
    'FontSize', 10, 'Position', [0.12 0.62, 0.35 0.35]);
legend({' calcium transient'}, 'Box', 'off', 'FontSize', 8, ...
    'Location', 'NorthEast');

% (B) the synchronous and asynchronous release rates.

subplot(223); hold all;
plot(time{1}, mean(S_rel{1}(:,3:4), 2), ...
    'LineWidth', 3, 'Color', Sr_color);
plot(time{1}, mean(A_rel{1}(:,3:4), 2), ...
    'LineWidth', 3, 'Color', Ar_color);
axis([0 100, 1e-9 0.1]);
xlabel('Time (ms)'); ylabel('Release Rate (ms^{-1})');
set(gca, 'YScale', 'log', 'YTick', 10.^(-8:2:0), 'XTick', 0:20:100, ...
    'FontSize', 10, 'Position', [0.12 0.12, 0.35 0.35]);
legend({' synchronous', ' asynchronous'}, 'Box', 'off', 'FontSize', 8, ...
    'Location', 'NorthEast');

% (C) Instantaneous time constants for Ca2+, synchronous, and asynchronous
% curves, calculated from the well-mixed model (see Eq (2)). Long release
% rate time constants (around 80 ms and 1000 ms; dashed lines) follow Ca2+
% curve due to slow un-buffering of latent Ca2+. Asynchronous starts high
% because fast and slow components have comparable magnitude and become
% conflated; it goes up to infinity where additive effects cause the curve
% to flatten.

dlX = @(t, X, X0) (log(X(3:end-0)-X0) - log(X(1:end-2)-X0)) ...
    ./ (t(3:end-0) - t(1:end-2));
TX = @(t, X, X0) -1 ./ dlX(t, X, X0);

warning('off', 'MATLAB:plot:IgnoreImaginaryXYPart');
warning('off', 'MATLAB:Axes:NegativeDataInLogAxis');

tmh = subplot(122); hold all;
plot(t(2:end-1), TX(t, Ca_Cb45, Ca_0), ...
    'Color', Ca_color, 'LineWidth', 2);
plot(t(2:end-1), TX(t, Sr_Cb45, sr0), ...
    'Color', Sr_color, 'LineWidth', 2);
plot(t(2:end-1), TX(t, Ar_Cb45, ar0), ...
    'Color', Ar_color, 'LineWidth', 2);

text(800, 1.2*80, '\tau = 80 ms', 'FontSize', 10, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
plot([1 1e4], [80 80], '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);

text(800, 1.2*1000, '\tau = 1000 ms', 'FontSize', 10, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
plot([1 1e4], [1000 1000], '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);

xlabel('Time (ms)');
ylabel('(\tau) Time Constant (ms)');
axis([1 1e4, 1e-1 1e4]);
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XTick', 10.^(0:4), 'YTick', 10.^(-1:4), 'FontSize', 10, ...
    'Position', [0.61 0.12, 0.35 0.86]);
legend({'[Ca^{2+}]_i time constant', ...
    'synchronous release', 'asynchronous release'}, ...
    'Box', 'off', 'Position', [0.75 0.25, 0.2 0.15], 'FontSize', 8);


%% Fig 6. Synchronous and Asynchronous Release Rates in Response to Ca2+
% Impulse at Different Resting Concentrations.
% Instantaneous impulse of Ca2+ delivered at 10 ms. Solid lines represent
% true release rate; dotted lines have spontaneous rates subtracted off to
% show secondary exponential components. Black lines show release rate
% decaying with a single exponential component with no baseline [Ca2+]i.
% For other curves, [Ca2+]i0 ranges from 0.001 ?M to 10 ?M.
% (A) Synchronous release rate over time: S(t).
% (B) Asynchronous release rate over time: A(t).
% (C) Instantaneous release rate decay time constants for synchronous and
% asynchronous mechanisms. Fast components (lower blue and red lines)
% determined from profiles with [Ca2+]i0=0 (black lines in A and B).
% Slower components (upper blue and red curves) determined from cases with
% small [Ca2+]i0.

load('data_files/impulse_response.mat');
load('data_files/SNARE_0.mat');

% (A) Synchronous release rate over time: S(t).
% (B) Asynchronous release rate over time: A(t).

% Plot the impulse response of release mechanisms for different Ca levels.
frac = sqrt(linspace(0, 0.75, length(Ca0)-2));

figure;
set(gcf, 'Units', 'inches', 'PaperUnits', 'inches');
set(gcf, 'Position', [1 1, 6 2.5], 'PaperPosition', [1 1, 6 2.5]);

% Synchronous Release
subplot(121); hold all;
for j = length(Ca0)-2:-2:2
    color = Sr_color * (1-frac(j)) + [1 1 1] * frac(j);  % white -> blue
    p1 = plot(tm, Sr(j,:), '-', 'Color', color, 'LineWidth', 2);
    pd = plot(tm, Sr(j,:)-Sr0(j), '-.', 'Color', color, 'LineWidth', 2);
end
p0 = plot(tm, Sr(1,:), 'k-', 'LineWidth', 2);
axis([0 250, 1e-25, 1]);
xlabel('Time (ms)'); ylabel('S(t) (ms^{-1})');
set(gca, 'XScale', 'lin', 'YScale', 'log', 'FontSize', 10);
set(gca, 'XTick', 0:50:250, 'YTick', 10.^(-25:5:0));
legend([p1, pd, p0], {' [Ca^{2+}]_{i0} > 0', ...
    ' S(t) - S(0)', ' [Ca^{2+}]_{i0} = 0'}, ...
    'Units', 'normalized', 'Position', [0.14 0.24, 0.3 0.15], ...
    'Box', 'off', 'FontSize', 8);

% Asynchronous Release
subplot(122); hold all;
for j = length(Ca0)-2:-2:2
    color = Ar_color * (1-frac(j)) + [1 1 1] * frac(j);  % white -> red
    p1 = plot(tm, Ar(j,:), '-', 'Color', color, 'LineWidth', 2);
    pd = plot(tm, Ar(j,:)-Ar0(j), '-.', 'Color', color, 'LineWidth', 2);
end
p0 = plot(tm, Ar(1,:), 'k-', 'LineWidth', 2);
axis([0 500, 1e-25, 1]);
xlabel('Time (ms)'); ylabel('A(t) (ms^{-1})');
set(gca, 'XScale', 'lin', 'YScale', 'log', 'FontSize', 10);
set(gca, 'XTick', 0:100:500, 'YTick', 10.^(-25:5:0));
legend([p1, pd, p0], {' [Ca^{2+}]_{i0} > 0', ...
	' A(t) - A(0)', ' [Ca^{2+}]_{i0} = 0'}, ...
    'Units', 'normalized', 'Position', [0.58 0.24, 0.3 0.15], ...
    'Box', 'off', 'FontSize', 8);

% (C) Instantaneous release rate decay time constants for synchronous and
% asynchronous mechanisms. Fast components (lower blue and red lines)
% determined from profiles with [Ca2+]i0=0 (black lines in A and B).
% Slower components (upper blue and red curves) determined from cases with
% small [Ca2+]i0.

% Uncover release rate decay time constants with and without baseline Ca.

% Show instantaneous time constants.
dlX = @(t, X, X0) (log(X(3:end-0)-X0) - log(X(1:end-2)-X0)) ...
    ./ (t(3:end-0) - t(1:end-2));
TX = @(t, X, X0) -1 ./ dlX(t, X, X0);

figure; hold all;
set(gcf, 'Units', 'inches', 'PaperUnits', 'inches');
set(gcf, 'Position', [1 1, 6 2.5], 'PaperPosition', [1 1, 6 2.5]);

% Gather data for all cases.
ix = round([50 150 300 450]/dt); id = round([125 350 499 499]/dt);
tx = {'\\tau = %.3f ms', '\\tau = %.2f ms', ...
      '\\tau = %.1f ms', '\\tau = %.1f ms'};
ha = {'left', 'left', 'right', 'right'};
va = {'bottom', 'top', 'bottom', 'top'}; vs = [1.25, 0.8, 1.25, 0.8];
TT = [TX(tm, Sr(1,:), Sr0(1)); TX(tm, Sr(2,:), Sr0(2)); ...
      TX(tm, Ar(1,:), Ar0(1)); TX(tm, Ar(2,:), Ar0(2))];
cl = 'bbrr';

% Plot all four cases (synchronous/asynchronous with no/some baseline Ca).
for j = 1:2
    plot(tm(2:id(j)+1), TT(j,1:id(j)), 'Color', Sr_color, 'LineWidth', 2);
    plot([id(j)*dt 500], TT(j,ix(j))*[1 1], [cl(j) ':'], 'LineWidth', 2);
    plot(ix(j)*dt, TT(j,ix(j)), 'k.', 'MarkerSize', 25);
    text(ix(j)*dt, TT(j,ix(j))*vs(j), sprintf(tx{j}, TT(j,ix(j))), ...
        'HorizontalAlignment', ha{j}, 'VerticalAlignment', va{j}, ...
        'FontSize', 10);
end
for j = 3:4
    plot(tm(2:id(j)+1), TT(j,1:id(j)), 'Color', Ar_color, 'LineWidth', 2);
    % plot([id(j)*dt 500], TT(j,ix(j))*[1 1], [cl(j) ':'], 'LineWidth', 2);
    plot(ix(j)*dt, TT(j,ix(j)), 'k.', 'MarkerSize', 25);
    text(ix(j)*dt, TT(j,ix(j))*vs(j), sprintf(tx{j}, TT(j,ix(j))), ...
        'HorizontalAlignment', ha{j}, 'VerticalAlignment', va{j}, ...
        'FontSize', 10);
end
axis([0 500, 0.1 100]);
xlabel('Time (ms)'); ylabel('Time Constant (ms)');
set(gca, 'XScale', 'lin', 'YScale', 'log', 'FontSize', 10);
set(gca, 'XTick', 0:100:500, 'YTick', 10.^(-1:2));

text(100, 1.2, 'synchronous', 'Color', Sr_color, 'FontSize', 10);
text(10,  50, 'asynchronous', 'Color', Ar_color, 'FontSize', 10);


%% Fig 7. Fitted Release Rate Histogram Profiles for a Single Spike.
% Parameter values given in Table 2.
% (A) Synchronous release rate: true histogram (blue) with estimated
% histogram (black).
% (B) Asynchronous release rate: true histogram (red) with estimated
% histogram (black).

load('data_files/RNP_Ca.mat');
load('data_files/RNP_runs.mat');
load('data_files/RNP_Sv_Av.mat');
load('data_files/SNARE_0.mat');
sr0 = Sr0(6); ar0 = Ar0(6);

[~, S_Fs] = fit_conv_profiles(time{1,1}{1}(:), S_rel{1,1}{1}(:), ...
    sr0, 0.4, S400, [], false);
[~, A_Fs] = fit_conv_profiles(time{1,1}{1}(:), A_rel{1,1}{1}(:), ...
    ar0, 0.4, A400, [], false);

figure; set(gcf, 'Units', 'inches', 'Position', [1 1, 6 3]);

% (A) Synchronous release rate: true histogram (blue) with estimated
% histogram (black).

subplot(121); hold all;
hsr = plot(time{1,1}{1}, S_rel{1,1}{1}, ...
    'LineWidth', 3, 'Color', 0.5*Sr_color+0.5*[1 1 1]);
hsf = plot(time{1,1}{1}, S_Fs, '-', 'LineWidth', 1, 'Color', [0 0 0]);
axis([0 100, 1e-9 2]);
xlabel('Time (ms)'); ylabel('Release Rate (ms^{-1})');
set(gca, 'YScale', 'log', 'YTick', 10.^(-8:2:0), 'XTick', 0:20:100, ...
    'FontSize', 10);
legend([hsr(1), hsf(1)], {' synchronous release', ' fitted profiles'}, ...
    'FontSize', 8, 'Box', 'off', 'Position', [0.2 0.7, 0.27 0.2]);

% (B) Asynchronous release rate: true histogram (red) with estimated
% histogram (black).

subplot(122); hold all;
har = plot(time{1,1}{1}, A_rel{1,1}{1}, ...
    'LineWidth', 3, 'Color', 0.5*Ar_color+0.5*[1 1 1]);
haf = plot(time{1,1}{1}, A_Fs, '-', 'LineWidth', 1, 'Color', [0 0 0]);
axis([0 100, 5e-6 1e-2]);
xlabel('Time (ms)'); ylabel('Release Rate (ms^{-1})');
set(gca, 'YScale', 'log', 'YTick', 10.^(-5:-1), 'XTick', 0:20:100, ...
    'FontSize', 10);
legend([har(1), haf(1)], {' asynchronous release', ' fitted profiles'}, ...
    'FontSize', 8, 'Box', 'off', 'Position', [0.65 0.7, 0.27 0.2]);


%% Fig 8. Convolutional Filter Applied to a Component of a Release Rate Function.
% Toy model with P=5, ?=10 ms, k=0.5 ms–1, ?=5 ms, and ?=1 ms.
% (A) Unfiltered release rate component.
% (B) MCMC ex-Gaussian filter shape.
% (C) Filtered release profile produced by convolving the release rate
% profile with the temporal delay filter.
% (D-F) Release rates in response to spike trains without applying delay filter.
% (G-I) Release rates in response to spike trains with delay filter applied.
% (D,G) Response to one spike.
% (E,H) Response to two spikes.
% (F,I) Response to multiple spikes. Dotted lines show how the histogram of
% response to one spike falls off with interference from the response to
% the following spike. Spike times at 0, 15, 20, 30, and 50 ms.

% exponential rate of release convolved with
% ex-Gaussian distribution of start times
XXG = @(t, mu, sg, K, P, T) P*K/(K*T-1) * ...
    (exp(-(t-mu-sg^2/2/T)/T) .* normcdf(t, mu+sg^2/T, sg) - ...
     exp(-(t-mu-sg^2/2*K)*K) .* normcdf(t, mu+sg^2*K, sg));

% ex-Gaussian distribution kernel
XG = @(t, mu, sg, K) ...
    K * exp(-(t-mu-sg^2/2*K)*K) .* normcdf(t, mu+sg^2*K, sg);

% probability of start time having occurred already
CXG = @(t, mu, sg, K) normcdf(t, mu, sg) - ...
    exp(-(t-mu-sg^2/2*K)*K) .* normcdf(t, mu+sg^2*K, sg);

% % probability of no release having occurred yet from a given profile
% DXXG = @(t, mu, sg, K, P, T) exp(P * ...
%     (K*T/(K*T-1) * exp(-(t-mu-sg^2/2/T)/T) .* normcdf(t, mu+sg^2/T, sg) ...
%      - 1/(K*T-1) * exp(-(t-mu-sg^2/2*K)*K) .* normcdf(t, mu+sg^2*K, sg) ...
%      - normcdf(t, mu, sg)));

% Set up parameters for simulation.
P = 5; T = 10; mu = 5; sg = 1; K = 0.5;
spks = cumsum([0 15 5 10 20]);  % spike times

% Build histograms for measuring release rates.
dt = 0.1; tt = -20:dt:60;

figure;
set(gcf, 'Units', 'inches', 'Position', [1 1, 6 1.5]);

% (A) Unfiltered release rate component.

a1 = subplot(131);
plot(tt, P/T * exp(-tt/T) .* (tt>0), 'LineWidth', 2, 'Color', [0.7 0 0]);
axis([-20 60, 0 P/T+0.1]);
xlabel('Time (ms)');
ylabel(sprintf('Release (ms^{-1})'));
text(80, 0.3, '\ast', 'FontSize', 24, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

% (B) MCMC ex-Gaussian filter shape.

a2 = subplot(132);
plot(tt, XG(tt, mu, sg, K), 'LineWidth', 2, 'Color', [0 0 0.7]);
axis([-20 60, 0 P/T+0.1]);
xlabel('Time (ms)');
text(80, 0.3, '=', 'FontSize', 24, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

% (C) Filtered release profile produced by convolving the release rate
% profile with the temporal delay filter.

a3 = subplot(133);
plot(tt, XXG(tt, mu, sg, K, P, T), 'LineWidth', 2, 'Color', [0.5 0 0.5]);
axis([-20 60, 0 P/T+0.1]);
xlabel('Time (ms)');

set([a1 a2 a3], 'XTick', -20:20:60, 'YTick', 0:0.2:0.6, ...
    'Box', 'off', 'FontSize', 10);
set([a2 a3], 'YTick', []);
drawnow;

x0 = get(a1, 'Position');
set(a2, 'Position', get(a2, 'Position') + [0.03 0, 0 0]);
set(a3, 'Position', get(a3, 'Position') + [0.06 0, 0 0]);


% (D-F) Release rates in response to spike trains without applying delay filter.
% (G-I) Release rates in response to spike trains with delay filter applied.
% (D,G) Response to one spike.
% (E,H) Response to two spikes.
% (F,I) Response to multiple spikes. Dotted lines show how the histogram of
% response to one spike falls off with interference from the response to
% the following spike. Spike times at 0, 15, 20, 30, and 50 ms.

% Build histograms for measuring release rates.
dt = 0.1; TC = 100; tc = -20:dt:TC;

figure;
set(gcf, 'Units', 'inches', 'Position', [1 1, 6 4]);

ss = [1,2,length(spks)];
for c = 0:1
for j = 1:length(ss)
    subplot(length(ss), 2, 2*j-1+c); s = ss(j);
    hold all;
    
    % Draw prediction curves, based on XXG and DXXG.
    Ys = zeros(s, length(tc));
    for s0 = 1:s
        if c == 0
            comp = XXG(tc, spks(s0), 0, 1e5, P, T);
        else
            comp = XXG(tc, mu+spks(s0), sg, K, P, T);
        end
        comp(isnan(comp)) = 0;
        Ys(s0,:) = comp;
        for sf = s0+1:s
            if c == 0
                filt = (1-CXG(tc, spks(sf), 0, 1e5));
            else
                filt = (1-CXG(tc, mu+spks(sf), sg, K));
            end
            filt(isnan(filt)) = 1;
            Ys(s0,:) = Ys(s0,:) .* filt;
        end
    end
    for s0 = 1:s-1
        plot(tc, sum(Ys(1:s0,:), 1), ':', 'LineWidth', 2);
    end
    plot(tc, sum(Ys(1:s,:), 1), '-', 'LineWidth', 2);
    
    axis([-20 TC, 0 P/T+0.1]);
    ylabel(sprintf('Release (ms^{-1})'));
    set(gca, 'XTick', -20:20:TC, 'YTick', 0:0.2:0.6, 'FontSize', 10);
    if c == 1
        set(gca, 'Position', get(gca, 'Position') + [0.06 0, 0 0]);
    end
end
xlabel('Time (ms)');
end


%% Fig 9. Empirical Facilitation in Synchronous and Asynchronous Release Rates.
% Synchronous release rate shown in blue; asynchronous shown in red. Dark
% colors represent initial spike ramp (common to all traces on a plot);
% light colors represent single probe spikes from different simulations.
% Associated Ca2+ traces omitted for clarity.
% (A) PPF decays with increasing ISI (probe ISIs of 2, 5, 10, 20, 50,
% 100, 200 ms).
% (B) 5-spike ramp with a 5-ms ISI shows strong facilitation.
% (C) 5-spike ramp with a 20-ms ISI shows weaker facilitation.
% Note the change in scale from A to B and C.

load('data_files/RNP_runs.mat');
load('data_files/RNP_fits.mat');

figure; set(gcf, 'Units', 'inches', 'Position', [1 1, 6 6]);

% (A) PPF decays with increasing ISI (probe ISIs of 2, 5, 10, 20, 50,
% 100, 200 ms).

subplot(321); hold all;
for P = 2:length(dat.prob)
    plot(time{1,1}{P}, S_rel{1,1}{P}, ...
        'LineWidth', 2, 'Color', (Sr_color + [1 1 1]) / 2);
end
plot(time{1,1}{1}, S_rel{1,1}{1}, 'LineWidth', 3, 'Color', Sr_color / 2);
text(150, 0.6*0.25, 'synchronous', 'Color', Sr_color, 'FontSize', 10);
xlabel('Time (ms)'); ylabel('Release Rate (ms^{-1})');
axis([0 400, 0 0.25]);
set(gca, 'XTick', 0:100:400, 'YTick', 0:0.05:0.25);

subplot(322); hold all;
for P = 2:length(dat.prob)
    plot(time{1,1}{P}, A_rel{1,1}{P}, ...
        'LineWidth', 2, 'Color', (Ar_color + [1 1 1]) / 2);
end
plot(time{1,1}{1}, A_rel{1,1}{1}, 'LineWidth', 3, 'Color', Ar_color / 2);
text(150, 0.6*1.2e-3, 'asynchronous', 'Color', Ar_color, 'FontSize', 10);
xlabel('Time (ms)'); ylabel('Release Rate (ms^{-1})');
axis([0 400, 0 1.2e-3]);
set(gca, 'XTick', 0:100:400, 'YTick', 0:0.2e-3:1.2e-3);

% (B) 5-spike ramp with a 5-ms ISI shows strong facilitation.

subplot(323); hold all;
for P = 2:length(dat.prob)
    plot(time{5,2}{P}, S_rel{5,2}{P}, ...
        'LineWidth', 2, 'Color', (Sr_color + [1 1 1]) / 2);
end
plot(time{5,2}{1}, S_rel{5,2}{1}, 'LineWidth', 3, 'Color', Sr_color / 2);
text(150, 0.6*1.4, 'synchronous', 'Color', Sr_color, 'FontSize', 10);
xlabel('Time (ms)'); ylabel('Release Rate (ms^{-1})');
axis([0 400, 0 1.4]);
set(gca, 'XTick', 0:100:400, 'YTick', 0:0.2:1.4);

subplot(324); hold all;
for P = 2:length(dat.prob)
    plot(time{5,2}{P}, A_rel{5,2}{P}, ...
        'LineWidth', 2, 'Color', (Ar_color + [1 1 1]) / 2);
end
plot(time{5,2}{1}, A_rel{5,2}{1}, 'LineWidth', 3, 'Color', Ar_color / 2);
text(150, 0.6*0.01, 'asynchronous', 'Color', Ar_color, 'FontSize', 10);
xlabel('Time (ms)'); ylabel('Release Rate (ms^{-1})');
axis([0 400, 0 1e-2]);
set(gca, 'XTick', 0:100:400, 'YTick', 0:0.2e-2:1e-2);

% (C) 5-spike ramp with a 20-ms ISI shows weaker facilitation.

subplot(325); hold all;
for P = 2:length(dat.prob)
    plot(time{5,4}{P}, S_rel{5,4}{P}, ...
        'LineWidth', 2, 'Color', (Sr_color + [1 1 1]) / 2);
end
plot(time{5,4}{1}, S_rel{5,4}{1}, 'LineWidth', 3, 'Color', Sr_color / 2);
text(150, 0.6*1.4, 'synchronous', 'Color', Sr_color, 'FontSize', 10);
xlabel('Time (ms)'); ylabel('Release Rate (ms^{-1})');
axis([0 400, 0 1.4]);
set(gca, 'XTick', 0:100:400, 'YTick', 0:0.2:1.4);

subplot(326); hold all;
for P = 2:length(dat.prob)
    plot(time{5,4}{P}, A_rel{5,4}{P}, ...
        'LineWidth', 2, 'Color', (Ar_color + [1 1 1]) / 2);
end
plot(time{5,4}{1}, A_rel{5,4}{1}, 'LineWidth', 3, 'Color', Ar_color / 2);
text(150, 0.6*0.01, 'asynchronous', 'Color', Ar_color, 'FontSize', 10);
xlabel('Time (ms)'); ylabel('Release Rate (ms^{-1})');
axis([0 400, 0 1e-2]);
set(gca, 'XTick', 0:100:400, 'YTick', 0:0.2e-2:1e-2);


%% Fig 10. Release Rate Parameters and Facilitation Metaparameters Fitted
% to Empirical Histogram Profiles.
% (A) Synchronous and asynchronous profiles fitted for baseline
% (un-facilitated) case, and for highly facilitated case (probe spike 5 ms
% after 5-spike ramp of 5-ms ISIs).
% (B, C) Release fidelity values fitted case-by-case (dark colors) overlaid
% with predictions from best-fit facilitation functions (light colors) for
% synchronous (B) and asynchronous (C) components.

load('data_files/RNP_runs.mat');
load('data_files/RNP_Ca.mat');
load('data_files/RNP_fits.mat');
load('data_files/SNARE_0.mat');
sr0 = Sr0(6); ar0 = Ar0(6);
dat = struct('rISI', rISI, 'nrAP', nrAP, 'prob', prob);

% (A) Synchronous and asynchronous profiles fitted for baseline
% (un-facilitated) case, and for highly facilitated case (probe spike 5 ms
% after 5-spike ramp of 5-ms ISIs).

tm_lo =  time{1,1}{1};        tm_hi =  time{5,2}{3} - 25;
ix_lo = (tm_lo <= 250);       ix_hi = (tm_hi >= 0) & (tm_hi <= 250);
Sr_lo = S_rel{1,1}{1}(ix_lo); Sr_hi = S_rel{5,2}{3}(ix_hi);
SF_lo =   S_F{1,1}{1};        SF_hi =   S_F{5,2}{3};
Ar_lo = A_rel{1,1}{1}(ix_lo); Ar_hi = A_rel{5,2}{3}(ix_hi);
AF_lo =   A_F{1,1}{1};        AF_hi =   A_F{5,2}{3};
tm_lo = tm_lo(ix_lo);         tm_hi = tm_hi(ix_hi);

figure; set(gcf, 'Units', 'inches', 'Position', [1 1, 6 2.5]);

subplot(121); hold all;
plot(tm_hi, Sr_hi, 'LineWidth', 3, 'Color', [0.6 0.5 0]);
plot(tm_hi, SF_hi+sr0, 'LineWidth', 2, 'Color', [0.3 0.1 0]);
plot(tm_lo, Sr_lo, 'LineWidth', 3, 'Color', [0 0.5 0.6]);
plot(tm_lo, SF_lo+sr0, 'LineWidth', 2, 'Color', [0 0.1 0.3]);
axis([0 250, 1e-9 1]);
xlabel('Time (ms)');
ylabel(sprintf('Synchronous\nRelease Rate (ms^{-1})'));
set(gca, 'XTick', 0:50:250, 'YScale', 'log', ...
    'YTick', 10.^(-9:3:0), 'FontSize', 10);
legend({' facilitated case', ' fit to facilitated', ...
    ' baseline case', ' fit to baseline'}, 'FontSize', 8, ...
    'Position', [0.26 0.6, 0.2 0.27], 'Box', 'off');

subplot(122); hold all;
plot(tm_hi, Ar_hi, 'LineWidth', 3, 'Color', [0.6 0.5 0]);
plot(tm_hi, AF_hi+ar0, 'LineWidth', 2, 'Color', [0.3 0.1 0]);
plot(tm_lo, Ar_lo, 'LineWidth', 3, 'Color', [0 0.5 0.6]);
plot(tm_lo, AF_lo+ar0, 'LineWidth', 2, 'Color', [0 0.1 0.3]);
axis([0 250, 1e-6 1]);
xlabel('Time (ms)');
ylabel(sprintf('Asynchronous\nRelease Rate (ms^{-1})'));
set(gca, 'XTick', 0:50:250, 'YScale', 'log', ...
    'YTick', 10.^(-6:2:0), 'FontSize', 10);
legend({' facilitated case', ' fit to facilitated', ...
    ' baseline case', ' fit to baseline'}, 'FontSize', 8, ...
    'Position', [0.7 0.6, 0.2 0.27], 'Box', 'off');

% (B, C) Release fidelity values fitted case-by-case (dark colors) overlaid
% with predictions from best-fit facilitation functions (light colors) for
% synchronous (B) and asynchronous (C) components.
% For clarity, only P_1 and P_3 are shown for synchronous facilitation and
% only P_1 for asynchronous. Values for spikes along ramps shown in dark
% solid colors, and those for probe spikes shown in light colors coming off
% their respective ramp spikes. Most facilitation seen in ramp and probe
% spikes of short ISI, while longer probe ISIs return P_c close to baseline
% (horizontal dashed lines).

% Plot full facilitation functions across all stimulus cases.
bases = {S400, A400};
F_Pvs = {S_Pvs, A_Pvs};
cases = {[1 3], 1};
ylims = {[1e-6 1e2], [1e-3 1e0]};
ytcks = {10.^(-6:2:2), 10.^(-3:0)};
ylbls = {'Synchronous P_c(n)', 'Asynchronous P_c(n)'};
yLhts = {10, 0.5};
figys = {3, 2};
txt = {{'P_{S1}', 'P_{S2}', 'P_{S3}'}, {'P_{A1}', 'P_{A2}'}};

nc = 8; cmap = parula(length(dat.rISI) * nc);
cmap = cmap(end+1-nc:-nc:1,:);
nR = length(dat.rISI); nN = length(dat.nrAP); nP = length(dat.prob);
for fg = 1:2
    
figure; hold all;
set(gcf, 'Units', 'inches', 'Position', [1 1, 6 figys{fg}]);

for j = cases{fg}
    
    F_map = bases{fg}.Ps(j) * facil_map(F_Pvs{fg}(j,:), dat);
    
    for R = 1:nR
        R0 = (R-1)*nN*nP;
    for N = 1:nN
        N0 = (N-1)*nP;
        plot(R0+N0+(0:nP-1), F_map(:,N,R), '.:', ...
            'Color', cmap(R,:), 'LineWidth', 2, 'MarkerSize', 15);
    end
        plot(R0+(0:nN-1)*nP, F_map(1,:,R), '.-', ...
            'Color', 0.7*cmap(R,:), 'LineWidth', 2, 'MarkerSize', 15);
    end
    
    plot([1 160], bases{fg}.Ps(j)*[1 1], '--', ...
        'Color', [0.3 0.3 0.3], 'LineWidth', 1);
    
    text(165, exp(mean(log([min(F_map(:)), max(F_map(:))]))), ...
        txt{fg}{j}, 'FontSize', 10, ... 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end

for R = 1:nR
    text(20 + 40*(R-1), yLhts{fg}, ...
        sprintf('%d-ms ISI ramp', dat.rISI(R)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 10, 'Color', 0.5*cmap(R,:));
end

xlabel('Stimulus Case #'); ylabel(ylbls{fg});
axis([0 161, ylims{fg}]);
set(gca, 'YScale', 'log', 'YTick', ytcks{fg}, ...
    'XTick', [1 20:20:160], 'FontSize', 10);

end


%% S1 Fig. State Diagrams for VDCC, Calbindin, and PMCA.
% All diagrams reproduced with permission from Nadkarni et al. [35].
% (A) VDCC state transition model adapted from Bischofberger et al. [68].
% Transition rates ?ij and ?ji depend on membrane potential v.
% (B) State transitions for calbindin (CB) at high-affinity (H) and
% medium-affinity (M) Ca2+-binding sites. On rates (kh+ and km+) are
% proportional to [Ca2+]i.
% (C) PMCA pump state diagram with Ca2+ interactions depicted on the
% relative side of the membrane. Ca2+ leakage occurs only in state PMCA0.
% Association rate kpm1 is proportional to [Ca2+]i.


%% S2 Fig. Action-Potential-Evoked Ca2+ Current.
% (A) Action-potential-like waveform applied to axon.
% (B) Probability of a single VDCC being in the open state in response to
% the action potential in panel (A) increases from about 10–5 to around 96%
% during the spike before quickly shutting off; computed from deterministic
% simulation of state probabilities.
% (C) Rate of Ca2+ influx through a single, pathologically open channel
% (red) and through a typical channel (blue), whose probability of being
% open follows (B).

% Plot the VDCC current relative to a single action potential.

load('data_files/mixed_sim_VDCC.mat');

figure; fig_sz = [3.5 4];
set(gcf, 'Units', 'inches', 'PaperUnits', 'inches');
set(gcf, 'Position', [1 1, fig_sz], 'PaperPosition', [1 1, fig_sz]);

% (A) Action-potential-like waveform applied to axon.

% single action potential waveform
h2 = subplot(311);
plot(t, v, 'Color', [0 0 0], 'LineWidth', 2);
axis([0 10, -80 80]);
ylabel('V_m (mV)'); 
set(gca, 'YTick', -80:40:80, 'XTick', 0:2:10);

% (B) Probability of a single VDCC being in the open state in response to
% the action potential in panel (A) increases from about 10–5 to around 96%
% during the spike before quickly shutting off; computed from deterministic
% simulation of state probabilities.

% probability of being in the open state as a function of time
h3 = subplot(312); hold all;
plot(t, VDCC(5,:), 'Color', [0 0 0.7], 'LineWidth', 2);
axis([0 10, 0 1]);
ylabel('P (open state)');
set(gca, 'YTick', 0:0.2:1, 'XTick', 0:2:10);

% (C) Rate of Ca2+ influx through a single, pathologically open channel
% (red) and through a typical channel (blue), whose probability of being
% open follows (B).

% inward calcium current response to action potential
h4 = subplot(313); hold all;
plot(t, CaI, 'Color', [0.7 0 0], 'LineWidth', 2);
plot(t, VDCC(5,:) .* CaI, 'Color', [0 0 0.7], 'LineWidth', 2);
axis([0 10, 1e-3 1e3]);
ylabel('I_{Ca} (ions / ms)'); xlabel('Time (ms)');

set(gca, 'YScale', 'log', 'YTick', 10.^(-3:2:3), 'XTick', 0:2:10);
set([h1 h2 h3 h4], 'Box', 'off', 'FontSize', 10);
legend({' per open channel', ' per average channel'}, 'Box', 'off', ...
    'Position', [0.65 0.16, 0.2 0.12], 'FontSize', 8);


%% S3 Fig. Spontaneous Rates of Vesicle Fusion Increase with [Ca2+]i0.
% For small [Ca2+]i0, S0=kS?([Ca2+]i0)5 and A0=kA?([Ca2+]i0)2, where
% kS?6×10–4 ms–1??M–5 and kA?2×10–3 ms–1??M–2. As [Ca2+]i0??, S0??S and
% A0??A. Values for S0 and A0 at [Ca2+]i0=100 nM, which is used throughout
% most of this paper, are pointed out for reference.

% Plot the rates of spontaneous release as a function of baseline Ca.

load('data_files/SNARE_0.mat');

relS =   2*3;       %  msec^-1  - synchronous-activated release rate
relA =   2*0.025;   %  msec^-1  - asynchronous-activated release rate

% Plot equilibrium ('mini') release rates as a function of baseline Ca.
figure; hold all;
set(gcf, 'Units', 'inches', 'PaperUnits', 'inches');
set(gcf, 'Position', [1 1, 4 3], 'PaperPosition', [1 1, 4 3]);

% Synchronous Release
hS = plot(Ca0, Sr0, 'b.-', 'LineWidth', 2, 'MarkerSize', 25);
plot([3e-4 3e2], relS*[1 1], 'b--', 'LineWidth', 1);
text(0.03, 2*relS, sprintf('\\gamma_S = %.2f ms^{-1}', relS), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
    'FontSize', 10, 'Color', [0 0 1]);

% Asynchronous Release
hA = plot(Ca0, Ar0, 'r.-', 'LineWidth', 2, 'MarkerSize', 25);
plot([3e-4 3e2], relA*[1 1], 'r--', 'LineWidth', 1);
text(0.03, 0.5*relA, sprintf('\\gamma_A = %.2f ms^{-1}', relA), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
    'FontSize', 10, 'Color', [1 0 0]);

% Ca_0 used in the simulations (0.1 uM)
plot(0.1*[1 1], [1e-20 1e4], 'k:', 'LineWidth', 1.5);
plot([0.1 2], Sr0(6)*[1 1], 'b:', 'LineWidth', 1.5);
text(2.5, Sr0(6), strrep(sprintf('S_0 = %.2e} ms^{-1}', Sr0(6)), ...
    'e-0', '\times10^{-'), ...
    'Color', [0 0 1], 'FontSize', 8);
plot([0.1 2], Ar0(6)*[1 1], 'r:', 'LineWidth', 1.5);
text(2.5, Ar0(6), strrep(sprintf('A_0 = %.2e} ms^{-1}', Ar0(6)), ...
    'e-0', '\times10^{-'), ...
    'Color', [1 0 0], 'FontSize', 8);

axis([3e-4 3e2, 1e-20 1e4]);
xlabel('[Ca^{2+}]_{i0} (\muM)'); ylabel('Release Rate (ms^{-1})');
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 10);
set(gca, 'XTick', 10.^(-3:2), 'YTick', 10.^(-20:4:4));

legend([hS, hA], {' synchronous', ' asynchronous'}, ...
    'Units', 'normalized', 'Position', [0.55 0.25, 0.3 0.15], ...
    'Box', 'off', 'FontSize', 8);

% Calculate exponents of fits to Ca-dependence of spontaneous release rate.
S_ab = log(Sr0(2:3)) / [log(Ca0(2:3)); ones(1,2)];
A_ab = log(Ar0(2:3)) / [log(Ca0(2:3)); ones(1,2)];
fprintf('Sr0 = %.2e * Ca0^%d\n', exp(S_ab(2)), round(S_ab(1)));
fprintf('Ar0 = %.2e * Ca0^%d\n', exp(A_ab(2)), round(A_ab(1)));


%% S4 Fig. Spatial Modeling Important for Capturing Fine-Grain Features of
% Ca2+ Transients.
% Color scheme identical to that used in Fig 2: yellow to blue represent
% proximal to distal Ca2+ sensors.
% (A) [Ca2+]i measured at increasing distance from VDCC source (yellow to
% blue), with well-mixed approximation overlaid for comparison (maroon).
% Inset focuses on shorter time scale.
% (B) Profiles with calbindin removed from MCell (yellow to blue) and
% well-mixed model (maroon). Note that peak [Ca2+]i for the most proximal
% case extends up to 81 ?M, but is cut off for clarity.

load('data_files/mixed_Ca.mat');
load('data_files/flux_array_Ca.mat');
load('data_files/flux_array_Ca_no_Cb.mat');

figure; set(gcf, 'Units', 'inches', 'Position', [1 1, 4 4.5]);
cmap = parula(2*length(locs));

% (A) [Ca2+]i measured at increasing distance from VDCC source (yellow to
% blue), with well-mixed approximation overlaid for comparison (maroon).
% Inset focuses on shorter time scale.

subplot(211); hold all;
plts = zeros(length(locs), 1);
lbls = cell(length(locs), 1);
for j = 1:length(locs)
    plts(j) = plot(time{1}, Ca_avg{1}(:,j), ...
        'LineWidth', 1, 'Color', 0.7*cmap(2*(length(locs)-j)+1,:));
    lbls{j} = sprintf(' %.3f \\mum', 1+locs(j));
end
wm = plot(t, Ca_Cb45, 'LineWidth', 2, 'Color', [0.5 0 0]);
axis([0 120, 3e-2 1e2]);
xlabel('Time (ms)'); ylabel('[Ca^{2+}]_i (\muM)');
set(gca, 'YScale', 'log', 'XTick', 0:20:120, 'YTick', 10.^(-2:2), ...
    'FontSize', 10);
legend([plts([1,9,17]); wm], [lbls([1,9,17]); ' well-mixed'], ...
    'FontSize', 8, 'Box', 'off', 'Position', [0.22 0.75, 0.3 0.15]);

axes('Position', [0.6 0.75, 0.3 0.18]); hold all;
for j = 1:length(locs)
    plot(time{1}, Ca_avg{1}(:,j), ...
        'LineWidth', 1, 'Color', 0.7*cmap(2*(length(locs)-j)+1,:));
end
plot(t, Ca_Cb45, 'LineWidth', 2, 'Color', [0.5 0 0]);
axis([0 10, 3e-2 1e2]);
set(gca, 'YScale', 'log', 'XTick', 0:2:10, 'YTick', 10.^(-2:2), ...
    'FontSize', 8);

% (B) Profiles with calbindin removed from MCell (yellow to blue) and
% well-mixed model (maroon). Note that peak [Ca2+]i for the most proximal
% case extends up to 81 ?M, but is cut off for clarity.

subplot(212); hold all;
plts = zeros(length(locs), 1);
lbls = cell(length(locs), 1);
for j = 1:length(locs)
    plts(j) = plot(tm_noCb, Ca_noCb(:,j), ...
        'LineWidth', 1, 'Color', 0.7*cmap(2*(length(locs)-j)+1,:));
    lbls{j} = sprintf(' %.3f \\mum', 1+locs(j));
end
wm = plot(t, Ca_Cb00, 'LineWidth', 2, 'Color', [0.5 0 0]);
axis([0 120, 0 30]);
xlabel('Time (ms)'); ylabel('[Ca^{2+}]_i (\muM)');
set(gca, 'YScale', 'lin', 'XTick', 0:20:120, 'YTick', 0:10:50, ...
    'FontSize', 10);
legend([plts([1,9,17]); wm], [lbls([1,9,17]); ' well-mixed'], ...
    'FontSize', 8, 'Box', 'off', 'Position', [0.22 0.25, 0.3 0.15]);

axes('Position', [0.6 0.25, 0.3 0.18]); hold all;
for j = 1:length(locs)
    plot(tm_noCb, Ca_noCb(:,j), ...
        'LineWidth', 1, 'Color', 0.7*cmap(2*(length(locs)-j)+1,:));
end
plot(t, Ca_Cb00, 'LineWidth', 2, 'Color', [0.5 0 0]);
axis([0 10, 3e-2 1e2]);
set(gca, 'YScale', 'log', 'XTick', 0:2:10, 'YTick', 10.^(-2:2), ...
    'FontSize', 8);


%% S5 Fig. Effect of Calbindin Buffer on Spike-Evoked Ca2+ Profile and Release Rates.
% Action-potential-like stimulus delivered to model axon starting at 0 ms.
% Diffusion is assumed to be instantaneous, and molecular state
% probabilities are tracked deterministically over time.
% (A) Free [Ca2+]i with no calbindin buffer decays linearly with time due
% to saturation of PMCA pumps.
% (B) Syt-1/7-mediated release rates are large but short-lived in response
% to unbuffered Ca2+.
% (C) Free [Ca2+]i with calbindin added to the axon has much smaller
% magnitude and much narrower peak but has much longer tail.
% (D) Vesicle release in response to buffered Ca2+ is much less pronounced.
% The calbindin buffer reduces the rate of synchronous transmission but
% extends the window for pronounced asynchronous transmission.

load('data_files/mixed_Ca.mat');
load('data_files/mixed_runs.mat');
% dt = 0.1;       % ms
Ca_0 = 0.100;   % uM

figure; set(gcf, 'Units', 'inches', 'Position', [1 1, 6 4]);

% (A) Free [Ca2+]i with no calbindin buffer decays linearly with time due
% to saturation of PMCA pumps.

sp1 = subplot(221); %axes('Position', [0.19 0.82, 0.75 0.16]);
plot(sp1, t, Ca_Cb00, 'Color', Ca_color, 'LineWidth', 2);
ylabel(sp1, '[Ca^{2+}]_i (\muM)');
axis(sp1, [0 120, 0 12.1]);
set(sp1, 'XTick', 0:20:120, 'YTick', 0:4:12, ...
    'Box', 'off', 'FontSize', 10);
% legend({' [CB] = 0 \muM'}, 'Location', 'NorthEast', ...'Position', [0.6 0.88, 0.2 0.04], ...[0.31 0.95, 0.21 0.03], ...
%     'FontSize', 8, 'Box', 'off');
text(mean(xlim), mean(ylim), '[CB] = 0 \muM', 'FontSize', 10);

% in1 = axes('Position', [0.65 0.89, 0.29 0.09]);
% plot(in1, t, dt*cumsum(Ca_Cb00-Ca_0), ...
%     'Color', [0.1 0.3 0.7], 'LineWidth', 2);
% axis(in1, [0.5 5e3, 3e-2 3e1]);
% set(in1, 'XScale', 'lin', 'YScale', 'log', ...
%     'XTick', 10.^(0:4), 'YTick', 10.^(-2:2));

% (B) Syt-1/7-mediated release rates are large but short-lived in response
% to unbuffered Ca2+.

sp2 = subplot(222); %axes('Position', [0.19 0.58, 0.75 0.16]);
plot(sp2, t, Sr_Cb00+Ar_Cb00, 'Color', Cr_color, 'LineWidth', 2);
xlabel(sp2, 'Time (ms)');
ylabel(sp2, sprintf('Release\nRate (ms^{-1})'));
axis(sp2, [0 120, 0 0.5]);
set(sp2, 'XTick', 0:20:120, 'YTick', 0:0.1:0.5, ...
    'Box', 'off', 'FontSize', 10);
% legend({' [CB] = 0 \muM'}, 'Location', 'NorthEast', ...'Position', [0.6 0.64, 0.2 0.04], ...[0.31 0.71, 0.21 0.03], ...
%     'FontSize', 8, 'Box', 'off');
text(mean(xlim), mean(ylim), '[CB] = 0 \muM', 'FontSize', 10);

% in2 = axes('Position', [0.65 0.65, 0.29 0.09]);
% plot(in2, t, dt*cumsum(Sr_Cb00-sr0+Ar_Cb00-ar0), ...
%     'Color', [0.8 0.5 0.1], 'LineWidth', 2);
% axis(in2, [0.5 5e3, 1e-5 1e-4]);
% set(in2, 'XScale', 'lin', 'YScale', 'log', ...
%     'XTick', 10.^(0:4), 'YTick', 10.^(-6:2:0));

% (C) Free [Ca2+]i with calbindin added to the axon has much smaller
% magnitude and much narrower peak but has much longer tail.

sp3 = subplot(223); %axes('Position', [0.19 0.32, 0.75 0.16]);
plot(sp3, t, Ca_Cb45, 'Color', Ca_color, 'LineWidth', 2);
ylabel(sp3, '[Ca^{2+}]_i (\muM)');
axis(sp3, [0 120, 0 5]);
set(sp3, 'XTick', 0:20:120, 'YTick', 0:5, ...
    'Box', 'off', 'FontSize', 10);
% legend({' [CB] = 45 \muM'}, 'Location', 'NorthEast', ...'Position', [0.6 0.38, 0.2 0.04], ...[0.31 0.45, 0.21 0.03], ...
%     'FontSize', 8, 'Box', 'off');
text(mean(xlim), mean(ylim), '[CB] = 45 \muM', 'FontSize', 10);

% in3 = axes('Position', [0.65 0.39, 0.29 0.09]);
% plot(in3, t, dt*cumsum(Ca_Cb45-Ca_0), ...
%     'Color', [0.1 0.3 0.7], 'LineWidth', 2);
% axis(in3, [0.5 5e3, 3e-2 3e1]);
% set(in3, 'XScale', 'lin', 'YScale', 'log', ...
%     'XTick', 10.^(0:4), 'YTick', 10.^(-2:2));

% (D) Vesicle release in response to buffered Ca2+ is much less pronounced.
% The calbindin buffer reduces the rate of synchronous transmission but
% extends the window for pronounced asynchronous transmission.

sp4 = subplot(224); %axes('Position', [0.19 0.08, 0.75 0.16]);
plot(sp4, t, Sr_Cb45+Ar_Cb45, 'Color', Cr_color, 'LineWidth', 2);
xlabel(sp4, 'Time (ms)');
ylabel(sp4, sprintf('Release\nRate (ms^{-1})'));
axis(sp4, [0 120, 0 3e-4]);
set(sp4, 'XTick', 0:20:120, 'YTick', 0:1e-4:3e-4, ...
    'Box', 'off', 'FontSize', 10);
% legend({' [CB] = 45 \muM'}, 'Location', 'NorthEast', ...'Position', [0.6 0.14, 0.2 0.04], ...[0.31 0.21, 0.21 0.03], ...
%     'FontSize', 8, 'Box', 'off');
text(mean(xlim), mean(ylim), '[CB] = 45 \muM', 'FontSize', 10);

% in4 = axes('Position', [0.65 0.15, 0.29 0.09]);
% plot(in4, t, dt*cumsum(Sr_Cb45-sr0+Ar_Cb45-ar0), ...
%     'Color', [0.8 0.5 0.1], 'LineWidth', 2);
% axis(in4, [0.5 5e3, 1e-5 1e-4]);
% set(in4, 'XScale', 'lin', 'YScale', 'log', ...
%     'XTick', 10.^(0:4), 'YTick', 10.^(-6:-3));


%% S6 Fig. Synchronous and Asynchronous Release Rates Decrease with
% Distance from the Ca2+ Source.
% Color scheme identical to that used in Fig 2 and S4 Fig: yellow to blue
% represent proximal to distal Ca2+ sensors.
% (A) Synchronous release rate.
% (B) Integrated probability of synchronous release falls off nearly
% exponentially with distance to a baseline level.
% (C) Asynchronous release rates.
% (D) Integrated probability of asynchronous release also decays with
% distance to some baseline, but not exponentially.

load('data_files/flux_array_runs.mat');
dt = 0.1;

figure; set(gcf, 'Units', 'inches', 'Position', [1 1, 6 4.5]);
cmap = parula(2*length(locs));

% (A) Synchronous release rate.

subplot(221); hold all;
plts = zeros(length(locs), 1);
lbls = cell(length(locs), 1);
for j = 1:length(locs)
    plts(j) = plot(time{1}, S_rel{1}(:,j), ...
        'LineWidth', 1.5, 'Color', 0.7*cmap(2*(length(locs)-j)+1,:));
    lbls{j} = sprintf(' %.3f \\mum', 1+locs(j));
end
axis([0 120, 1e-9 1e1]);
xlabel('Time (ms)'); ylabel('Release Rate (ms^{-1})');
set(gca, 'YScale', 'log', 'XTick', 0:20:120, 'YTick', 10.^(-9:2:1), ...
    'FontSize', 10);
legend(plts([1,9,17]), lbls([1,9,17]), 'FontSize', 8, ...
    'Box', 'off', 'Position', [0.24 0.75, 0.3 0.15]);

% (B) Integrated probability of synchronous release falls off nearly
% exponentially with distance to a baseline level.

subplot(222);
plot(1+locs, sum(S_rel{1}, 1) * dt, 'k', 'LineWidth', 2);
axis([0 2, 5e-6 2]);
xlabel('Distance (\mum)');
ylabel(sprintf('Synchronous\nRelease Fidelity'));
set(gca, 'YScale', 'log', 'FontSize', 10, 'Box', 'off', ...
    'XTick', 1+locs(1:4:17), 'YTick', 10.^(-5:1));

% (C) Asynchronous release rates.

subplot(223); hold all;
plts = zeros(length(locs), 1);
lbls = cell(length(locs), 1);
for j = 1:length(locs)
    plts(j) = plot(time{1}, A_rel{1}(:,j), ...
        'LineWidth', 1.5, 'Color', 0.7*cmap(2*(length(locs)-j)+1,:));
    lbls{j} = sprintf(' %.3f \\mum', 1+locs(j));
end
axis([0 120, 3e-6 3e-2]);
xlabel('Time (ms)'); ylabel('Release Rate (ms^{-1})');
set(gca, 'YScale', 'log', 'XTick', 0:20:120, 'YTick', 10.^(-6:0), ...
    'FontSize', 10);
legend(plts([1,9,17]), lbls([1,9,17]), 'FontSize', 8, ...
    'Box', 'off', 'Position', [0.24 0.3, 0.3 0.15]);

% (D) Integrated probability of asynchronous release also decays with
% distance to some baseline, but not exponentially.

subplot(224);
plot(1+locs, sum(A_rel{1}, 1) * dt, 'k', 'LineWidth', 2);
axis([0 2, 0 0.1]);
xlabel('Distance (\mum)');
ylabel(sprintf('Asynchronous\nRelease Fidelity'));
set(gca, 'YScale', 'lin', 'FontSize', 10, 'Box', 'off', ...
    'XTick', 1+locs(1:4:17), 'YTick', 0:0.02:0.1);


%% S7 Fig. Parametric Fits to Release Histogram Profiles at Increasing
% Distance from the Ca2+ Source.
% (A,B) Fitted release profiles (black) imposed over the true histograms
% for synchronous (A, blue) and asynchronous (B, red).
% (C) Parameter values as a function of distance for synchronous release.
% (D) The same for asynchronous release.

load('data_files/flux_array_runs.mat');
load('data_files/flux_array_Sv_Av.mat');
load('data_files/SNARE_0.mat');
sr0 = Sr0(6); ar0 = Ar0(6);

% (A,B) Fitted release profiles (black) imposed over the true histograms
% for synchronous (A, blue) and asynchronous (B, red).

[~, S_Fs] = fit_conv_profiles(time{1}, S_rel{1}, ...
    sr0, 1+locs, Sv, [], false);
[~, A_Fs] = fit_conv_profiles(time{1}, A_rel{1}, ...
    ar0, 1+locs, Av, [], false);

figure; set(gcf, 'Units', 'inches', 'Position', [1 1, 6 2.5]);

subplot(121); hold all;
hsr = plot(time{1}, S_rel{1}, 'LineWidth', 3, 'Color', (Sr_color+1)/2);
hsf = plot(time{1}, S_Fs, 'LineWidth', 1, 'Color', [0 0 0]);
axis([0 100, 1e-9 2]);
xlabel('Time (ms)'); ylabel('Release Rate (ms^{-1})');
set(gca, 'YScale', 'log', 'YTick', 10.^(-8:2:0), 'XTick', 0:20:100, ...
    'FontSize', 10);
legend([hsr(1), hsf(1)], {' synchronous release', ' fitted profiles'}, ...
    'FontSize', 8, 'Box', 'off', 'Position', [0.2 0.7, 0.27 0.2]);

subplot(122); hold all;
har = plot(time{1}, A_rel{1}, 'LineWidth', 3, 'Color', (Ar_color+1)/2);
haf = plot(time{1}, A_Fs, 'LineWidth', 1, 'Color', [0 0 0]);
axis([0 100, 5e-6 1e-2]);
xlabel('Time (ms)'); ylabel('Release Rate (ms^{-1})');
set(gca, 'YScale', 'log', 'YTick', 10.^(-5:-1), 'XTick', 0:20:100, ...
    'FontSize', 10);
legend([har(1), haf(1)], {' asynchronous release', ' fitted profiles'}, ...
    'FontSize', 8, 'Box', 'off', 'Position', [0.65 0.7, 0.27 0.2]);

% (C) Parameter values as a function of distance for synchronous release.

figure; set(gcf, 'Units', 'inches', 'Position', [1 1, 3 2]);

subplot(221); plot(1+locs, Sv.Ps, 'LineWidth', 1);
axis([0 2, 1e-14 10]);
% xlabel('Distance (\mum)');
ylabel('P_c');
set(gca, 'YScale', 'log', 'XTick', [0 1 2], ...
    'YTick', 10.^([-16 -8 0]), 'FontSize', 8, 'Box', 'on');

subplot(222); plot(1+locs, Sv.Ks, 'LineWidth', 1);
axis([0 2, 1e-5 1e3]);
% xlabel('Distance (\mum)');
ylabel('k_c (ms^{-1})');
set(gca, 'YScale', 'log', 'XTick', [0 1 2], ...
    'YTick', 10.^([-6 -1 4]), 'FontSize', 8, 'Box', 'on');

subplot(223); plot(1+locs, Sv.Ms, 'LineWidth', 1);
axis([0 2, 1 1e3]);
xlabel('Distance (\mum)');
ylabel('\mu_c (ms)');
set(gca, 'YScale', 'log', 'XTick', [0 1 2], ...
    'YTick', 10.^(0:3), 'FontSize', 8, 'Box', 'on');

subplot(224); plot(1+locs, Sv.Ss, 'LineWidth', 1);
axis([0 2, 1e-2 1e2]);
xlabel('Distance (\mum)');
ylabel('\sigma_c (ms)');
set(gca, 'YScale', 'log', 'XTick', [0 1 2], ...
    'YTick', 10.^([-2 0 2]), 'FontSize', 8, 'Box', 'on');

% (D) The same for asynchronous release.

figure; set(gcf, 'Units', 'inches', 'Position', [1 1, 3 2]);

subplot(221); plot(1+locs, Av.Ps, 'LineWidth', 1);
axis([0 2, 1e-6 1]);
% xlabel('Distance (\mum)');
ylabel('P_c');
set(gca, 'YScale', 'log', 'XTick', [0 1 2], ...
    'YTick', 10.^([-6 -3 0]), 'FontSize', 8, 'Box', 'on');

subplot(222); plot(1+locs, Av.Ks, 'LineWidth', 1);
axis([0 2, 1e-3 1e1]);
% xlabel('Distance (\mum)');
ylabel('k_c (ms^{-1})');
set(gca, 'YScale', 'log', 'XTick', [0 1 2], ...
    'YTick', 10.^([-3 -1 1]), 'FontSize', 8, 'Box', 'on');

subplot(223); plot(1+locs, Av.Ms, 'LineWidth', 1);
axis([0 2, 1 1e3]);
xlabel('Distance (\mum)');
ylabel('\mu_c (ms)');
set(gca, 'YScale', 'log', 'XTick', [0 1 2], ...
    'YTick', 10.^(0:3), 'FontSize', 8, 'Box', 'on');

subplot(224); plot(1+locs, Av.Ss, 'LineWidth', 1);
axis([0 2, 1e-2 1e2]);
xlabel('Distance (\mum)');
ylabel('\sigma_c (ms)');
set(gca, 'YScale', 'log', 'XTick', [0 1 2], ...
    'YTick', 10.^([-2 0 2]), 'FontSize', 8, 'Box', 'on');


%% S8 Fig. Change in the Balance of Binding Kinetics and Internal State
% Distribution of Ca2+ Sensor with Spike History.
% State diagrams the same as shown in Fig 3.
% (A) Synchronous state diagrams. At baseline [Ca2+]i (first red dot),
% unbinding kinetics (left arrows) overpower binding (right arrows),
% biasing Syt-1 toward unbound state (S0; top diagram), with almost no
% probability of having any Ca2+ ions bound before an action potential
% (left pie chart). During peak Ca2+ influx (second red dot), binding rates
% (thicker right arrows) overpower unbinding, biasing Syt-1 toward its
% fully-bound releasable state (S5; lower diagram), with much greater
% probability of having at least some Ca2+ bound (right pie chart).
% (B) The same for asynchronous release with Syt-7, whose releasable state
% requires two Ca2+ ions bound (A2). Slower kinetics lead to only slight
% bias in favor of binding during an action potential (slightly thicker
% right arrows in lower diagram), leading to miniscule increase in
% probability of being in the releasable state on later spikes (right pie
% chart). Release becomes more probable on subsequent spikes because
% previous activity has pushed synaptotagmin into higher-bound states,
% making reaching the releasable state easier.

% SNARE State Transitions

% Load SNARE binding/unbinding rates and state probabilities,
% and convert them into normalized values.
load('data_files/SNARE_states.mat');
to_range = @(data, mxD, p, mx) round(10*(data / mxD).^p * mx) / 10;

% normalizing linewidths for state transitions
mnDT = min([S_n, S0p, S1p, A_n, A0p, A1p]);
mxDT = max([S_n, S0p, S1p, A_n, A0p, A1p]);
mxT  =  8;  pT = log(0.3/ mxT) / log(mnDT/mxDT);
mxAT = 16; pAT = log(3.0/mxAT) / log(mnDT/mxDT);

% normalizing linewidths for release transitions
mnDR = min([S0rel, S1rel, A0rel, A1rel]);
mxDR = max([S0rel, S1rel, A0rel, A1rel]);
mxR  =  4;  pR = log(0.3/ mxR) / log(mnDR/mxDR);
mxAR =  8; pAR = log(3.0/mxAR) / log(mnDR/mxDR);

font_max = 12;
font_med = 10;
font_min =  8;

% data for drawing state circles
radius = 0.17;                  % inches
rel_xr = 0.35; rel_yr = 0.15;   % inches
Drel   = 0.35;                  % inches

tt = linspace(0, 2*pi, 1+4*20);
xs = radius*cos(tt); ys = radius*sin(tt);   % state circle
xr = rel_xr*cos(tt); yr = rel_yr*sin(tt);   % release ellipse

% Determine where things go within the figure.
AB_sz = [5.75, 0.75];                       % inches
CD_sz = [AB_sz(1)/2, 2.2];                  % inches
fig_sz = [AB_sz(1), CD_sz(2)+2*AB_sz(2)];   % inches

L_buf = 0.2; R_buf = 0.2;
xcA = L_buf + radius;                               % start
xcB = AB_sz(1) - 2*rel_xr - Drel - radius - R_buf;  % end
xcR = AB_sz(1) - rel_xr - R_buf; yc = AB_sz(2) / 2; % positions

pie_rad = 0.6;       % inches - pie chart radius

% Draw figures with steady-state and peak state-transitions and pie charts,
% with [Ca++] insets, for both synchronous and asynchronous mechanisms.
for fig = 1:2

% Drawing parameters based on which figure.
switch fig
    % (A) Synchronous state diagrams. At baseline [Ca2+]i (first red dot),
    % unbinding kinetics (left arrows) overpower binding (right arrows),
    % biasing Syt-1 toward unbound state (S0; top diagram), with almost no
    % probability of having any Ca2+ ions bound before an action potential
    % (left pie chart). During peak Ca2+ influx (second red dot), binding rates
    % (thicker right arrows) overpower unbinding, biasing Syt-1 toward its
    % fully-bound releasable state (S5; lower diagram), with much greater
    % probability of having at least some Ca2+ bound (right pie chart).
    case 1
        xc = linspace(xcA, xcB, 6); D_xc = xc(2) - xc(1) - 2*radius;
        color = color_grad(S0_color, Cr_color, 6);
        type = 'S'; 
        neg_L = to_range(S_n, mxDT, pT, mxT);   % unbinding linewidths
        neg_A = to_range(S_n, mxDT, pAT, mxAT); % unbinding arrowheads
        
        pos0L = to_range(S0p, mxDT, pT, mxT);   % binding linewidths 0
        pos0A = to_range(S0p, mxDT, pAT, mxAT); % binding arrowheads 0
        
        pos1L = to_range(S1p, mxDT, pT, mxT);   % binding linewidths 1
        pos1A = to_range(S1p, mxDT, pAT, mxAT); % binding arrowheads 1
        
        rel0L = to_range(S0rel, mxDR, pR, mxR);     % release linewidth 0
        rel0A = to_range(S0rel, mxDR, pAR, mxAR);   % release arrowhead 0
        
        rel1L = to_range(S1rel, mxDR, pR, mxR);     % release linewidth 1
        rel1A = to_range(S1rel, mxDR, pAR, mxAR);   % release arrowhead 1
        
        pie_0 = S0pie;
        pie0lbl = {'S_0', 'S_1', '', '', '', 'S_5'};
        pie0Dxy = {[0 0], [0.2 0], [0 0], [0 0], [0 0], [-0.2 0.1]};
        pie_1 = S1pie;
        pie1lbl = {'S_0', 'S_1', 'S_2', 'S_3', 'S_4', 'S_5'};
        pie1Dxy = {[0 0], [0 0], [0 0], [0 0], [0 0], [-0.2 0.1]};
        
        sf = figure;
        
    % (B) The same for asynchronous release with Syt-7, whose releasable state
    % requires two Ca2+ ions bound (A2). Slower kinetics lead to only slight
    % bias in favor of binding during an action potential (slightly thicker
    % right arrows in lower diagram), leading to miniscule increase in
    % probability of being in the releasable state on later spikes (right pie
    % chart). Release becomes more probable on subsequent spikes because
    % previous activity has pushed synaptotagmin into higher-bound states,
    % making reaching the releasable state easier.
    case 2
        xc = linspace(xcA, xcB, 3); D_xc = xc(2) - xc(1) - 2*radius;
        color = color_grad(A0_color, Cr_color, 3);
        type = 'A';
        neg_L = to_range(A_n, mxDT, pT, mxT);   % unbinding linewidths
        neg_A = to_range(A_n, mxDT, pAT, mxAT); % unbinding arrowheads
        
        pos0L = to_range(A0p, mxDT, pT, mxT);   % binding linewidths 0
        pos0A = to_range(A0p, mxDT, pAT, mxAT); % binding arrowheads 0
        
        pos1L = to_range(A1p, mxDT, pT, mxT);   % binding linewidths 1
        pos1A = to_range(A1p, mxDT, pAT, mxAT); % binding arrowheads 1
        
        rel0L = to_range(A0rel, mxDR, pR, mxR);     % release linewidth 0
        rel0A = to_range(A0rel, mxDR, pAR, mxAR);   % release arrowhead 0
        
        rel1L = to_range(A1rel, mxDR, pR, mxR);     % release linewidth 1
        rel1A = to_range(A1rel, mxDR, pAR, mxAR);   % release arrowhead 1
        
        pie_0 = A0pie;
        pie0lbl = {'A_0', 'A_1', 'A_2'};
        pie0Dxy = {[0 0], [0.1 0], [-0.2 0.1]};
        pie_1 = A1pie;
        pie1lbl = {'A_0', 'A_1', 'A_2'};
        pie1Dxy = {[0 0], [0 0], [-0.2 0.1]};
        
        af = figure;
end

% Create figure.
set(gcf, 'Units', 'inches', 'PaperUnits', 'inches');
set(gcf, 'Position', [1 1, fig_sz], 'PaperPosition', [1 1, fig_sz]);

for diag = 1:2

% Set up subplot axis for this state transitio diagram.
base_y = CD_sz(2)+(2-diag)*AB_sz(2);
ax_pos = [0 base_y, AB_sz] ./ [fig_sz, fig_sz];
ax = axes('Position', ax_pos); hold all;

% Drawing parameters depend on which diagram.
switch diag
    case 1  % steady-state
        pos_L = pos0L; pos_A = pos0A;
        rel_L = rel0L; rel_A = rel0A;
    case 2  % peak
        pos_L = pos1L; pos_A = pos1A;
        rel_L = rel1L; rel_A = rel1A;
end

% Draw the circles along the state diagram.
for j = 1:length(xc)
    fill(xc(j) + xs, yc + ys, color(j,:));
    ht = text(xc(j), yc, sprintf('%s_%d', type, j-1));
    set(ht, 'FontSize', font_max, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
end

% Draw the arrows between the states.
for j = 1:length(pos_L)
    % Make positive and negative arrows.
    arP = annotation('arrow'); arN = annotation('arrow');
    set([arP arN], 'HeadStyle', 'vback1', 'Units', 'inches');
    
    % Ca-binding (positive) arrow size
    arP.LineWidth  = pos_L(j);
    arP.HeadWidth  = pos_A(j);
    arP.HeadLength = pos_A(j);
    
    % Ca-unbinding (negative) arrow size
    arN.LineWidth  = neg_L(j);
    arN.HeadWidth  = neg_A(j);
    arN.HeadLength = neg_A(j);
    
    xL = xc( j ) + radius;
    xR = xc(j+1) - radius;
    
    arP.Position = [xL yc+0.4*radius+base_y,  D_xc 0];
    arN.Position = [xR yc-0.4*radius+base_y, -D_xc 0];
end

% Draw release state.
LR = 1;
xrc = xc(end) + (1+1.5)*radius + LR;
fill(xcR + xr, yc + yr, [0.9 0.9 0.9]);
ht = text(xcR, yc, 'release');
set(ht, 'FontSize', font_med, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

% Draw release arrow.
arR = annotation('arrow');
set(arR, 'HeadStyle', 'vback1', 'Units', 'inches');
arR.LineWidth  = rel_L;
arR.HeadWidth  = rel_A;
arR.HeadLength = rel_A;
arR.Position = [xc(end)+radius yc+base_y, Drel 0];

% Add label.
r_lbl = text(xc(end)+radius+Drel/2, yc+0.6*radius, ...
    sprintf('\\gamma_%s', type), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
    'FontSize', font_min, 'FontWeight', 'bold');

% After plots, size the axes and make them invisible.
axis(ax, [0 AB_sz(1), 0 AB_sz(2)]); axis off;

end

% Draw the state pie charts.
for chart = 1:2
    base_x = (chart-1)*CD_sz(1);
    ax = axes('Position', [base_x 0, CD_sz] ./ [fig_sz, fig_sz]);
    hold all;
    
    % Different pie charts use different data.
    switch chart
        case 1
            states = pie_0;
            pie_lbl = pie0lbl;
            pie_Dxy = pie0Dxy;
            xL = [-0.35*CD_sz(1), 0.65*CD_sz(1)];
        case 2
            states = pie_1;
            pie_lbl = pie1lbl;
            pie_Dxy = pie1Dxy;
            xL = [-0.65*CD_sz(1), 0.35*CD_sz(1)];
    end
    hp = pie(states, pie_lbl);
    
    % Change colors of wedges and sizes of labels.
    for j = 1:length(states)
        hp(2*j-1).FaceColor = color(j,:);
        hp(2*j-1).Vertices = pie_rad * hp(2*j-1).Vertices;
        set(hp(2*j), 'FontSize', font_med, 'FontWeight', 'bold');
        hp(2*j).Position = pie_rad * (hp(2*j).Position + [pie_Dxy{j}, 0]);
    end
    
    % Set the releasable state apart.
    verts = hp(end-1).Vertices; NV = size(verts, 1);
    verts = mean(verts, 1); verts = verts / sqrt(sum(verts.^2));
    hp(end-1).Vertices = hp(end-1).Vertices ...
        + 0.2*pie_rad*(ones(NV,1)*verts);
    
    % After plots, size the axes and make them invisible.
    axis(ax, [xL, -0.5*CD_sz(2) 0.5*CD_sz(2)]); axis off;
end

% Draw inset to show where on the [Ca++] trace data are drawn from.
ax = axes('Position', [0.4 0.15, 0.25 0.35]); hold all;
plot(ax, Ca_t, Ca_c, 'LineWidth', 2, 'Color', Ca_color);
plot(ax, Ca_x, Ca_y, '.', 'MarkerSize', 25, 'Color', [0.8 0 0]);
ylabel('[Ca^{2+}]_i (\muM)'); xlabel('Time (ms)');
axis(ax, [0 40, 0 20]);
set(ax, 'FontSize', font_med);
set(ax, 'XTick', 0:10:40, 'YTick', 0:5:20);

end


%% S9 Fig. Empirical Facilitation in Release Probability is a Nonlinear
% Function of Spike History and Ca2+ Buildup.
% (A) [Ca2+]i and (B) release rate in response to a 5-spike ramp stimulus
% with a 10-ms ISI (black and dark green), followed by a single probe spike
% at increasing delay from the end of the ramp (gray and light green;
% multiple cases overlaid on the same plot). Release rate grows much faster
% than Ca2+ buildup can account for.

load('data_files/RNP_Ca.mat');
load('data_files/RNP_runs.mat');

% Show how facilitation of Ca2 is insufficient to explain
% the highly nonlinear facilitation in release rate.

figure; set(gcf, 'Units', 'inches', 'Position', [1 1, 5 4]);

% Plot Ca and (combined) release rate with a
% 5-spike ramp at 10-ms ISI for all probe spikes.
N = 5; R = 3;
for P = 2:8
    subplot(211); hold all;
    plot(time{N,R}{P}, Ca_avg{N,R}{P}, ...
        'Color', (Ca_color+1)/2, 'LineWidth', 2);
    subplot(212); hold all;
    plot(time{N,R}{P}, S_rel{N,R}{P}+A_rel{N,R}{P}, ...
        'Color', (Cr_color+1)/2, 'LineWidth', 2);
end

% Plot ramp spikes over probe spikes in darker colors.
subplot(211);
plot(time{N,R}{1}, Ca_avg{N,R}{1}, ...
    'Color', (Ca_color+0)/2, 'LineWidth', 2);
axis([0 250, 0 30]);
ylabel('[Ca^{2+}]_i (\muM)'); xlabel('Time (ms)');
set(gca, 'YTick', 0:10:30, 'XTick', 0:50:250, 'FontSize', 10);

subplot(212);
plot(time{N,R}{1}, S_rel{N,R}{1}+A_rel{N,R}{1}, ...
    'Color', (Cr_color+0)/2, 'LineWidth', 2);
axis([0 250, 0 1]);
ylabel('Release Rate (ms^{-1})'); xlabel('Time (ms)');
set(gca, 'YTick', 0:0.2:1, 'XTick', 0:50:250, 'FontSize', 10);


%% S10 Fig. Empirical Facilitation in Release Probability is a Nonlinear
% Function of Spike History.
% Integrated release fidelity (P(n)) relative to baseline (P(0)) for the
% various stimulus cases explored. Ramp # indicates the number of spikes in
% the ramp preceding the probe spike, and ?t represents the ISI between the
% last ramp spike and the probe spike. Spike history noticeably affects the
% growth of facilitation, as seen for ramps with 2-ms ISIs (A), 5-ms ISIs (B),
% 10-ms ISIs (C), and 20-ms ISIs (D). Different colors distinguish
% facilitation functions with different spike histories. Dark lines follow
% relative release fidelity for spikes along spike ramps, and dotted lines
% follow relative release fidelity for probe spikes.

load('data_files/true_facil.mat');
nR = length(dat.rISI); nN = length(dat.nrAP); nP = length(dat.prob);

P_rel = P_Srel + P_Arel;
% P_rel = P_Srel;
% P_rel = P_Arel;
FR = P_rel / P_rel(1,1,1);
DT = dat.prob;
DT(1) = 1;

nc = 8; cmap = parula(nR * nc);
cmap = cmap(end+1-nc:-nc:1,:);
subpos = [ ...
    0.10 0.60, 0.36 0.35; ...
    0.60 0.60, 0.36 0.35; ...
    0.10 0.08, 0.36 0.35; ...
    0.60 0.08, 0.36 0.35];
% subpos = [ ...
%     0.1323 0.5838, 0.3324 0.3412; ...
%     0.5726 0.5838, 0.3324 0.3412; ...
%     0.1323 0.1100, 0.3324 0.3412; ...
%     0.5726 0.1100, 0.3324 0.3412];

figure; set(gcf, 'Units', 'inches', 'Position', [1 1, 6 5]);
for R = 1:nR
    subplot(2, 2, R);
    hold all; ax = gca;
    for N = 1:nN
        plot3(DT, N+0*(1:nP), FR(:,N,R), '.:', ...
            'Color', cmap(R,:), 'LineWidth', 2, 'MarkerSize', 15);
    end
    plot3(1+0*(1:nN), 1:nN, FR(1,:,R), '.-', ...
        'Color', 0.7*cmap(R,:), 'LineWidth', 2, 'MarkerSize', 15);
    
    title(sprintf('%d-ms ISI ramp', dat.rISI(R)), 'Color', 0.7*cmap(R,:));
    xlabel('\Deltat (ms)')
    ylabel('Ramp #');
    zlabel('P(n) / P(0)');
    axis([1 250, 0.5 5.5, 1 100]);
    set(ax, 'XScale', 'log', 'YScale', 'lin', 'ZScale', 'log', ...
        'XTick', [2 20 200], 'YTick', dat.nrAP, ...
        'Position', subpos(R,:));
%     set(ax.XAxis.Label, 'Position', [10, 0, 0.4]);
    view([45 20]);
end


%% S11 Fig. Saturation of Facilitation Parameters.
% (A) Facilitation parameter f(?) increases almost linearly from one spike
% (f(n-1)) to the next (f(n)), until it approaches some limit N?1.
% (B) Curves represent the unseen change in f(?) between spikes. Dots
% represent actual values observed at spike times, values determined by the
% Ca2+-triggered increment in release fidelity at each spike. Steady-state
% value for facilitation parameter limited by stimulus frequency and by
% value of N. No facilitation above baseline occurs for N=1.

% (A) Facilitation parameter f(?) increases almost linearly from one spike
% (f(n-1)) to the next (f(n)), until it approaches some limit N?1.

% one-component positive facilitation:
%
% g = f(n-1) * exp(-t / T)
% f(n) = g + 1 - (g/N)^N
% F(n) = f(n)^X
%
% N = number of facilitation steps to reach saturation
% T = facilitation decay time constant
% X = nonlinear scaling of facilitation

% Describe facilitation on the next spike, following a decay and a jump.
fn = @(fn1, T, N) (fn1.*exp(-T)) + 1 - ((fn1.*exp(-T))./N).^N;

maxN = 6;   % the maximum number of steps to saturation (for this figure)
fg_sz = 2.5/maxN * [maxN maxN+1];   % sizing proportional axes

figure; hold all;
set(gcf, 'Units', 'inches', 'PaperUnits', 'inches');
set(gcf, 'Position', [1 1, fg_sz], 'PaperPosition', [1 1, fg_sz]);

% Draw six facilitation plots.
plts = zeros(1, maxN);
lbls =  cell(1, maxN);
k = 3; cmap = parula(maxN*k);
for N = maxN:-1:1
    plts(N) = plot(0:0.01:N, fn(0:0.01:N, 0, N), ...
        'Color', cmap((N-1)*k+1,:), 'LineWidth', 3);
    lbls{N} = sprintf(' N = %d', N);
    plot([N N], [0 N], '--', ...
        'Color', cmap((N-1)*k+1,:), 'LineWidth', 2);
end

% Set axis labels and properties.
plot([0 maxN], [0 maxN]  , 'k--', 'LineWidth', 1);
plot([0 maxN], [0 maxN]+1, 'k--', 'LineWidth', 1);
axis([0 maxN, 0 maxN+1]);
xlabel('f (n - 1) \cdot exp(-\Deltat /\tau)');
ylabel('f (n)');
set(gca, 'XTick', 0:maxN, 'YTick', 0:maxN+1, 'FontSize', 10);
legend(plts(end:-1:1), lbls(end:-1:1), ...
    'Box', 'off', 'Position', [0.32 0.55, 0 0.4], 'FontSize', 8);

% (B) Curves represent the unseen change in f(?) between spikes. Dots
% represent actual values observed at spike times, values determined by the
% Ca2+-triggered increment in release fidelity at each spike. Steady-state
% value for facilitation parameter limited by stimulus frequency and by
% value of N. No facilitation above baseline occurs for N=1.

% Show the decay-and-jump dynamic working out.

fg_sz = [4.6, 4];

figure;
set(gcf, 'Units', 'inches', 'PaperUnits', 'inches');
set(gcf, 'Position', [1 1, fg_sz], 'PaperPosition', [1 1, fg_sz]);

dt = 0.01; T = 10; tm = 0:dt:T;
NpT = [5 2 0.5];    % number of spikes per time constant
N = [4, 2, 1];      % number of spikes to reach saturation

cmap = lines(length(N));
plts = zeros(size(N));  % plot handles
lbls = cell(size(N));   % legend labels
yL = [0.9 0.53 0.2];    % y-positions of legends
for p = 1:length(N)
    subplot(length(N), 1, p); hold all;
for j = 1:length(NpT)
    ff = 0*tm;
    spks = (1:(NpT(j)*T-1)) / NpT(j);
    fv = 0*spks;
    lbls{j} = sprintf('\\Deltat = %.3g \\tau', 1/NpT(j));
for s = 1:length(spks)
    % Facilitate at the specified spike time.
    ix = find(tm >= spks(s), 1, 'first');
    df = 1 - (ff(ix)/N(p))^N(p);
    fv(s) = ff(ix) + df;
    
    % Facilitation decays after this time.
    ff(ix:end) = (ff(ix) + df) * exp(-(tm(ix:end)-spks(s)));
end
    % Plot decay curves and points at which facilitation is measured.
    plts(j) = plot(tm, ff, 'LineWidth', 1.5, 'Color', cmap(j,:));
    plot(spks, fv, '.', 'MarkerSize', 16, 'Color', cmap(j,:));
end
    ylabel(sprintf('f (t)\nN = %d', N(p)));
    axis([0 T, 0 max(N)]);
    set(gca, 'YTick', 0:4, 'XTick', 0:T, 'FontSize', 10);
    legend(plts, lbls, 'Box', 'off', ...
        'Position', [0.12 yL(p), 0.8 0.1], ...
        'Orientation', 'horizontal', 'FontSize', 8);
end
xlabel('t / \tau');


%% S12 Fig. Release Rate Parameters and Facilitation Metaparameters Fitted
% to Empirical Histogram Profiles.
% Errors across all cases in linear and logarithmic space for the
% predictive model.

load('data_files/RNP_runs.mat');
load('data_files/RNP_Ca.mat');
load('data_files/RNP_fits.mat');
load('data_files/SNARE_0.mat');
sr0 = Sr0(6); ar0 = Ar0(6);
dat = struct('rISI', rISI, 'nrAP', nrAP, 'prob', prob);

S_err = reshape(sum(S_FVU, 1)', ...
    [length(dat.prob), length(dat.nrAP)*length(dat.rISI)]);
A_err = reshape(sum(A_FVU, 1)', ...
    [length(dat.prob), length(dat.nrAP)*length(dat.rISI)]);

figure; set(gcf, 'Units', 'inches', 'Position', [1 1, 6 4]);
subplot(211); imagesc(S_err); caxis([0 0.2]);
hc = colorbar; ylabel(hc, 'FVU Error');
title('Error of Synchronous Facilitation Fits');
xlabel('Ramp Spikes (1-5)'); ylabel('Probe Spikes (ms)');
set(gca, 'YTick', 1:length(dat.prob), 'YTickLabel', int2str(dat.prob'), ...
    'FontSize', 10, 'XTick', 3+5*(0:3), 'XTickLabel', ...
    {'2-ms ISI', '5-ms ISI', '10-ms ISI', '20-ms ISI'});
subplot(212); imagesc(A_err); caxis([0 0.2]);
hc = colorbar; ylabel(hc, 'FVU Error');
title('Error of Asynchronous Facilitation Fits');
xlabel('Ramp Spikes (1-5)'); ylabel('Probe Spikes (ms)');
set(gca, 'YTick', 1:length(dat.prob), 'YTickLabel', int2str(dat.prob'), ...
    'FontSize', 10, 'XTick', 3+5*(0:3), 'XTickLabel', ...
    {'2-ms ISI', '5-ms ISI', '10-ms ISI', '20-ms ISI'});

