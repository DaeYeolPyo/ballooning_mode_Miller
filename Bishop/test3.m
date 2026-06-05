clc
clear
close all
%% Read GEQDSK and select flux surfaces
%gfile = './curavture_analysis/geqdsk/scan_B2.5_C16_G09_H04.geqdsk';
gfile = './geqdsk_PT0.6';
psiN = 0.3;

eq = read_geqdsk(gfile);

figure;
contour(eq.rgrid, eq.zgrid, eq.psirz.', 40); hold on;

psins = [0.2 0.4 0.6 0.8 1.0];
for k = 1:numel(psins)
    s = extract_flux_surface(eq, psins(k));
    plot(s.R, s.Z, 'LineWidth', 2);
end

plot(eq.rmaxis, eq.zmaxis, 'ko', 'MarkerFaceColor', 'y');
axis equal;
xlabel('R [m]');
ylabel('Z [m]');
title('Selected flux surfaces');

%% 
opts = struct();
opts.N = 1024;

% Rc(l)이 아직 삐죽거리면 30~40으로 낮추고,
% surface shape을 너무 많이 깎는 느낌이면 60~100으로 올리면 됩니다.
opts.smoothCutoff = 100;

% local rho range. GEQDSK length unit이 m이면 rho도 m입니다.
opts.rho = linspace(-0.03, 0.03, 41);

opts.plot = true;

surf = build_local_flux_surface_from_geqdsk(gfile, psiN, opts);

%% 
% 1. 먼저 이전에 만든 surface 생성
optsSurf = struct();
optsSurf.N = 1024;
optsSurf.smoothCutoff = 50;
optsSurf.rho = linspace(-0.03, 0.03, 41);
optsSurf.plot = true;

surf = build_local_flux_surface_from_geqdsk(gfile, psiN, optsSurf);

% 2. Chapter 3 ballooning coefficients 계산
optsBln = struct();
optsBln.plot = true;

% Bp나 dB2/dl이 아직 조금 튀면 cutoff를 낮추기
optsBln.BpSmoothCutoff = 60;
optsBln.derivSmoothCutoff = 60;

bln = compute_bishop_ch3_ballooning_terms(surf, optsBln);