clc
clear
close all

ntheta = 300;
ngeom = 600;

p = DshapeMillerParams();
p.F = p.A*p.r*p.B0;

theta = linspace(0, 2*pi, ntheta+1).';
theta(end) = [];

[Bp, raw] = Bpol_Dshape(p, ngeom);
u = raw.u(:);
R = raw.R(:);
Z = raw.Z(:);
Ru = raw.Ru(:);
Zu = raw.Zu(:);
jac = raw.jac(:);
Bp = Bp(:);

pitch = raw.F.*abs(jac)./(raw.dpdr.*R);
theta_raw = cumtrapz(u, pitch)./p.q;
theta_raw = theta_raw - theta_raw(1);
theta_raw = 2*pi*theta_raw./theta_raw(end);

keep = [true; diff(theta_raw) > 1e-12];
theta_u = theta_raw(keep);
keep_idx = find(keep);
if numel(theta_u) > 1 && abs(theta_u(end) - theta_u(1) - 2*pi) < 1e-10
    keep_idx(end) = [];
    theta_u(end) = [];
end
keep = false(size(theta_raw));
keep(keep_idx) = true;

R = R(keep);
Z = Z(keep);
Ru = Ru(keep);
Zu = Zu(keep);
jac = jac(keep);
Bp = Bp(keep);

theta_ext = [theta_u-2*pi; theta_u; theta_u+2*pi];
R_ext = [R; R; R];
Z_ext = [Z; Z; Z];
Bp_ext = [Bp; Bp; Bp];
Ru_ext = [Ru; Ru; Ru];
Zu_ext = [Zu; Zu; Zu];
jac_ext = [jac; jac; jac];

tq = mod(theta(:), 2*pi);
surf = struct();
surf.theta = theta(:);
surf.R = interp1(theta_ext, R_ext, tq, 'pchip');
surf.Z = interp1(theta_ext, Z_ext, tq, 'pchip');
surf.Bp = interp1(theta_ext, Bp_ext, tq, 'pchip');
surf.Ru = interp1(theta_ext, Ru_ext, tq, 'pchip');
surf.Zu = interp1(theta_ext, Zu_ext, tq, 'pchip');
surf.jac = interp1(theta_ext, jac_ext, tq, 'pchip');
surf.F = raw.F;
surf.dpdr = raw.dpdr;
surf.q = p.q;
surf.Bphi = raw.F./surf.R;
surf.B2 = surf.Bp.^2 + surf.Bphi.^2;
surf.raw = raw;