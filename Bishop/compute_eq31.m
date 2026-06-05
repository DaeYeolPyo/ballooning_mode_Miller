function out = compute_eq31(geom, eq29, varargin)
%COMPUTE_EQ31 Compute Bishop Eq. (31) curvature term
%
%   out = COMPUTE_EQ31(geom, eq29)
%   out = COMPUTE_EQ31(geom, eq29, 'Name', value, ...)
%
% Inputs
%   geom : output from build_bishop_geometry
%   eq29 : output from compute_eq29
%
% Name-value options
%   'UseMu0PprimeDrive' : true/false, default false
%       If true, also returns 2*(mu0 p')*Eq31, which is the quantity
%       appearing directly in Eq. (32) under your convention.
%
% Required fields in geom
%   geom.l
%   geom.R
%   geom.h0
%   geom.Bp
%   geom.u
%   geom.R0
%   geom.I0
%   geom.mu0_pprime0
%
% Required fields in eq29
%   eq29.integral_term
%
% Output fields
%   out.l
%   out.B2
%   out.dSdl
%   out.dSdrho
%   out.dQdl
%   out.dQdrho
%   out.curvatureTerm      = ((B x grad S)·K)/B^2
%   out.driveTerm          = 2*(mu0 p')*curvatureTerm   [optional use]
%
% Notes
%   Here Q = 2p + B^2, so that
%       dQ/dl    = d(B^2)/dl     on a flux surface
%       dQ/drho  = 2*Bp^2/R - 2*I^2*sin(u)/R^3
%
%   The implemented formula is
%       ((B x grad S)·K)/B^2
%         = Bp/(2*B^4) * [ dSdrho*dQdl - dSdl*dQdrho ]

    p = inputParser;
    addParameter(p, 'UseMu0PprimeDrive', false, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    %-----------------------------------
    % 1) Pull data
    %-----------------------------------
    l  = geom.l(:);
    R  = geom.R(:);
    h0 = geom.h0(:);
    Bp = geom.Bp(:);
    u  = geom.u(:);

    X0 = geom.R0;
    I0 = geom.I0;

    int29 = eq29.integral_term(:);

    if numel(l) ~= numel(int29)
        error('compute_eq31:SizeMismatch', ...
            'geom.l and eq29.integral_term must have the same length.');
    end

    %-----------------------------------
    % 2) Surface quantities
    %-----------------------------------
    % B^2 on the surface
    B2 = Bp.^2 + (I0^2) ./ (R.^2);

    % dS/dl from B·grad S = 0 on the surface
    dSdl = - I0 ./ (R.^2 .* Bp);

    % dS/drho from Eq. (29) integral term
    dSdrho = (h0 .* Bp .* I0 ./ X0) .* int29;

    %-----------------------------------
    % 3) Derivatives of Q = 2p + B^2
    %-----------------------------------
    % Along the surface p is constant, so dQ/dl = d(B^2)/dl
    dQdl = gradient(B2, l);

    % Normal derivative on the surface
    dQdrho = 2 .* (Bp.^2) ./ R - 2 .* (I0^2) .* sin(u) ./ (R.^3);

    %-----------------------------------
    % 4) Eq. (31)
    %-----------------------------------
    curvatureTerm = (Bp ./ (2 .* B2.^2)) .* ...
        ( dSdrho .* dQdl - dSdl .* dQdrho );

    %-----------------------------------
    % 5) Optional drive term for Eq. (32)
    %-----------------------------------
    if opt.UseMu0PprimeDrive
        driveTerm = 2 .* geom.mu0_pprime0 .* curvatureTerm;
    else
        driveTerm = [];
    end

    %-----------------------------------
    % 6) Output
    %-----------------------------------
    out = struct();
    out.l = l;

    out.B2 = B2;

    out.dSdl = dSdl;
    out.dSdrho = dSdrho;

    out.dQdl = dQdl;
    out.dQdrho = dQdrho;

    out.curvatureTerm = curvatureTerm;
    out.driveTerm = driveTerm;
end