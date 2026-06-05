function out = compute_eq29(geom, varargin)
%COMPUTE_EQ29 Compute |grad S|^2 from Bishop Eq. (29)
%
%   out = COMPUTE_EQ29(geom)
%   out = COMPUTE_EQ29(geom, 'Name', value, ...)
%
% Inputs
%   geom : output from build_bishop_geometry
%
% Name-value options
%   'L0Index' : starting index corresponding to l0, default 1
%   'Closed'  : true/false, default true
%
% Required fields in geom
%   geom.l
%   geom.R
%   geom.h0
%   geom.Bp
%   geom.u
%   geom.R0
%   geom.I0
%   geom.Iprime0
%   geom.mu0_pprime0
%
% Output fields
%   out.l
%   out.integrand
%   out.integral_term
%   out.term1
%   out.term2
%   out.absGradS2
%
% Notes
% - Uses mu0*p' through geom.mu0_pprime0.
% - For a closed surface, the arrays are circularly shifted so that l0 is first.
% - The returned out.l starts at 0 at the chosen l0.

    p = inputParser;
    addParameter(p, 'L0Index', 1, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'Closed', true, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    %-----------------------------------
    % 1) Pull required fields
    %-----------------------------------
    l  = geom.l(:);
    R  = geom.R(:);
    h0 = geom.h0(:);
    Bp = geom.Bp(:);
    u  = geom.u(:);

    X0 = geom.R0;
    I0 = geom.I0;
    Iprime0 = geom.Iprime0;
    mu0_pprime0 = geom.mu0_pprime0;

    n = numel(l);
    i0 = round(opt.L0Index);

    if i0 < 1 || i0 > n
        error('compute_eq29:BadL0Index', ...
            'L0Index must be between 1 and %d.', n);
    end

    %-----------------------------------
    % 2) Reorder so that l0 is first
    %-----------------------------------
    if opt.Closed
        idx = [i0:n, 1:i0-1];
        l  = l(idx);
        R  = R(idx);
        h0 = h0(idx);
        Bp = Bp(idx);
        u  = u(idx);

        % rebuild l starting from 0
        dl = hypot(diff(R), diff(geom.Z(idx)));
        l = [0; cumsum(dl)];
    else
        l = l - l(i0);
        if any(diff(l) <= 0)
            error('compute_eq29:OpenCurveOrdering', ...
                'For open curves, l must increase monotonically after shifting.');
        end
    end

    %-----------------------------------
    % 3) Common pieces
    %-----------------------------------
    Bp2 = Bp.^2;
    h02 = h0.^2;

    A = 1 + (I0^2) ./ (X0^2 .* h02 .* Bp2);

    % integrand inside the curly braces
    integrand = (1 ./ (h02 .* Bp2)) .* ( ...
          (Iprime0 / I0) .* A ...
        + mu0_pprime0 ./ Bp2 ...
        + (2 ./ R) .* (1 + (R .* sin(u)) ./ (X0 .* h0)) ...
          .* (1 ./ (X0 .* h0 .* Bp)) );

    %-----------------------------------
    % 4) Integral from l0 to l
    %-----------------------------------
    integral_term = cumtrapz(l, integrand);

    %-----------------------------------
    % 5) Eq. (29)
    %-----------------------------------
    term1 = (1 ./ (X0^2 .* h02)) .* A;

    term2 = ((h02 .* Bp2 .* I0^2) ./ X0^2) .* (integral_term.^2);

    absGradS2 = term1 + term2;

    %-----------------------------------
    % 6) Output
    %-----------------------------------
    out = struct();
    out.l = l;

    out.integrand = integrand;
    out.integral_term = integral_term;

    out.term1 = term1;
    out.term2 = term2;

    out.absGradS2 = absGradS2;
end