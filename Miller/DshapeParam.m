function [R, Z, u, p] = DshapeParam(p, ntheta)
    p = fill_Dshape_defaults(p);

    u = linspace(0, 2*pi, ntheta);
    r = p.r;
    R0 = p.A*r;
    x = asin(p.delta);

    R = R0 + r*cos(u + x*sin(u));
    Z = p.kappa*r*sin(u);
end

function p = fill_Dshape_defaults(p)
    if ~isfield(p, 'r'), p.r = 1.0; end
    if ~isfield(p, 'B0'), p.B0 = 1.0; end
end
