function shape = normalize_indented_cshape_boundary(aspectRatio, kappa, a, delta, indent, opts)
%NORMALIZE_INDENTED_CSHAPE_BOUNDARY Preserve actual box dimensions.
%   a, kappa and aspectRatio are measured geometric values, not the
%   construction parameters of Indented_Miller_boundary. indent is I0.
%   The returned boundary has no duplicate periodic endpoint and is CCW.
    arguments
        aspectRatio (1,1) double {mustBeFinite,mustBeGreaterThan(aspectRatio,1)}
        kappa (1,1) double {mustBeFinite,mustBePositive}
        a (1,1) double {mustBeFinite,mustBePositive}
        delta (1,1) double {mustBeFinite}
        indent (1,1) double {mustBeFinite,mustBeNonnegative}
        opts.NBoundary (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(opts.NBoundary,16)} = 128
        opts.ExtremumSamples (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(opts.ExtremumSamples,128)} = 4096
    end
    setup_indented_cshape_scan;
    if abs(delta)>=1
        error('cshape:Boundary:Delta','Require abs(delta) < 1.');
    end
    x = asin(delta);
    raw = @(t) cos(t+x*sin(t))-indent*((1+cos(t))/2).^3;
    % For each Z, the two branches must stay distinct. Their radial gap
    % is a0*cos(t)*[2*cos(x*sin(t))-I*(3+cos(t)^2)/4]. A zero at the
    % tips also gives a degenerate tangent, so it is excluded.
    limit = @(t) 8*cos(x*sin(t))./(3+cos(t).^2);
    [indentLimit,~] = extrema_periodic(limit,opts.ExtremumSamples);
    if indent >= indentLimit-1.e-10*max(1,indentLimit)
        error('cshape:Boundary:Topology', ...
            'Indentation %.6g reaches the self-contact/tangent limit %.6g.',indent,indentLimit);
    end
    [rmin,rmax] = extrema_periodic(raw,opts.ExtremumSamples);
    s = (rmax-rmin)/2;
    m = (rmax+rmin)/2;
    a0 = a/s;
    Rc = aspectRatio*a-a0*m;
    kappa0 = kappa*s;
    bnd = Indented_Miller_boundary(Rc,a0,delta,kappa0,indent, ...
        NTheta=int32(opts.NBoundary+1));
    names = fieldnames(bnd);
    for j=1:numel(names)
        value = bnd.(names{j});
        if isvector(value) && numel(value)==opts.NBoundary+1
            value=value(1:end-1);
            bnd.(names{j})=value(:);
        end
    end
    if min(bnd.R)<=0 || min(hypot(bnd.Rt,bnd.Zt))<=1.e-12*a
        error('cshape:Boundary:InvalidCurve','Require R > 0 and a nonzero tangent.');
    end
    area = sum(bnd.R.*circshift(bnd.Z,-1)-circshift(bnd.R,-1).*bnd.Z)/2;
    if area<=0
        error('cshape:Boundary:Orientation','Expected a CCW simple boundary.');
    end
    shape = struct('boundary',bnd,'aspectRatio',aspectRatio,'kappa',kappa, ...
        'aMinor',a,'Rgeo',aspectRatio*a,'delta',delta,'indent',indent, ...
        'Rmin',(aspectRatio-1)*a,'Rmax',(aspectRatio+1)*a, ...
        'bMinor',kappa*a,'indentationSide','outboard', ...
        'iOut',(rmax-(1-indent))/(2*s), ...
        'dOut',a0*(rmax-(1-indent)), ...
        'deltaGeometric',(m+delta+indent/8)/s, ...
        'midplaneConcave',3*indent/4>(1+x)^2, ...
        'construction',struct('R0',Rc,'a',a0,'kappa',kappa0), ...
        'rawHalfWidth',s,'rawMidpoint',m,'indentationLimit',indentLimit, ...
        'polygonArea',area);
end

function [lo,hi] = extrema_periodic(fun,n)
    step = 2*pi/n;
    t = (0:n-1)'*step;
    y = fun(t);
    lower = find(y<=circshift(y,1) & y<=circshift(y,-1));
    upper = find(y>=circshift(y,1) & y>=circshift(y,-1));
    lo = min(y); hi = max(y);
    options = optimset('TolX',1.e-12,'Display','off');
    for j=lower(:).'
        [~,v] = fminbnd(fun,t(j)-step,t(j)+step,options);
        lo = min(lo,v);
    end
    for j=upper(:).'
        [~,v] = fminbnd(@(z)-fun(z),t(j)-step,t(j)+step,options);
        hi = max(hi,-v);
    end
end
