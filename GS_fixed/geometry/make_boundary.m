function geom = make_boundary(params)
%MAKE_BOUNDARY Build a sampled fixed-boundary LCFS in the (R,Z) plane.
%
%   GEOM = MAKE_BOUNDARY() returns a circular verification boundary.
%
%   GEOM = MAKE_BOUNDARY(PARAMS) uses an analytic or tabulated LCFS.
%
%   Harmonic (default, or PARAMS.type='harmonic'):
%       R(theta) = A + B*sin(theta) + C*cos(2*theta)
%       Z(theta) = G*cos(theta) - H*sin(2*theta).
%
%   Localized indentation (PARAMS.type='localized-c'):
%       R(theta) = A + B*sin(theta)
%                - D*exp(beta*(sin(theta)-1))
%       Z(theta) = G*cos(theta) - H*sin(2*theta).
%
%   Pointed boomerang (PARAMS.type='boomerang-c'):
%       R_inner(u) = Rtip +(Rend-Rtip )*|u|^pInner
%       R_back (u) = Rback+(Rend-Rback)*|u|^pOuter
%       Z(u) = Zend*u,                         -1 <= u <= 1.
%
%   The inner and back arcs meet exactly at (Rend,+/-Zend), giving pointed
%   upper and lower ends. At Z=0 the two crossings are Rtip and Rback.
%
%   Third-harmonic cosine boundary (PARAMS.type='cosine-c'):
%       R(theta) = R0 + A*cos(theta) + B*cos(2*theta)
%                  + C*cos(3*theta)
%       Z(theta) = kappa*A*sin(theta).
%
%   Miller D boundary (PARAMS.type='miller-d'):
%       R(theta) = R0 + a*cos(theta + asin(delta)*sin(theta)),
%       Z(theta) = kappa*a*sin(theta).
%
%   Tabulated boundary (PARAMS.type='tabulated'):
%       PARAMS.R and PARAMS.Z contain ordered points on a closed LCFS.
%       They are resampled by arclength to PARAMS.nBoundary points.
%
%   Here B<0 moves the upper/lower points outward while moving the average
%   midplane radius inward. The vertical half-height is exactly kappa*A.
%
%   The localized term equals D at the outboard midplane theta=pi/2 and
%   decays exponentially away from it. Thus R_tip=A+B-D, and D>B places the
%   indentation tip inside the nominal major radius A without forcing the
%   inboard boundary inward. The endpoint at 2*pi is not duplicated.
%
%   This function defines only the LCFS. It deliberately does not prescribe
%   any interior flux surfaces; those are determined later by the
%   Grad-Shafranov solve.

narginchk(0, 1);

if nargin == 0 || isempty(params)
    params = struct();
end

if ~isstruct(params) || ~isscalar(params)
    error('GS:geometry:InvalidParameters', ...
        'params must be a scalar struct.');
end

type = boundary_type(params);
defaults = parameter_defaults(type);
unknown = setdiff(fieldnames(params), fieldnames(defaults));
if ~isempty(unknown)
    error('GS:geometry:UnknownParameter', ...
        'Unknown boundary parameter: %s', strjoin(unknown, ', '));
end

names = fieldnames(defaults);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(params, name)
        params.(name) = defaults.(name);
    end
end
params.type = type;

if strcmp(type,'tabulated')
    validateattributes(params.R,{'numeric'}, ...
        {'real','finite','vector'},mfilename,'params.R');
    validateattributes(params.Z,{'numeric'}, ...
        {'real','finite','vector','numel',numel(params.R)}, ...
        mfilename,'params.Z');
    if numel(params.R)<3
        error('GS:geometry:TooFewTabulatedPoints', ...
            'A tabulated boundary requires at least three distinct points.');
    end
elseif strcmp(type, 'miller-d')
    validateattributes(params.R0, {'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'params.R0');
    validateattributes(params.a, {'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'params.a');
    validateattributes(params.kappa, {'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'params.kappa');
    validateattributes(params.delta, {'numeric'}, ...
        {'real','finite','scalar','>',-1,'<',1},mfilename,'params.delta');
elseif strcmp(type, 'boomerang-c')
    validateattributes(params.A, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'params.A');
    validateattributes(params.Rback, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'params.Rback');
    validateattributes(params.Rtip, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'params.Rtip');
    validateattributes(params.Rend, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'params.Rend');
    validateattributes(params.Zend, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'params.Zend');
    validateattributes(params.pInner, {'numeric'}, ...
        {'real', 'finite', 'scalar', '>', 1}, mfilename, 'params.pInner');
    validateattributes(params.pOuter, {'numeric'}, ...
        {'real', 'finite', 'scalar', '>', 1}, mfilename, 'params.pOuter');
    if ~(params.Rback < params.Rtip && params.Rtip < params.Rend)
        error('GS:geometry:InvalidBoomerangRadii', ...
            'boomerang-c requires Rback < Rtip < Rend.');
    end
    % The exponents independently control the indentation and back arcs.
    % Their combination is accepted whenever the two arcs stay separated.
    uProbe = linspace(0, 1, 10001).';
    gapProbe = params.Rtip-params.Rback ...
        +(params.Rend-params.Rtip)*uProbe.^params.pInner ...
        -(params.Rend-params.Rback)*uProbe.^params.pOuter;
    gapTol = 1e-12*max(1, params.Rend-params.Rback);
    if any(gapProbe(1:end-1) <= gapTol)
        error('GS:geometry:CrossingBoomerangArcs', ...
            ['The boomerang inner and back arcs cross before the pointed ' ...
             'ends. Change the radii or exponents.']);
    end
elseif strcmp(type, 'cosine-c')
    validateattributes(params.A, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'params.A');
    validateattributes(params.R0, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'params.R0');
    validateattributes(params.B, {'numeric'}, ...
        {'real', 'finite', 'scalar'}, mfilename, 'params.B');
    validateattributes(params.C, {'numeric'}, ...
        {'real', 'finite', 'scalar'}, mfilename, 'params.C');
    validateattributes(params.kappa, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'params.kappa');
else
    validateattributes(params.A, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'params.A');
    validateattributes(params.B, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'params.B');
    validateattributes(params.G, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, 'params.G');
    validateattributes(params.H, {'numeric'}, ...
        {'real', 'finite', 'scalar'}, mfilename, 'params.H');
    if strcmp(type, 'harmonic')
        validateattributes(params.C, {'numeric'}, ...
            {'real', 'finite', 'scalar'}, mfilename, 'params.C');
    else
        validateattributes(params.D, {'numeric'}, ...
            {'real', 'finite', 'scalar', 'nonnegative'}, mfilename, 'params.D');
        validateattributes(params.beta, {'numeric'}, ...
            {'real', 'finite', 'scalar', 'positive'}, mfilename, 'params.beta');
    end
end
validateattributes(params.nBoundary, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', '>=', 16}, ...
    mfilename, 'params.nBoundary');
if strcmp(type, 'boomerang-c') && mod(params.nBoundary, 2) ~= 0
    error('GS:geometry:OddBoomerangResolution', ...
        'boomerang-c requires an even nBoundary.');
end

if strcmp(type,'tabulated')
    [R,Z] = resample_closed_curve(params.R,params.Z,params.nBoundary);
    theta = linspace(0,2*pi,params.nBoundary+1).';
    theta(end) = [];
    outboardMidplaneR = max(R);
    inboardMidplaneR = min(R);
    nominalMajorRadius = 0.5*(outboardMidplaneR+inboardMidplaneR);
    [endPointZ,topIndex] = max(Z);
    endPointR = R(topIndex);
    indentationDepth = 0;
    hasPointedEnds = false;
    outboardMidplaneCurvature = NaN;
    outboardMidplaneCurvatureRadius = NaN;
elseif strcmp(type, 'miller-d')
    theta = linspace(0,2*pi,params.nBoundary+1).';
    theta(end) = [];
    triangularAngle = asin(params.delta);
    R = params.R0+params.a*cos( ...
        theta+triangularAngle*sin(theta));
    Z = params.kappa*params.a*sin(theta);
    indentationDepth = params.a*params.delta;
    outboardMidplaneR = params.R0+params.a;
    inboardMidplaneR = params.R0-params.a;
    endPointR = params.R0-params.a*params.delta;
    endPointZ = params.kappa*params.a;
    hasPointedEnds = false;
    outboardMidplaneCurvature = NaN;
    outboardMidplaneCurvatureRadius = NaN;
elseif strcmp(type, 'boomerang-c')
    % Traverse the indentation arc from the upper point to the lower point,
    % then the back arc in the opposite direction. Endpoints occur once.
    nArc = params.nBoundary/2 + 1;
    uInner = linspace(1, -1, nArc).';
    RInner = params.Rtip + (params.Rend-params.Rtip) ...
        .*abs(uInner).^params.pInner;
    ZInner = params.Zend*uInner;

    uBack = linspace(-1, 1, nArc).';
    RBack = params.Rback + (params.Rend-params.Rback) ...
        .*abs(uBack).^params.pOuter;
    ZBack = params.Zend*uBack;

    R = [RInner; RBack(2:end-1)];
    Z = [ZInner; ZBack(2:end-1)];
    theta = linspace(0, 2*pi, params.nBoundary + 1).';
    theta(end) = [];
    indentationDepth = params.Rend-params.Rtip;
    outboardMidplaneR = params.Rtip;
    inboardMidplaneR = params.Rback;
    endPointR = params.Rend;
    endPointZ = params.Zend;
    hasPointedEnds = true;
    outboardMidplaneCurvature = NaN;
    outboardMidplaneCurvatureRadius = NaN;
elseif strcmp(type, 'cosine-c')
    theta = linspace(0, 2*pi, params.nBoundary + 1).';
    theta(end) = [];
    R = params.R0 + params.A*cos(theta) + params.B*cos(2*theta) ...
        + params.C*cos(3*theta);
    Z = params.kappa*params.A*sin(theta);
    midplaneR = [params.R0+params.A+params.B+params.C, ...
                 params.R0-params.A+params.B-params.C];
    outboardMidplaneR = max(midplaneR);
    inboardMidplaneR = min(midplaneR);
    endPointR = params.R0-params.B;
    endPointZ = params.kappa*params.A;
    indentationDepth = endPointR-outboardMidplaneR;
    hasPointedEnds = abs(params.C-params.A/3) ...
        <= 1e-10*max(1, params.A);
    if params.A+params.C >= 0
        radialSecondDerivative = -params.A-4*params.B-9*params.C;
    else
        radialSecondDerivative = params.A-4*params.B+9*params.C;
    end
    outboardMidplaneCurvature = abs(radialSecondDerivative) ...
        /(params.kappa*params.A)^2;
    if outboardMidplaneCurvature > 0
        outboardMidplaneCurvatureRadius = ...
            1/outboardMidplaneCurvature;
    else
        outboardMidplaneCurvatureRadius = Inf;
    end
else
    theta = linspace(0, 2*pi, params.nBoundary + 1).';
    theta(end) = [];
    if strcmp(type, 'harmonic')
        R = params.A + params.B*sin(theta) + params.C*cos(2*theta);
        indentationDepth = params.C;
        outboardMidplaneR = params.A+params.B-params.C;
        inboardMidplaneR = params.A-params.B-params.C;
    else
        localizedWeight = exp(params.beta*(sin(theta)-1));
        R = params.A + params.B*sin(theta)-params.D*localizedWeight;
        indentationDepth = params.D;
        outboardMidplaneR = params.A+params.B-params.D;
        inboardMidplaneR = params.A-params.B ...
                          -params.D*exp(-2*params.beta);
    end
    Z = params.G*cos(theta) - params.H*sin(2*theta);
    endPointR = NaN;
    endPointZ = NaN;
    hasPointedEnds = false;
    outboardMidplaneCurvature = NaN;
    outboardMidplaneCurvatureRadius = NaN;
end
points = [R, Z];

if any(~isfinite(points(:)))
    error('GS:geometry:NonFiniteBoundary', ...
        'The boundary contains a non-finite coordinate.');
end
if min(R) <= 0
    error('GS:geometry:NonPositiveRadius', ...
        ['The LCFS reaches R <= 0. The Grad-Shafranov weak form contains ' ...
         '1/R, so the complete domain must remain at positive R.']);
end

scale = max([1; max(R)-min(R); max(Z)-min(Z)]);
coordTol = 1e-11*scale;
if has_self_intersection(points, coordTol)
    error('GS:geometry:SelfIntersection', ...
        ['The sampled LCFS intersects itself. Reduce the shaping or use a ' ...
         'different boundary parameter set.']);
end

closed = [points; points(1,:)];
edgeCross = closed(1:end-1,1).*closed(2:end,2) ...
          - closed(2:end,1).*closed(1:end-1,2);
signedArea = 0.5*sum(edgeCross);
areaTol = 100*eps(scale^2)*params.nBoundary;
if abs(signedArea) <= areaTol
    error('GS:geometry:DegenerateBoundary', ...
        'The sampled LCFS encloses negligible area.');
end

centroidR = sum((closed(1:end-1,1) + closed(2:end,1)).*edgeCross) ...
          / (6*signedArea);
centroidZ = sum((closed(1:end-1,2) + closed(2:end,2)).*edgeCross) ...
          / (6*signedArea);
edgeDelta = diff(closed, 1, 1);
edgeLengths = sqrt(sum(edgeDelta.^2, 2));

turn = polygon_turns(points);
turnTol = 1e-10*scale^2;
hasPositiveTurn = any(turn > turnTol);
hasNegativeTurn = any(turn < -turnTol);

geom = struct();
geom.type = type;
geom.params = params;
geom.theta = theta;
geom.R = R;
geom.Z = Z;
geom.points = points;
geom.closedPoints = closed;
geom.nBoundary = params.nBoundary;
geom.signedArea = signedArea;
geom.area = abs(signedArea);
geom.perimeter = sum(edgeLengths);
geom.centroid = [centroidR, centroidZ];
geom.characteristicLength = sqrt(4*geom.area/pi);
geom.RLim = [min(R), max(R)];
geom.ZLim = [min(Z), max(Z)];
geom.radialWidth = geom.RLim(2)-geom.RLim(1);
geom.verticalHeight = geom.ZLim(2)-geom.ZLim(1);
geom.boundingBoxElongation = geom.verticalHeight/geom.radialWidth;
if strcmp(type,'tabulated')
    geom.nominalMajorRadius = nominalMajorRadius;
    geom.nominalElongation = geom.boundingBoxElongation;
elseif strcmp(type, 'cosine-c')
    geom.nominalMajorRadius = params.R0;
    geom.nominalElongation = params.kappa;
elseif strcmp(type,'miller-d')
    geom.nominalMajorRadius = params.R0;
    geom.nominalElongation = params.kappa;
else
    geom.nominalMajorRadius = params.A;
    geom.nominalElongation = NaN;
end
geom.indentationDepth = indentationDepth;
geom.outboardMidplaneR = outboardMidplaneR;
geom.inboardMidplaneR = inboardMidplaneR;
geom.tipInsideNominalMajorRadius = ...
    outboardMidplaneR < geom.nominalMajorRadius;
geom.endPointR = endPointR;
geom.endPointZ = endPointZ;
geom.hasPointedEnds = hasPointedEnds;
geom.outboardMidplaneCurvature = outboardMidplaneCurvature;
geom.outboardMidplaneCurvatureRadius = ...
    outboardMidplaneCurvatureRadius;
geom.orientation = orientation_name(signedArea);
geom.isConcave = hasPositiveTurn && hasNegativeTurn;

end


function type = boundary_type(params)
if isfield(params, 'type')
    value = params.type;
else
    value = 'harmonic';
end
if isstring(value) && isscalar(value)
    value = char(value);
end
if ~ischar(value) || ~isrow(value)
    error('GS:geometry:InvalidBoundaryType', ...
        ['params.type must be harmonic, localized-c, boomerang-c, ' ...
        'cosine-c, miller-d, or tabulated.']);
end
type = lower(strrep(strtrim(value), '_', '-'));
if ismember(type, {'localized', 'local-c'})
    type = 'localized-c';
end
if ismember(type, {'boomerang', 'crescent'})
    type = 'boomerang-c';
end
if ismember(type, {'cosine', 'fourier-c', 'harmonic3', 'third-harmonic'})
    type = 'cosine-c';
end
if ismember(type, {'miller','d-shape','dshape','miller-dshape'})
    type = 'miller-d';
end
if ismember(type,{'points','geqdsk','sampled'})
    type = 'tabulated';
end
if ~ismember(type, ...
        {'harmonic','localized-c','boomerang-c','cosine-c','miller-d', ...
         'tabulated'})
    error('GS:geometry:InvalidBoundaryType', ...
        ['params.type must be harmonic, localized-c, boomerang-c, ' ...
        'cosine-c, miller-d, or tabulated.']);
end
end


function defaults = parameter_defaults(type)
if strcmp(type,'tabulated')
    angle = linspace(0,2*pi,257).';
    angle(end) = [];
    defaults = struct('type',type,'R',3+cos(angle), ...
        'Z',sin(angle),'nBoundary',256);
elseif strcmp(type, 'miller-d')
    defaults = struct('type',type,'R0',3.0,'a',0.9, ...
        'kappa',1.66,'delta',0.416,'nBoundary',256);
elseif strcmp(type, 'boomerang-c')
    defaults = struct('type', type, 'A', 3.0, ...
        'Rback', 2.0, 'Rtip', 2.65, 'Rend', 3.8, 'Zend', 2.7, ...
        'pInner', 1.6, 'pOuter', 1.4, 'nBoundary', 256);
elseif strcmp(type, 'cosine-c')
    defaults = struct('type', type, 'R0', 3.0, 'A', 0.9, ...
        'B', -0.6, 'C', -0.2, 'kappa', 3.0, 'nBoundary', 256);
else
    defaults = struct('type', type, 'A', 3.0, 'B', 1.0, ...
        'G', 1.0, 'H', 0.0, 'nBoundary', 256);
    if strcmp(type, 'harmonic')
        defaults.C = 0.0;
    else
        defaults.D = 0.0;
        defaults.beta = 5.0;
    end
end
end


function [Rout,Zout] = resample_closed_curve(R,Z,nBoundary)
points = [R(:),Z(:)];
scale = max([1;max(points,[],1).'-min(points,[],1).']);
tolerance = 1e-12*scale;
if norm(points(end,:)-points(1,:))<=tolerance
    points(end,:) = [];
end
if size(points,1)<3
    error('GS:geometry:TooFewTabulatedPoints', ...
        'A tabulated boundary requires at least three distinct points.');
end
if size(points,1)==nBoundary
    Rout = points(:,1);
    Zout = points(:,2);
    return
end
closed = [points;points(1,:)];
edgeLength = hypot(diff(closed(:,1)),diff(closed(:,2)));
if any(edgeLength<=tolerance)
    error('GS:geometry:DuplicateTabulatedPoints', ...
        'Tabulated boundary points must not contain adjacent duplicates.');
end
arc = [0;cumsum(edgeLength)];
query = (0:nBoundary-1).'*arc(end)/nBoundary;
Rout = interp1(arc,closed(:,1),query,'linear');
Zout = interp1(arc,closed(:,2),query,'linear');
end


function name = orientation_name(signedArea)
if signedArea > 0
    name = 'counterclockwise';
else
    name = 'clockwise';
end
end


function turn = polygon_turns(points)
n = size(points, 1);
turn = zeros(n, 1);
for i = 1:n
    iPrev = mod(i - 2, n) + 1;
    iNext = mod(i, n) + 1;
    incoming = points(i,:) - points(iPrev,:);
    outgoing = points(iNext,:) - points(i,:);
    turn(i) = cross2(incoming, outgoing);
end
end


function tf = has_self_intersection(points, coordTol)
n = size(points, 1);
tf = false;
for i = 1:n
    iNext = mod(i, n) + 1;
    p1 = points(i,:);
    p2 = points(iNext,:);

    for j = i+1:n
        jNext = mod(j, n) + 1;

        % Adjacent polygon edges share one endpoint by construction.
        if j == iNext || (i == 1 && jNext == 1)
            continue
        end

        q1 = points(j,:);
        q2 = points(jNext,:);
        if segments_intersect(p1, p2, q1, q2, coordTol)
            tf = true;
            return
        end
    end
end
end


function tf = segments_intersect(p1, p2, q1, q2, coordTol)
crossTol = coordTol*max([1, norm(p2-p1), norm(q2-q1)]);
o1 = cross2(p2-p1, q1-p1);
o2 = cross2(p2-p1, q2-p1);
o3 = cross2(q2-q1, p1-q1);
o4 = cross2(q2-q1, p2-q1);

proper = ((o1 > crossTol && o2 < -crossTol) || ...
          (o1 < -crossTol && o2 > crossTol)) && ...
         ((o3 > crossTol && o4 < -crossTol) || ...
          (o3 < -crossTol && o4 > crossTol));

touching = (abs(o1) <= crossTol && on_segment(p1, p2, q1, coordTol)) || ...
           (abs(o2) <= crossTol && on_segment(p1, p2, q2, coordTol)) || ...
           (abs(o3) <= crossTol && on_segment(q1, q2, p1, coordTol)) || ...
           (abs(o4) <= crossTol && on_segment(q1, q2, p2, coordTol));

tf = proper || touching;
end


function tf = on_segment(a, b, p, tol)
tf = p(1) >= min(a(1), b(1)) - tol && ...
     p(1) <= max(a(1), b(1)) + tol && ...
     p(2) >= min(a(2), b(2)) - tol && ...
     p(2) <= max(a(2), b(2)) + tol;
end


function value = cross2(a, b)
value = a(1)*b(2) - a(2)*b(1);
end
