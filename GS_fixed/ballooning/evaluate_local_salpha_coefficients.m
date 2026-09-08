function coeff = evaluate_local_salpha_coefficients( ...
        model,sHat,alpha,theta0)
%EVALUATE_LOCAL_SALPHA_COEFFICIENTS Evaluate frozen-geometry PEST g,c,f.

narginchk(4,4);
validateattributes(sHat,{'numeric'}, ...
    {'real','finite','scalar'},mfilename,'sHat');
validateattributes(alpha,{'numeric'}, ...
    {'real','finite','scalar'},mfilename,'alpha');
validateattributes(theta0,{'numeric'}, ...
    {'real','finite','scalar'},mfilename,'theta0');
if ~isstruct(model) || ~isscalar(model) || ...
        ~all(isfield(model,{'base','qprimePerSHat', ...
                           'pprimePerAlpha','targetIndex'}))
    error('GS:ballooning:InvalidSAlphaModel', ...
        'model must be returned by prepare_local_salpha_model.');
end

map = model.base.normalizedMap;
k = model.targetIndex;
qprime = sHat*model.qprimePerSHat;
pprime = alpha*model.pprimePerAlpha;

% A local linear q profile supplies exactly the requested derivative at
% the target while leaving q itself fixed. Only target pprime is used by
% the coefficient formula, so a constant radial profile is sufficient.
map.q_profile = model.q+qprime*(map.psi-map.psi(k));
map.pprime = pprime*ones(size(map.pprime));
thetaExtent = model.base.coefficients.theta;

coeff = calculate_ballooning_coeffs( ...
    map,model.base.metrics,model.targetPsiN, ...
    Theta0=theta0,ThetaExt=thetaExtent);
coeff.sHat = sHat;
coeff.alpha = alpha;
coeff.scanQprime = qprime;
coeff.scanPprime = pprime;
coeff.geometryAssumption = model.assumption;
end
