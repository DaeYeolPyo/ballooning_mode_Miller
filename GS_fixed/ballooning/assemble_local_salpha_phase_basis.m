function basis = assemble_local_salpha_phase_basis( ...
        model,theta0,elementsPerPeriod)
%ASSEMBLE_LOCAL_SALPHA_PHASE_BASIS Assemble polynomial P2 matrix bases.
%
% The frozen-geometry coefficients have the exact dependence
%
%   g = g0 + s*gS + s^2*gSS,
%   f = f0 + s*fS + s^2*fSS,
%   c = alpha*cA + alpha^2*cAA + s*alpha*cSA.
%
% Three calls to the generic P2 assembler therefore replace one assembly
% at every point of a large s-alpha scan.

narginchk(3,3);
validateattributes(theta0,{'numeric'}, ...
    {'real','finite','scalar'},mfilename,'theta0');
validateattributes(elementsPerPeriod,{'numeric'}, ...
    {'real','finite','scalar','integer','positive'}, ...
    mfilename,'elementsPerPeriod');
if elementsPerPeriod<2*model.options.poloidalModeCutoff+2
    error('GS:ballooning:AngularGridTooCoarse', ...
        'elementsPerPeriod does not resolve the retained poloidal modes.');
end
add_sturm_liouville_functions();

zero = evaluate_local_salpha_coefficients(model,0,0,theta0);
sPlus = evaluate_local_salpha_coefficients(model,1,0,theta0);
sMinus = evaluate_local_salpha_coefficients(model,-1,0,theta0);
aPlus = evaluate_local_salpha_coefficients(model,0,1,theta0);
aMinus = evaluate_local_salpha_coefficients(model,0,-1,theta0);
crossPoint = evaluate_local_salpha_coefficients(model,1,1,theta0);

nElements = model.options.nPeriods*elementsPerPeriod;
thetaGrid = linspace(-model.options.nPeriods*pi, ...
    model.options.nPeriods*pi,nElements+1).';
sample = @(values,theta) interp1(theta,values,thetaGrid,'pchip');

g0 = sample(zero.g,zero.theta);
gPlus = sample(sPlus.g,sPlus.theta);
gMinus = sample(sMinus.g,sMinus.theta);
gS = 0.5*(gPlus-gMinus);
gSS = 0.5*(gPlus+gMinus)-g0;

f0 = sample(zero.f,zero.theta);
fPlus = sample(sPlus.f,sPlus.theta);
fMinus = sample(sMinus.f,sMinus.theta);
fS = 0.5*(fPlus-fMinus);
fSS = 0.5*(fPlus+fMinus)-f0;

cPlus = sample(aPlus.c,aPlus.theta);
cMinus = sample(aMinus.c,aMinus.theta);
cA = 0.5*(cPlus-cMinus);
cAA = 0.5*(cPlus+cMinus);
cCross = sample(crossPoint.c,crossPoint.theta);
cSA = cCross-cA-cAA;

[G0,CA,F0,thetaDof] = construct_matrix(thetaGrid,g0,cA,f0);
[GS,CAA,FS,thetaDof2] = construct_matrix(thetaGrid,gS,cAA,fS);
[GSS,CSA,FSS,thetaDof3] = construct_matrix(thetaGrid,gSS,cSA,fSS);
if max(abs(thetaDof-thetaDof2))>eps || ...
        max(abs(thetaDof-thetaDof3))>eps
    error('GS:ballooning:InconsistentPhaseBasisGrid', ...
        'Polynomial basis assemblies produced inconsistent DOF grids.');
end

basis = struct();
basis.theta0 = theta0;
basis.thetaGrid = thetaGrid;
basis.thetaDof = thetaDof;
basis.coefficients = struct('g0',g0,'gS',gS,'gSS',gSS, ...
    'f0',f0,'fS',fS,'fSS',fSS, ...
    'cA',cA,'cAA',cAA,'cSA',cSA);
basis.matrices = struct('G0',G0,'GS',GS,'GSS',GSS, ...
    'F0',F0,'FS',FS,'FSS',FSS, ...
    'CA',CA,'CAA',CAA,'CSA',CSA);
end


function add_sturm_liouville_functions()
gsRoot = fileparts(fileparts(mfilename('fullpath')));
folder = fullfile(fileparts(gsRoot),'sturm_liouville_solver');
if ~isfolder(folder)
    error('GS:ballooning:MissingSturmLiouvilleFolder', ...
        'The sibling sturm_liouville_solver folder was not found.');
end
addpath(folder);
end
