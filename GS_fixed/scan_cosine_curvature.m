function result = scan_cosine_curvature(opts)
%SCAN_COSINE_CURVATURE Search cosine-C shapes for reduced bad curvature.
%
%   RESULT = SCAN_COSINE_CURVATURE() scans a nonlinear fixed-boundary GS
%   equilibrium over A, B, C/A while choosing kappa so every candidate has
%   the requested actual bounding-box elongation. Invalid, low-quality,
%   asymmetric, multi-axis, X-point, and non-nested candidates are rejected.
%
%   With the main_cos convention p'(psi)<0 and psi_axis<psi_LCFS,
%   kappa dot grad(p)>0 is classified as bad curvature. The primary score is
%   the positive fraction among sign-active toroidal volume, excluding the
%   near-zero pressure-curvature region. RESULT.ranking is sorted ascending
%   by this metric.
%
%   Supported OPTS fields:
%       AValues                    default 0.28:0.04:0.48
%       BValues                    default -(0.50:0.05:0.75)
%       COverAValues               default [-0.1,0,0.1]
%       R0                         default 3
%       targetElongation           actual bounding-box value; default 3
%       maximumRtip                default 2.85
%       minimumR                   default 0.5
%       maximumIndentationRadius   default 3.5
%       nBoundary                  default 120
%       targetH                    default 0.17
%       minimumMeshQuality         default 0.1
%       maximumAxisAbsZ            default 0.05
%       verbose                    default true
%
%   This is a screening scan. Re-run the leading candidates with main_cos
%   resolution, strict Picard tolerances, and continuation before selection.

narginchk(0,1);
if nargin<1 || isempty(opts)
    opts = struct();
end
opts = validate_scan_options(opts);
add_project_paths();

profiles = make_profiles(struct( ...
    'pprime',struct('type','normalized-power', ...
                    'axisValue',-8e5,'exponent',2), ...
    'FFprime',struct('type','normalized-power', ...
                     'axisValue',-0.8,'exponent',1)));
picard = struct('tolerance',3e-6,'residualTolerance',1e-5, ...
    'maxIterations',800,'omega',0.22,'boundaryValue',0, ...
    'axisMode','min','updateFluxAxis',true, ...
    'failOnNonconvergence',false,'verbose',false);

variableNames = {'A','B','C','kappa','COverA','Rtip','Rback', ...
    'Rend','Zend','actualElongation','indentationCurvatureRadius', ...
    'Rmin','Rmax','minimumMeshQuality','nNodes','nIterations', ...
    'finalResidual','axisR','axisAbsZ','badActiveVolumeFraction', ...
    'goodActiveVolumeFraction','badVolumeFraction', ...
    'goodVolumeFraction','nearZeroVolumeFraction', ...
    'volumeWeightedMeanKappaDotGradP'};
nCombinations = numel(opts.AValues)*numel(opts.BValues) ...
    *numel(opts.COverAValues);
rows = NaN(nCombinations,numel(variableNames));
nAccepted = 0;
counters = struct('total',0,'validGeometry',0,'converged',0, ...
    'symmetric',0,'topologyPassed',0);

for A = opts.AValues(:).'
    for B = opts.BValues(:).'
        for cOverA = opts.COverAValues(:).'
            counters.total = counters.total+1;
            C = cOverA*A;
            try
                unitParams = struct('type','cosine-c','R0',opts.R0, ...
                    'A',A,'B',B,'C',C,'kappa',1, ...
                    'nBoundary',opts.nBoundary);
                unitGeom = make_boundary(unitParams);
                kappa = opts.targetElongation*unitGeom.radialWidth/(2*A);
                params = unitParams;
                params.kappa = kappa;
                geom = make_boundary(params);
                if ~geom.isConcave || ~geom.tipInsideNominalMajorRadius || ...
                        geom.outboardMidplaneR>opts.maximumRtip || ...
                        geom.RLim(1)<opts.minimumR || ...
                        geom.outboardMidplaneCurvatureRadius ...
                            >opts.maximumIndentationRadius
                    continue
                end
                mesh = generate_mesh(geom,struct( ...
                    'targetH',opts.targetH,'boundaryClearance',0.12));
                if mesh.quality.minimum<opts.minimumMeshQuality
                    continue
                end
                counters.validGeometry = counters.validGeometry+1;

                K = assemble_stiffness(mesh);
                psi0 = initial_guess(mesh,K,struct( ...
                    'boundaryValue',0,'axisValue',-1));
                [psi,solver] = solve_picard( ...
                    mesh,profiles,psi0,picard);
                if ~solver.converged
                    continue
                end
                counters.converged = counters.converged+1;
                axisData = find_axis(mesh,psi,struct( ...
                    'mode','min','boundaryValue',0));
                if abs(axisData.Z)>opts.maximumAxisAbsZ
                    continue
                end
                counters.symmetric = counters.symmetric+1;
                topology = check_topology(mesh,psi,axisData,struct( ...
                    'normalizedLevels',(0.1:0.1:0.9).', ...
                    'boundaryValue',0,'failOnFailure',false));
                if ~topology.passed
                    continue
                end
                counters.topologyPassed = counters.topologyPassed+1;

                span = -axisData.psi;
                psiN = min(max((psi-axisData.psi)/span,0),1);
                F = sqrt(3^2+0.8*span.*(1-psiN).^2);
                field = compute_B(mesh,psi,F);
                curvature = compute_magnetic_curvature( ...
                    mesh,psi,field,profiles,struct( ...
                        'psiAxis',axisData.psi,'boundaryValue',0));
                stats = curvature.statistics;
                badActive = stats.positiveActiveVolumeFraction;
                goodActive = stats.negativeActiveVolumeFraction;

                nAccepted = nAccepted+1;
                rows(nAccepted,:) = [A,B,C,kappa,cOverA, ...
                    geom.outboardMidplaneR,geom.inboardMidplaneR, ...
                    geom.endPointR,geom.endPointZ, ...
                    geom.boundingBoxElongation, ...
                    geom.outboardMidplaneCurvatureRadius,geom.RLim, ...
                    mesh.quality.minimum,mesh.nNodes,solver.nIterations, ...
                    solver.finalResidual,axisData.R,abs(axisData.Z), ...
                    badActive,goodActive,stats.positiveVolumeFraction, ...
                    stats.negativeVolumeFraction, ...
                    stats.nearZeroVolumeFraction, ...
                    stats.volumeWeightedMeanKappaDotGradP];
            catch exception
                if opts.verbose
                    fprintf('Rejected A=%.4g B=%.4g C=%.4g: %s\n', ...
                        A,B,C,exception.identifier)
                end
            end
        end
    end
end

rows = rows(1:nAccepted,:);
ranking = array2table(rows,'VariableNames',variableNames);
if ~isempty(ranking)
    ranking = sortrows(ranking,'badActiveVolumeFraction','ascend');
end
result = struct('ranking',ranking,'counters',counters, ...
    'options',opts,'profiles',profiles,'picardOptions',picard, ...
    'signConvention',[ ...
        'For the main_cos orientation, positive kappa dot grad(p) is ' ...
        'classified as bad curvature.']);

if opts.verbose
    fprintf(['Curvature scan: %d total, %d valid geometry/mesh, ' ...
             '%d converged, %d symmetric, %d topology-passed.\n'], ...
        counters.total,counters.validGeometry,counters.converged, ...
        counters.symmetric,counters.topologyPassed)
    disp(ranking(1:min(20,height(ranking)),:))
end
end


function opts = validate_scan_options(opts)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:scan:InvalidOptions','opts must be a scalar struct.');
end
defaults = struct('AValues',0.28:0.04:0.48, ...
    'BValues',-(0.50:0.05:0.75), ...
    'COverAValues',[-0.1,0,0.1], ...
    'R0',3,'targetElongation',3,'maximumRtip',2.85, ...
    'minimumR',0.5,'maximumIndentationRadius',3.5, ...
    'nBoundary',120,'targetH',0.17,'minimumMeshQuality',0.1, ...
    'maximumAxisAbsZ',0.05,'verbose',true);
unknown = setdiff(fieldnames(opts),fieldnames(defaults));
if ~isempty(unknown)
    error('GS:scan:UnknownOption', ...
        'Unknown scan option: %s',strjoin(unknown,', '));
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts,names{k})
        opts.(names{k}) = defaults.(names{k});
    end
end
for name = {'AValues','BValues','COverAValues'}
    validateattributes(opts.(name{1}),{'numeric'}, ...
        {'real','finite','vector','nonempty'},mfilename,['opts.' name{1}]);
end
for name = {'R0','targetElongation','maximumRtip','minimumR', ...
        'maximumIndentationRadius','targetH','minimumMeshQuality', ...
        'maximumAxisAbsZ'}
    validateattributes(opts.(name{1}),{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,['opts.' name{1}]);
end
validateattributes(opts.nBoundary,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',32},mfilename,'opts.nBoundary');
validateattributes(opts.verbose,{'logical','numeric'}, ...
    {'real','finite','scalar'},mfilename,'opts.verbose');
opts.verbose = logical(opts.verbose);
end


function add_project_paths()
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'geometry'),fullfile(root,'fem'), ...
    fullfile(root,'physics'),fullfile(root,'solver'),fullfile(root,'post'));
end
