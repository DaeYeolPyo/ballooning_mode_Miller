function config = make_indented_cshape_scan_config(opts)
%MAKE_INDENTED_CSHAPE_SCAN_CONFIG Small prescribed-p/q production setup.
%   Defaults preserve the actual dimensions of the existing fixed-q test
%   (R0=3,a0=.8,kappa0=3,delta=-.4,I=.8). CentralBeta controls p0,
%   not the volume beta displayed by the scan.
    arguments
        opts.AspectRatio double = []
        opts.Kappa double = []
        opts.MinorRadius double = []
        opts.Delta (1,1) double = -0.4
        opts.Indentation double = [0.6,0.8,1.0]
        opts.BReference (1,1) double {mustBeFinite,mustBePositive} = 3
        opts.CentralBeta double = [0.02,0.1,0.3]
        opts.PressureAmplitude double = []
        opts.GridSize (1,2) double {mustBeInteger,mustBeGreaterThanOrEqual(opts.GridSize,8)} = [48,64]
        opts.NBoundary (1,1) double = 128
        opts.NProfile (1,1) double = 101
        opts.NBasis (1,1) double = 8
    end
    setup_indented_cshape_scan;
    reference=normalize_indented_cshape_boundary(4,3,1,opts.Delta,0.8);
    a=.8*reference.rawHalfWidth;
    R=3+.8*reference.rawMidpoint;
    config=struct();
    config.AspectRatio=opts.AspectRatio;
    config.Kappa=opts.Kappa;
    config.MinorRadius=opts.MinorRadius;
    if isempty(config.AspectRatio), config.AspectRatio=R/a; end
    if isempty(config.Kappa), config.Kappa=2.4/a; end
    if isempty(config.MinorRadius), config.MinorRadius=a; end
    config.Delta=opts.Delta;
    config.Indentation=opts.Indentation;
    config.BReference=opts.BReference;
    config.PressureAmplitude=opts.PressureAmplitude;
    if isempty(config.PressureAmplitude)
        validateattributes(opts.CentralBeta,{'double'},{'vector','nonempty','finite','nonnegative'});
        config.PressureAmplitude=opts.CentralBeta*opts.BReference^2/(2*(4*pi*1.e-7));
    end
    config.QCoefficients=[1.03,0.84375,-1.1775,3.50375];
    config.Grid=struct('nBoundary',opts.NBoundary,'NR',opts.GridSize(1), ...
        'NZ',opts.GridSize(2),'quadratureDegree',7,'nProfile',opts.NProfile);
    config.Picard=struct('omega',.5,'maxIterations',250, ...
        'updateTolerance',2.e-7,'residualTolerance',5.e-7, ...
        'axisMode','max','symmetryTolerance',1.e-12,'fluxSpanTolerance',1.e-12);
    config.PQConstraint=struct();
    config.PQConstraint.outer=struct('maxIterations',30, ...
        'qRelativeTolerance',2.e-3,'qRmsLogTolerance',1.e-3, ...
        'minimumIterations',2,'minimumQPoints',3,'errorOnFailure',false,'verbose',false);
    config.PQConstraint.FFprime=struct('nBasis',opts.NBasis, ...
        'enforceZeroEdge',false,'regularization',1.e-10,'jacobianRidge',1.e-12, ...
        'initialDamping',1.e-3,'maxCoefficientStep',5, ...
        'initialScales',[1,.5,0],'fitRange',[0,.95],'axisWeight',1);
    config.QEvaluation=struct('interiorPsiN',[.05:.05:.90,.95].','contourNR',200);
    config.Continuation=struct('Enabled',true, ...
        'BootstrapPressureAmplitude',.002*opts.BReference^2/(2*(4*pi*1.e-7)), ...
        'MinimumAbsoluteStep',10,'MinimumRelativeStep',1.e-3, ...
        'StepGrowth',1.5,'MaxAttempts',24,'MaximumQIterationRestarts',1);
    config.Ballooning=struct('RadialLabels',[.05:.1:.85,.90,.95], ...
        'Theta0',2*pi*((0:7)-4)/8,'ThetaBoundary',5*pi,'CellsPerTurn',24, ...
        'NSurfacePoints',1024,'FourierModes',48,'ProjectionIterations',2, ...
        'AdapterNR',257,'EnforceDirectQ',true,'EnforcePeriodicity',true, ...
        'LambdaTolerance',1.e-6,'MaximumFluxRms',2.e-3, ...
        'MaximumFluxError',1.e-2,'MinimumNormalAlignment',.95, ...
        'MaximumQCorrection',.05,'MaximumInputClosureDefect',.05, ...
        'StoreWorstEigenfunction',true,'Verbose',false);
    config.BetaQuadratureDegree=10;
    config.StoreEquilibria=false;
    config.CheckpointFile='';
    config.Resume=false;
    config.RetryUnknown=false;
    config.Verbose=true;
end
