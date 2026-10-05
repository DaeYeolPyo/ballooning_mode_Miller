function config = make_chance_prl_figure3_config(opts)
%MAKE_CHANCE_PRL_FIGURE3_CONFIG Paper geometry/profiles; explicit test resolution.
    arguments
        opts.Indentation double {mustBeFinite,mustBeNonnegative} = [.11,.185,.25,.304,.325]
        opts.TargetBetaPercent double {mustBeFinite,mustBePositive} = [.8,5,15,24]
        opts.GridSize (1,2) double = [48,64]
        opts.NBasis (1,1) double = 12
        opts.ContourNR (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(opts.ContourNR,200)} = 800
    end
    config=make_indented_cshape_scan_config(AspectRatio=4,Kappa=1.386, ...
        MinorRadius=1,Indentation=opts.Indentation,GridSize=opts.GridSize, ...
        NBoundary=128,NProfile=201,NBasis=opts.NBasis);
    config.GeometryModel='Chance1983_Eq2';
    config.BenchmarkSchemaVersion=2;
    config.BeanA=1; config.BeanB=.6;
    config.TargetBetaPercent=opts.TargetBetaPercent;
    config.BetaRelativeTolerance=.015;
    config.MaximumBetaIterations=7;
    config.PQConstraint.outer.verbose=false;
    config.Continuation.MaxAttempts=14;
    % Stability at .95 needs q on both sides of that surface. Fitting only
    % through .95 leaves its radial derivative unconstrained by outer data.
    config.QEvaluation.interiorPsiN=[.05:.05:.9 .925 .95 .9625 .975].';
    config.QEvaluation.contourNR=opts.ContourNR;
    config.PQConstraint.FFprime.fitRange=[0 .975];
    config.Ballooning.MaximumInputClosureDefect=.01;
    config.Ballooning.AdapterNR=opts.ContourNR;
    config.Ballooning.CellsPerTurn=96;
    config.Ballooning.NSurfacePoints=2048;
    config.PQConstraint.FFprime.regularization=1.e-10;
    config.PQConstraint.FFprime.finiteDifferenceRelativeStep=.002;
    config.PQConstraint.FFprime.finiteDifferenceAbsoluteStep=.005;
    config.PressureCompensatedWarmStart=true;
    config.PaperReference=struct('DOI','10.1103/PhysRevLett.51.1963', ...
        'figure',3,'boundaryEquation',2,'aspectRatio',4,'kappa',1.386, ...
        'A',1,'B',.6,'qAxis',1.03,'qEdge',4.2,'qSlopeAxis',.84375, ...
        'qSlopeEdge',9,'pressure','p0*(1-psiN^2)^2');
end
