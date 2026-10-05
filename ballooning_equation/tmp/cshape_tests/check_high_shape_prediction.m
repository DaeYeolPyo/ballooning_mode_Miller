root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_indented_cshape_scan;
saved=load(fullfile(root,'tmp','cshape_tests','full_shape_10.mat'),'scan');
config=make_indented_cshape_scan_config(Indentation=.8);
shape=normalize_indented_cshape_boundary(config.AspectRatio,config.Kappa, ...
    config.MinorRadius,config.Delta,.8,NBoundary=config.Grid.nBoundary);
p0=config.PressureAmplitude(3);
input=build_indented_cshape_scan_input(config,shape,p0);
prepared=prepare_GS_problem(input);
% Diagnostic FF' predictor from the adjacent, fully converged shape.
% No other-shape flux field is injected into the target FEM mesh.
warm=saved.scan.equilibria{3};
[equilibrium,diagnostics,predictor]=solve_indented_cshape_equilibrium( ...
    config,shape,p0,prepared,warm);
save(fullfile(root,'tmp','cshape_tests','high_shape_prediction_08.mat'), ...
    'equilibrium','diagnostics','predictor','config');
fprintf('High-pressure adjacent-shape predictor succeeded: %d\n',~isempty(equilibrium));
if ~isempty(equilibrium)
    fprintf('GS residual %.6e, max q %.6e, RMS q %.6e\n', ...
        equilibrium.solverResult.finalRelativeResidual, ...
        equilibrium.constraintResult.finalMaxAbsRelativeQ, ...
        equilibrium.constraintResult.finalRmsLogQ);
end
