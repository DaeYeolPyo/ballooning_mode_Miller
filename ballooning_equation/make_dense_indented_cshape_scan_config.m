function config = make_dense_indented_cshape_scan_config()
%MAKE_DENSE_INDENTED_CSHAPE_SCAN_CONFIG Original outboard C-shape dense map.
%   Central beta is only a pressure control. The plotted ordinate is the
%   actual 2*mu0*integral(p dV)/integral(B^2 dV) of each converged equilibrium.
    config=make_indented_cshape_scan_config(Indentation=.6:.025:1, ...
        CentralBeta=[.002,.01:.01:.12,.14:.02:.30],NBasis=12);
    config.PressureCompensatedWarmStart=true;
    config.PQConstraint.outer.maxIterations=60;
    config.Ballooning.CellsPerTurn=96;
    config.Ballooning.AdapterNR=513;
    config.Ballooning.NSurfacePoints=2048;
    config.StoreEquilibria=true;
end
