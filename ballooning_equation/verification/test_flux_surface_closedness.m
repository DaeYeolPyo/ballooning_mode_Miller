function test_flux_surface_closedness
%TEST_FLUX_SURFACE_CLOSEDNESS A masked open arc is not a flux surface.
    setup_indented_cshape_scan;
    eq=struct('rgrid',linspace(2,4,201),'zgrid',linspace(-1,1,201), ...
        'simag',0,'sibry',1,'rmaxis',3,'zmaxis',0);
    [R,Z]=ndgrid(eq.rgrid,eq.zgrid); eq.psirz=(R-3).^2+Z.^2;
    func=struct('F',@(u)1+0*u,'q',@(u)1+0*u, ...
        'pprime',@(u)0*u,'qprime',@(u)0*u);
    surf=extract_flux_surface(eq,func,.25);
    assert(abs(polyarea(surf.R,surf.Z)-pi/4)<5.e-4);
    assert(max(abs((surf.R-3).^2+surf.Z.^2-.25))<1.e-4);
    % The remaining open arc surrounds the axis if artificially connected
    % by a chord; the old unconditional closure wrongly accepted it.
    eq.psirz(R<2.7)=NaN;
    rejected=false;
    try
        extract_flux_surface(eq,func,.25);
    catch problem
        rejected=strcmp(problem.identifier,'equilibrium:FluxSurface:NoClosedContour');
    end
    assert(rejected,'An open contour must not enter local geometry.');
    fprintf('Closed flux contour accepted; masked open contour rejected.\n');
end
