function report = diagnose_chance_q_derivative(equilibrium)
%DIAGNOSE_CHANCE_Q_DERIVATIVE Independent contour derivative vs profile/target.
    setup_indented_cshape_scan;
    [eq,func]=build_ballooning_adapter(equilibrium);
    coefficients=[1.03 .84375 -1.1775 3.50375];
    q=@(u)polyval(fliplr(coefficients),u);
    dq=@(u)polyval(fliplr((1:3).*coefficients(2:4)),u);
    labels=[.85 .90 .925 .94 .945 .95 .955 .96 .97 .975].';
    mesh=equilibrium.mesh;
    cache=build_q_contour_cache(mesh,400);
    native=equilibrium.flux;
    geom=evaluate_q_geometry(mesh,equilibrium.psi,native.axis,native.nativeSpan, ...
        labels,400,equilibrium.axisNode,equilibrium.solverResult.axisMode,cache);
    values=evaluate_q_from_G(geom,equilibrium.profile.psiN, ...
        equilibrium.profile.Gpol,sign(equilibrium.profile.Fpol(end)));
    ua=values.psiN(:); qa=values.q(:);
    a=find(abs(ua-.945)<1.e-10); b=find(abs(ua-.955)<1.e-10);
    contourDerivative=(qa(b)-qa(a))/.01;
    fprintf('qprimeN at .95: adapter %.8g, paper %.8g, contour central difference %.8g\n', ...
        func.qprimeN(.95),dq(.95),contourDerivative);
    prescribed=func;
    prescribed.q=q; prescribed.qprimeN=dq;
    prescribed.qprime=@(u)dq(u)/abs(native.nativeSpan);
    u=[.85,.9,.95]; states=cell(numel(u),2);
    for k=1:numel(u)
        states{k,1}=prepare_salpha_equilibrium(eq,func,u(k),Model='Direct', ...
            NSurfacePoints=1024,FourierModes=48);
        states{k,2}=prepare_salpha_equilibrium(eq,prescribed,u(k),Model='Direct', ...
            NSurfacePoints=1024,FourierModes=48);
        fprintf('u %.3f closure: computed-profile %.8g, prescribed-profile %.8g\n', ...
            u(k),states{k,1}.thetaDiag.input_closure_defect, ...
            states{k,2}.thetaDiag.input_closure_defect);
    end
    report=struct('psiN',ua,'qContour',qa,'qTarget',q(ua), ...
        'contourDerivative',contourDerivative,'adapterDerivative',func.qprimeN(.95), ...
        'targetDerivative',dq(.95),'states',{states});
end
