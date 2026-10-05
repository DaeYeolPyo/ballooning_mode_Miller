function report = test_chance_prl_geometry
    root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
    depths=[.11,.185,.25,.304,.325]; shapes=cell(numel(depths),1);
    for k=1:numel(depths)
        s=chance_prl_bean_boundary(depths(k),NBoundary=8192); shapes{k}=s;
        b=s.boundary;
        a=(max(b.R)-min(b.R))/2; R=(max(b.R)+min(b.R))/2;
        measured=(s.construction.xbar+.4-min(b.R))/(2*a);
        assert(abs(measured-depths(k))<2.e-7);
        assert(abs(R/a-4)<2.e-6 && abs(max(b.Z)/a-1.386)<2.e-6);
        assert(s.polygonArea>0 && min(hypot(b.Rt,b.Zt))>0);
        h=1.e-5; t=b.theta; p=s.construction;
        f=@(x)p.xbar+p.A*(1+p.B*cos(x)).*cos(p.C*sin(x));
        z=@(x)p.E*p.A*(1+p.B*cos(x)).*sin(p.C*sin(x));
        assert(max(abs((f(t+h)-f(t-h))/(2*h)-b.Rt))<1.e-8);
        assert(max(abs((z(t+h)-z(t-h))/(2*h)-b.Zt))<1.e-8);
        fprintf('Bean d/(2a) %.6f: C %.8f, E %.8f, a %.8f\n', ...
            s.iIn,p.C,p.E,s.aMinor);
    end
    cfg=make_chance_prl_figure3_config;
    u=[0;1]; q=polyval(fliplr(cfg.QCoefficients),u);
    qp=polyval(fliplr((1:3).*cfg.QCoefficients(2:4)),u);
    assert(norm(q-[1.03;4.2],inf)<1.e-12 && norm(qp-[.84375;9],inf)<1.e-12);
    report=struct('shapes',{shapes},'source',cfg.PaperReference);
    fprintf('Chance Eq. (2) geometry and prescribed-q endpoint checks passed.\n');
end
