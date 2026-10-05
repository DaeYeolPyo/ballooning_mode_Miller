function result = compute_volume_beta(equilibrium, opts)
%COMPUTE_VOLUME_BETA 2*mu0*integral(p dV)/integral(|B|^2 dV).
%   Integrates the native P3 field directly. Squared gradients make this
%   independent of the native flux orientation. No Cartesian cache is used.
    arguments
        equilibrium (1,1) struct
        opts.QuadratureDegree (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(opts.QuadratureDegree,1),mustBeLessThanOrEqual(opts.QuadratureDegree,20)} = 10
    end
    setup_indented_cshape_scan;
    mesh = equilibrium.mesh;
    profile = equilibrium.profile;
    psi = equilibrium.psi(:);
    axis = equilibrium.global.simag;
    span = equilibrium.global.sibry-axis;
    if ~isfinite(span) || span==0 || numel(psi)~=size(mesh.points,1)
        error('cshape:Beta:Flux','Invalid dimensional P3 flux or flux span.');
    end
    u = profile.psiN(:);
    p = profile.pressure(:);
    F = profile.Fpol(:);
    if numel(u)~=numel(p) || numel(u)~=numel(F) || ...
            any(~isfinite([u;p;F;psi])) || any(diff(u)<=0) || u(1)>0 || u(end)<1
        error('cshape:Beta:Profiles','Finite pressure and F profiles must cover [0,1].');
    end
    if min(p)<-1.e-10*max(1,max(abs(p)))
        error('cshape:Beta:Pressure','Negative pressure.');
    end
    [xi,eta,w] = Gauss_Dunavant_quadrature(opts.QuadratureDegree);
    xi=xi(:); eta=eta(:); w=w(:);
    if abs(sum(w)-0.5)>1.e-12
        error('cshape:Beta:Weights','Expected reference triangle weights summing to 1/2.');
    end
    nq=numel(w); N=zeros(nq,10); Nx=N; Ny=N;
    for k=1:nq
        [N(k,:),dN]=construct_P3_shape_functions(xi(k),eta(k));
        Nx(k,:)=dN(1,:); Ny(k,:)=dN(2,:);
    end
    volume=0; pIntegral=0; bpIntegral=0; btIntegral=0;
    maximumOvershoot=0;
    for e=1:size(mesh.elements,1)
        V=mesh.P1points(mesh.P1elements(e,:),:);
        J=[V(2,:)-V(1,:); V(3,:)-V(1,:)].';
        detJ=det(J);
        if detJ<=0
            error('cshape:Beta:Orientation','Inverted element %d.',e);
        end
        R=(1-xi-eta)*V(1,1)+xi*V(2,1)+eta*V(3,1);
        if any(R<=0)
            error('cshape:Beta:Radius','Quadrature point R <= 0.');
        end
        pe=psi(mesh.elements(e,:));
        uq=(N*pe-axis)/span;
        overshoot=max([0;-uq;uq-1]);
        maximumOvershoot=max(maximumOvershoot,overshoot);
        if overshoot>1.e-3
            error('cshape:Beta:FluxRange','psiN overshoot %.3g in element %d.',overshoot,e);
        end
        uq=min(max(uq,0),1);
        pressure=interp1(u,p,uq,'linear');
        fpol=interp1(u,F,uq,'pchip');
        dR=(J(2,2)*Nx-J(2,1)*Ny)/detJ;
        dZ=(-J(1,2)*Nx+J(1,1)*Ny)/detJ;
        bp2=((dR*pe).^2+(dZ*pe).^2)./R.^2;
        bt2=fpol.^2./R.^2;
        dV=2*pi*R.*w*detJ;
        volume=volume+sum(dV);
        pIntegral=pIntegral+sum(pressure.*dV);
        bpIntegral=bpIntegral+sum(bp2.*dV);
        btIntegral=btIntegral+sum(bt2.*dV);
    end
    bIntegral=bpIntegral+btIntegral;
    if volume<=0 || bIntegral<=0 || any(~isfinite([volume,pIntegral,bIntegral]))
        error('cshape:Beta:Integral','Invalid volume or total magnetic-field integral.');
    end
    beta=2*(4*pi*1.e-7)*pIntegral/bIntegral;
    result=struct('beta',beta,'betaPercent',100*beta,'volume',volume, ...
        'pressureIntegral',pIntegral,'Bp2VolumeIntegral',bpIntegral, ...
        'Bphi2VolumeIntegral',btIntegral,'B2VolumeIntegral',bIntegral, ...
        'pressureVolumeAverage',pIntegral/volume, ...
        'quadratureDegree',opts.QuadratureDegree,'maximumPsiNOvershoot',maximumOvershoot);
end
