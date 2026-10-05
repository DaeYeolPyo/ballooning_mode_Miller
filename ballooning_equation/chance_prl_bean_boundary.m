function shape = chance_prl_bean_boundary(indentation, opts)
%CHANCE_PRL_BEAN_BOUNDARY Chance et al., PRL 51, 1963 (1983), Eq. (2).
%   R=xbar+rho*cos(gamma), Z=E*rho*sin(gamma),
%   rho=A*(1+B*cos(t)), gamma=C*sin(t). C fixes d/(2a).
%   E and xbar fix the measured b/a and R/a for every C.
    arguments
        indentation (1,1) double {mustBeFinite,mustBeNonnegative}
        opts.A (1,1) double {mustBeFinite,mustBePositive} = 1
        opts.B (1,1) double {mustBeFinite,mustBePositive} = .6
        opts.AspectRatio (1,1) double {mustBeFinite,mustBeGreaterThan(opts.AspectRatio,1)} = 4
        opts.Kappa (1,1) double {mustBeFinite,mustBePositive} = 1.386
        opts.NBoundary (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(opts.NBoundary,16)} = 128
    end
    if opts.B>=1 || indentation>=.5
        error('chance:Bean:Parameters','Require 0<B<1 and 0<=d/(2a)<0.5.');
    end
    A=opts.A; B=opts.B;
    lower=sqrt(B/(1-B)); upper=pi-1.e-5;
    if lower>=upper
        error('chance:Bean:Branch','The simple bean branch does not include these parameters.');
    end
    if indentation==0
        C=lower;
    else
        maximum=depth(upper);
        if indentation>=maximum
            error('chance:Bean:Depth','Requested indentation is beyond the simple C<pi branch.');
        end
        C=fzero(@(c)depth(c)-indentation,[lower,upper],optimset('TolX',1.e-12));
    end
    radial=@(t)A*(1+B*cos(t)).*cos(C*sin(t));
    vertical=@(t)A*(1+B*cos(t)).*sin(C*sin(t));
    [rmin,rmax]=extrema(radial); [~,zmax]=extrema(vertical);
    a=(rmax-rmin)/2; midpoint=(rmax+rmin)/2;
    E=opts.Kappa*a/zmax; xbar=opts.AspectRatio*a-midpoint;
    t=(0:opts.NBoundary-1).'*2*pi/opts.NBoundary;
    rho=A*(1+B*cos(t)); rt=-A*B*sin(t); rtt=-A*B*cos(t);
    gamma=C*sin(t); gt=C*cos(t); gtt=-C*sin(t);
    cg=cos(gamma); sg=sin(gamma);
    R=xbar+rho.*cg; Z=E*rho.*sg;
    Rt=rt.*cg-rho.*sg.*gt;
    Zt=E*(rt.*sg+rho.*cg.*gt);
    Rtt=rtt.*cg-2*rt.*sg.*gt-rho.*cg.*gt.^2-rho.*sg.*gtt;
    Ztt=E*(rtt.*sg+2*rt.*cg.*gt-rho.*sg.*gt.^2+rho.*cg.*gtt);
    area=sum(R.*circshift(Z,-1)-circshift(R,-1).*Z)/2;
    if min(R)<=0 || min(hypot(Rt,Zt))<1.e-12*a || area<=0
        error('chance:Bean:Geometry','Require positive R, a nondegenerate tangent and CCW orientation.');
    end
    measured=(A*(1-B)-rmin)/(2*a);
    shape=struct('boundary',struct('theta',t,'R',R,'Z',Z,'Rt',Rt,'Zt',Zt, ...
        'Rtt',Rtt,'Ztt',Ztt),'aspectRatio',opts.AspectRatio,'kappa',opts.Kappa, ...
        'aMinor',a,'Rgeo',opts.AspectRatio*a,'bMinor',opts.Kappa*a, ...
        'Rmin',xbar+rmin,'Rmax',xbar+rmax,'indent',measured,'iIn',measured, ...
        'iOut',NaN,'dIn',A*(1-B)-rmin,'indentationSide','inboard', ...
        'polygonArea',area,'construction',struct('A',A,'B',B,'C',C,'E',E,'xbar',xbar), ...
        'source','Chance et al. PRL 51, 1963 (1983), Eq. (2)');

    function value=depth(c)
        f=@(t)A*(1+B*cos(t)).*cos(c*sin(t));
        [lo,hi]=extrema(f);
        value=(A*(1-B)-lo)/(hi-lo);
    end
end

function [lo,hi]=extrema(fun)
    n=2048; dt=2*pi/n; t=(0:n-1).'*dt; v=fun(t);
    small=find(v<=circshift(v,1) & v<=circshift(v,-1));
    large=find(v>=circshift(v,1) & v>=circshift(v,-1));
    lo=min(v); hi=max(v); opts=optimset('TolX',1.e-12,'Display','off');
    for k=small(:).'
        [~,value]=fminbnd(fun,t(k)-dt,t(k)+dt,opts); lo=min(lo,value);
    end
    for k=large(:).'
        [~,value]=fminbnd(@(x)-fun(x),t(k)-dt,t(k)+dt,opts); hi=max(hi,-value);
    end
end
