function cache = build_q_contour_cache(P3, nR)
%BUILD_Q_CONTOUR_CACHE Prepare the fixed q-contour grid once per mesh.
    validateattributes(nR,{'double'},{'scalar','integer','>=',16});
    Rlim=[min(P3.P1points(:,1)),max(P3.P1points(:,1))];
    Zlim=[min(P3.P1points(:,2)),max(P3.P1points(:,2))];
    nZ=max(200,round(nR*diff(Zlim)/diff(Rlim)));
    Rvec=linspace(Rlim(1),Rlim(2),nR);
    Zvec=linspace(Zlim(1),Zlim(2),nZ);
    [Rgrid,Zgrid]=meshgrid(Rvec,Zvec);
    cache=build_P3_value_cache(P3,[Rgrid(:),Zgrid(:)]);
    cache.RGrid=Rvec; cache.ZGrid=Zvec; cache.gridSize=size(Rgrid);
    cache.nR=nR;
end
