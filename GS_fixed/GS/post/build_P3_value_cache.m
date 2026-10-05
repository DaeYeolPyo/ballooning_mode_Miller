function cache = build_P3_value_cache(P3, Xq)
%BUILD_P3_VALUE_CACHE Fixed-mesh interpolation at fixed Cartesian points.
%   Reuses the canonical scalar P3 basis; no field derivative or numerical
%   formulation changes. Outside points stay NaN when applying the cache.
    validateattributes(Xq,{'double'},{'real','finite','2d','ncols',2});
    P1=triangulation(P3.P1elements,P3.P1points);
    [elementID,barycentric]=pointLocation(P1,Xq);
    inside=find(~isnan(elementID));
    weights=zeros(numel(inside),10);
    for k=1:numel(inside)
        lambda=barycentric(inside(k),:);
        weights(k,:)=construct_P3_shape_functions(lambda(2),lambda(3));
    end
    columns=P3.elements(elementID(inside),:);
    rows=repmat(inside,1,10);
    operator=sparse(rows(:),columns(:),weights(:),size(Xq,1),size(P3.points,1));
    cache=struct('operator',operator,'inside',~isnan(elementID), ...
        'P1points',P3.P1points,'P1elements',P3.P1elements,'elements',P3.elements);
end
