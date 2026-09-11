function mesh = generate_P1_mesh(vertices, boundary, plot_mesh)
    arguments
        vertices (:, 2) double
        boundary (:, 2) double
        plot_mesh (1, 1) logical = false
    end

    % Constrained Delaunay triangulation covers the convex hull.
    DT = delaunayTriangulation(vertices, boundary);

    % Retain only triangles inside the closed constrained domain.
    interiorTriangle = isInterior(DT);

    if ~any(interiorTriangle)
        error('GENERATE_P1_MESH:EmptyDomain', ...
            'No triangles were found inside the prescribed boundary.');
    end

    connectivity = DT.ConnectivityList(interiorTriangle, :);

    % A general triangulation is required because a subset of the
    % Delaunay connectivity is being retained.
    mesh = triangulation(connectivity, DT.Points);

    if plot_mesh
        triplot(mesh);
    end
end