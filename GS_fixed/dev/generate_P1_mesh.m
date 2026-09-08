function mesh = generate_P1_mesh(vertices, boundary, plot_mesh)
    arguments
        % Triangular vertices
        vertices (:, 2) double
        % Plasma boundary
        boundary (:, 2) double
        % Show the mesh
        plot_mesh (1, 1) logical = false
    end

    mesh = delaunayTriangulation(vertices, boundary);

    if plot_mesh
        triplot(mesh);
    end
end