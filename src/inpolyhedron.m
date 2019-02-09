function in = inpolyhedron(polytope, points)
% inpolyhedron - Perform point in polyhedron test and return a boolean array
%
% in = inpolyhedron(polytope, vertices)
%
% polytope = a set of vertices of a polytope, each ROW of which is one vertex
% points = points to be tested (nx3 array). each row is a single point
% in = a boolean vector of length n, indicating inside or outside result

[A, b] = vert2con(polytope);
C = points*A' - repmat(b', size(points, 1), 1); % Ax - b
in = all(C < 0, 2);
