% Main routine for 3D reconstruction from unorganized planar cross-sections
% 
% Author: Ojaswa Sharma, <ojaswa@iiitd.ac.in>
%
%% 3D reconstruction from arbitrary cutting planes
close all;
clear;
clearvars;
clc;

use_hermite_interpolation = true; % Normalize barycentric coordinates and convert to higher order variant
save_mesh = false;
verbose = false;

%% Step 1: Import model and cross sections
% Select a model
 g_model_name = 'torus';
% g_model_name = 'dragon';
% g_model_name = 'hand';
% g_model_name = 'rook';
% g_model_name = 'femur';
% g_model_name = 'duck';

% Call for loading the model
[g_model, l_section_names] = load_3D(g_model_name);

%% Expand the bounding box
centre = g_model.bbox_origin;
scaling_factor = 1.1;
tansf_1 = [1 0 0 centre(1); 0 1 0 centre(2); 0 0 1 centre(3); 0 0 0 1];
tansf_2 = eye(4);
tansf_2(1:3, 1:3) = scaling_factor * eye(3);
tansf_3 = [1 0 0 -centre(1); 0 1 0 -centre(2); 0 0 1 -centre(3); 0 0 0 1];

l_trans = tansf_1 * tansf_2 * tansf_3;
l_tmp = [g_model.bbox, [1; 1]]* l_trans';
g_model.bbox = l_tmp(:,1:3);

% Get the cross sections
for i=1:size(l_section_names,1)
    l_section = read_wobj(l_section_names{i});
    g_sections{i}.vertices = l_section.vertices;
    l_lines = [l_section.objects.data];
    l_lines = [l_lines.vertices];
    g_sections{i}.raw_lines = [l_lines(1:2:(end-1))' l_lines(2:2:end)'];
    
    % Find plane equation
    [~, l_idx1] = max(dot(l_section.vertices, l_section.vertices, 2));
    l_verts_temp = l_section.vertices;
    l_verts_temp(l_idx1,:) = [];
    [~, l_idx2] = max(dot(l_verts_temp, l_verts_temp, 2));
    
    if(l_idx1 == l_idx2)
        l_idx2 = l_idx2 + 1;
    end
    
    l_a = l_section.vertices(1,:);
    l_b = l_section.vertices(l_idx1,:);
    l_c = l_section.vertices(l_idx2,:);
    ab = [l_b(1)-l_a(1) l_b(2)-l_a(2) l_b(3)-l_a(3)];
    ac = [l_c(1)-l_a(1) l_c(2)-l_a(2) l_c(3)-l_a(3)];
    pts = cross(ab,ac);
    d = -(pts(1)*l_a(1) + pts(2)*l_a(2) + pts(3)*l_a(3));
    g_sections{i}.plane = [pts(1) pts(2) pts(3) d];
end

%% Find Plane coordiantes (A,B,C,D)
for i = 1 : size(g_sections,2)
    plane_param = g_sections{i}.plane;
    if(any(plane_param))
        new_coords = plane_param/norm(plane_param(1:3));       % [A B C]
        g_sections{i}.plane = [new_coords(1:3) -new_coords(4)];
    else
        g_sections{i}.plane = plane_param;
    end
end

% Preprocess cross sections - store disconnected chain of indices into vertex list
for i=1:size(g_sections,2)
    l_section = g_sections{i};
    % Compute adjacency matrix
    l_adj = sparse(size(l_section.vertices,1), size(l_section.vertices,1));
    for j=1:size(l_section.raw_lines,1)
        l_adj(l_section.raw_lines(j,1), l_section.raw_lines(j,2)) = 1;
        l_adj(l_section.raw_lines(j,2), l_section.raw_lines(j,1)) = 1;
    end
    
    % Convert into a linear chain.
    l_chain_id = 1;
    l_chains = {};
    while(nnz(l_adj))
        [l_first, l_tmp] = find(l_adj,1);
        if(l_first == [])
            break;
        end
        l_idx0 = l_first;
        l_idx1 = 0;
        l_chain = [l_idx0];
        while(true)
            l_idx1 = find(l_adj(l_idx0,:), 1 );
            l_adj(l_idx0, l_idx1) = 0;
            l_adj(l_idx1, l_idx0) = 0;
            if(l_idx1 == l_first)
                break; % This completes the chain.
            end
            l_chain = [l_chain; l_idx1];
            l_idx0 = l_idx1;
        end
        l_chains{l_chain_id} = l_chain;
        l_chain_id  = l_chain_id + 1;
    end
    
    if(verbose), fprintf('Section %d: %d chains detected.\n', i, size(l_chains, 2)); end
    g_sections{i}.chains = l_chains;        
end
clear l_section l_chain l_chains l_chain_id l_idx0 l_idx1 l_adj l_first l_tmp;

%% Step 1.1 Process chains / build topology
tic
g_TransformationMatrix = {};
for i=1:size(g_sections,2)
    l_chains = g_sections{i}.chains;
    l_plane = g_sections{i}.plane;
    l_vertices = g_sections{i}.vertices;
    
    l_origin = l_vertices(1, :); % Set the first point on polytope face to origin
    l_nodes = l_vertices - repmat(l_origin, size(l_vertices, 1), 1); % Translate
    
    % Transform vertices to XY plane - compute the matrix
    l_w = l_plane(1:3);
    [~, l_idx] = max(dot(l_nodes, l_nodes, 2));
    l_bbox_center = (g_model.bbox(1,:) + g_model.bbox(2,:))/2;
    l_u = l_nodes(l_idx, :);
    
    l_u = l_u / sqrt(dot(l_u, l_u));
    l_v = cross(l_w, l_u);
    l_M = [l_u; l_v; l_w];
    g_TransformationMatrix{i} = l_M;
    
    l_section.vertices = [];
    l_section.chains = {};
    for j=1:size(l_chains, 2)
        l_chain = l_chains{j};
        l_poly = l_vertices(l_chain, :) - repmat(l_origin, size(l_chain, 1), 1);
        l_poly = l_poly * l_M';                 % Third column will now be zero!
        
        % Rectify chain orientation if not already counterclockwise
        if polygonArea(l_poly(:,1:2)) < 0
            l_poly = flipud(l_poly);
            if(verbose), fprintf('Reversed section %d, chain %d to CCW orientation.\n', i, j); end
        end
        
        % Simplify polygon
        l_poly_simplified = simplifyPolygon(l_poly(:,1:2), g_model.bbox_diag/2000);
        l_poly_simplified(:,3) = zeros(size(l_poly_simplified, 1), 1);
        l_idx0 = size(l_section.vertices,1) + 1;
        l_section.vertices = [l_section.vertices; l_poly_simplified * l_M + repmat(l_origin, size(l_poly_simplified, 1), 1)];
        l_idx1 = size(l_section.vertices,1);
        l_section.chains{j} = [l_idx0:l_idx1]';
        if(verbose), fprintf('Simplified section %d, chain %d (%d -> %d verts)\n', i, j, size(l_poly, 1), size(l_poly_simplified, 1)); end
    end
    
    % Build topology (classify contours as holes or contours
    l_verts_2d = l_section.vertices - repmat(l_origin, size(l_section.vertices, 1), 1);
    l_verts_2d = l_verts_2d * l_M';
    [l_chain_parents, l_contours, l_holes] = classifyChains(l_verts_2d, l_section.chains);
    
    % Store results back to g_sections
    l_section.plane = l_plane;
    l_section.chain_parents = l_chain_parents;
    l_section.contours = l_contours;
    l_section.holes = l_holes;
    % Holes must be clockwise oriented.
    if ~isempty(l_section.holes)
        for j=1:size(l_section.holes, 1)
            l_holeID = l_section.holes(j, 1);
            l_section.chains{l_holeID} = flipud(l_section.chains{l_holeID});
        end
    end
    g_sections{i} = l_section;
end

l_timeElapsed = toc;
fprintf('Chain processing took: %f sec.\n', l_timeElapsed);

clear l_chains l_chain l_vertices l_max l_axis
clear l_section l_idx0 l_idx1 l_poly l_poly_simplified
clear l_u l_v l_w l_bbox l_bbox_center l_M l_offset
clear l_holeID l_holes l_contours l_chain_parents l_verts_2d

%% Find average chain length for all the cross-sections
sum_l = 0;
num = 0;
for i=1:size(g_sections, 2)
    chains = g_sections{i}.chains;
    verts = g_sections{i}.vertices;
    if isempty(chains)
        continue
    end
    for j=1:size(chains, 2)
        vecs = verts(chains{j}, :) - circshift(verts(chains{j}, :), -1, 1);
        lengths = sqrt(dot(vecs, vecs, 2));
        sum_l = sum_l + sum(lengths, 1);
        num = num + size(vecs, 1);
    end
end

avg_length = sum_l/num;

%% Step 2(a): Partition model bounding box into polytopes
tic
%Create bounding box
g_bbox = createCube;
l_bbox = g_model.bbox;
g_bbox.vertices = horzcat(repmat(l_bbox(:,1), 4, 1), ...
    repmat([l_bbox(1,2); l_bbox(1,2); l_bbox(2,2); l_bbox(2,2)], 2, 1), ...
    [repmat(l_bbox(1, 3), 4, 1); repmat(l_bbox(2, 3), 4, 1)]);

g_partitions = {};
l_polytope.vertices = g_bbox.vertices;
l_polytope.faces = g_bbox.faces;
l_polytope.sections = cell(6, 1);
g_partitions{1} = l_polytope;
l_partitions = {};

for i = 1 : size(g_sections,2)        % For all cutting planes, partition
    l_ABCD = g_sections{i}.plane;
    l_npartitions = size(g_partitions,2);
    k = 1;
    for j = 1 : l_npartitions
        l_polytope = g_partitions{j};
        
        l_plane = createPlane(l_ABCD(4)*l_ABCD(1:3), l_ABCD(1:3));
        [l_n1, l_f1] = clipConvexPolyhedronHPNew(l_polytope.vertices, l_polytope.faces, l_plane);
        
        % Polytope 1
        if(~isempty(l_n1))
            l_partitions{k}.vertices = l_n1;
            l_partitions{k}.faces = l_f1;
            k = k + 1;
        end
        
        % Polytope 2
        l_plane = createPlane(l_ABCD(4)*l_ABCD(1:3), -l_ABCD(1:3)); % Reverse the normal
        [l_n2, l_f2] = clipConvexPolyhedronHPNew(l_polytope.vertices, l_polytope.faces, l_plane);
        if(~isempty(l_n2))
            l_partitions{k}.vertices = l_n2;
            l_partitions{k}.faces = l_f2;
            k = k + 1;
        end
    end
    % Remove duplicate nodes, faces and null polytopes
    l_partitions = remove_duplicates(l_partitions);
    
    g_partitions = l_partitions;
end
l_timeElapsed = toc;
fprintf('Polytope space partitioning took: %f sec.\n', l_timeElapsed);

%% Correctly reoient polytope faces  (face normal pointing out)
for i=1: size(g_partitions, 2)
    l_polytope = g_partitions{i};
    l_normals = faceNormal(l_polytope.vertices, l_polytope.faces);
    l_polytope_centre = sum(l_polytope.vertices, 1)/size(l_polytope.vertices, 1);
    for j=1:size(l_polytope.faces, 2)
        l_face_vertices  = l_polytope.vertices(l_polytope.faces{j}, :);
        l_face_center = sum(l_face_vertices, 1)/size(l_face_vertices, 1);
        l_out_dir = l_face_center - l_polytope_centre;
        if dot(l_normals(j,:), l_out_dir) < 0 % reorient face
            l_polytope.faces{j} = fliplr(l_polytope.faces{j});
            if(verbose), fprintf('Reoriented face %d on polytope %d\n', j, i); end
        end
    end
    l_polytope.normals = faceNormal(l_polytope.vertices, l_polytope.faces); % Store updated normals
    g_partitions{i} = l_polytope;
end

clear l_polytope l_normals l_face_center;
clear i j l_polytope_centre l_face_vertices l_out_dir;

%% Step 2(b): Assign global cutting plane IDs to each polytope face
for i = 1 : size(g_partitions,2)                % For all Polytopes
    l_polytope = g_partitions{i};
    l_nodes = l_polytope.vertices;
    
    l_planeId = zeros(size(l_polytope.faces,2),1);
    for j = 1 : size(l_polytope.faces,2)        % For each face of the polytope
        l_face = l_polytope.faces{j};
        
        g_eqn = [];
        for k = 1 : size(g_sections,2)          % For all Cutting Planes
            plane = g_sections{k}.plane;
            
            % Check, face cordinates satisfies which plane equation
            l_eqn = [];
            for l = 1 : size(l_face, 2)         % For each vertex of face
                point = l_nodes(l_face(l),:);
                x = point(1);
                y = point(2);
                z = point(3);
                l_eqn(l) = plane(1)*x + plane(2)*y + plane(3)*z - plane(4);
            end
            g_eqn(k,:) = l_eqn(:);      % 3xn - 3 cutting planes & n - no. of vertices in a polytope face
        end
        
        % Rounding off to 1 decimal place
        eqn = roundToDP(g_eqn,6);
        
        % Find number of 0s obtined(after putting each vertex to plane eq) for the current face
        no_of_zeros = [];
        for m = 1 : size(g_sections,2)
            no_of_zeros(m) =  sum(eqn(m,:) == 0);
        end
        
        val = max(no_of_zeros);
        id = find(no_of_zeros == val);
        
        % Assign plane id to the Current face
        if(val < size(l_face,2))
            l_planeId(j) = 0;
        else
            if(size(id) == 1)
                l_planeId(j) = id;
            end
        end
    end
    g_partitions{i}.planeId = l_planeId;
end

clear l_polytope l_nodes l_planeId l_face l_eqn val id;
clear g_eqn plane point no_of_zeros i j k l m x y z;

%% Step 2(c):Add Cross sections on each face of Polytopes generated
tic
flag = 0;                       % For Holes Orientation
for i = 1 : size(g_partitions,2)                % For all Polytopes
    l_polytope = g_partitions{i};
    l_polytope.sections = {};
    for j = 1 : size(l_polytope.faces,2)        % For each face of the polytope
        l_polytope.identification{j} = [];
        l_face = l_polytope.faces{j};
        face_vertices3d = l_polytope.vertices(l_face(:),:);
        plane_id = l_polytope.planeId(j);
        
        if(plane_id == 0)
            l_polytope.sections{j} = {};
            l_polytope.identification{j} = [];
        else
            % Clipping Window
            face_vertices2d = face_vertices3d * g_TransformationMatrix{plane_id}';
            face_vertices2d = face_vertices2d(:,1:2);
            if polygonArea(face_vertices2d) < 0
                face_vertices2d = flipud(face_vertices2d);
            end
            
            % IDs of contour and holes
            contour_id = g_sections{plane_id}.contours;
            hole_id = g_sections{plane_id}.holes;
            
            l_polygon = g_sections{plane_id};
            for k = 1 : size(g_sections{plane_id}.chains, 2)                % For each chain
                chain = l_polygon.chains{k};
                polygon_vertices3d = l_polygon.vertices(chain(:),:);
                
                % Target Polygon, that is to be clipped against Clipping window
                polygon_vertices2d = polygon_vertices3d * g_TransformationMatrix{plane_id}';
                z_val = polygon_vertices2d(1,3);
                polygon_vertices2d = polygon_vertices2d(:,1:2);
                
                % Check if this polygon is contour or hole
                if(ismember(k, contour_id))
                    l_polytope.identification{j} = [l_polytope.identification{j}; 1];
                elseif(ismember(k, hole_id))
                    l_polytope.identification{j} = [l_polytope.identification{j}; 0];
                    polygon_vertices2d = flipud(polygon_vertices2d);        % Required in counter-clockwise for clipping
                    flag = 1;
                end
                
                % Check if polygon is completely inside clipping window
                in = inpolygon(polygon_vertices2d(:,1), polygon_vertices2d(:,2), face_vertices2d(:,1), face_vertices2d(:,2));
                if(all(in))
                    l_polytope.sections{j}{k} = polygon_vertices3d;
                else
                    u_face_vertices2d = roundToDP(face_vertices2d,10);
                    u_face_vertices2d = unique(u_face_vertices2d,'rows','stable');
                    
                    % Find clipped Polygon
                    clippedPolygon2d = [];
                    clippedPolygon2d = sutherlandHodgman(polygon_vertices2d, u_face_vertices2d);
                    if(isempty(clippedPolygon2d))
                        l_polytope.sections{j}{k} = [];
                    else
                        clippedPolygon2d(:,3) = z_val;
                        clippedPolygon3d = clippedPolygon2d * g_TransformationMatrix{plane_id};
                        if(flag == 1)
                            clippedPolygon3d = flipud(clippedPolygon3d);
                            flag = 0;
                        end
                        l_polytope.sections{j}{k} = clippedPolygon3d;
                    end
                end
            end
        end
    end
    
    g_partitions{i}.sections = l_polytope.sections;             % Update g_partitions
    g_partitions{i}.identification = l_polytope.identification; % 1-Contour; 0-Hole
end
l_timeElapsed = toc;
fprintf('Adding cross-section on each face took: %f sec.\n', l_timeElapsed);

clear l_polytope l_face l_polygon polygon_vertices3d polygon_vertices2d clippedPolygon3d;
clear face_vertices3d plane_id face_vertices face_vertices2d clippedPolygon2d;
clear contour_id hole_id chain flag i j k in x y;

%% Step 2(d): Find out original index of cross section points (w.r.t the cutting plane)
% This is better done at the time of clipping. Register a zero if the point
% does not belong to the original set
tic
for i=1:size(g_partitions, 2)
    l_polytope = g_partitions{i};
    l_sections_idx = cell(1, size(l_polytope.faces,2));
    for j = 1:size(l_polytope.faces,2)
        l_sections = l_polytope.sections{j}; % sectional curves for this face, possibily multiple.
        if isempty(l_sections)
            continue;
        end
        l_original_sections = g_sections{l_polytope.planeId(j)}.vertices;
        
        l_original_sections = roundToDP(l_original_sections,6);
        
        l_face_sections_idx = cell(1, size(l_sections, 2));
        for k=1:size(l_sections, 2)
            l_section = l_sections{k};
            if isempty(l_section)
                continue;
            end
            l_section_idx = zeros(size(l_section, 1), 1);
            l_dist_mat = pdist2(l_original_sections, l_section, 'cityblock'); % For our purpose 'cityblock' is a good metric.
            [l_r, l_c] = find(l_dist_mat < 0.099999); % l_c will contain index of l_section
            
            l_section_idx(l_c) = l_r;
            l_face_sections_idx{k} = l_section_idx;
        end
        l_sections_idx{j} = l_face_sections_idx;
    end
    g_partitions{i}.sections_idx = l_sections_idx;
end
l_timeElapsed = toc;

fprintf('X-section index mapping took: %f sec.\n', l_timeElapsed);

clear l_polytope l_sections_idx l_sections l_original_sections l_face_sections_idx l_section
clear l_section_idx l_timeElapsed i j k l_dist_mat l_r l_c;

%% Step 3: Triangulate cutting planes and BBox faces (with creases and intersection boundary as constraints)
% Calculate minimum area based on edge length heuristic
l_timeElapsed = 0;
l_area = getMaximumTriangleAreaHeuristic(g_sections);
for i = 1 : size(g_partitions,2)            % For each Polytope
    l_polytope = g_partitions{i};
    l_normals = l_polytope.normals;
    for j = 1:size(l_polytope.faces,2) % Triangulate each face
        % Accumulate vertices and boundary for triangulation
        l_face = l_polytope.faces{j};
        l_verts = l_polytope.vertices(l_face, :);
        l_origin = l_verts(1, :); % Set the first point on polytope face to origin
        l_nodes = l_verts - repmat(l_origin, size(l_verts, 1), 1); % Translate
        
        verts = l_nodes;
        
        l_boundary = [1:size(l_nodes, 1); [2:size(l_nodes, 1) 1]]';
        l_markers = zeros(size(l_nodes, 1), 1);
        l_sections = l_polytope.sections{j};
        l_face_sections_idx = l_polytope.sections_idx{j};
        if ~isempty(l_sections)
            for k=1:size(l_sections, 2)
                flag = 0;                       % cross section is polyline or polygon; 0-polyline
                l_section = l_sections{k};
                l_section_idx = l_face_sections_idx{k};
                if ~isempty(l_section)
                    k1 = size(l_nodes, 1)+1;
                    l_nodes = [l_nodes; (l_section - repmat(l_origin, size(l_section, 1), 1))];
                    k2 = size(l_nodes, 1);
                    l_boundary = [l_boundary; [k1:k2; [(k1+1):k2 k1]]'];
                    l_marker = l_section_idx + (l_section_idx > 0); % increment non-zero values by 1. {0, 1} are reserved markers in Triangle
                    l_markers = [l_markers; l_marker];
                end
            end
        end
        
        % Transform 3D vertices to 2D
        l_w = l_normals(j, :);
        % Select a far-away point on the
        [~, l_idx] = max(dot(l_nodes, l_nodes, 2));
        l_u = l_nodes(l_idx, :);
        l_u = l_u / sqrt(dot(l_u, l_u));
        l_v = cross(l_w, l_u);
        l_M = [l_u; l_v; l_w];
        l_nodes = l_nodes * l_M';
        z1 = l_nodes(1,3);
        
        % Compute 2D CDT        
        l_args = num2str(5*l_area);
        [l_CDT_nodes, l_CDT_triangles, l_CDT_node_markers, l_time] = triangle(l_nodes(:, 1:2), l_boundary, [], l_markers, ['-QXa', l_args]);
        l_timeElapsed = l_timeElapsed + l_time;
        % Transform 2D nodes to 3D
        l_CDT_nodes = [l_CDT_nodes', z1*zeros(size(l_CDT_nodes, 2), 1)];
        l_CDT_nodes = l_CDT_nodes * l_M + repmat(l_origin, size(l_CDT_nodes, 1), 1);
        
        % Save 'nodes' & 'triangles' in g_partitions corresponding to each face
        g_partitions{i}.triangulated_nodes{j} = l_CDT_nodes;
        g_partitions{i}.triangles{j} = l_CDT_triangles';
        g_partitions{i}.triangulated_node_markers{j} = l_CDT_node_markers - (l_CDT_node_markers > 0); % Convert back to section node IDs
    end
end


for i = 1 : size(g_partitions,2)
    tris = g_partitions{i}.triangles;
    for j = 1 : size(tris,2)
        l_tri = tris{j};
        if(isempty(l_tri))
            g_partitions{i}.triangles{j} = [1 2 3];
        end
    end
end
fprintf('Face triangulation took: %f sec.\n', l_timeElapsed);

clear l_area l_polytope l_normals l_face l_verts flag_vert eps dist col row l_origin
clear l_nodes l_boundary l_markers l_sections l_face_sections_idx flag l_section l_section_idx
clear k1 k2 l_marker l_w l_u l_v l_tmp l_idx l_M l_args l_CDT_nodes l_CDT_triangles l_CDT_node_markers

%% Step 4: Compute SDT at all vertices on all polytopes
% Compute array of line segmentsg_sections
tic
l_all_lines = [];
for i = 1 : size(g_sections,2)
    l_chains = g_sections{i}.chains;
    l_verts = g_sections{i}.vertices;
    for j = 1 : size(l_chains,2)
        l_chain = g_sections{i}.chains{j};
        l_all_lines = [l_all_lines; [l_verts(l_chain, :) l_verts(circshift(l_chain, [-1, 1]), :)]];
    end
end

for i = 1 : size(g_partitions,2)            % For all polytopes
    l_polytope = g_partitions{i};
    
    for j = 1 : size(l_polytope.faces,2)    % For each Face of Polytope
        l_triangulatedNodes = l_polytope.triangulated_nodes{j};
        l_num_nodes = size(l_triangulatedNodes, 1);
        l_min_dist = zeros(l_num_nodes, 1);
        for k=1:l_num_nodes
            l_min_dist(k, 1) = min(distancePointEdge3d(l_triangulatedNodes(k, :), l_all_lines), [], 1);
        end
        g_partitions{i}.min_dist{j} = l_min_dist;
        
        % SDT Caclculation
        l_face = l_polytope.faces{j};
        l_vertices = l_polytope.vertices(l_face,:);
        l_plane_id = l_polytope.planeId(j);
        l_triangulatedNodes = l_polytope.triangulated_nodes{j};
        
        % Transformation Matrix for conversion to 2D
        transf_matrix = [];
        if(l_plane_id == 0)                 % No cross sections on this face%
            % Find equation of face
            
            % Normal to the face
            points = l_vertices(1:3,:);     % 1st 3 points of face
            v1 = points(2,:)-points(1,:);
            v2 = points(3,:)-points(1,:);
            normal = cross(v1, v2);         % Normal to the face
            d = dot(normal, points(1,:));   % Offset from origin
            
            % Transformation matrix
            l_w = normal;
            l_w = l_w/sqrt(dot(l_w, l_w));
            l_u = l_vertices(1,:);
            l_u = l_u / sqrt(dot(l_u, l_u));
            l_v = cross(l_w, l_u);
            l_M = [l_u; l_v; l_w];
            transf_matrix = l_M;
            sections = [];
            
        else
            transf_matrix = g_TransformationMatrix{l_plane_id};
            
            % Combine all cross sections on the cutting plane separated by NaN
            sections = [];
            % For all chains on the cutting plane
            for k = 1 : size(g_sections{l_plane_id}.chains,2)
                sections3d = [];
                l_chain = g_sections{l_plane_id}.chains{k};
                if(isempty(sections))
                    sections3d = g_sections{l_plane_id}.vertices(l_chain,:);
                    sections3d = sections3d * transf_matrix';
                    sections = sections3d(:,1:2);
                else
                    sections = [sections; NaN NaN];
                    sections3d = g_sections{l_plane_id}.vertices(l_chain,:);
                    sections3d = sections3d * transf_matrix';
                    sections = [sections; sections3d(:,1:2)];
                end
            end
            
        end
        
        l_triangulatedNodes = l_triangulatedNodes * transf_matrix';
        triangulatedNodes = l_triangulatedNodes(:,1:2);     % Points of which sign has to be determined
        
        % Find Sign value
        if(isempty(sections))
            in = zeros(size(triangulatedNodes,1),1);
            on = in;
        else
            [in, on] = inpolygon(triangulatedNodes(:,1), triangulatedNodes(:,2), sections(:,1), sections(:,2));
        end
        l_sign = (in)*2 - ones(size(in, 1), 1);
        g_partitions{i}.sdt{j} = l_sign.*g_partitions{i}.min_dist{j};
    end
end

l_timeElapsed = toc;
fprintf('SDT computation took: %f sec.\n', l_timeElapsed);

clear all_sections l_chain l_chains l_verts l_polytope l_triangulatedNodes l_all_lines;
clear l_sections l_face l_vertices l_plane_id sections3d sections l_num_nodes ;
clear id min_dist transf_matrix points triangulatedNodes sign l_min_dist;
clear id i j k l in v1 v2 normal d l_w l_u l_v l_M;

%% Step 4(b) Compute and store gradients with every face of the polytope
%Also flip the plane if orientated incorrectly
tic
for i=1:size(g_partitions,2)
    l_polytope = g_partitions{i};
    l_gradients = cell(1, size(l_polytope.faces, 2));
    for temp_f = 1 : size(l_polytope.faces, 2)
        if l_polytope.planeId(temp_f) > 0
            l_normal = l_polytope.normals(temp_f, :);
            l_chains = g_sections{l_polytope.planeId(temp_f)}.chains;
            l_plane = g_sections{l_polytope.planeId(temp_f)}.plane;
            if dot(l_normal, l_plane(1:3)) < 0 % Then reverse the order of chains in l_chains
                for j=1: size(l_chains, 2)
                    l_chains{j} = flipud(l_chains{j});
                end
            end
            l_face_gradients = computeGradients(l_polytope.triangulated_node_markers{temp_f}, l_normal, l_chains, g_sections{l_polytope.planeId(temp_f)}.vertices);
        else
            l_face_gradients = zeros(size(l_polytope.sdt{temp_f}, 1), 3);
        end
        l_gradients{temp_f} = l_face_gradients;
    end
    g_partitions{i}.gradients = l_gradients;
end

l_timeElapsed = toc;
fprintf('Gradient computation took: %f sec.\n', l_timeElapsed);

clear l_polytope l_gradients temp_f l_normal l_chains l_plane l_face_gradients epsilon
clear l_normals l_vertices l_triangles l_func l_origin l_vertices2d l_tmp l_idx l_w l_u l_v l_M
clear grad_x grad_y grad_mag l_nan i;

%% Step 5 (a): Interpolate distance values inside polytopes (at regular grid points) using MVC for meshes
% Create a sample grid of points
l_max_dim = max(g_model.bbox(2,:) - g_model.bbox(1,:));
g_max_res = 100;   % Resolution along the maximum dimension
l_grid_delta = l_max_dim / (g_max_res - 1 + 2);

g_x_vec = g_model.bbox(1,1):l_grid_delta:g_model.bbox(2,1);
g_x_vec([1, end]) = [];

g_y_vec = g_model.bbox(1,2):l_grid_delta:g_model.bbox(2,2);
g_y_vec([1, end]) = [];

g_z_vec = g_model.bbox(1,3):l_grid_delta:g_model.bbox(2,3);
g_z_vec([1, end]) = [];

[g_X, g_Y, g_Z] = meshgrid(g_x_vec, g_y_vec, g_z_vec);

% Create label array - every element is a polyhedron ID
g_L = zeros(numel(g_X), 1);
l_all_sample_points = [g_X(:) g_Y(:) g_Z(:)];
tic
for i = 1 : size(g_partitions, 2)
    l_poly_verts  = g_partitions{i}.vertices;
    l_in = inpolyhedron(l_poly_verts, l_all_sample_points);
    g_L = max(g_L, i * l_in);
end
g_L = reshape(g_L, size(g_X, 1), size(g_X, 2), size(g_X, 3));
l_timeElapsed = toc;

fprintf('Inpolyhedron took: %f sec to perform %d tests.\n', l_timeElapsed, size(l_all_sample_points, 1)*size(g_partitions, 2));
clear l_in l_poly_verts l_grid_delta l_max_dim l_timeElapsed

%% Step 5 (b):Perform MVC interpolation of distance values at these sample points
tic
epsilon = 10^-8;
fprintf('Total # of polytopes: %d\n', size(g_partitions, 2));
g_F = zeros(size(g_X));
fprintf('Computing MVC: ');
tic
% For all partitions
for n = 1:size(g_partitions,2)
    fprintf('%d \n',n);
    % Prepare mesh data for MVC computation
    l_polytope = g_partitions{n};
    l_nnodes =  sum(cell2mat(cellfun(@(x) size(x, 1), l_polytope.triangulated_nodes, 'UniformOutput', false)));
    l_ntriangles  = sum(cell2mat(cellfun(@(x) size(x, 1), l_polytope.triangles, 'UniformOutput', false)));
    l_meshVertices = zeros(3, l_nnodes);
    l_meshTriangles = zeros(3, l_ntriangles, 'uint32');
    l_meshSDT = zeros(1, l_nnodes);
    l_meshGradients = zeros(3, l_nnodes);
    l_idp0 = 1; l_idp1 = 1; l_idt0 = 1; l_idt1 = 1;

    for f = 1:size(l_polytope.faces, 2)
        l_idt1 = l_idt0 + size(l_polytope.triangles{f}, 1) - 1;
        l_idp1 = l_idp0 + size(l_polytope.triangulated_nodes{f}, 1) - 1;
        
        l_meshTriangles(:, l_idt0:l_idt1) = l_polytope.triangles{f}' + l_idp0 - 1 - 1;
        
        l_meshVertices(:,l_idp0:l_idp1) = l_polytope.triangulated_nodes{f}';
        l_meshSDT(:,l_idp0:l_idp1) = l_polytope.sdt{f}';
        l_meshGradients(:,l_idp0:l_idp1) = l_polytope.gradients{f}';
        
        l_idp0 = l_idp1 + 1;
        l_idt0 = l_idt1 + 1;
    end
    
    % Find out sampled voxels inside the polytope
    l_voxels = find(g_L == n);
    
    % Compute weights and interpolate
    l_weights = zeros(1, size(l_meshVertices, 2));
    l_uvec = zeros(3, size(l_meshVertices, 2));
    l_dist = zeros(1, size(l_meshVertices, 2));
    for v = 1:size(l_voxels, 1)
        l_idx = l_voxels(v);
        l_x = [g_X(l_idx); g_Y(l_idx); g_Z(l_idx)];
        [l_weights, l_uvec, l_dist] = mvc(l_x, l_meshVertices, l_meshTriangles);
        
        if use_hermite_interpolation
            l_howeights = higherOrderBary(l_weights);
            l_howeights = l_howeights / sum(l_howeights); % Normalize
            l_func = dot(l_howeights, l_meshSDT + dot(l_meshGradients, l_uvec).*l_dist);
        else
            l_func = dot(l_weights, l_meshSDT);
        end
        g_F(l_idx) = l_func; 
    end
end
l_timeElapsed = toc;
fprintf('MVC interpolation (CPU) took: %f sec.\n', l_timeElapsed);

%% Step 6: Extract reconstructed object using Marching cubes
figure; 
isosurface(g_X,g_Y,g_Z,g_F,0);
axis equal

%% Extract iso-surface and save to obj file
% Compute sosurface 
if save_mesh
    FV = isosurface(g_X,g_Y,g_Z,g_F,0);
    
    % Calculate Iso-Normals of the surface
    N=isonormals(g_X,g_Y,g_Z,g_F,FV.vertices);
    
    L=sqrt(N(:,1).^2+N(:,2).^2+N(:,3).^2)+eps;
    N(:,1)=N(:,1)./L; N(:,2)=N(:,2)./L; N(:,3)=N(:,3)./L;
    % Display the iso-surface
    figure, patch(FV,'facecolor',[1 0 0],'edgecolor','none'); view(3);camlight; axis equal
    % Invert Face rotation
    FV.faces=[FV.faces(:,3) FV.faces(:,2) FV.faces(:,1)];
    
    % Make a material structure
    material(1).type='newmtl';
    material(1).data='skin';
    material(2).type='Ka';
    material(2).data=[0.8 0.4 0.4];
    material(3).type='Kd';
    material(3).data=[0.8 0.4 0.4];
    material(4).type='Ks';
    material(4).data=[1 1 1];
    material(5).type='illum';
    material(5).data=2;
    material(6).type='Ns';
    material(6).data=27;
    
    clear OBJ
    OBJ.vertices = FV.vertices;
    OBJ.vertices_normal = N;
    OBJ.material = material;
    OBJ.objects(1).type='g';
    OBJ.objects(1).data='skin';
    OBJ.objects(2).type='usemtl';
    OBJ.objects(2).data='skin';
    OBJ.objects(3).type='f';
    OBJ.objects(3).data.vertices=FV.faces;
    OBJ.objects(3).data.normal=FV.faces;
    write_wobj(OBJ,'../Results/Torus/torus_nohermite250.obj');
end