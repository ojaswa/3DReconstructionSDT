
% triangulate - interface to 'triangle' library. 
%
%   [vertex, face] = compute_constrained_delaunay(vertex, edges, holes);

function [nodes, triangles, output_node_markers, timeSpent] = triangle(nodes, boundary, holes, input_node_markers, args)
% Markers is an array of same length as nodes, containing a marker for
% certain nodes to be identified in the output. A zero marker is not
% appended to the node in inout file

if size(nodes,1)<size(nodes,2)
    nodes = nodes';
end
if size(boundary,1)<size(boundary,2)
    boundary = boundary';
end
if size(holes,1)<size(holes,2)
    holes = holes';
end

nodes = [(1:size(nodes,1))', nodes];
boundary = [(1:size(boundary,1))', boundary];
holes = [(1:size(holes,1))', holes];

n = size(nodes, 1);
m = size(boundary, 1);
h = size(holes, 1);

% nodes = roundToDP(nodes,4);
% Write out domain in PSLG format.
fid = fopen('mesh.poly', 'w');
if fid<0
    error('Unable to create file mesh.poly.');
end
fprintf(fid, '%d 2 0 1\n', n);
for i = 1:n
    fprintf(fid, '%d %f %f', nodes(i,1), nodes(i,2), nodes(i,3));
    fprintf(fid, ' %d', input_node_markers(i));  
    fprintf(fid, '\n');
end
fprintf(fid, '%d 0\n', m);
for j = 1:m
    fprintf(fid, '%d %d %d\n', boundary(j,1), ...
        boundary(j,2), ...
        boundary(j,3));
end
fprintf(fid, '%d\n', h);
for k = 1:h
    fprintf(fid, '%f %f %f\n', holes(k,1), holes(k,2), holes(k,3));
end
fclose(fid);

% Triangulate domain.
executable = [mfilename('fullpath') ' '];
cmd = [executable args ' mesh.poly'];
% cmd = ['triangle ' args ' mesh.poly'];
tic
status = system(cmd);
if status ~= 0
    fprintf('Exe crashed!\n');
end
timeSpent = toc;

% Read in new nodes.
fid = fopen('mesh.1.node', 'r');
if fid<0
    cmd = [executable ' -i' args ' mesh.poly'];          % No steiner opints will be added to boundary
    tic
    status = system(cmd, '-echo');
    if status ~= 0
        fprintf('Exe crashed!\n');
    end
    timeSpent = toc;
    fid = fopen('mesh.1.node', 'r');
    if fid < 0
        cmd = [executable ' ' args ' -Y mesh.poly'];          % No steiner opints will be added to boundary
        tic
        status = system(cmd, '-echo');
        if status ~= 0
            fprintf('Exe crashed!\n');
        end
        timeSpent = toc;
        fid = fopen('mesh.1.node', 'r');
        if fid<0
            cmd = [executable ' ' args ' -YY mesh.poly'];     % No steiner opints will be added to inside & boundary
            tic
            status = system(cmd, '-echo');
            if status ~= 0
                fprintf('Exe crashed!\n');
            end
            timeSpent = toc;
            fid = fopen('mesh.1.node', 'r');
            if fid<0
                cmd = [executable ' -Q mesh.poly'];     % No options
                tic
                status = system(cmd, '-echo');
                if status ~= 0
                    fprintf('Exe crashed!\n');
                end
                timeSpent = toc;
                fid = fopen('mesh.1.node', 'r');
            end
        end
    end
end
n = fscanf(fid, '%d', 1);
dummy = fscanf(fid, '%d', 3);
nodes = [];
output_node_markers = zeros(n, 1); % Markers array for output vertices
for i = 1:n
    point = [ fscanf(fid, '%d', 1)' fscanf(fid, '%f', 2)' ];
    nodes = [ nodes; point ];
    output_node_markers(i) = fscanf(fid, '%d', 1);
end
fclose(fid);

% Read in triangle mesh.
fid = fopen('mesh.1.ele', 'r');
m = fscanf(fid, '%d', 1);
dummy = fscanf(fid, '%d', 2);
triangles = [];
for j = 1:m
    triangles = [ triangles; fscanf(fid, '%d', 4)' ];
end
fclose(fid);

triangles = triangles(:,2:end)';
nodes = nodes(:,2:end)';

% Check for nan or infinity nodes
[x_1 x_2]= find(isnan(nodes) | isinf(nodes));

if(~isempty(x_1))
    disp('nan');
end

% Delete temporary files.
!rm -f mesh.poly mesh.1.poly mesh.1.node mesh.1.ele



