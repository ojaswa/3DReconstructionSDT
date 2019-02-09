function [g_model, l_section_names] = load_3D (model_name)
%
% [g_model, l_section_names] = load_3D (model_name)
% 
% Loader function for 3D model and its cross-sections
% 
% Author: Ojaswa Sharma, <ojaswa@iiitd.ac.in>
%
%% Torus model (BBox origin - (0,0,0))
if strcmp(model_name, 'torus')
    l_model = read_wobj('../data/Torus/torus.obj');
    g_model.bbox  = [-4, -1, -4; 4, 1, 4]; % torus [min_x min_y min_z max_x max_y max_z]
    g_model.bbox_diag = 11.489125; % torus
    g_model.bbox_origin = (g_model.bbox(1,:) + g_model.bbox(2,:))/2;
    g_model.vertices = l_model.vertices;
    g_model.faces = l_model.objects(4).data.vertices; % torus
    g_sections = {};
    l_section_names = {
        '../data/Torus/torus_section0.obj'
        '../data/Torus/torus_section1.obj'
        '../data/Torus/torus_section2.obj'
        '../data/Torus/torus_section3.obj'
        };
end
