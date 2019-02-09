function gradients = computeGradients(nodeIDs, plane_normal, chains, chain_verts)
%
% gradients = computeGradients(nodeIDs, plane_normal, chains, chain_verts)
% 
% Gradient computation
% 
% Author: Ojaswa Sharma, <ojaswa@iiitd.ac.in>
%
gradients = zeros(size(nodeIDs, 1), 3);
q_nz_nodeIDs = find(nodeIDs > 0);
nz_node_IDs  = nodeIDs(q_nz_nodeIDs);

% Perform a membsership test nodeID <-> curveID
membership = zeros(size(q_nz_nodeIDs, 1), 1);
for i=1:size(chains, 2)
    chain = chains{i};
    membership = membership + i*ismember(nz_node_IDs, chain);
end
%assert isequal(unique(membership), [1:size(curves, 2)]');

% Compute normal for each point
for i=1:size(q_nz_nodeIDs, 1)
    chain = chains{membership(i)};
    n = size(chain, 1);
    nodeID = nz_node_IDs(i);
    
    % Find out index of nodeID in its chain
    idx_cur = find(chain == nodeID);
    idx_nxt = mod(idx_cur,n)+1;
    idx_prv = mod(idx_cur-2,n)+1;
    
    % Compute gradient in the plane -- outgoing
    v01 = normalizeVector3d(chain_verts(chain(idx_cur), :) - chain_verts(chain(idx_prv), :));
    v12 = normalizeVector3d(chain_verts(chain(idx_nxt), :) - chain_verts(chain(idx_cur), :));
    grad01 = cross(v01, plane_normal);
    grad12 = cross(v12, plane_normal);
    
    val = grad01 + grad12;
    if any(isnan(val))
        gradients(q_nz_nodeIDs(i), :) = [0 0 0];
    elseif(val == 0)
        gradients(q_nz_nodeIDs(i), :) = [0 0 0];
    else
        gradients(q_nz_nodeIDs(i), :) = normalizeVector3d(val);
    end
    
end