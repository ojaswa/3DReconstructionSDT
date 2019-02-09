function [parents, contours, holes] = classifyChains(vertices, chains)
%
% [parents, contours, holes] = classifyChains(vertices, chains)
% 
% Classify input chain
% 
% Author: Ojaswa Sharma, <ojaswa@iiitd.ac.in>
%
% Input: vertices
%        chains
% Output: parents
%         contours
%         holes
%
n_chains = size(chains, 2);
parents = cell(n_chains, 1);
parent_sums = zeros(n_chains, 1);

for j=1:n_chains %Compute ownership
    for k=1:n_chains
        j_chain = vertices(chains{j},1:2);
        if k == j
            continue
        end
        k_chain = vertices(chains{k}, 1:2);
        if all(inpolygon(k_chain(:,1), k_chain(:,2), j_chain(:,1), j_chain(:,2))) % True => v_chain is inside q_chain
            parents{k} = [parents{k} j]; % Indicate j as a potential parent of k
        end
    end
end

% Count parent sums of chains with only one parent in the parent-list
n_parents = 1;
l_idx = find(cellfun('length', parents) == n_parents);
if ~isempty(l_idx)
    for j=1:size(l_idx, 1)
        l_parent_sum  = 1;
        l_idx_next = parents{l_idx(j,1)}(1, 1);
        while (~isempty(parents{l_idx_next}))
            l_idx_next = parents{l_idx_next}(1, 1); % This will have only one parent!
            l_parent_sum = l_parent_sum + 1;
        end
        parent_sums(l_idx(j,1)) = l_parent_sum;
    end
end

% Eliminate multiple parents in the parents table. Keep only the immidiate parent
% Corollary: Any chain with n parents will have all parents (in its list) referring to less than n parents.
%            Reducing parent list from low to high, we make sure that any chain will have parents referring to no more than 1 parent
n_parents = 2; % Start reducing elements with 2 parebts, increment and process longer lists over iterations
while  max(cellfun('length', parents)) > 1
    l_idx = find(cellfun('length', parents) == n_parents);
    if ~isempty(l_idx)
        for j=1:size(l_idx, 1)
            l_parent_sum = ones(n_parents,1); % keeps track of number of parents in heirarchy for each candidate parent in list
            for k=1:n_parents % For every parent in the candidate list, add up the parents in heirarchy
                l_idx_next = parents{l_idx(j,1)}(1, k);
                while (~isempty(parents{l_idx_next}))
                    l_idx_next = parents{l_idx_next}(1, 1); % This will have only one parent!
                    l_parent_sum(k,1) = l_parent_sum(k,1) + 1;
                end
            end
            % Select the longest chain of the parents
            [l_max, l_selected] = max(l_parent_sum, [], 1);
            parents{l_idx(j,1)} = parents{l_idx(j,1)}(1, l_selected);
            parent_sums(l_idx(j,1)) = l_max;
        end
    end
    n_parents = n_parents + 1;
end

% Deduce contour/hole classification from the reduced parent table
contours = find(~mod(parent_sums,2)); % Even # of parents
holes = find(mod(parent_sums,2)); % Odd # of parents
