function area = getMaximumTriangleAreaHeuristic(sections)
%
% [parents, contours, holes] = classifyChains(vertices, chains)
% 
% Calculate an estimate of maximum triangle area for 2D triangulation
% using average chain length
% 
% Author: Ojaswa Sharma, <ojaswa@iiitd.ac.in>
%
% Input: section
%
% Output: area value
% 
sum_lengths = 0;
num_segments = 0;
for i=1:size(sections, 2)
    chains = sections{i}.chains;
    verts = sections{i}.vertices;
    if isempty(chains)
        continue
    end
    for j=1:size(chains, 2)
        vecs = verts(chains{j}, :) - circshift(verts(chains{j}, :), -1, 1);
        lengths = sqrt(dot(vecs, vecs, 2));
        sum_lengths = sum_lengths + sum(lengths, 1);
        num_segments = num_segments + size(vecs, 1);
    end
end

avg_length = sum_lengths/num_segments;
area = avg_length*avg_length/2;