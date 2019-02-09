function x = roundToDP(x,n)
%
% [parents, contours, holes] = classifyChains(vertices, chains)
% 
% Round the matrix x to n decimal places
% 
% Author: Ojaswa Sharma, <ojaswa@iiitd.ac.in>
%
for i = 1 : size(x,1)
    for j = 1 : size(x,2)
        x(i,j) = round(x(i,j) * 10^n) / 10^n;
    end
end
