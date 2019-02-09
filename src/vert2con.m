function [A,b] = vert2con(V)
% VERT2CON - convert a set of points to the set of inequality constraints
%            which most tightly contain the points; i.e., create
%            constraints to bound the convex hull of the given points
%
% [A,b] = vert2con(V)
%
% V = a set of points, each ROW of which is one point
% A,b = a set of constraints such that A*x <= b defines
%       the region of space enclosing the convex hull of
%       the given points
%
% For n dimensions:
% V = p x n matrix (p vertices, n dimensions)
% A = m x n matrix (m constraints, n dimensions)
% b = m x 1 vector (m constraints)


k = convhulln(V);
c = mean(V(unique(k),:));
V=V-repmat(c,[size(V,1) 1]);
A  = NaN*zeros(size(k,1),size(V,2));
rc=0;
for ix = 1:size(k,1)
    F = V(k(ix,:),:);
    if rank(F,1e-5) == size(F,1)
        rc=rc+1;
        A(rc,:)=F\ones(size(F,1),1);
    end
end
A=A(1:rc,:);
b=ones(size(A,1),1);
b=b+A*c';
% eliminate dumplicate constraints:
[null,I]=unique(num2str([A b],6),'rows');
A=A(I,:); % rounding is NOT done for actual returned results
b=b(I);
return

