function u_polytopes = remove_duplicates(polytopes)
%
% u_polytopes = remove_duplicates(polytopes)
% 
% Remove duplicate enteries from polytope faces
% 
% Author: Ojaswa Sharma, <ojaswa@iiitd.ac.in>
%
%
m = 1;
for i = 1 : size(polytopes,2)
    l_poly = polytopes{i};
    verts = l_poly.vertices;
    faces = l_poly.faces;
    
    if(length(faces) < 3)
        continue;
    end
    
    % Round vertices upto 10 decimal places
    verts = roundToDP(verts,10);
    [u_verts, ia, ic] = unique(verts,'rows','stable');
    n2 = size(u_verts,1);
    n1 = size(verts,1);
    
    if(n1 == n2)
        u_polytopes{m}.vertices = u_verts;
        u_polytopes{m}.faces = faces;
        m = m + 1;
    else
        counts = arrayfun( @(x)sum(ic == x), unique(ic));
        nums = find(counts > 1);
        
        % Remove duplicate entries from faces
        for l = 1 : length(nums)
            dup_ids = find(ic == nums(l));
            for j = 1 : length(faces)
                f_ids = [];
                for k = 1 : length(dup_ids)
                    f_ids = [f_ids; find(faces{j} == dup_ids(k))'];
                end
                if(~isempty(f_ids))
                    faces{j}(f_ids) = dup_ids(1);
                    faces{j} = unique(faces{j},'stable');
                    if(length(faces{j}) < 3)
                        faces{j} = [];
                    end
                end
            end
            faces(cellfun('isempty',faces)) = [];
        end
        
        % Map face ids to new unique node ids
        for j = 1 : length(faces)
            f_idx = arrayfun(@(x)find(ia==x,1),faces{j});
            faces{j} = f_idx;
        end
        if(~isempty(faces))
            u_polytopes{m}.vertices = u_verts;
            u_polytopes{m}.faces = faces;
            m = m + 1;
        end
    end
end

end