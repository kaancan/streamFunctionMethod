%% Symmetric coil calculation
% 23.04.2018 - M. Kaan Can
% Function to detect the surface boundaries and nodes associated with
% these. Works fine with most coil geometries, but always good to check.
% Inputs: 
% Coil geometry information:
% surface = 
%     x: X coordinates of coil, ncoilx1 array
%     y: Y coordinates of coil, ncoilx1  array
%     z: Z coordinates of coil, ncoilx1  array
%     n: normal vectors of each vertice, ncoilx3 matrix
%     faces: nfacesx3 matrix connecting vertices to form faces
% Outputs:
% edgeNodes : A vector that has the boundary nodes. If the
%   geometry has more than one surface edge, gives a matrix.

%%
function edgeNodes = detectEdgev2(surface)
tic
% Find the triangular face edges that is related to a single face. These
% are the connecting the nodes on the surface edge. Then, start by the
% first node and keep going until all the points on the edge is found. If
% there are edge nodes not related to his edge, start again and go on until
% all the edges and nodes on those are found.
trep = triangulation(surface.faces, surface.x, surface.y, surface.z);
e = edges(trep);
nedges = length(e);
nfaces = length(surface.faces);
ncoil = length(surface.x);

ee = zeros(nedges,1);
for i = 1:nfaces
    eface = combnk(surface.faces(i,:),2);
    for j = 1:3
        ee = ee + ismember(e,sort(eface(j,:)),'rows');
    end
end

eee = 2 - ee;
edgeNodesd = e(eee ==1,:);
edgeNodesd = unique(edgeNodesd(:));

dummy = edgeNodesd;
i = edgeNodesd(1);
exi = 1;
noe = 1;
n = 1;
 while(dummy)
        neighbours = 0;
        cnct = 0;
        for j = 1:nfaces
            if(sum(surface.faces(j,:) == i))
                cnct = [cnct j];
            end
        end
        neighbours = surface.faces(cnct(2:end),:);
        neighbours = unique(neighbours(:));
        neighbours(neighbours == i) = [];
        
        C = intersect(neighbours, edgeNodesd);
        dumb = i;
        if(i == dummy(1))
            i = C(2);
            exi = dumb;
        else
            if(length(C) > 2)
                dist(:) = ((surface.x(i)-surface.x(C)).^2 + (surface.y(i)-surface.y(C)).^2 + (surface.z(i)-surface.z(C)).^2);
                [a,b] = max(dist);
                C(b) = [];
            end 
            i = C(C ~= exi);
            exi = dumb;
        end
        edgeNodes(noe,n) = i;
        n = n+1;
        dummy(dummy == i) = [];
        
        if(isempty(dummy) | i == edgeNodesd(noe,1) & exi ~= edgeNodesd(noe,1))
            if(isempty(dummy))
                text1 = ['Edge detection completed.\nNumber of edges = ',int2str(noe),'\n'];
                fprintf(text1)
                
                for k = 1:noe
                    text = ['Number of points on edge ', int2str(k),' = ', int2str(length(unique(edgeNodes(k,:)))),'\n'];
                    fprintf(text)
                end
                break
            else
                i = dummy(1);
                noe = noe + 1;
                n = 1;
            end
        end  
 end
toc
end