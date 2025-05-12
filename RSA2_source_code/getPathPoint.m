function [p,N,R] = getPathPoint(A, Slen, nodes, edges, path, R)
% getPathPoints: returns coordinate of the path's edges at distance R
%
% A                 adjacency matrix (edges are edge indices)
% Slen              length of the edges 
% nodes             node points
% edges             list of points corresponding to the edge (E,P)
% path              edges path
% R                 length of path
%
% p                 coordinate of the along the edges at distance R
% N                 index of the node in the path
% R                 approx R
%
% Copyright 2012-2013 Daniel Leitner. See license.txt for details.
%

lind = sub2ind(size(A),path(1:end-1),path(2:end)); % linear indices
lc_= cumsum(Slen(lind)); % cumulative lengths

if lc_(end)<=R % return end point
    
    N = length(path);
    ep=edges{A(lind(end))}; 
    if isequal(ep(1,:),nodes(path(end-1),:))         
        p = ep(end,:);
    else
        p = ep(1,:);
    end
    R = lc_(end);
    
else
    
    N = find(lc_>R,1,'first');
    
    % the wanted point is within this edge
    ep=edges{A(lind(N))};
    if ~isequal(ep(1,:),nodes(path(N),:)) % edges points shall connect p(i) and p(i+1)
        ep(:,1) = ep(end:-1:1,1);
        ep(:,2) = ep(end:-1:1,2);
    end
    
    if N>1
        l = lc_(N-1);
    else
        l = 0;
        p(1,:) = ep(1,:);
    end
    
    i = 2;
    while l < R && i<=size(ep,1)          
        l = l + sqrt(sum((ep(i,:)-ep(i-1,:)).^2));
        i = i+1;
    end
    p = ep(i-1,:);
    R = l;
    
end


