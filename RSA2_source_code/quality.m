function q = quality(A, nodes, edges, np, Sva, Slen, Sa, x1, R1, Rs, Ra)
% quality: returns the quality of path
%
% A                 adjacency matrix (entries are edge indices)
% nodes             points representing the nodes (|N|,2)
% edges             list of points corresponding to the edges {|E|}(|P|)
% np                path (vector of nodes)
% Sva               visited area
% Slen              length of the edges
% Sa                area of the edges
% x1
% R1
%
% q                 quality of the path
%
%
% Copyright 2012-2013 Daniel Leitner. See license.txt for details.
%

E=2;

if length(np)>1
    
    lind = sub2ind(size(A),np(1:end-1),np(2:end));
    npl = sum(Slen(lind));  % new path length
    nplR = min(npl,Rs); 
    theta = 1- (npl-nplR)/Slen(lind(end)); % percentage of last edge
    
    % new path area    
    npaR = sum(Sa(lind(1:end-1))) + theta*Sa(lind(end)) - ...
        1.5 *(sum(Sva(lind(1:end-1))) + theta*Sva(lind(end))); % px^2    
    mina = min([Sa(lind(1:end-1)), Sa(lind(end))]);
        
    % straightness    
    [x2,~,R2] = getPathPoint(A, Slen, nodes, edges, np, Ra);
    straight = sqrt(sum((x2-x1).^2))/(R1+R2);            
        
    q = npaR/((10+nplR)*2) * straight^E * (mina>1);  % units is pixel 
    
    % npaR/(2*nplR) = radius
    
else
    
    q=0;
    
end


