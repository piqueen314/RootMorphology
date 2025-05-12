function tip = getTip(path,A,edges,D)
% getTip: creates a tip structure with length and area
%
% path              node path of the new lateral
% A                 adjacency matrix
% edges             list of coordinates corresponding to the edge (E,P)
% D                 distance matrix
%
% tip               root tip structure with correct length and area
%
% See also: RSA_GUI
%
% Copyright 2012-2013 Daniel Leitner. See license.txt for details.
%

l = 0;
a = 0;
for i = 1 : length(path)-1
    
    p_ = edges{A(path(i),path(i+1))};   
    p = p_(2:end,:) - p_(1:end-1,:);
    l = l +sum(sqrt( p(:,1).^2 + p(:,2).^2 ));    
    pmid = 0.5*(p_(1:end-1,:) + p_(2:end,:)); % segment mid points
    d = 2*D(sub2ind(size(D), round(pmid(:,2)),round(pmid(:,1))));  % diameters at midpoints
    a= a + sum(double(d).*sqrt( p(:,1).^2 + p(:,2).^2 ));
        
end

tip = struct('number',0,'order',0,'ct',0,'predecessor',0,...
    'prelength',0,'length',l,'area',a,'path',path);
