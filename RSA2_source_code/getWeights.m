function W = getWeights(A,edges,D,type)
% getWeights: calculates edge weights
% 
% A                 adjacency matrix (entries are edge indices)
% edges             list of points corresponding to the edge 
% D                 distance map 
% type              type of weigth   
%
% W                 sparse matrix containing edge weights
%
% Copyright 2012-2013 Daniel Leitner. See license.txt for details.
% 

ne = length(edges); % number of edges
ew = zeros(1,ne); % edges weights

if strcmp(type,'radius')
    for i = 1 : ne
        p = edges{i};
        r = D(sub2ind(size(D), round(p(:,2)),round(p(:,1))));
        ew(i)=mean(r);
    end    
elseif strcmp(type,'rradius') % reciprocal radius
    for i = 1 : ne
        p = edges{i};
        r = D(sub2ind(size(D), round(p(:,2)),round(p(:,1))));
        ew(i)=1/(mean(r)+1e-11);
    end        
elseif strcmp(type,'area') 
    for i = 1 : ne
        p = edges{i};
        pmid = 0.5*(p(1:end-1,:) + p(2:end,:)); % segment mid points
        d = 2*D(sub2ind(size(D), round(pmid(:,2)),round(pmid(:,1))));  % diameter at midpoints      
        v = p(2:end,:) - p(1:end-1,:);        
        l = sqrt( v(:,1).^2 + v(:,2).^2 );
        ew(i)=sum(d.*l);
    end        
elseif strcmp(type,'length')
    for i = 1 : ne
        p_ = edges{i};        
        p = p_(2:end,:) - p_(1:end-1,:);
        ew(i) =sum(sqrt( p(:,1).^2 + p(:,2).^2 ));        
    end
elseif strcmp(type,'rvolume') 
    for i = 1 : ne
        p_ = edges{i};
        if size(p_,1)>1
            r = D(sub2ind(size(D), round(p_(:,2)),round(p_(:,1))));
            p = p_(2:end,:) - p_(1:end-1,:);
            l =sum(sqrt( p(:,1).^2 + p(:,2).^2 )) + 1e-11;
        end
        ew(i) = l / (mean(r)^2+0.5);
    end
else
    disp('invalid option');    
end

[r,c,s] = find(A);
s=ew(s);
W=sparse(r,c,s,size(A,1),size(A,2));
