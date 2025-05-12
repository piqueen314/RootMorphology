function paths = getPaths(A, Slen, path, r)
% getPaths: returns all possible node paths in the graph 
%
% stops if the edges of the node path are longer than the search radius
% 
% A                 adjacency matrix 
% Slen              length of edges 
% path              initial node path
% r                 search radius
%
% paths             a list of node paths
%
% See also: getWeights
% 
% Copyright 2012-2013 Daniel Leitner. See license.txt for details.
% 

global pcount; % total path counter

pcount = 1;
paths{1} = path; % paths

c = 1; % depth counter
newpaths = getPathsR(A,Slen,path, r, c);
c=c+1;

while ~isempty(newpaths)
   
    paths = [paths, newpaths];    
    newpaths = getPathsR(A, Slen, path, r, c);
    c=c+1;
    
end

clear pcount;


% recursive function
function paths = getPathsR(A,Slen,path,rad,depth)

global pcount;

n = path(end);
n_ = find(A(n,:)); % find neighbours

paths = {};

c=1;
for i = 1 : length(n_)
    
    if ~ismember(n_(i),path)
        
        l=Slen(n,n_(i));
        if l<rad && depth>1 && pcount<500
            
            npaths = getPathsR(A,Slen,[path,n_(i)], rad-l, depth-1);
            paths(c:c+length(npaths)-1) = npaths;
            c=c+length(npaths);
            pcount = pcount+length(npaths);
            
        elseif depth==1 % only return the path if depth is reached
            
            paths{c} = [path,n_(i)];
            c=c+1;
            pcount = pcount+1;
            
        end
        
    end
    
end
