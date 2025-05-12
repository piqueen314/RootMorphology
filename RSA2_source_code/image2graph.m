function [A,nodes,edges] = image2graph(imb, imcol, sigma, erosion, small, wb)
% imgae2graph: creates a graph representation of a binary image
%
% imb               binary image
% (imcol)           images containing dynamic information (default = [])
% (sigma)           for gauss filtering (default = 0)
% (erosion)         size of erosion before skeletonization (default = 0)
% (small)           lenght of small edges, that will be neglected (default
%                   = 0.3% of image width)
% (wb)              if a waitbar shall be shown (default = false)
%
% A                 adjacency matrix (entries are edge indices)
% nodes             points representing the nodes (|N|,2)
% edges             list of points corresponding to the edges {|E|}(|P|)
%
% Example:
% im = imread('artificial.tif');
% imb = mean(im,3) > 50; 
% [A,nodes,edges] = image2graph(imb,[],4);
% plotGraph(A,edges,A>0,'b');
%
%
% Copyright 2013 Daniel Leitner. See license.txt for details.
% 

if nargin<2
    imcol = [];
end

if nargin<3
    sigma = 0;
end

if nargin<4
    erosion = 0;
end

if nargin<5
    small = 3e-3*size(imb,1); 
end

if nargin<6
    wb = 0;
end

% skeleton
if wb, h=waitbar(0.25, 'Create skeleton...'); end
ims = skeleton(imb,erosion);

% graph
if wb, waitbar(0.5,h, 'Create graph...'); end
if isempty(imcol)
    [A, nodes, edges] = createGraph(double(ims));
else 
    [A, nodes, edges] = createGraph(double(ims), imcol);
end

%smooth edges
if wb, waitbar(0.75,h, 'Smooth edge points...'); end
if sigma>0
    edges = smoothEdgePoints(edges,sigma);
end

% remove small edges
if wb, waitbar(1,h, 'Remove small edges...'); end
[A, nodes, edges] = removeSmall(A, nodes, edges, small); 
[A, nodes, edges] = removeSmall(A, nodes, edges, small); 
[A, nodes, edges] = removeSmall(A, nodes, edges, small); 

if wb, delete(h); end



function ims = skeleton(imb, erosion)
% skeleton: create a skeleton from a binary image
%
% imb           binary image
% erosion       size of initial erosion
%
% ims          skeleton
%
% Copyright 2012 Daniel Leitner. See license.txt for details.
% 

if erosion>0
    ims = bwmorph(imb,'erode',erosion);
else 
    ims = imb;
end

for i = 1 : 20
    ims =bwmorph(ims,'thin',1);
    ims=bwmorph(ims,'close');
end
ims = bwmorph(ims,'thin',Inf);

% ims = bwlabel(ims); % consider only the largest connected component
% lind = find(ims);
% aC = accumarray(ims(lind),ones(size(lind)));
% [~,mi] = max(aC);
% ims = ims==mi;



function [A, nodes, edges] = createGraph(ims, imcol)
% createGraph: creates a graph from a skeleton
%
% ims              skeleton
% (imcol)          dynamic information (create a node if colour changes)
%
% A                adjacency matrix (entries are edge indices)
% nodes            points representing the nodes 
% edges            list of points corresponding to the edges 
%
% Copyright 2012-2013 Daniel Leitner. See license.txt for details.
% 

hp = [2:size(ims,1),size(ims,1)-1]; % indices: vertical  plus 1 
hm = [2,1:size(ims,1)-1]; 
vp = [2:size(ims,2),size(ims,2)-1];
vm = [2,1:size(ims,2)-1];

if nargin<2
    cc = 0; % color change
else
    cs = ims.*imcol;  % colored skelton
    csc = ( (cs(1:end, vp)==cs)+(cs(hm, vp)==cs)+(cs(hm, 1:end)==cs)+...
        (cs(hm, vm)==cs)+(cs(1:end, vm)==cs)+(cs(hp, vm)==cs)+...
        (cs(hp,1:end)==cs)+(cs(hp, vp)==cs) ).*ims;
    cc = csc<2 & ims; % color change
end

% count pixels in 8-neighbours
imc = (ims(1:end, vp)+ims(hm, vp)+ims(hm, 1:end)+ims(hm, vm)+...
    ims(1:end, vm)+ims(hp, vm)+ims(hp,1:end)+ims(hp, vp)).*ims;

% get nodes
im_nodes = bwlabel(imc==1 | imc>2 | cc); % leaf or branch or colorchange

% final node positions are at the center of the connected object
[I,J] = find(im_nodes);
lind = sub2ind(size(im_nodes),I,J);
aI = accumarray(im_nodes(lind),I);
aJ = accumarray(im_nodes(lind),J);
aC = accumarray(im_nodes(lind),ones(size(J)));
nodes = [aJ./aC, aI./aC]; 

% get edges
im_edges = bwlabel((~(im_nodes>0)) & (imc~=0));

[I,J] = find(im_edges);  % coordinates of the points
lind = sub2ind(size(im_edges),I,J);
edgeindex = im_edges(lind); % index of the edge the points belong to
noe = max(im_edges(:)); % number of edges
edges = cell(noe,1); % edge list 
AI = zeros(noe,1); % indices for adjacency matrix
AJ = zeros(noe,1);

N = [0 1; 1 1; -1 1; 1 0; -1 0; 0 -1; 1 -1;  -1 -1]; % neighbours 

for i = 1 : noe
    
    li = find(edgeindex==i); % indices of the points of edge i 
    X = [I(li),J(li)]; % coordinates of the points
    
    iN = kron(X,ones(8,1)) + repmat(N,size(X,1),1); % coordinates of neighbours    
    liN = sub2ind(size(im_nodes),iN(:,1),iN(:,2));    
    
    nv = im_nodes(liN); % values of neighbours (node indices)
    liN2 = find(nv); % start and end node      
    [n1, ni1] = min(nv(liN2));
    [n2, ni2] = max(nv(liN2));        
            
    % add the two node coordinates to the edge points
    [i1,j1] = ind2sub(size(im_nodes),liN(liN2(ni1)));
    [i2,j2] = ind2sub(size(im_nodes),liN(liN2(ni2)));
    X(end+1,:) = [i1,j1];
    X(end+1,:) = [i2,j2];
        
    % copy edge points into a list in the order from node 1 to node 2
    if n1~=n2 % remove self circles at a later point                
        P = zeros(size(X,1),1);  % order for the permutation
        ind = size(X,1)-1; % current index (start at node 1)
        c = 1;
        P(ind) = c; % first node
        c = c+1;
        P(end) = size(P,1); % last node
        for j = 1 : size(P,1)-2
            ind = find( P==0 & X(ind,1)<=X(:,1)+1 & X(ind,1)>=X(:,1)-1 & ...
                X(ind,2)<=X(:,2)+1 & X(ind,2)>=X(:,2)-1, 1, 'first' ); % find next point
            P(ind) = c;
            c = c+1;
        end        
    else
        P = 1 : size(X,1);
    end    
    X(end-1,:) = nodes(nv(liN2(ni1)),end:-1:1); % replace by real node position
    X(end,:) = nodes(nv(liN2(ni2)),end:-1:1);  
    X(P,1) = X(:,1); % order points
    X(P,2) = X(:,2); % order points
    
    if ~ismember([n1,n2], [AI, AJ],'rows') && ~ismember([n2,n1], [AI, AJ],'rows')
        AI(i) = n1; % indices for the adjacency matrix
        AJ(i) = n2;  
        edges{i} = [X(:,2),X(:,1)]; % store in list
    else
        mi = round(size(X,1)/2);
        nodes(end+1,:) = X(mi,[2,1]); % add a new node
        AI(i) = n1; 
        AJ(i) = size(nodes,1);
        AI(end+1) = size(nodes,1);
        AJ(end+1) = n2;        
        edges{i} = [X(1:mi,2),X(1:mi,1)];
        edges{end+1} = [X(mi:end,2),X(mi:end,1)];
    end      
    
end

% remove self circles in graph
ind = AI(:) == AJ(:);
edges(ind) = [];
AI(ind)=[];
AJ(ind)=[];

% build adjacency matrix
s = 1:length(edges);
A1 = sparse(AI, AJ, s, size(nodes,1),size(nodes,1)); 
A2 = sparse(AJ, AI, s, size(nodes,1),size(nodes,1)); 
A= A1 + A2.*(A1==0); % avoid n1->n2 & n2->n1 

function edges = smoothEdgePoints(edges,sigma)
% smoothEdges: gauss smooths the edge points
%
% edges         list of points corresponding to the edges {|E|}(|P|)             
% sigma         standard deviation for gauss filter
%
% Copyright 2013 Daniel Leitner. See license.txt for details.
% 

fg = fspecial('gauss',5,sigma);

for i = 1 : length(edges)
    X = edges{i};
    x1 = X(1,:);
    x2 = X(end,:);    
    X(:,1) = imfilter(X(:,1),fg,'replicate');
    X(:,2) = imfilter(X(:,2),fg,'replicate'); 
    X(1,:) = x1;
    X(end,:) = x2;
    edges{i}=X;
end



function [A,nodes,edges] = removeSmall(A,nodes,edges,small)
% removeSmall: removes terminal edges with small length
%
% A                 adjacency matrix (entries are edge indices)
% nodes             points representing the nodes (|N|,2)
% edges             list of points corresponding to the edges {|E|}(|P|)
% small             lenght of small edges, that will be neglected 
%
% Copyright 2013 Daniel Leitner. See license.txt for details.
% 

Slen = getWeights(A,edges,0,'length');

[ei,ej] = find(Slen<small & Slen>0);

c=0;
for i = 1 : length(ei)
    n1=ei(i);
    n2=ej(i);
    e = A(n1,n2); % edge number
    if e>0
        if size(edges{e},1) < small
            if nnz(A(n1,:))==1 || nnz(A(n2,:))==1 % only if it is an ending
                A(n1,n2)=0;
                A(n2,n1)=0;
                A(A>e) = A(A>e) - 1; % update edge indices
                edges(e) = [];
                c=c+1;
            end
        end
    end
end

% remove lonely nodes
c=1;
for i = 1 : size(A,1)
    if nnz(A(c,:))==0
        A(c,:) = [];
        A(:,c) = [];
        nodes(c,:) = [];
    else
        c=c+1;
    end
end
