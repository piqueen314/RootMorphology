function tips = trackRoots(A,nodes,edges,D,p,tips,finaltips)
% dynamicFill: detects coherent roots in a graph representation
%
% A                 adjacency matrix (edges are edge indices)
% nodes             node coordinates
% edges             list of coordinates corresponding to the edge
% D                 distance map
% (p)               algorithm parameters (p.v: growth speed (px/day),
%                   p.ldt: lateral delay time (days), p.dt: time step
%                   (days) p.Rs: search scale (px), p.Ra: variational scale
%                   (px), p.ft: final threshold (px).
% (tips)            initial growing root tips
% (finaltips)       initial static roots
%
% tips              structure describing the individual roots
%
% Example:
% imb = mean(imread('artificial.tif'),3) > 50;
% [A,nodes,edges] = image2graph(imb);
% tips = trackRoots(A,nodes,edges,bwdist(~imb));
% plotRootSystem(A,edges,nodes,tips);
%
% Copyright 2012-2013 Daniel Leitner. See license.txt for details.
%

if nargin<5 % default parameters
    p=struct('v',75.8519,'ldt',3,'dt',0.0417,'Rs',227.5556,'Ra',37.9259,'ft',1);
end

if nargin<6 % default initial path
    [~,n1] = min(nodes(:,2)); % smallest y-coordinate
    n2 = find(A(n1,:),1,'first');
    tips(1) = struct('number',1,'order',0,'ct',0,'predecessor',0,...
        'prelength',0,'length',0,'area',0,'path',[n1 n2]);
end

if nargin<7 % create empty struct
    finaltips = struct('number',{},'order',{},'ct',{},'predecessor',{},...
        'prelength',{},'length',{},'area',{},'path',{});
end

%
% initialize
%
n = size(A,1); % number of nodes (constant)
Slen = getWeights(A,edges,D,'length'); % the length of each edge (constant)
Sa = getWeights(A,edges,D,'area'); % area of each edge (constant)

tips = rebuildTips(Sa,tips,finaltips);
branched = setBranched([tips,finaltips],n); % indicates if ndoes have already branched
rc = max([finaltips.number tips.number])+1; % root counter (to calculate new id)
time = 0; % simulation time

h = waitbar(0, 'Detection...'); % progress bar

%
% 'grow' until there are no active tips left
%
while ~isempty(tips)
    
    Sva = rebuildSva([tips([tips.length]>0), finaltips],Slen); % rebuild visited area matrix Sva
    Sv = rebuildSv([tips([tips.length]>0), finaltips], size(Slen,1));
    
    time = time+p.dt; % increase time step
    
    waitbar(sum(branched(:))/length(branched),h, 'Detection...');
    
    %
    % start iteration
    %
    i=1;
    while i<=length(tips) % for each growing root tip
        
        ind = []; % index of optimal path
        targetlength = (time-tips(i).ct)*p.v;
        
        % update root tip i
        while tips(i).length<targetlength  % grow until target length is reached
            
            %            
            % 1) find all possible paths
            %
            if tips(i).length==0
                if Sv(tips(i).path(1),tips(i).path(2))==0
                    x1 = nodes(tips(i).path(1),:);
                    R1 = 0;
                    newpaths = getPaths(A,Slen,tips(i).path,p.Rs);
                    ip = 1; % index in newpath where the newpath begins
                else
                    newpaths=[];
                end
            else 
                [x1,~,R1] = getPathPoint(A, Slen, nodes, edges, tips(i).path(end:-1:1), p.Ra);
                [~,ip] = getPathPoint(A, Slen, nodes, edges, tips(i).path(end:-1:1), p.Rs/2);
                newpaths = getPaths(A,Slen,tips(i).path(end-ip+1:end),p.Rs);
            end
            
            % determine quality of all paths
            if ~isempty(newpaths)
                q=zeros(length(newpaths),1);
                for j = 1 : length(newpaths)
                    q(j) = quality(A, nodes,edges, newpaths{j}(ip:end), Sva,...
                        Slen, Sa, x1, R1, p.Rs, p.Ra);
                end
                % pick maximum
                [m,ind] = max(q);
            else
                m=0;
            end
            
            %
            % 2) follow best path or finish
            %
            if m>p.ft && ~(ind==1 && tips(i).length~=0)
                
                if tips(i).length~=0 % growing tip
                    on = newpaths{ind}(ip); % last node of the path
                    nn = newpaths{ind}(ip+1); % first new node of winning path
                    tips(i).path = [tips(i).path, nn]; % add node
                else % emerging tip
                    on = newpaths{ind}(1);
                    nn = newpaths{ind}(2);
                end
                tips(i).length = tips(i).length+Slen(on,nn); % update tip'A length
                tips(i).area = tips(i).area+Sa(on,nn); % update tip'A area
                Sva(on,nn) = Sva(on,nn) + Sa(on,nn); % update visits
                Sva(nn,on) = Sva(nn,on) + Sa(on,nn);
                
                %
                % 3) create new laterals
                %
                if branched(on)==0
                    
                    n_ = find(A(on,:));
                    n_(n_==tips(i).path(end-2))=[]; % node before lateral
                    n_(n_==tips(i).path(end))=[]; % node after lateral
                    
                    for j = 1 : length(n_) % for all valid neighbours create laterals
                        pl = tips(i).length-Slen(on,nn); % prelength
                        ct = pl/p.v+tips(i).ct+p.ldt; % creation time
                        tips(end+1) = struct('number',rc, 'order',tips(i).order+1, ...
                            'ct',ct,'predecessor',tips(i).number,'prelength',pl,...
                            'length',0,'area',0,'path',[on, n_(j)]);
                        rc = rc+1;
                    end
                    branched(on) = 1; % this node has laterals and is not allowed to branch anymore
                    
                end
                
            else
                
                %
                % 4) finish root tip
                %
                if tips(i).length>0
                    finaltips(end+1) = tips(i);
                end
                tips(i) = [];
                ind = [];
                i=i-1;
                break
                
            end
                                    
        end % while length<targetlength
        
        %
        % reserve prospective path in Sva regarding tips with higher index
        %
        if ~isempty(ind) && length(newpaths{ind}) > ip+1 % future path of tip i
            fp = newpaths{ind}(ip+1:end); % future path (reservation)
            md = tips(i).area/tips(i).length; % current mean diameter of tip i
            lind = [sub2ind([n,n],fp(1:end-1),fp(2:end)),...
                sub2ind([n,n],fp(2:end),fp(1:end-1))];
            Sva(lind) = Sva(lind) + Slen(lind)*md;
        end
        
        i=i+1; % next root tip
        
    end % while i <= length(tips)
    
    
end % while ~isempty(tips) && time<finaltime

close(h); % close waitbar;

tips = [tips, finaltips];






function Sva = rebuildSva(tips,Slen)
% rebuildSva: rebuilds the sparse visited area matrix
%
n=size(Slen,1);
AI = []; AJ = []; A = [];
for i = 1 : length(tips) % active roots
    md = tips(i).area/tips(i).length; % mean diameter of tip i
    p = tips(i).path;
    AI = [AI, p(1:end-1)];
    AJ = [AJ, p(2:end)];
    A = [A, ones(1,length(p)-1)*md];
end
Sva=sparse([AI,AJ],[AJ,AI],[A,A],n,n).*Slen;



function Sv = rebuildSv(tips,n)
% rebuildSv: counts how often the edges are visited
%
AI = []; AJ = []; A = [];
for i = 1 : length(tips)
    p = tips(i).path;
    AI = [AI, p(1:end-1)];
    AJ = [AJ, p(2:end)];
    A = [A, ones(1,length(p)-1)];
end
Sv=sparse([AI,AJ],[AJ,AI],[A,A],n,n);



function b = setBranched(tips,n)
% setBranched: set the start nodes of the tips to branched
%
b = zeros(n,1);
for i = 1 : length(tips)
    b(tips(i).path(1))=1;
end



function tips = rebuildTips(Sa,tips,finaltips)
% rebuildTips: rebuilds areas of the roots with current Sa
%
n=size(Sa,1);
Sv = rebuildSv([tips,finaltips],n);
for i = 1 : length(tips)
    p = tips(i).path;
    lind = sub2ind([n,n],p(1:end-1),p(2:end));
    tips(i).area = full(sum(Sa(lind)./Sv(lind)));
end

