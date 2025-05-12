function tips = emergingTips(itips, v, ldt, A, nodes, edges)
% emergingTips: create new root tips along existing root axis
%
% itips             initial lateral roots
% v                 growth rate (px/day)
% ldt               lateral delay time (days)
% nodes             node points
% A                 adjacency matrix 
% edges             list of coordinates corresponding to the edge (E,P)
%
% tips              emerging tips
%
% See also: RSA_GUI
%
% Copyright 2012-2013 Daniel Leitner. See license.txt for details.
%

tips = struct('number',{},'order',{},'ct',{},'predecessor',{},...
    'prelength',{},'length',{},'area',{},'path',{});

etc = 1; % emerging tip counter
number = max([itips.number])+1;

Sv = rebuildSv(itips,size(nodes,1));

for i = 1 : length(itips)
    
    tip = itips(i);
    
    % ct and prelength for every node along root 
    nct = zeros(length(tip.path),1);
    npl = zeros(length(tip.path),1);
    nct(1) = tip.ct+ldt;
    npl(1) = 0;
    l = 0;
    for j = 1 : length(tip.path)-1
        e = A(tip.path(j),tip.path(j+1));
        p_ = edges{e};
        p = p_(2:end,:) - p_(1:end-1,:);
        l = l +sum(sqrt( p(:,1).^2 + p(:,2).^2 ));
        nct(j+1) = l/v+tip.ct+ldt;
        npl(j+1) = l;
    end
    
    % create emerging root tips
    for j = 1 : length(tip.path)-1
        
        n = tip.path(j);
        n_ = find(A(n,:));
        
        for k = 1 :length(n_)
            
            if Sv(n,n_(k))==0
                tips(etc) = struct('number',number,'order',tip.order+1,...
                    'ct',nct(j),'predecessor',tip.number,'prelength',npl(j),...
                    'length',0,'area',0,'path',[n,n_(k)]);
                number=number+1;
                etc = etc +1;                
            end
            
        end
        
    end
    
end

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

