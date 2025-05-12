function h = plotRootSystem(A, edges, nodes, tips, maxorder)
% plotRootSystem: plots detected root axis
%
% A                 adjacency matrix 
% nodes             coordinates representing the nodes (|N|,2)
% edges             list of coordinates corresponding to the edges {|E|}(|P|)
% tips              roots
% (maxorder)        plots the roots up to branching order maxorder
%                   (default = Inf)
%
% h                 handles 
%
% Copyright 2012-2013 Daniel Leitner. See license.txt for details.
%

if nargin<5
    maxorder = Inf;
end

hold on;
set(gca,'YDir','reverse');
opts = {'r','g','c','m','y','k'};
h = [];

for i = 1 : length(tips)
    
    p = tips(i).path;
    n = tips(i).order;

    if n<=maxorder
        
        n = mod(n,6)+1;        
        lind = sub2ind(size(A),p(1:end-1),p(2:end));
        
        for j = 1 : length(lind)
            ep=edges{A(lind(j))};
            if ~isequal(ep(1,:),nodes(p(j),:)) % edge points shall connect p(i-1) and p(i)
                ep(:,1) = ep(end:-1:1,1);
                ep(:,2) = ep(end:-1:1,2);
            end
            h = [h; plot(ep(:,1),ep(:,2),opts{n},'LineWidth',2)];   
            set(h(end),'UserData',num2str(i));
        end
        
        h = [h; plot(ep(end,1),ep(end,2),[opts{n},'*'])]; % marker for root tip
        set(h(end),'UserData',num2str(i));
        
    end
    
end
