function plotGraph(A,edges,C,range)
% plotGraph: plots the edges of graph A
%
% A                  adjacency matrix
% edges              list of coordinates per edge
% C                  scalar for coloring
% range              (optional) color range or line options (default:
%                    range = 'w')
%
% Copyright 2013 Daniel Leitner. See license.txt for details.
% 

hold on;
set(gca,'YDir','reverse');
[r,c] = find(A);

if nargin<4
    range = 'w';
end

if ~ischar(range)
    map = colormap();
    x = linspace(0,1,64)';
    mins = range(1);
    maxs = range(2);
    f = @(i,j) (C(i,j)-mins)/(maxs-mins);
end


for i = 1 : size(r,1)
    
    ei = A(r(i),c(i));
    p = edges{ei};
    
    if ischar(range)
        if C(r(i),c(i))>0
            plot(p(:,1),p(:,2),range,'LineWidth',1);
        end
    else
        v = f(r(i),c(i));
        v = min(v,1);
        v = max(v,0);
        c_ = interp1(x,map,v);
        plot(p(:,1),p(:,2),'Color',c_,'LineWidth',1);
    end
end
