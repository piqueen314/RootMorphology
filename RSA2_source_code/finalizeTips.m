function tips = finalizeTips(tips, v, ldt, A, edges)
% finalTips: finalizes tips structure (order, prelength, ct)
%
% tips              root tips
% v                 growth rate (px/day)
% ldt               lateral delay time (days)
% ft                final time (days) 
% A                 adjacency matrix 
% edges             list of coordinates corresponding to the edge (E,P)
%
% tips              root tips
%
% See also: RSA_GUI
%
% Copyright 2012-2013 Daniel Leitner. See license.txt for details.
%

c = max([tips.number]);

% number roots without id number
for i = 1 : length(tips)
    if tips(i).number==0
        c=c+1;
        tips(i).number = c;
    end
end

% find predecessors for roots without predecessor
i=1;
while i <= length(tips);
    if tips(i).predecessor==0
        for j = 1 : length(tips)            
            if i~=j && tips(i).path(1)==tips(j).path(end) % append
                tips(j).path = [tips(j).path, tips(i).path(2:end)];
                tips(j).length = tips(j).length+tips(i).length;
                tips(j).area = tips(j).area+tips(i).area;
                ind = find([tips.predecessor]==tips(i).number);
                for k = 1 : length(ind)
                    tips(ind(k)).predecessor = tips(j).number;
                end
                tips(i)=[];
                i=i-1;
                break; % j loop
            elseif i~=j && ismember(tips(i).path(1), tips(j).path(2:end-1)) % lateral
                tips(i).predecessor = tips(j).number;
                break; % j loop
            end
        end
    end
    i=i+1;
end

% set all root orders
for j = 1 : length(tips)
    order = 0;
    t = tips(j);
    while t.predecessor~=0        
        t = tips(t.predecessor==[tips.number]);        
        order = order +1;
    end    
    tips(j).order = order;    
end

% set prelength and approximated creation time
[~,idx] = sort([tips.order]); 
for i =  idx
    if tips(i).predecessor>0
        
        % prelength
        j = find(tips(i).predecessor==[tips.number]);
        path = tips(j).path(1:find(tips(i).path(1) == tips(j).path));        
        l=0;
        for j = 1 : length(path)-1
            p_ = edges{A(path(j),path(j+1))};   
            p = p_(2:end,:) - p_(1:end-1,:);
            l = l +sum(sqrt( p(:,1).^2 + p(:,2).^2 ));   
        end
        tips(i).prelength = l;        
        
        % approximated creation time            
        tips(i).ct = tips(tips(i).predecessor==[tips.number]).ct + ...
            tips(i).prelength/v + ldt;    
                
    else
        tips(i).prelength = 0;
        tips(i).ct = 0;
    end        
end
