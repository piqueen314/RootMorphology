function bw=enhanceBW(bw, small)
% enhanceBW: repairs gaps by connecting all compenents 
% 
% bw            binary image
% small         remove small components based on the major axis length (px)
%
% bw            enhancedBW
% 
% Copyright 2014 Daniel Leitner. See license.txt for details.
%

maxlength = small*5; % maximal gap length (px)

if nargin  == 0 % for testing only
    bw = imread('C:\Users\Andrea\Desktop\images\BW_fotos_8K\MT.tiff');
    small = 10; 
end

% remove small components
disp('remove small');
CC = bwconncomp(bw);
s = regionprops(CC, 'MajorAxisLength', 'Area', 'PixelList');
p = [s.MajorAxisLength];
for i = 1 : length(p)
    if p(i)<=small
        fprintf('*');
        bw(sub2ind(size(bw),s(i).PixelList(:,2),s(i).PixelList(:,1)))=0;
    end
end
fprintf('\n');
% update s
CC = bwconncomp(bw);
s = regionprops(CC, 'MajorAxisLength', 'Area', 'PixelList');

% connect all parts to the biggest part
disp('connecting components');

[~,IX] = sort([s.Area],'descend');
boundaries = bwmorph(bw,'remove'); 
N = max(round(small/4),1); % very small

idx = CC.PixelIdxList{IX(1)}; 
idx_ = idx(boundaries(idx)==1); % amazingly clever 
[I,J] = ind2sub(size(bw),idx_);
X0=[I(1:N:end),J(1:N:end)]; % every Nth pixel of the boundary of component 1

for i = 2 : length(IX)

        idx = CC.PixelIdxList{IX(i)};
        idx_ = idx(boundaries(idx)==1); % amazingly clever
        [I,J] = ind2sub(size(bw),idx_);
        X=[I(1:N:end),J(1:N:end)]; % every Nth pixel of the  boundary of component IX(i)
        
        n = size(X0,1); 
        m = size(X,1);        
        %disp(n*m);
        
        if n*m<1e8 % safety first
           
            D = pdist2(X0,X);                        
            [~,midx] = min(D(:));            
            [ii,jj] = ind2sub(size(D),midx);
            
            if D(ii,jj) < maxlength
                
                fprintf('*'); % connect
                
                % line from X0(ii,:) to X(jj,:)                
                ml = max(abs(X0(ii,1)-X(jj,1)),abs(X0(ii,2)-X(jj,2)));
                y_ = round(linspace(X0(ii,1), X(jj,1), ml+2));
                x_ = round(linspace(X0(ii,2), X(jj,2), ml+2));                
                bw(sub2ind(size(bw),y_,x_))=2;                                
                
                % update X0
                X0 = [X0; X];
                                        
            else 
                fprintf('F'); % too far away                
                bw(idx)=0; % delelte X                                
            end
            
        else    
            fprintf('X'); % too big
        end

end
fprintf('\n');
