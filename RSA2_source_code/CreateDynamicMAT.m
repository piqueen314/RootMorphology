%
% Helps to create a workspace that can be loaded into Root System Analyser,
% from an image sequence of individual images located in a folder (sorted
% in ascending order)
%
% please add measurement times in line 45. 
%
% Copyright 2013 Daniel Leitner. See license.txt for details.
% 

clear all;

folder_name = uigetdir;
files = dir(folder_name);
[~,I] = sort({files.name});

c=0;

disp('read files');
for j = 1 : length(files)
    
    i=I(j);
    
    if ~files(i).isdir
                
         try
            im = imread([folder_name, '/', files(i).name]); % open image
            c=c+1; % counter
            imb = mean(im,3) > 50; % create binary image
            D{c} = bwdist(bwmorph(~imb,'close')); % create distance function
            
            if exist('imcol','var')
                imcol = imcol + imb;
            else % initialize
                imcol = imb; 
            end
        catch
            disp([files(i).name ' is not an image']);
        end
        
    end
    
end

disp('create graph');
[A,nodes,edges] = image2graph(imcol>0, imcol);

times(1:c) = 1:c; % <----- MANUALLY set measurement times

tips = []; % nothing detected yet

save('worspace.mat','A','edges','nodes','D','times','tips'); % save relevant
