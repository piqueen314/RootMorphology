function writeRSML(filename,tslabel,N,nFiles, A, edges, nodes, tips, p, Slen, D)%(name, A, edges, nodes, tips, p, Slen, D)
% writes the RSML file
name = [filename, '(', num2str(N), ')'];
%
% Define Meta Data
%
m.version = 1;
m.unit = 'cm';
m.resolution = p.pxcm;
m.last_DASH_modified=date;
m.software = 'Root System Analyzer';
%m.user = '';
%m.file_DASH_key = ''
im.label = name;
m.image = im;
m.property_DASH_definitions.property_DASH_definition(1).label = 'length';
m.property_DASH_definitions.property_DASH_definition(1).type = 'real';
m.property_DASH_definitions.property_DASH_definition(1).unit = 'pixel';
%
m.property_DASH_definitions.property_DASH_definition(2).label = 'insertion_DASH_angle';
m.property_DASH_definitions.property_DASH_definition(2).type = 'real';
m.property_DASH_definitions.property_DASH_definition(2).unit = 'degree';
%
m.property_DASH_definitions.property_DASH_definition(3).label = 'parent_DASH_node';
m.property_DASH_definitions.property_DASH_definition(3).type = 'integer';
m.property_DASH_definitions.property_DASH_definition(3).unit = '1';
%
%if N>1
    m.time_DASH_sequence.label=tslabel;
    m.time_DASH_sequence.index=[num2str(N),'/',num2str(nFiles)];
    m.time_DASH_sequence.unified='true';
%end

%
tree.metadata = m;

%
% Define Roots
%
if ~isempty(tips)
    ind0 = find([tips.order]==0);
    for i = 1 : length(ind0)
        root0(i) = getRoot(A,edges,nodes,tips,ind0(i),Slen,p,D);
    end
else
    root0 = [];
end
tree.scene.plant.root = root0;

%
% Write file
%
rsmlname = [name, '.rsml'];
wPref.StructItem = false;
xml_write(rsmlname,tree,'rsml',wPref);
    


function X = getPoints(A,edges,nodes,tip)
% returns the geometry of a single root

p = tip.path;

if length(p)==1 % should not be 
    disp('warning, a root has only 1 node');
    X = nodes(p,:); 
else
    lind = sub2ind(size(A),p(1:end-1),p(2:end));
    X=[];
    for j = 1 : length(lind)
        ep=edges{A(lind(j))};
        if ~isequal(ep(1,:),nodes(p(j),:)) % edge points shall connect p(i-1) and p(i)
            ep(:,1) = ep(end:-1:1,1);
            ep(:,2) = ep(end:-1:1,2);
        end
        ep = ep(round(linspace(1,size(ep,1),max(size(ep,1)/4,2))),:); %reduce number of point
        X = [X; ep];
    end
end


function root = getRoot(A,edges,nodes,tips,i,Slen,p,D)
% converts a single root
id = tips(i).number; % find all successors
suc = find([tips.predecessor] == id);
% attribute
root.ATTRIBUTE.ID = id;

% geometry 
X = getPoints(A,edges,nodes,tips(i));
for j = 1 : size(X,1)
    root.geometry.polyline.point(j).ATTRIBUTE.x = X(j,1);
    root.geometry.polyline.point(j).ATTRIBUTE.y = X(j,2);    
end                                           

% properties
root.properties.length.ATTRIBUTE.value = tips(i).length;
%%%%%%%

if ~isempty(suc)
    for j = 1 : length(suc) % parse successors
        rootS(j) = getRoot(A,edges,nodes,tips,suc(j),Slen,p,D);
        x0 = rootS(j).geometry.polyline.point(1).ATTRIBUTE;  
        % add parent-node attribute 
        pn = find(x0.x==X(:,1) & x0.y==X(:,2),1,'first');          
        rootS(j).properties.parent_DASH_node.ATTRIBUTE.value = pn;      
    end    
    root.root = rootS;
else
    root.root = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ip = find(tips(i).predecessor==[tips.number]); % index predecessor
if ip>0
    p1 = tips(i).path; % node path of selected tip
    p2 = tips(ip).path; % node path of predecessor tip
    ei = find(p2==p1(1),1,'first');
    c0 = nodes(p1(1),:);
    if length(p1)==1 % should not be
        c1 = nodes(p1,:);
    else
        c1 = getPathPoint(A, Slen, nodes, edges, p1, p.Ra);
    end
    if ei==length(p2)
        c2 = nodes(p2(end),:);
    else
        c2 = getPathPoint(A, Slen, nodes, edges, p2(ei:end), p.Ra);
    end
    v1 = c1-c0;
    v2 = c2-c0;
    v1 = v1/norm(v1);
    v2 = v2/norm(v2);
    theta = acos(v1*v2')/pi*180;
else
    theta = [];
end
root.properties.insertion_DASH_angle.ATTRIBUTE.value=theta;

% functions
root.functions.function.ATTRIBUTE.domain='polyline';
root.functions.function.ATTRIBUTE.name='diameter';
for j = 1 : size(X,1)
    root.functions.function.sample(j).ATTRIBUTE.value=2*D(round(X(j,2)),round(X(j,1)));      
end
