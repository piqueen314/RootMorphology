function varargout = RSA_GUI(varargin)
% RSA_GUI MATLAB code for RSA_GUI.fig
%      RSA_GUI, by itself, creates a new RSA_GUI or raises the existing
%      singleton*.
%
%      H = RSxA_GUI returns the handle to a new RSA_GUI or the handle to
%      the existing singleton*.
%
%      RSA_GUI('CALLBACK',hObject,eventData,h,...) calls the local
%      function named CALLBACK in RSA_GUI.M with the given input arguments.
%
%      RSA_GUI('Property','Value',...) creates a new RSA_GUI or raises theh
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RSA_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RSA_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIh

% Edit the above text to modify the response to help RSA_GUI
% Last Modified by GUIDE v2.5 25-Feb-2014 14:40:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @RSA_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @RSA_GUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before RSA_GUI is made visible.
function RSA_GUI_OpeningFcn(hObject, eventdata, h, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% h          structure with h and user data (see GUIDATA)
% varargin   command line arguments to RSA_GUI (see VARARGIN)

% Choose default command line output for RSA_GUI
h.output = hObject;

% Update h structure
guidata(hObject, h);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using RSA_GUI.
% if strcmp(get(hObject,'Visible'),'off')
%     plot(rand(5));
% end

% UIWAIT makes RSA_GUI wait for user response (see UIRESUME)
% uiwait(h.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RSA_GUI_OutputFcn(hObject, eventdata, h)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% h          structure with h and user data (see GUIDATA)

% Get default command line output from h structure
varargout{1} = h.output;

function edit1_Callback(hObject, eventdata, h)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h          structure with h and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, h)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h          empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, h)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, h)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, h)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, h)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, h)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h          structure with h and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, h)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h          empty - h not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
updateView(hObject,1);

% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4

h=guidata(hObject);
h.sroot = 0;
h.sn1 = 0;
guidata(hObject,h);
updateView(hObject,1);

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end









%
% --- OPEN
%
function uipushtool1_ClickedCallback(hObject, eventdata, h)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

[file,path,fi] = uigetfile({'*.jpg;*.tif;*.png;*.gif','Image files';'*.mat','Matlab workspace';'*.*', 'Image sequence' });

if ~isequal(file, 0)
  
    if fi==2 % Matlab file
        
        hl = load([path file]); % load
        % store userdata
        h.edges = hl.edges;
        h.nodes = hl.nodes;
        h.D = hl.D;
        h.tips = hl.tips;
        h.A = hl.A;
        h.times = hl.times;
        % set drop down menu
        str={};
        for i = 1 : length(h.times)
            str{i}  = sprintf('Measurement %i (Day %i)',i,h.times(i));
        end
        set(h.popupmenu4,'String',str);
        % update gui
        if isfield(hl, 'imwidth')
            set(h.edit1,'String',num2str(hl.imwidth));
        end
        if isfield(hl, 'searchradius')
            set(h.edit2,'String',num2str(hl.searchradius));
        end
        if isfield(hl, 'angleradius')
            set(h.edit3,'String',num2str(hl.angleradius));
        end
        h.nFiles=length(h.times);
        h.label = file;
        
    elseif fi==1 % image file
        
        % create graph
        disp('read image');
        im = imread([path file]); % read binary image
        if size(im,3)>1 %            
%             disp('tresholding image');
%             im = rgb2gray(im);
%             level = graythresh(im);
%             imb = im2bw(im,level);            
            disp('clustering & threholding image'); % slow but better
            im1 = im(:,:,1);
            im2 = im(:,:,2);
            im3 = im(:,:,3);
            im1 = im1(:);
            im2 = im2(:);
            im3 = im3(:);
            im_ = double([im1 im2 im3]); % there might be a faster way            
            [idx,~] = kmeans(im_,2);            
            imb = reshape(idx,size(im,1),size(im,2))==2;   
        else
            disp('tresholding image'); % for gray scale
            level = graythresh(im);
            imb = im2bw(im,level);   c
        end
        % inverse in case background is white
        iw = sum(imb(:)>=1)/(size(imb,1)*size(imb,2));
        if iw>0.5 
            imb = ~imb;
            disp('inverting image');
        end
        %enhance the original image
        small = 5e-3*size(imb,1);
        fprintf('enhance image: diameter of small components is %g px\n',small);
        imb = enhanceBW(imb,small);
        % create graph (roots = 1, background = 0)        
        [A,nodes,edges] = image2graph(imb, [], small, 0, small, 1);
        % store userdata
        h.A = A;
        h.edges = edges;
        h.nodes = nodes;
        h.times= 28; % default value
        h.D{1} = bwdist(~imb);
        h.tips{1} = [];
        set(h.popupmenu4,'String','Single Measurement');
        h.nFiles=1;
        h.label = file;        
        
    elseif fi==3 % image sequence ('folder')
        
        files = dir(path);
        [~,I] = sort({files.name});
        c=0;
        disp('read files');
        for j = 1 : length(files)
            i=I(j);
            if ~files(i).isdir
                try
                    imb = imread([path, '/', files(i).name]); % open image
                    disp(['reading ' path, '/', files(i).name]);                                        
                    if size(imb,3)>1
                        imb = mean(imb,3)>128;
                    end
                    iw = sum(imb(:))/(size(imb,1)*size(imb,2));
                    if iw>0.5 % inverse in case background is white
                        imb = ~imb;
                    end
                    small = 1e-3*size(imb,1);
                    imb = enhanceBW(imb,small); % enhance
                    if exist('imcol','var')
                        imcol = imcol + imb;
                    else % initialize
                        imcol = imb;
                    end
                    c=c+1; % counter
                    h.D{c} = bwdist(~imb); % create distance function                    
                catch
                    disp([files(i).name ' is not a suitable image']);
                end
                
            end
        end % for
%         n=equalStrs(files(end).name,files(end-1).name); % woher kommt denn die funktion?
%         str1=files(end).name;
%         file = [str1(1:n-1),str1(end-3:end)];
         h.label = 'sequence'; %file;
        [A,nodes,edges] = image2graph(imcol>0, imcol, small, 0, small, 1);
        % store userdata
        h.A = A;
        h.edges = edges;
        h.nodes = nodes;
        h.times= 1:c; % TODO
        h.nFiles=length(h.times);
       
        for i = 1:c 
            h.tips{i} = [];
        end
        % set drop down menu
        str={};
        for i = 1 : length(h.times)
            str{i}  = sprintf('Measurement %i (Day %i)',i,h.times(i));
        end
        set(h.popupmenu4,'String',str);
    end
    
    h.Slen = getWeights(h.A,h.edges,h.D,'length');
    
    h.sn1 = 0; % selected node 1
    h.sroot = 0; % selected roots
    h.oldsel = []; % old selection
    h.scatter = []; % scatter plot handle
    h.filename = [path file]; % save the name
    
    % update h structure
    guidata(hObject, h);
    
    % plot results
    updateView(hObject,2);
end



%
% --- SAVE
%
function uipushtool2_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=guidata(hObject);

if isfield(h,'filename')
    
    name = h.filename;
    
    % save relevant
    s = struct;
    s.A = h.A;
    s.edges = h.edges;
    s.nodes = h.nodes;
    s.times = h.times;
    s.D = h.D;
    s.tips = h.tips;
    s.imwidth = str2double(get(h.edit1,'String'));
    s.searchradius = str2double(get(h.edit2,'String'));
    s.angleradius = str2double(get(h.edit3,'String'));
    save([name(1:end-4), '.mat'], '-struct', 's');
    
    % copy to workspace
    assignin('base','A',h.A);
    assignin('base','edges',h.edges);
    assignin('base','nodes',h.nodes);
    assignin('base','D',h.D);
    assignin('base','tips',h.tips);
    assignin('base','times',h.times);
    
    msgbox(['File saved: ', name(1:end-4), '.mat'] );
    
end


%
% --- DETECT (this measurement)
%
function pushbutton1_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

h=guidata(hObject);

if isfield(h,'nodes')
    
    N = get(h.popupmenu4,'Value'); % current measurement
    
    if isempty(h.tips{N}) % use node with smallest y-coordinate
        [~,n1] = min(h.nodes(:,2));
        n2 = find(h.A(n1,:),1,'first');
        h.tips{N} = struct('number',1,'order',0,'ct',0,'predecessor',0,...
            'prelength',0,'length',0,'area',0,'path',[n1 n2]);
    end
    
    p = getParameters(h);
    h.tips{N} = finalizeTips(h.tips{N}, p.v, p.ldt, h.A, h.edges); % update root tips
    tips = emergingTips(h.tips{N}, p.v, p.ldt, h.A, h.nodes, h.edges); % create emerging tips
    h.tips{N} = trackRoots(h.A,h.nodes,h.edges,h.D{N},p,[h.tips{N},tips]); % detect
    % deselect all
    h.sn1 = 0;
    h.sroot = 0;
    % update h structure
    guidata(hObject, h);
    updateView(hObject,1);
    
end

%
% --- DETECT (next measurement)
%
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=guidata(hObject);

if isfield(h,'nodes')
    
    N = get(h.popupmenu4,'Value');
    
    if ~isempty(h.tips{N}) && N+1<=length(h.times) % growing initial roots
        
        p = getParameters(h); % Parameters
        h.tips{N} = finalizeTips(h.tips{N}, p.v, p.ldt, h.A, h.edges); % update root tips
        tips = emergingTips(h.tips{N}, p.v, p.ldt, h.A, h.nodes, h.edges); % create emerging tips
        h.tips{N+1} = trackRoots(h.A,h.nodes,h.edges,h.D{N+1},p,[h.tips{N},tips]); % detect
        % deselect all
        h.sn1 = 0;
        h.sroot = 0;
        % check
        set(h.popupmenu4,'Value',N+1);
        % update h structure
        guidata(hObject, h);
        updateView(hObject,1);
    end
    
end



function p = getParameters(h)
% returns parameters for the detection algorithm
%
p.pxcm = size(h.D{1},2)/str2double(get(h.edit1,'String'));
p.Rs = str2double(get(h.edit2,'String'))*p.pxcm; % search scale (px)
p.Ra = str2double(get(h.edit3,'String'))*p.pxcm; % angular scale (px)
p.v = 2*p.pxcm; % growth speed = 2 (cm/day) = p.v (px/day)
p.ldt = 3; % lateral delay time (days) -> apical zone of 6 cm
p.dt = 1/24; % time step (day)
p.ft = 0.5; % final threshold




%
% --- SELECT
%
function axes3_ButtonDownFcn(hObject, eventdata, h)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

h=guidata(hObject);

if isfield(h,'nodes')
    
    p = get(hObject, 'Parent'); % Obtain click coordinates
    c = get(p, 'CurrentPoint');
    l = c(1,1:2);
    N = get(h.popupmenu4,'Value');
    
    if get(h.radiobutton2,'Value') % select root
        
        if ~isempty(h.tips{N})
            i = get(hObject, 'UserData');
            if ~isempty(i)
                h.sroot = str2double(i);
            else
                h.sroot=0;
            end
            h.sn1=0;
            guidata(hObject,h);
            updateView(hObject,0);
        end
        
    elseif get(h.radiobutton3,'Value') % add root
        
        k = dsearchn(h.nodes, l);
        if h.sn1==0
            h.sn1=k;
            h.sroot=0;
            guidata(hObject,h);
            updateView(hObject,0);
        else
            h.sn2=k;
            if get(h.popupmenu2,'Value')==2 % maximal diameter
                type = 'rvolume'; % minimize reciprocal radius
            elseif get(h.popupmenu2,'Value')==1 % minimal length
                type = 'length';
            end
            % find shortest path
            G = getWeights(h.A,h.edges,h.D{N},type);
            
            [path, ~] = dijkstra_sparse(G, h.sn1, h.sn2); % thanks to Y Simson
            % no need for bioninformatics toolbox anymore
            %[~, path, ~]=graphshortestpath(G,h.sn1,h.sn2); 
            
            % add root tip
            nt = getTip(path,h.A,h.edges,h.D{N});
            if isempty(h.tips{N})
                nt.number = 1;
                h.tips{N} = nt;
            else
                nt.number = max([h.tips{N}.number])+1;
                h.tips{N}(end+1) = nt;
            end
            p = getParameters(h);
            h.tips{N} = finalizeTips(h.tips{N}, p.v, p.ldt, h.A, h.edges);
            h.sn1=0;
            h.sroot = 0;
            guidata(hObject,h);
            updateView(hObject,1);
        end
        
    elseif get(h.radiobutton6,'Value') % remove
        
        i = get(hObject, 'UserData');
        if ~isempty(i)
            oldtip = [];
            h.sroot = str2double(i);
            if N>1
                oldtip = h.tips{N-1}([h.tips{N-1}.number]==h.tips{N}(h.sroot).number);
            end
            ind = find([h.tips{N}.predecessor]==h.tips{N}(h.sroot).number);
            for i = 1 : length(ind)
                h.tips{N}(ind(i)).predecessor=0;
            end
            if isempty(oldtip)
                h.tips{N}(h.sroot) = [];
            else
                h.tips{N}(h.sroot) = oldtip;
            end
        end
        p = getParameters(h);
        h.tips{N} = finalizeTips(h.tips{N}, p.v, p.ldt, h.A, h.edges);
        h.sn1=0;
        h.sroot = 0;
        guidata(hObject,h);
        updateView(hObject,1);
        
    end
    
end


%
% --- RESET
%
function pushbutton5_Callback(hObject, eventdata, h)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

h=guidata(hObject);
if isfield(h,'times')
    h.tips=cell(length(h.times),1);
    h.sn1=0;
    h.sroot=0;
    guidata(hObject,h);
    updateView(hObject,1);
end


%
% --- change radio button selection (DESELECT)
%
function uipanel1_SelectionChangeFcn(hObject, eventdata, h)
% hObject    handle to the selected object in uipanel1
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% h    structure with h and user data (see GUIDATA)

h=guidata(hObject);
h.sn1=0;
h.sroot=0;
guidata(hObject,h);
updateView(hObject,0);
zoom off;
pan off;


%
% --- UPDATE AXIS AND LABELS
%
function updateView(hObject,mode)

h=guidata(hObject);
hObject = gcbf; % since the axis might be deleted

% which measurement
N = get(h.popupmenu4,'Value');

if isfield(h,'nodes')
    
    if mode>0 % redraw everything
        
        if mode==1
            v = axis;
        end
        cla; % clear axis
        h.oldsel = []; % gone
        h.scatter = []; % gone
        
        % draw image
        colormap('gray');
        ih = imagesc(h.D{N}>0);
        hold on;
        set(ih,'ButtonDownFcn',{@axes3_ButtonDownFcn, h});
        if mode==1
            axis(v);
        else
            axis tight;
        end
        axis equal;
        
        % draw graph
        [r,c] = find(h.A);
        for i = 1 : size(r,1)
            ei = h.A(r(i),c(i));
            p = h.edges{ei};
            line('XData',p(:,1),'YData',p(:,2),'Color','b','LineWidth',1);
        end
        
        % draw axis
        pxcm = size(h.D{1},2)/str2double(get(h.edit1,'String'));
        for i = 1 : 5
            ticksx(i) = (i-1)/4*size(h.D{1},2);
            strx{i} = sprintf('%2.1f cm',(i-1)/4*size(h.D{1},2)/pxcm);
            ticksy(i) = (i-1)/4*size(h.D{1},1);
            stry{i} = sprintf('%2.1f cm',(i-1)/4*size(h.D{1},1)/pxcm);
        end
        set(gca,'XTick',ticksx);
        set(gca,'XTickLabel',strx);
        set(gca,'YTick',ticksy);
        set(gca,'YTickLabel',stry);
        
        % update root system
        i = get(h.popupmenu3,'Value');
        i(i==1) = Inf;
        i(i==2) = 2;
        i(i==3) = 1;
        if i~=4
            lh = plotRootSystem(h.A,h.edges,h.nodes,h.tips{N},i);
            set(lh,'ButtonDownFcn',{@axes3_ButtonDownFcn, h});
        end
        
        set(gcf,'Name',h.filename);
        
    end % end mode 1
    
    % draw nodes
    if ~isempty(h.scatter)
        delete(h.scatter);
        h.scatter = [];
    end
    h.scatter =  plot(h.nodes(:,1),h.nodes(:,2),'bo');
    set(h.scatter,'ButtonDownFcn',{@axes3_ButtonDownFcn, h});
    
    % delete old selection
    if ~isempty(h.oldsel)
        for i = 1 : length(h.oldsel)
            delete(h.oldsel(i));
        end
        h.oldsel = [];
    end
    
    % draw selection: ROOT
    k = h.sroot;
    % Cece
    tipSTATS = {};
    if k>0
        p = h.tips{N}(k).path;
        h.oldsel = zeros(length(p)-1,1);
        for i = 1 : length(p)-1
            e = h.A(p(i), p(i+1));
            eps = h.edges{e};
            h.oldsel(i) = plot(eps(:,1),eps(:,2),'y*-');
        end
        set(h.oldsel,'ButtonDownFcn',{@axes3_ButtonDownFcn, h});
        
        % update root text
        tip = h.tips{N}(k);
        
        pxcm = str2double(get(h.edit1,'String')) / size(h.D{1},2);
        l = tip.length * pxcm ;
        l2 = tip.prelength * pxcm ;
        d = (tip.area/tip.length)*pxcm;
        
        % calculate theta
        ip = find(tip.predecessor==[h.tips{N}.number]); % index predecessor
        if ip>0
            p = getParameters(h);
            p1 = tip.path; % node path of selected tip
            p2 = h.tips{N}(ip).path; % node path of predecessor tip
            ei = find(p2==p1(1),1,'first');
            c0 = h.nodes(p1(1),:);
            c1 = getPathPoint(h.A, h.Slen, h.nodes, h.edges, p1, p.Ra);
            c2 = getPathPoint(h.A, h.Slen, h.nodes, h.edges, p2(ei:end), p.Ra);
            v1 = c1-c0;
            v2 = c2-c0;
            v1 = v1/norm(v1);
            v2 = v2/norm(v2);
            theta = acos(v1*v2')/pi*180;
        else
            theta = nan;
        end
        % Cece
        tipSTATS = [tipSTATS,l];
        
        str = sprintf('Root #%i (#%i): \nLength: %3.4g cm\nDiameter: %3.4g cm\nPrelength: %3.4g cm\nOrder: %i\nApprox time: %2.2g days \nBranching angle: %3.4g', ...
            tip.number,tip.predecessor, l, d, l2, tip.order, tip.ct,theta);
        set(h.text3,'String',str);

       
        
    end
    
    
    % draw selection: NODES
    k = h.sn1;
    if k>0
        h.oldsel = plot(h.nodes(k,1),h.nodes(k,2),'ro',...
            h.nodes(k,1),h.nodes(k,2),'r*');
        str=['Node #', num2str(k)]; % update label
        set(h.text3,'String',str);
    end
    
    % update label
    str=sprintf('Graph: %i Nodes, %i Edges. Root system: %i roots', ...
        size(h.nodes,1),length(h.edges), length(h.tips{N}));
    set(h.text1,'String',str);
    % Cece Code Try to make an array of roots for stats
    rootArray = h.tips{N};
    assignin('base',"GraphRepresentation",rootArray);
    % Cece Code Try loop through tips and create a tip length array
    tipNumber = {};
    tipPredecessor ={};
    tipLen= {};
    tipDiameter = {};
    tipPrelength = {};
    tipOrder ={};
    tipTime ={};
    tipTheta ={};

    for i =1 : length(rootArray)
        tip = h.tips{N}(i);
        pxcm = str2double(get(h.edit1,'String')) / size(h.D{1},2);
        tipNumber = [tipNumber, tip.number];
        tipPredecessor = [tipPredecessor, tip.predecessor];
        l = tip.length * pxcm ;
        tipLen = [tipLen, l];
        d = (tip.area/tip.length)*pxcm;
        tipDiameter = [tipDiameter, d];
        l2 = tip.prelength * pxcm ;
        tipPrelength = [tipPrelength,l2];
        tipOrder = [tipOrder, tip.order];
        tipTime = [tipTime, tip.ct];
         % calculate theta
        ip = find(tip.predecessor==[h.tips{N}.number]); % index predecessor
        if ip>0
            p = getParameters(h);
            p1 = tip.path; % node path of selected tip
            p2 = h.tips{N}(ip).path; % node path of predecessor tip
            ei = find(p2==p1(1),1,'first');
            c0 = h.nodes(p1(1),:);
            c1 = getPathPoint(h.A, h.Slen, h.nodes, h.edges, p1, p.Ra);
            c2 = getPathPoint(h.A, h.Slen, h.nodes, h.edges, p2(ei:end), p.Ra);
            v1 = c1-c0;
            v2 = c2-c0;
            v1 = v1/norm(v1);
            v2 = v2/norm(v2);
            theta = acos(v1*v2')/pi*180;
        else
            theta = nan;
        end

        tipTheta = [tipTheta,theta]
    end
    % Create a structure array
    
    rootStats = struct("Number", tipNumber, "Predecessor", tipPredecessor, "Length", tipLen, "Prelength", tipPrelength,"Diameter",tipDiameter, "Order", tipOrder, "Time", tipTime, "Theta", tipTheta)
    assignin('base',"RootStats",rootStats);

    % % update root text
    %     tip = h.tips{N}(k);
    % 
    %     pxcm = str2double(get(h.edit1,'String')) / size(h.D{1},2);
    %     l = tip.length * pxcm ;
    %     l2 = tip.prelength * pxcm ;
    %     d = (tip.area/tip.length)*pxcm;
    % 
    %     % calculate theta
    %     ip = find(tip.predecessor==[h.tips{N}.number]); % index predecessor
    %     if ip>0
    %         p = getParameters(h);
    %         p1 = tip.path; % node path of selected tip
    %         p2 = h.tips{N}(ip).path; % node path of predecessor tip
    %         ei = find(p2==p1(1),1,'first');
    %         c0 = h.nodes(p1(1),:);
    %         c1 = getPathPoint(h.A, h.Slen, h.nodes, h.edges, p1, p.Ra);
    %         c2 = getPathPoint(h.A, h.Slen, h.nodes, h.edges, p2(ei:end), p.Ra);
    %         v1 = c1-c0;
    %         v2 = c2-c0;
    %         v1 = v1/norm(v1);
    %         v2 = v2/norm(v2);
    %         theta = acos(v1*v2')/pi*180;
    %     else
    %         theta = nan;
    %     end
    


    % update h
    guidata(hObject,h);
    
    
end


%
% save to RSML
%
function uipushtool3_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=guidata(hObject);

if isfield(h,'nodes')
    N = get(h.popupmenu4,'Value');  % which measurement    
    name = [ h.filename(1:end-4), '(', num2str(N), ')'];
  
    p = getParameters(h);
    writeRSML(h.filename(1:end-4),h.label(1:end-4),N,h.nFiles, h.A, h.edges, h.nodes, h.tips{N}, p, h.Slen, h.D{N});    
    msgbox(['File saved: ', name, '.rsml']);
end



