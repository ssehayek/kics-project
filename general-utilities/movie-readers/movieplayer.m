% Created by H.B.
% last edited 2012/02/24
% usage: movieplayer(imageSeries)
% modified to use imagesc not imshow

function varargout = movieplayer(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @movieplayer_OpeningFcn, ...
                   'gui_OutputFcn',  @movieplayer_OutputFcn, ...
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


% --- Executes just before movieplayer is made visible.
function movieplayer_OpeningFcn(hObject, eventdata, handles, varargin)
global imageSeriesDiff;
imageSeriesDiff = varargin{1};
%load('imageSeriesDiff.mat');
max = size(imageSeriesDiff,3); min = 1;
set(handles.figure1,'resize','on');
set(handles.axes1,'units','normalized');
set(handles.slider1,'Max',[max],'Min',[min],'value',1,'units','normalized');
set(handles.text1,'string',num2str(1),'units','normalized');
minstep = 1/size(imageSeriesDiff,3);
set(handles.slider1,'sliderstep',[minstep,0.1]);
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to movieplayer (see VARARGIN)

% Choose default command line output for movieplayer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes movieplayer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = movieplayer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
global imageSeriesDiff;
f = round(get(handles.slider1,'value'));
set(handles.text1,'string',num2str(f));
colormap(pink)
%imshow(imageSeriesDiff(:,:,f));
imagesc(imageSeriesDiff(:,:,f));

%%%%%%%%%%%%%% for ROI visualization
% numY = 16;
% numX = 16;
% stepsize = 32;
% for j=1:numY
%     for k=1:numX 
%         minY = (j-1)*stepsize + 1;
%         minX = (k-1)*stepsize + 1;
%         rectangle('Position',[minY,minX,stepsize,stepsize],'linewidth',1);
%     end
% end
%%%%%%%%%%%%%%%


% colormap(jet)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end






