function varargout = ManualCurateGUI(varargin)
% MANUALCURATEGUI MATLAB code for ManualCurateGUI.fig
%      MANUALCURATEGUI, by itself, creates a new MANUALCURATEGUI or raises the existing
%      singleton*.
%
%      H = MANUALCURATEGUI returns the handle to a new MANUALCURATEGUI or the handle to
%      the existing singleton*.
%
%      MANUALCURATEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANUALCURATEGUI.M with the given input arguments.
%
%      MANUALCURATEGUI('Property','Value',...) creates a new MANUALCURATEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ManualCurateGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ManualCurateGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ManualCurateGUI

% Last Modified by GUIDE v2.5 23-Sep-2022 19:28:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ManualCurateGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ManualCurateGUI_OutputFcn, ...
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


% --- Executes just before ManualCurateGUI is made visible.
function ManualCurateGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ManualCurateGUI (see VARARGIN)

% Choose default command line output for ManualCurateGUI
handles.output = hObject;

initialData = NaN(600,1);
set(handles.uitable1,'data',initialData)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ManualCurateGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ManualCurateGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- LOAD SEGMENTATION DATA: Executes on button press
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile('*.mat','Select Calcium .mat File');
load([path file])
handles.C = C;
handles.ROI = ROI;
handles.DeltaFoverF = DeltaFoverF;

%show avg data image
imagesc(handles.axes1,handles.C);hold on
cl = caxis;
set(handles.axes1,'YTickLabel','','XTickLabel','','YTick','','XTick','')


%set the slider value to 1 and set the max/min values depending on numbre
%of ROIs
sliderMin = 1;
sliderMax = length(handles.ROI);
set(handles.slider1,'SliderStep',[1, 1] / (sliderMax - sliderMin),'max',length(handles.ROI),'min',1,'Value',1)
x=get(handles.slider1,'Value');

%set the current component string textbox to 1
set(handles.edit1,'string',num2str(x));

%plot outline of current component
for ii = 1:length(handles.ROI{x})
    plot(handles.axes1,handles.ROI{x}{ii}(1),handles.ROI{x}{ii}(2),'r.')
end

x1=handles.ROI{x}{1}(1);
y1=handles.ROI{x}{1}(2);

if(x1>50 && x1<462)
    xleft=x1-50;
    xright=x1+50;
elseif(x1<50)
    xleft=1;
    xright=x1+50;
elseif(x1>462)
    xleft=x1-50;
    xright=512;
end
if(y1>50 && y1<462)
    ytop=y1-50;
    ybottom=y1+50;
elseif(y1<50)
    ytop=1;
    ybottom=y1+50;
elseif(y1>462)
    ytop=y1-50;
    ybottom=512;
end
hold(handles.axes3,'off')
imagesc(handles.axes3,handles.C(ytop:ybottom,xleft:xright))
set(handles.axes3,'YTickLabel','','XTickLabel','','YTick','','XTick','')
caxis(cl)
hold(handles.axes3,'on')
for ii = 1:length(handles.ROI{x})
    plot(handles.axes3,handles.ROI{x}{ii}(1)-(xleft),handles.ROI{x}{ii}(2)-(ytop),'r.')
end
% Plots fluorescent values
plot(handles.axes4,handles.DeltaFoverF(x,:));

%if the variable badComponents exists, update these Bad Components table
if(exist('badComponents','var'))
    initialData = NaN(600,1);
    initialData(1:length(badComponents)) = badComponents;
    set(handles.uitable1,'data',initialData); 
end

%update handles
guidata(hObject, handles);


% --- SLIDER: Executes on slider movement
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%get new slider value
x=get(handles.slider1,'Value');

%update the text box to new slider value
set(handles.edit1,'string',num2str(x));

%replot the avg image and plot the ROI of new component
hold(handles.axes1,'off')
imagesc(handles.axes1,handles.C);hold on
cl=caxis;
set(handles.axes1,'YTickLabel','','XTickLabel','','YTick','','XTick','')
for ii = 1:length(handles.ROI{x})
    plot(handles.axes1,handles.ROI{x}{ii}(1),handles.ROI{x}{ii}(2),'r.')
end

% Plots fluorescent values
plot(handles.axes4,handles.DeltaFoverF(x,:));

x1=handles.ROI{x}{1}(1);
y1=handles.ROI{x}{1}(2);

if(x1>50 && x1<462)
    xleft=x1-50;
    xright=x1+50;
elseif(x1<50)
    xleft=1;
    xright=x1+50;
elseif(x1>462)
    xleft=x1-50;
    xright=512;
end
if(y1>50 && y1<462)
    ytop=y1-50;
    ybottom=y1+50;
elseif(y1<50)
    ytop=1;
    ybottom=y1+50;
elseif(y1>462)
    ytop=y1-50;
    ybottom=512;
end
hold(handles.axes3,'off')
imagesc(handles.axes3,handles.C(ytop:ybottom,xleft:xright))
set(handles.axes3,'YTickLabel','','XTickLabel','','YTick','','XTick','')
caxis(cl)
hold(handles.axes3,'on')
for ii = 1:length(handles.ROI{x})
    plot(handles.axes3,handles.ROI{x}{ii}(1)-(xleft),handles.ROI{x}{ii}(2)-(ytop),'r.')
end


guidata(hObject, handles);


% --- MARK AS BAD PUSHBUTTON
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get the current component selection
currentComponent = str2num(get(handles.edit1,'string'));
%get the current data in the table
tableData=get(handles.uitable1,'data');
%find all NaN
nans = find(isnan(tableData));
%set the first NaN to the current component selection
tableData(nans(1)) = currentComponent;
%reset the table data
set(handles.uitable1,'data',tableData);
guidata(hObject, handles);
% Increment slider
%set(handles.slider1,'Value',currentComponent+1);


% --- SEND TO WORKSPACE BUTTON: Executes on button press.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get the table data with all bad components
tableData=get(handles.uitable1,'data');
%find all NaN
nans = find(isnan(tableData));
%set all NaNs to empty
tableData(nans) = [];
%send table data to workspace
assignin('base','badComponents',tableData)


% --- CURRENT COMPONENT #: Executes during number input
function edit1_Callback(hObject,eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

%get the input value
x=str2num(get(handles.edit1,'string'));
set(handles.slider1,'Value',x);
imagesc(handles.axes1,handles.C);hold on
set(handles.axes1,'YTickLabel','','XTickLabel','','YTick','','XTick','')
for ii = 1:length(handles.ROI{x})
    plot(handles.axes1,handles.ROI{x}{ii}(1),handles.ROI{x}{ii}(2),'r.')
end

x1=handles.ROI{x}{1}(1);
y1=handles.ROI{x}{1}(2);

if(x1>50 && x1<462)
    xleft=x1-50;
    xright=x1+50;
elseif(x1<50)
    xleft=1;
    xright=x1+50;
elseif(x1>462)
    xleft=x1-50;
    xright=512;
end
if(y1>50 && y1<462)
    ytop=y1-50;
    ybottom=y1+50;
elseif(y1<50)
    ytop=1;
    ybottom=y1+50;
elseif(y1>462)
    ytop=y1-50;
    ybottom=512;
end
hold(handles.axes3,'off')
imagesc(handles.axes3,handles.C(ytop:ybottom,xleft:xright))
set(handles.axes3,'YTickLabel','','XTickLabel','','YTick','','XTick','')
hold(handles.axes3,'on')
for ii = 1:length(handles.ROI{x})
    plot(handles.axes3,handles.ROI{x}{ii}(1)-(xleft),handles.ROI{x}{ii}(2)-(ytop),'r.')
end


guidata(hObject, handles);



% --- CURRENT COMPONENT #: executes during creation
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes4
