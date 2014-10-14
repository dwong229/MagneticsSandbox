function varargout = Fonaxis_2coilsGUI(varargin)
% FONAXIS_2COILSGUI MATLAB code for Fonaxis_2coilsGUI.fig
%      FONAXIS_2COILSGUI, by itself, creates a new FONAXIS_2COILSGUI or raises the existing
%      singleton*.
%
%      H = FONAXIS_2COILSGUI returns the handle to a new FONAXIS_2COILSGUI or the handle to
%      the existing singleton*.
%
%      FONAXIS_2COILSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FONAXIS_2COILSGUI.M with the given input arguments.
%
%      FONAXIS_2COILSGUI('Property','Value',...) creates a new FONAXIS_2COILSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Fonaxis_2coilsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Fonaxis_2coilsGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Fonaxis_2coilsGUI

% Last Modified by GUIDE v2.5 11-Sep-2014 14:35:02
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Fonaxis_2coilsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Fonaxis_2coilsGUI_OutputFcn, ...
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


% --- Executes just before Fonaxis_2coilsGUI is made visible.
function Fonaxis_2coilsGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Fonaxis_2coilsGUI (see VARARGIN)

%% 9/11/14 from Fonaxis_2coils.m
% COIL 1
handles.centerLeft = 0; % 
handles.radiusLeft = 0.01;
handles.currentLeft = 10; % current going in the CCW direction when axis of coil is pointing toward you

%fprintf('Coil located at x = %4.3fm with radius, r = %4.3fm and current, I = %4.1fA \n',centerLeft,radiusLeft,currentLeft)

% COIL 2
handles.centerRight = 0.5;
handles.radiusRight = 0.01;
handles.currentRight = 10;

%fprintf('Coil located at x = %4.3fm with radius, r = %4.3fm and current, I = %4.1fA \n',centerRight,radiusRight,currentRight)

% magnet properties
v = 1;
handles.Mx = 10*1;
%fprintf('Magnet properties, Magnetization, m_x = %4.1fNm/T \n',v*Mx)
% workspace
handles.dx = 0.001;
handles.x = handles.centerLeft:handles.dx:handles.centerRight;

% Define magnet position
handles.magnets = [.1;.4]; % x position for each magnet in a vector
magnetIdx = (handles.magnets-handles.centerLeft)/handles.dx + 1;
%% compute B
gradBLeft = computegradBonaxis(handles.x,handles.currentLeft,handles.radiusLeft,handles.centerLeft);
gradBRight = computegradBonaxis(handles.x,handles.currentRight,handles.radiusRight,handles.centerRight);

% compute F
FLeft = computeFelectromagnet(gradBLeft,handles.Mx);
FRight = computeFelectromagnet(gradBRight,handles.Mx);
FTotal = FLeft + FRight;
%% end from Fonaxis_2coils

%h_Force(1) = 
axes(handles.axes1)
title('Force profile')
handles.FLeft = plot(handles.x,FLeft,'b');
hold on
handles.FRight = plot(handles.x,FRight,'r');
handles.FTotal = plot(handles.x,FTotal,'.k');
xlabel('Position')
ylabel('Force')
legend('Left','Right','Total')

%% subplot(2,1,2)
% print zero crossing
zeroCrossExists = checkZeroCross(FTotal);

%% Compute force vector at magnet location
axes(handles.axes3)
title('Coil and Magnet Geometry')

% draw ellipse for coils
hold on
handles.h_coil(1) = rectangle('Position',[handles.centerLeft-0.005,-handles.radiusLeft,0.01,2*handles.radiusLeft],'Curvature',[1,1]);
handles.h_coil(2) = rectangle('Position',[handles.centerRight-0.005,-handles.radiusRight,0.01,2*handles.radiusRight],'Curvature',[1,1]);
axis([-0.005 max(handles.x)+0.005 -handles.radiusLeft handles.radiusLeft]);

% draw circles for magnet
h_magnetposn = plot(handles.magnets,zeros(1,length(handles.magnets)),'ok','MarkerSize',15,'MarkerFaceColor','k');
% plot direction of current 

% draw Force arrows
magnetIdx = (handles.magnets-handles.centerLeft)/handles.dx + 1;
magnetForce = FTotal(magnetIdx);
qscale = 0.5;
handles.h_magnetforce = quiver(handles.magnets,zeros(1,length(handles.magnets)),magnetForce,zeros(1,length(handles.magnets)),qscale);
set(handles.edit3,'String',sprintf('%0.3g',magnetForce(1)));
set(handles.edit4,'String',sprintf('%0.3g',magnetForce(2)));

% print and plot zerocross
if zeroCrossExists
    [~,zeroCrossIdx] = min(FTotal(1:end-1).^2);
    zeroCross = handles.x(zeroCrossIdx);
    fprintf('Force balanced at x = %4.4f \n',zeroCross);
    handles.zerocross = plot(zeroCross,0,'xr','MarkerSize',10);
else 
    disp('No zero cross')
    handles.zerocross = plot(0,0,'xr','MarkerSize',10,'Visible','off');
end


% Choose default command line output for Fonaxis_2coilsGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Fonaxis_2coilsGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Fonaxis_2coilsGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderVal = get(hObject,'Value');
% change value of text box
set(handles.edit1,'String',sliderVal);
handles.currentLeft = sliderVal;
% update force plot
% compute B
gradBLeft = computegradBonaxis(handles.x,handles.currentLeft,handles.radiusLeft,handles.centerLeft);
% compute F
FLeft = computeFelectromagnet(gradBLeft,handles.Mx);
FTotal = FLeft + get(handles.FRight,'YData');
set(handles.FLeft,'YData',FLeft);
set(handles.FTotal,'YData',FTotal);


% update magnet plot
zeroCrossExists = checkZeroCross(FTotal);
% draw Force arrows
magnetIdx = (handles.magnets-handles.centerLeft)/handles.dx + 1;
magnetForce = FTotal(magnetIdx);
qscale = 0.5;
set(handles.h_magnetforce,'UData',magnetForce);
set(handles.edit3,'String',sprintf('%0.3g',magnetForce(1)));
set(handles.edit4,'String',sprintf('%0.3g',magnetForce(2)));


% print and plot zerocross
if zeroCrossExists
    [~,zeroCrossIdx] = min(FTotal(1:end-1).^2);
    zeroCross = handles.x(zeroCrossIdx);
    fprintf('Force balanced at x = %4.4f \n',zeroCross);
    set(handles.zerocross,'XData',zeroCross,'Visible','on');
else 
    disp('No zero cross')
    set(handles.zerocross,'Visible','off');
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderVal = get(hObject,'Value');
set(handles.edit2,'String',sliderVal);
handles.currentRight = sliderVal;
% update force plot
% compute B
gradBRight = computegradBonaxis(handles.x,handles.currentRight,handles.radiusRight,handles.centerRight);
% compute F
FRight = computeFelectromagnet(gradBRight,handles.Mx);
FTotal = FRight + get(handles.FLeft,'YData');
set(handles.FRight,'YData',FRight);
set(handles.FTotal,'YData',FTotal);

% update magnet plot
zeroCrossExists = checkZeroCross(FTotal);
% draw Force arrows
magnetIdx = (handles.magnets-handles.centerLeft)/handles.dx + 1;
magnetForce = FTotal(magnetIdx);
qscale = 0.5;
set(handles.h_magnetforce,'UData',magnetForce);
set(handles.edit3,'String',sprintf('%0.3g',magnetForce(1)));
set(handles.edit4,'String',sprintf('%0.3g',magnetForce(2)));

% print and plot zerocross
if zeroCrossExists
    [~,zeroCrossIdx] = min(FTotal(1:end-1).^2);
    zeroCross = handles.x(zeroCrossIdx);
    fprintf('Force balanced at x = %4.4f \n',zeroCross);
    set(handles.zerocross,'XData',zeroCross,'Visible','on');
else 
    disp('No zero cross')
    set(handles.zerocross,'Visible','off');
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
