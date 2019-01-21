% prepData 

function varargout = prepData_UI(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prepData_UI_OpeningFcn, ...
                   'gui_OutputFcn',  @prepData_UI_OutputFcn, ...
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

function prepData_UI_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
guidata(hObject, handles);

function varargout = prepData_UI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function checkbox1_Callback(hObject, eventdata, handles)
handles.properties.collapse = get(hObject,'Value');
guidata(hObject,handles)

function edit7_Callback(hObject, eventdata, handles)
handles.properties.minGroupSize = str2double(get(hObject,'String'));
guidata(hObject,handles)

function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, handles)
handles.properties.collectionName = get(hObject,'String');
guidata(hObject,handles)

function edit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox2_Callback(hObject, eventdata, handles)
handles.properties.saveSet = get(hObject, 'Value');
guidata(hObject,handles)

function edit1_Callback(hObject, eventdata, handles)
file_name = get(hObject,'String');
handles.properties.file = file_name{1};
guidata(hObject,handles)

function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit4_Callback(hObject, eventdata, handles)
handles.properties.sheet = get(hObject,'String');
guidata(hObject,handles)

function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton1_Callback(hObject, eventdata, handles)
try;
    [output] = prepData_GUI_170917(handles.properties);
    assignin('base', 'output', output);
catch e;
    assignin('base','e',e)
    disp(e)
    msgbox({'Error while running prepare data module',...
        '                                                              ',...
        'Ensure inputs are complete'},'Error', 'Error');
end




% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)

file_name = uigetfile('*.xlsx','Select excel file for preparation');
try; handles.properties; catch; handles.properties = []; end;
file_name = erase(file_name, '.xlsx');
set(handles.edit1, 'String', file_name);
handles.properties.file = file_name;
guidata(hObject,handles);
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
close all; clear all; 
msgbox({'Prepare data module has been closed'}, 'Cancelled');
error('Prepare data module has been closed');
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
