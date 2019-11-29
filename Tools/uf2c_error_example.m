function varargout = uf2c_error_example(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @uf2c_error_example_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_error_example_OutputFcn, ...
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



function uf2c_error_example_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

guidata(hObject, handles);

function varargout = uf2c_error_example_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function axes1_CreateFcn(hObject, eventdata, handles)
img = imread(fullfile(fileparts(which('error_example')),'reflex_error.png'));
imagesc(img)
set(gca,'Visible','off')
