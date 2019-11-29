function varargout = ChCompat(varargin)
% UF²C - User Friendly Functional Connectivity
% Brunno Machado de Campos
% University of Campinas, 2015
%
% Copyright (c) 2016, Brunno Machado de Campos
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ChCompat_OpeningFcn, ...
                   'gui_OutputFcn',  @ChCompat_OutputFcn, ...
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


function ChCompat_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

guidata(hObject, handles);

function varargout = ChCompat_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function chbtm_Callback(hObject, eventdata, handles)

if ispc
    set(handles.text25,'String','Current: Windows')
end

if ismac
    set(handles.text25,'String','Current: Mac OS')
end

if isunix
    set(handles.text25,'String','Current: Unix')
end

Tester = ver;
SPMex = 0;
IPTex = 0;
SPTex = 0;
STex = 0;

for i = 1:size(Tester,2)
    if isequal(Tester(i).Name,'MATLAB')
        set(handles.text10,'String',['Current: ' Tester(i).Version ' ' Tester(i).Release])
        if str2num(Tester(i).Version)>=7.10
            set(handles.text15,'Visible','on')
            set(handles.text20,'Visible','off')
        else
            set(handles.text15,'Visible','off')
            set(handles.text20,'Visible','on')
        end
    end
    if isequal(Tester(i).Name,'Statistical Parametric Mapping')
        SPMex = 1;
        set(handles.text11,'String',['Current: '  Tester(i).Release(2:end-1)])
        if isequal(Tester(i).Release,'(SPM12)')
            set(handles.text16,'Visible','on')
            set(handles.text21,'Visible','off')
        else
            set(handles.text16,'Visible','off')
            set(handles.text21,'Visible','on')
        end
    end

    if isequal(Tester(i).Name,'Image Processing Toolbox')
        IPTex = 1;
    end
    if isequal(Tester(i).Name,'Signal Processing Toolbox')
        SPTex = 1;
    end
    if isequal(Tester(i).Name,'Statistics and Machine Learning Toolbox')
        STex = 1;
    end
    if isequal(Tester(i).Name,'Statistics Toolbox')
        STex = 1;
    end
end

if ~SPMex
    set(handles.text26,'Visible','on')
end

if IPTex
    set(handles.text29,'Visible','on')
    set(handles.text24,'Visible','off')
else
    set(handles.text24,'Visible','on')
end

if SPTex
    set(handles.text28,'Visible','on')
    set(handles.text23,'Visible','off')
else
    set(handles.text23,'Visible','on')
end

if STex
    set(handles.text27,'Visible','on')
    set(handles.text22,'Visible','off')
else
    set(handles.text22,'Visible','on')
end

wc = which('uf2c','-all');
if length(wc) > 1
    set(handles.text30,'Visible','on')
else
    set(handles.text30,'ForegroundColor',[0 0.49 0])
    set(handles.text30,'String','No problems here   :-)')
    set(handles.text30,'Visible','on')
end

uf2cdir = fileparts(which('uf2c'));
[status,values] = fileattrib(uf2cdir);

if isequal(values.UserWrite,0)
    set(handles.text31,'Visible','on')
else
    set(handles.text31,'ForegroundColor',[0 0.49 0])
    set(handles.text31,'String','No problems here   :-)')
    set(handles.text31,'Visible','on')
end
    




    






   
