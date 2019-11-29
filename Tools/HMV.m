function varargout = HMV(varargin)
% UF²C M-file for HMV.fig
% UF²C - User Friendly Functional Connectivity
% Brunno Machado de Campos
% University of Campinas, 2015
%
% Copyright (c) 2015, Brunno Machado de Campos
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
                   'gui_OpeningFcn', @HMV_OpeningFcn, ...
                   'gui_OutputFcn',  @HMV_OutputFcn, ...
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

function HMV_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = HMV_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function Addtxt_Callback(hObject, eventdata, handles)
global filetxt pathtxt
cla(handles.XYZg,'reset')
cla(handles.RPYg,'reset')

[filetxt,pathtxt] = uigetfile({'rp_*.txt','TEXT file'},...
    'Add a Realigment re_*.txt file','MultiSelect','off');

vetore = importdata([pathtxt filesep filetxt]);
vetore(:,4:6) = (180/pi).*vetore(:,4:6);

set(handles.addstr,'String',[pathtxt filesep filetxt])
set(plot(vetore(:,1)),'Parent', handles.XYZg)

hold on
set(plot(vetore(:,2)),'Parent', handles.XYZg)

hold on
set(plot(vetore(:,3)),'Parent', handles.XYZg)

title(handles.XYZg,'Translation movement')
legend(handles.XYZg,'X','Y','Z','Location','southwest')
ylabel(handles.XYZg,'mm')
xlabel(handles.XYZg,'volume')


set(plot(vetore(:,4)),'Parent', handles.RPYg)
hold on
set(plot(vetore(:,5)),'Parent', handles.RPYg)
hold on
set(plot(vetore(:,6)),'Parent', handles.RPYg)
title(handles.RPYg,'Rotational movement')
legend(handles.RPYg,'Row','Pitch','Yaw','Location','southwest')
ylabel(handles.RPYg,'degrees')
xlabel(handles.RPYg,'volume')

MaxX = max(abs(vetore(:,1)));
MaxY = max(abs(vetore(:,2)));
MaxZ = max(abs(vetore(:,3)));

MaxR = max(abs(vetore(:,4)));
MaxP = max(abs(vetore(:,5)));
MaxYa = max(abs(vetore(:,6)));

MaxT = max([MaxX,MaxY,MaxZ]);
MaxRot = max([MaxR,MaxP,MaxYa]);

switch  MaxT
    case MaxX
        strHd = 'X';
        PosiDH = find(abs(vetore(:,1))==MaxX);
    case MaxY
        strHd = 'Y';
        PosiDH = find(abs(vetore(:,2))==MaxY);
    case MaxZ
        strHd = 'Z';
        PosiDH = find(abs(vetore(:,3))==MaxZ);
end
set(handles.HDp,'String',PosiDH)
MaxT = num2str(MaxT);
HDvV = [MaxT(1:5) ' on '   strHd];
set(handles.HDv,'String',HDvV)

switch  MaxRot
    case MaxR
        strHt = 'Row';
        PosiTH = find(abs(vetore(:,4))==MaxR);
    case MaxP
        strHt = 'Pitch';
        PosiTH = find(abs(vetore(:,5))==MaxP);
    case MaxYa
        strHt = 'Yaw';
        PosiTH = find(abs(vetore(:,6))==MaxYa);
end
set(handles.HTp,'String',PosiTH)
MaxRot = num2str(MaxRot);
HTvV = [MaxRot(1:5) ' on '   strHt];
set(handles.HTv,'String',HTvV)


function LTp_Callback(hObject, eventdata, handles)

function LTp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LTv_Callback(hObject, eventdata, handles)

function LTv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function HTp_Callback(hObject, eventdata, handles)

function HTp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function HTv_Callback(hObject, eventdata, handles)

function HTv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LDv_Callback(hObject, eventdata, handles)

function LDv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LDp_Callback(hObject, eventdata, handles)
function LDp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function HDv_Callback(hObject, eventdata, handles)
function HDv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function HDp_Callback(hObject, eventdata, handles)
function HDp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function uipushtool1_ClickedCallback(hObject, eventdata, handles)

function pushbutton2_Callback(hObject, eventdata, handles)
global filetxt pathtxt
[Mot_fig,HDvV,HTvV] = uf2c_plot_motion([pathtxt filetxt],'off');
imgRR = getframe(Mot_fig);
imwrite(imgRR.cdata, [pathtxt,filetxt(4:6),'_Realignment_Parameters_Plot.png']);
saveas(Mot_fig,[pathtxt,filetxt(4:6),'_Realignment_Parameters_Plot'],'fig');
close(Mot_fig)
