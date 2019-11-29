function varargout = checkcorres(varargin)
% UF²C - User Friendly Functional Connectivity
% Brunno Machado de Campos
% University of Campinas, 2017
%
% Copyright (c) 2017, Raphael Fernandes Casseb
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
%

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @checkcorres_OpeningFcn, ...
                   'gui_OutputFcn',  @checkcorres_OutputFcn, ...
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


function checkcorres_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = checkcorres_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function uitable1_CreateFcn(hObject, eventdata, handles)
global Image4D
set(hObject,'Data', []);

set(hObject,'Data', cell(size(Image4D,1),1));
set(hObject,'Data',Image4D')

function uitable2_CreateFcn(hObject, eventdata, handles)
global Image4D MarkerF

set(hObject,'Data', []);

if size(Image4D,2)>size(MarkerF,2)
    diff = size(Image4D,2)-size(MarkerF,2);
    gg = cell(1,diff);
    gg(1,:) = deal({'Empty'});
    MarkerF = [MarkerF gg];
end

set(hObject,'Data', cell(size(Image4D,2),2));
set(hObject,'Data',[MarkerF',num2cell(1:size(MarkerF,2))'])

function uitable3_CreateFcn(hObject, eventdata, handles)
global Image4D rpF
set(hObject,'Data', []);

set(hObject,'Data', cell(size(Image4D,1),2));
set(hObject,'Data',[rpF',num2cell(1:size(rpF,2))'])

function save_Callback(hObject, eventdata, handles)
global Image4D rpF MarkerF

mkd = get(handles.uitable2,'Data');
if ~isempty(mkd)
    mkd = mkd(:,2);
    tester1 = unique(cell2mat(mkd));
    if  size(tester1,1) < size(Image4D,2)
        warndlg('You probably attributed repetitive values into the Markers order','Attention')
        return
    end
    for y = 1:size(MarkerF,2)
        MarkerF2(y) = MarkerF(strfind(cell2mat(mkd)',y));
    end
     MarkerF = MarkerF2;
end

regd = get(handles.uitable3,'Data');
if ~isempty(regd)
    regd = regd(:,2);
    tester2 = unique(cell2mat(regd));
    if size(tester2,1) < size(Image4D,2)
        warndlg('You probably attributed repetitive values into the Regressors order','Attention')
        return
    end
    for y = 1:size(rpF,2)
        rpF2(y) = rpF(strfind(cell2mat(regd)',y));
    end
    rpF = rpF2;
end

close(checkcorres)
