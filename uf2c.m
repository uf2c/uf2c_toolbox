function varargout = uf2c(varargin)
% UF²C M-file for uf2c.fig
% UF²C - User Friendly Functional Connectivity
% Brunno Machado de Campos
% University of Campinas, 2019
%
% Copyright (c) 2019, Brunno Machado de Campos
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
                   'gui_OpeningFcn', @uf2c_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_OutputFcn, ...
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

function uf2c_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
a = fopen('uf2c_defaults.m');
a1 = fread(a);
b = fopen('uf2c_RDRF.m');
b1 = fread(b);
if isequal(a1,b1)
    set(handles.REBT,'Visible','off')
    set(handles.OnDefTxt,'Visible','on')
else
    set(handles.REBT,'Visible','on')
    set(handles.OnDefTxt,'Visible','off')
end

function varargout = uf2c_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function background_CreateFcn(hObject, eventdata, handles)
clc
disp('                                            _____                     ')
disp('     ____          ____     _____________  |___  |   ________________ ')
disp('    |    |        |    |   |             |  ___| |  |                |')
disp('    |    |        |    |   |     ________| |  ___|  |     ___________|')
disp('    |    |        |    |   |    |          | |___   |    |')
disp('    |    |        |    |   |    |          |_____|  |    |')
disp('    |    |        |    |   |    |______             |    |')
disp('    |    |        |    |   |           |            |    |')
disp('    |    |        |    |   |     ______|            |    |')      
disp('    |    |        |    |   |    |                   |    |')
disp('    |    |        |    |   |    |                   |    |')
disp('    |    |        |    |   |    |                   |    |')
disp('    |    |________|    |   |    |                   |    |____________ ')
disp('    |                  |   |    |                   |                 |')
disp('    |__________________|   |____|                   |_________________|')
disp(' ')
disp('                   USER FRIENDLY FUNCTIONAL CONNECTIVITY')
disp('                          (For Matlab with SPM12)')

uf2cdir = fileparts(which('uf2c'));
img = imread(fullfile(uf2cdir,'fIGcAPA.png'));
imagesc(img)
set(gca,'Visible','off')


delete([uf2cdir filesep 'Analysis' filesep 'FC_tmp' filesep '*.*'])

function figure1_CreateFcn(hObject, eventdata, handles)

function axes3_CreateFcn(hObject, eventdata, handles)

function Preprop_Callback(hObject, eventdata, handles)

function CrosCorre_Callback(hObject, eventdata, handles)

function analy_Callback(hObject, eventdata, handles)

function Tools_Callback(hObject, eventdata, handles)

function imEd_Callback(hObject, eventdata, handles)
uf2c_ErrorCorrection

function FilChanc_Callback(hObject, eventdata, handles)
uf2c_RenameFiles

function funcon_Callback(hObject, eventdata, handles)
uf2c_FuncCon

function Tseg_Callback(hObject, eventdata, handles)
uf2c_Segment

function Volum_Callback(hObject, eventdata, handles)
uf2c_Volumetry

function abto_Callback(hObject, eventdata, handles)

function about_Callback(hObject, eventdata, handles)
open uf2c_About

function FI_Callback(hObject, eventdata, handles)
uf2c_FuncInte

function SWFuncCn_Callback(hObject, eventdata, handles)
uf2c_FuncConSW

function interp2tmp_Callback(hObject, eventdata, handles)
uf2c_InterpGUI

function CrossCorr_Callback(hObject, eventdata, handles)
uf2c_TS_CrossCor

function CrossCor_Second_Callback(hObject, eventdata, handles)
uf2c_CrossCorSecLe

function Credits_Callback(hObject, eventdata, handles)
uf2c_Credits

function MovView_Callback(hObject, eventdata, handles)
uf2c_HMV

function MovAna_Callback(hObject, eventdata, handles)
uf2c_HMA

function r2z_Callback(hObject, eventdata, handles)
uf2c_r2z

function AnonyDCM_Callback(hObject, eventdata, handles)
uf2c_Anonymize_DCM_Header

function StdPreproc_Callback(hObject, eventdata, handles)
uf2c_JustPreproc

function GenericPreproc_Callback(hObject, eventdata, handles)
uf2c_GeneralPreProc

function uf2c2nbs_Callback(hObject, eventdata, handles)
uf2c_to_NBS

function compat_Callback(hObject, eventdata, handles)
uf2c_ChCompat

function BDA_Callback(hObject, eventdata, handles)
uf2c_FirstLevel_BD

function Untitled_1_Callback(hObject, eventdata, handles)

function Untitled_4_Callback(hObject, eventdata, handles)


function SecLevelCorr_Callback(hObject, eventdata, handles)
uf2c_CrossCorSecLe_CORRE

function MapsInter_Callback(hObject, eventdata, handles)
uf2c_Maps_Intersection

function MapQuanti_Callback(hObject, eventdata, handles)
uf2c_Maps_Quantific

function AbortedfMRI_Callback(hObject, eventdata, handles)
uf2c_AbortedHeader

function ConvPreproc_Callback(hObject, eventdata, handles)
uf2c_JustPreproc

function RMC_Preproc_Callback(hObject, eventdata, handles)
uf2c_Preproc_RMC

function Outlier_Callback(hObject, eventdata, handles)
uf2c_OutliersTool

function eegfmri_Callback(hObject, eventdata, handles)
uf2c_EEG_fMRI_Ana

function WokDir_Callback(hObject, eventdata, handles)

function DefWorkDir_Callback(hObject, eventdata, handles)
WDir = uigetdir('','Define a work Directory');
cd(WDir);

function REBT_Callback(hObject, eventdata, handles)
uf2cdir = fileparts(which('uf2c'));
copyfile([uf2cdir,filesep,'Utilities',filesep,'uf2c_RDRF.m'],...
    [uf2cdir,filesep,'Utilities',filesep,'uf2c_defaults.m'],'f')
set(handles.REBT,'Visible','off')
set(handles.OnDefTxt,'Visible','on')
fprintf('UF²C Default Parameters redefined to originals settings\r\n')

function OpnDeff_Callback(hObject, eventdata, handles)
open uf2c_defaults
