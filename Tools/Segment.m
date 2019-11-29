function varargout = Segment(varargin)
% SEGMENT MATLAB code for Segment.fig
%      SEGMENT, by itself, creates a new SEGMENT or raises the existing
%      singleton*.
%
%      H = SEGMENT returns the handle to a new SEGMENT or the handle to
%      the existing singleton*.
%
%      SEGMENT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEGMENT.M with the given input arguments.
%
%      SEGMENT('Property','Value',...) creates a new SEGMENT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Segment_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Segment_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Segment

% Last Modified by GUIDE v2.5 28-Sep-2011 13:59:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Segment_OpeningFcn, ...
                   'gui_OutputFcn',  @Segment_OutputFcn, ...
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


% --- Executes just before Segment is made visible.
function Segment_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

guidata(hObject, handles);

function varargout = Segment_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function addfiles_Callback(hObject, eventdata, handles)
global fileseg pathseg nfi

[fileseg,pathseg] = uigetfile({'*.nii','*.nii (NiFTI)'},'Select all files', 'MultiSelect', 'on');
nf = cellstr(fileseg);
nf = size(nf);
nfi = nf(1,2);

if nfi>1
set(handles.textfiles, 'String', 'files added')
else
    set(handles.textfiles, 'String', 'file added')
end

set(handles.nfiles, 'String', sprintf('%d',nfi))



function x_Callback(hObject, eventdata, handles)

function x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_Callback(hObject, eventdata, handles)

function y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function z_Callback(hObject, eventdata, handles)

function z_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
global fileseg pathseg nfi


aux3 = str2double(get(handles.x, 'String'));
aux1 = str2double(get(handles.y, 'String'));
aux2 = str2double(get(handles.z, 'String'));
aux = [aux3 aux1 aux2];

if nfi>1
fileseg2 = cell(1,nfi);
for i = 1:nfi
    fileseg{1,i} = sprintf('%s,1',fileseg{1,i});
fileseg2{1,i} = fullfile(pathseg,fileseg{1,i});
end

fileseg2=transpose(fileseg2);
else 
    fileseg = sprintf('%s,1',fileseg);
    fileseg2 = fullfile(pathseg,fileseg);
end

for i=1:nfi

    if nfi>1
matlabbatch{1}.spm.spatial.preproc.data = fileseg2(i,1);
    else
        matlabbatch{1}.spm.spatial.preproc.data = {fileseg2};
    end

matlabbatch{1}.spm.spatial.preproc.output.GM = [1 1 1];
matlabbatch{1}.spm.spatial.preproc.output.WM = [1 1 1];
matlabbatch{1}.spm.spatial.preproc.output.CSF = [1 1 1];
matlabbatch{1}.spm.spatial.preproc.output.biascor = 1;
matlabbatch{1}.spm.spatial.preproc.output.cleanup = 0;
% [pathw,filew,extw] = fileparts(which('spm'));
% matlabbatch{1}.spm.spatial.preproc.opts.tpm = {fullfile(pathw,'tpm','grey.nii'), fullfile(pathw,'tpm','white.nii'),fullfile(pathw,'tpm','csf.nii')};
matlabbatch{1}.spm.spatial.preproc.opts.ngaus = [2;2;2;4];
matlabbatch{1}.spm.spatial.preproc.opts.regtype = 'mni';
matlabbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{1}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.opts.samp = 3;
matlabbatch{1}.spm.spatial.preproc.opts.msk = {''};

if get(handles.NormY, 'Value')==1
    
matlabbatch{2}.spm.spatial.normalise.write.subj.matname(1) = cfg_dep;
matlabbatch{2}.spm.spatial.normalise.write.subj.matname(1).tname = 'Parameter File';
matlabbatch{2}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{2}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{2}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{2}.spm.spatial.normalise.write.subj.matname(1).tgt_spec{1}(2).value = 'e';
matlabbatch{2}.spm.spatial.normalise.write.subj.matname(1).sname = 'Segment: Norm Params Subj->MNI';
matlabbatch{2}.spm.spatial.normalise.write.subj.matname(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.spm.spatial.normalise.write.subj.matname(1).src_output = substruct('()',{1}, '.','snfile', '()',{':'});
matlabbatch{2}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep;
matlabbatch{2}.spm.spatial.normalise.write.subj.resample(1).tname = 'Images to Write';
matlabbatch{2}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{2}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(1).value = 'image';
matlabbatch{2}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{2}.spm.spatial.normalise.write.subj.resample(1).tgt_spec{1}(2).value = 'e';
matlabbatch{2}.spm.spatial.normalise.write.subj.resample(1).sname = 'Segment: Bias Corr Images';
matlabbatch{2}.spm.spatial.normalise.write.subj.resample(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.spm.spatial.normalise.write.subj.resample(1).src_output = substruct('()',{1}, '.','biascorr', '()',{':'});
matlabbatch{2}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{2}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50
                                                          78 76 85];
matlabbatch{2}.spm.spatial.normalise.write.roptions.vox = aux;
matlabbatch{2}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{2}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.normalise.write.roptions.prefix = 'w';

end
save('Segment.mat','matlabbatch');
spm_jobman('run',matlabbatch)

end

function NormY_Callback(hObject, eventdata, handles)
set(handles.Native, 'Value', 0)
set(handles.Nattext, 'Visible', 'off')
set(handles.Normtext, 'Visible', 'on')

if get(handles.NormY, 'Value') == 0
    set(handles.Nattext, 'Visible', 'on')
    set(handles.Native, 'Value', 1)
    set(handles.Normtext, 'Visible', 'off')

end

function Native_Callback(hObject, eventdata, handles)
set(handles.NormY, 'Value', 0)
set(handles.Nattext, 'Visible', 'on')
set(handles.Normtext, 'Visible', 'off')

if get(handles.Native, 'Value') == 0
      set(handles.Normtext, 'Visible', 'on')
      set(handles.Nattext, 'Visible', 'off')
      set(handles.NormY, 'Value', 1)
end
