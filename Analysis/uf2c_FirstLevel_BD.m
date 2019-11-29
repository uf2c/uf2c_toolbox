function varargout = uf2c_FirstLevel_BD(varargin)
% UF2C_FIRSTLEVEL_BD MATLAB code for uf2c_FirstLevel_BD.fig
%
% UF²C - User Friendly Functional Connectivity
% Raphael Fernandes Casseb
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
%
% LOGS:
% code to output thresholded SPM images (.nii) Brunno M. de Campos, 2016
% Code to output slice view figure (options added to the GUI) Brunno M. de
% Campos, 2017


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @uf2c_FirstLevel_BD_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_FirstLevel_BD_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

warning('off','all')

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function uf2c_FirstLevel_BD_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

axes(handles.FigDesMatrix)
set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])
set(gca, 'xcolor', 'w', 'ycolor', 'w');
% set(gca, 'visible', 'off');

axes(handles.FigParadigm)
set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])
set(gca, 'xcolor', 'w', 'ycolor', 'w');
% set(gca, 'visible', 'off');

% axes(handles.axes3)
% set(gca,'xtick',[],'ytick',[])
% set(gca,'xticklabel',[],'yticklabel',[])
% set(gca, 'xcolor', 'w', 'ycolor', 'w');
% set(gca, 'visible', 'off');

set(handles.IncludeRP,'Value',1)

guidata(hObject, handles);


function varargout = uf2c_FirstLevel_BD_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function SPMbb_Callback(hObject, eventdata, handles)

% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
set(handles.statusS,'String','Running...')
drawnow
rehash

% Checking parameters:
%---------------------
if isempty(get(handles.EPIpref,'String'))
    warndlg('You must provide a functional file prefix','Attention!')
    return
end
if or(isempty(get(handles.Onsets,'String')),...
        ~isequal(size(str2num(get(handles.Onsets,'String')),2),1))
    warndlg('You must provide a nx1 matrix for Onsets','Attention!')
    return
end
if isempty(get(handles.Dur,'String'))
    warndlg('You must provide a scalar for Duration','Attention!')
    return
end
if isempty(get(handles.TR,'String'))
    warndlg('You must provide a scalar for TR','Attention!')
    return
end
if get(handles.SVcheck,'Value')
    if isempty(get(handles.T1pref,'String'))
        warndlg('You must provide a Structural file prefix','Attention!')
        return
    end
end    
if isempty(get(handles.nSlices,'String'))
    set(handles.handles.nSlices,'String','-70:2:80')
    drawnow
end

% Parameters to provide if glass brain is to be saved:
glass_b = get(handles.SaveGlassBrain,'Value');
% if glass_b
    if isempty(get(handles.TaskName,'String'))
        warndlg('You must provide a Task Name','Attention!')
    end
    if isempty(get(handles.PValue,'String'))
        warndlg('You must provide a scalar for Alpha level','Attention!')
    end
    if isempty(get(handles.ClustThresh,'String'))
        warndlg('You must provide a scalar for Cluster Threshold','Attention!')
    end
    if (~get(handles.FWE,'Value') && ~get(handles.FDR,'Value') && ...
            ~get(handles.None,'Value'))
        warndlg('You must choose a correction method','Attention!')
    end
    glass_b  = get(handles.SaveGlassBrain,'Value');
    pthresh  = str2num(get(handles.PValue,'String'));
    
    clust    = str2num(get(handles.ClustThresh,'String'));
    % Getting mult. correlation correction method:
    if get(handles.FWE,'Value') == 1
        mult_cor = 'FWE'; 
    elseif get(handles.None,'Value') == 1
        mult_cor = 'none';
    elseif get(handles.FDR,'Value') == 1
        mult_cor = 'FDR';
    end

% User inputs:
%-------------
task_name = get(handles.TaskName,'String');
onsets    = str2num(get(handles.Onsets,'String'));
duration  = str2num(get(handles.Dur,'String'));
TR        = str2num(get(handles.TR,'String'));
contrast  = 1; % By default
RP        = get(handles.IncludeRP,'Value');
FDRcInfe = get(handles.FDRclc,'Value');

% Running:
%---------
if get(handles.SVcheck,'Value')
    sl_view = get(handles.SVcheck,'Value');
%     T1_nameF = handles.T1name;
    SlicesN = str2num(get(handles.nSlices,'String'));
    
    switch get(handles.menu,'Value')
        case 1
            Orie = 'axial';
        case 2
             Orie = 'sagittal';
        case 3
             Orie = 'sagittal';
    end
             
    for i=1:size(handles.folders,2)
        uf2c_FirstLevel_BDaux(handles.folders{1,i},handles.EPI_pref,task_name,...
            onsets,duration,TR,RP,contrast,pthresh,mult_cor,clust,FDRcInfe,glass_b,...
            sl_view,handles.T1_pref,SlicesN,Orie)
    end
    
else
    sl_view = 0;
    T1_pref = [];
    SlicesN = [];
    Orie = [];

    for i=1:size(handles.folders,2)
        uf2c_FirstLevel_BDaux(handles.folders{1,i},handles.EPI_pref,task_name,...
            onsets,duration,TR,RP,contrast,pthresh,mult_cor,clust,glass_b,...
            sl_view,T1_pref,SlicesN,Orie)
    end
    
end
set(handles.EPIpref,'String','')
set(handles.T1pref,'String','')
set(handles.statusS,'String','Done!')
drawnow
if isequal(mult_cor,'FDR')
    warndlg({'You should to change back the SPM defaults.';...
        'The UF²C will open the SPM Default file,';...
        'set the option "defaults.stats.topoFDR" to 1 and save:';...
        'E.g.: defaults.stats.topoFDR     = 1;'},'Attention')
    open spm_defaults
    rehash
end

guidata(hObject, handles);


function AddFolders_Callback(hObject, eventdata, handles)
handles.folders = uipickfiles('Output','cell','Prompt','Add all subjects directories','REFilter','\');

% Checking if every folder has something inside:
%-----------------------------------------------
for i=1:size(handles.folders,2)
    clear Content
    Content = dir(fullfile(handles.folders{1,i}));
    Content(1:1) = [];
    if isempty(Content)
        warndlg(sprintf('Please check directory %s. It seems to be nothing inside',handles.folders{1,i}),'Attention!')
    end
end

set(handles.textFolders,'String',sprintf('%d folder(s) added',...
    size(handles.folders,2)))

guidata(hObject, handles);


% --- Executes on button press in PlotDesMatrix.
function PlotDesMatrix_Callback(hObject, eventdata, handles)

%warndlg('Des. matrix example will be based on first subject''s folder','Attention!')

curr_dir = pwd;

% Checking input parameters:
%---------------------------
if or(isempty(get(handles.Onsets,'String')),...
        ~isequal(size(str2num(get(handles.Onsets,'String')),2),1))
    warndlg('You must provide a nx1 matrix for Onsets','Attention!')
end
if isempty(get(handles.Dur,'String'))
    warndlg('You must provide a scalar for Duration','Attention!')
end
if isempty(get(handles.TR,'String'))
    warndlg('You must provide a scalar for TR','Attention!')
end
if isempty(get(handles.TaskName,'String'))
    warndlg('You must provide a Task Name','Attention!')
end


if ~get(handles.IncludeRP,'Value')
    % If rp files are not going to be included:
    handles.rp_file = '';
else
    handles.rp_file = dir(fullfile(handles.folders{1,1},'rp_*.txt'));
end


spm_jobman('initcfg')
clear matlabbatch

matlabbatch{1}.spm.stats.fmri_design.dir = {fileparts(which('uf2c'))};
matlabbatch{1}.spm.stats.fmri_design.timing.units = 'scans';
matlabbatch{1}.spm.stats.fmri_design.timing.RT = str2num(get(handles.TR,'String'));
matlabbatch{1}.spm.stats.fmri_design.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_design.timing.fmri_t0 = 8;

if size(handles.EPI_file,1)==1
    matlabbatch{1}.spm.stats.fmri_design.sess.nscan = handles.nvol;
    matlabbatch{1}.spm.stats.fmri_design.sess.cond.name = get(handles.TaskName,'String');
    matlabbatch{1}.spm.stats.fmri_design.sess.cond.onset = str2num(get(handles.Onsets,'String'));
    matlabbatch{1}.spm.stats.fmri_design.sess.cond.duration = str2num(get(handles.Dur,'String'));
    matlabbatch{1}.spm.stats.fmri_design.sess.cond.tmod = 0;
    matlabbatch{1}.spm.stats.fmri_design.sess.cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_design.sess.cond.orth = 1;
    matlabbatch{1}.spm.stats.fmri_design.sess.multi = {''};
    matlabbatch{1}.spm.stats.fmri_design.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_design.sess.multi_reg = ...
        {fullfile(handles.folders{1,1},handles.rp_file(1,1).name)};
    matlabbatch{1}.spm.stats.fmri_design.sess.hpf = 128;
else
    for ss = 1:size(handles.EPI_file,1)
        matlabbatch{1}.spm.stats.fmri_design.sess(ss).nscan = handles.nvol;
        matlabbatch{1}.spm.stats.fmri_design.sess(ss).cond.name = get(handles.TaskName,'String');
        matlabbatch{1}.spm.stats.fmri_design.sess(ss).cond.onset = str2num(get(handles.Onsets,'String'));
        matlabbatch{1}.spm.stats.fmri_design.sess(ss).cond.duration = str2num(get(handles.Dur,'String'));
        matlabbatch{1}.spm.stats.fmri_design.sess(ss).cond.tmod = 0;
        matlabbatch{1}.spm.stats.fmri_design.sess(ss).cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_design.sess(ss).cond.orth = 1;
        matlabbatch{1}.spm.stats.fmri_design.sess(ss).multi = {''};
        matlabbatch{1}.spm.stats.fmri_design.sess(ss).regress = struct('name', {}, 'val', {});        
        if ~get(handles.IncludeRP,'Value')
            % If rp files are not going to be included:
            matlabbatch{1}.spm.stats.fmri_design.sess(ss).multi_reg = {''};
        else
            matlabbatch{1}.spm.stats.fmri_design.sess(ss).multi_reg = ...
                {fullfile(handles.folders{1,1},handles.rp_file(ss,1).name)};
        end
        matlabbatch{1}.spm.stats.fmri_design.sess(ss).hpf = 128;
    end
end

matlabbatch{1}.spm.stats.fmri_design.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_design.bases.hrf.derivs = [1 1];
matlabbatch{1}.spm.stats.fmri_design.volt = 1;
matlabbatch{1}.spm.stats.fmri_design.global = 'None';
matlabbatch{1}.spm.stats.fmri_design.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_design.cvi = 'AR(1)';

spm_jobman('run',matlabbatch)

axes(handles.FigDesMatrix)

load(fullfile(fileparts(which('uf2c')),'SPM.mat'))
image((spm_DesMtx('sca',SPM.xX.X,SPM.xX.name)+1)*32)
colormap gray
set(gca,'FontSize',7);

delete(fullfile(fileparts(which('uf2c')),'SPM.mat'))
clear SPM

cd(curr_dir)

guidata(hObject, handles);

function TaskName_Callback(hObject, eventdata, handles)

function TaskName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Onsets_Callback(hObject, eventdata, handles)

function Onsets_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Dur_Callback(hObject, eventdata, handles)

function Dur_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function IncludeRP_Callback(hObject, eventdata, handles)
% Checking inputs:
%-----------------
% for i=1:size(handles.folders,2)
%     clear rp_file
%     handles.rp_file = dir(fullfile(handles.folders{1,i},'rp_*.txt'));
%     if isempty(handles.rp_file)
%         warndlg(sprintf('Please check directory %s. It seems to be no movement parameter file (rp file) inside',handles.folders{1,i}),'Attention!')
%     end
% end


function PValue_Callback(hObject, eventdata, handles)

function PValue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SaveGlassBrain_Callback(hObject, eventdata, handles)
% axes(handles.axes3)
% 
% TB_Path = which('uf2c');
% TB_Path = fileparts(TB_Path);
% 
% img = imread(fullfile(TB_Path,'Utilities','Glass_brain.png'));
% imagesc(img)
% axis off
% t   = title('Example');
% Pos = get(t,'Position');
% set(t,'Position',[Pos(1) Pos(2)+5 Pos(3)])
% set(t,'FontSize',7)

function TR_Callback(hObject, eventdata, handles)

function ClustThresh_Callback(hObject, eventdata, handles)

function ClustThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TR_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PlotParadigm_Callback(hObject, eventdata, handles)
% Checking input parameters:
%---------------------------
if or(isempty(get(handles.Onsets,'String')),...
        ~isequal(size(str2num(get(handles.Onsets,'String')),2),1))
    warndlg('You must provide a nx1 matrix for Onsets','Attention!')
end
if isempty(get(handles.Dur,'String'))
    warndlg('You must provide a scalar for Duration','Attention!')
end

% Plotting:
%----------
nvol = handles.nvol;
axes(handles.FigParadigm)
onsets = str2num(get(handles.Onsets,'String'));
dur    = str2num(get(handles.Dur,'String'));

ceil  = 1;
floor = 0;
ends  = onsets+dur;
x_points  = sort([0;onsets;onsets;ends;ends;nvol]);

data = zeros(size(x_points,1),2);
data(:,1) = x_points;

for i=1:(size(data,1)/2)
    clear tmp
    if rem(i,2)==0
        tmp = ceil;
    else
        tmp = floor;
    end
    data(2*i-1,2) = tmp;
    data(2*i,2)   = tmp;
end

plot(data(:,1),data(:,2))
xlim([x_points(1)-5,x_points(end)+5])
ylim([-0.5,1.5])
set(gca,'xcolor','k','ycolor','k');
set(gca,'XTick',sort([1;onsets;ends;nvol])); % Which points will be labeled
set(gca,'XTickLabel',sort([1;onsets;ends;nvol])); % Labels for XTick
set(gca,'YTick',[0;1]);
set(gca,'FontSize',7);

guidata(hObject, handles);

function nVols_Callback(hObject, eventdata, handles)
function nVols_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Run_ButtonDownFcn(hObject, eventdata, handles)

function EPIpref_Callback(hObject, eventdata, handles)

handles.EPI_pref = get(handles.EPIpref,'String');

if isempty(handles.folders)
    warndlg('You must add subjects folders first','Attention!')
else
    for i=1:size(handles.folders,2)
        clear EPI_file
        handles.EPI_file = dir(fullfile(handles.folders{1,i},[handles.EPI_pref '*.nii']));
        if isempty(handles.EPI_file)
            warndlg(sprintf('Please check directory %s. It seems to be no functional image inside',handles.folders{1,i}),'Attention!')
        end
    end
end

% Getting number of volumes for PlotParadigm and PlotDesMatrix: 
%--------------------------------------------------------------
% Just using the first EPI image (assuming all of them are based on the
% same block design)
img = nifti(fullfile(handles.folders{1,size(handles.folders,2)},...
    handles.EPI_file(1,1).name));
handles.nvol = img.dat.dim(4);

set(handles.PlotParadigm,'Enable','On')
set(handles.PlotDesMatrix,'Enable','On')
% TRvalue = img.timing.tspace;
% set(handles.TR,'String',num2str(TRvalue))
% drawnow
guidata(hObject, handles);

function EPIpref_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mult_cor_panel_SelectionChangeFcn(hObject, eventdata, handles)

if get(handles.FWE,'Value')
    set(handles.PValue,'String','0.05')
elseif get(handles.None,'Value')
    set(handles.PValue,'String','0.001')
elseif get(handles.FDR,'Value')
    set(handles.PValue,'String','0.05')
    topoFDRt = spm_get_defaults('stats.topoFDR');
    if topoFDRt
        warndlg({'You should to edit the SPM default to use FDR.';...
            'The UF²C will open the SPM Default file,';...
            'set the option "defaults.stats.topoFDR" to 0 and save:';...
            'E.g.: defaults.stats.topoFDR     = 0;';...
            'If you want to continue you need to close all and start again'},'Attention')
        open spm_defaults
    end
end

function mult_cor_panel_CreateFcn(hObject, eventdata, handles)

function pushbutton6_Callback(hObject, eventdata, handles)

function T1pref_Callback(hObject, eventdata, handles)
handles.T1_pref = get(handles.T1pref,'String');

if isempty(handles.folders)
    warndlg('You must add subjects folders first','Attention!')
% else
%     for i=1:size(handles.folders,2)
%         clear T1_file
%         T1_file{i} = dir(fullfile(handles.folders{1,i},[handles.T1_pref '*.nii']));
%         if isempty(T1_file)
%             warndlg(sprintf('Please check directory %s. It seems to be no structural image inside',handles.folders{1,i}),'Attention!')
%         end
%     end
end
% img = nifti(fullfile(handles.folders{1,size(handles.folders,2)},...
%     EPI_file.name));
% handles.nvol = img.dat.dim(4);
% 
% imgT1 = nifti(fullfile(handles.folders{1,size(handles.folders,2)},...
%     T1_file.name));

% handles.T1name = T1_file.name;

guidata(hObject, handles);

function T1pref_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton7_Callback(hObject, eventdata, handles)

function checkbox6_Callback(hObject, eventdata, handles)

function SVcheck_Callback(hObject, eventdata, handles)
if get(handles.SVcheck,'Value')
    set(handles.text20,'Enable','on')
    set(handles.T1pref,'Enable','on')
    set(handles.pushbutton7,'Enable','on')
    set(handles.text21,'Enable','on')
    set(handles.nSlices,'Enable','on')
    set(handles.menu,'Enable','on')
else
    set(handles.text20,'Enable','off')
    set(handles.T1pref,'Enable','off')
    set(handles.pushbutton7,'Enable','off')
    set(handles.text21,'Enable','off')
    set(handles.nSlices,'Enable','off')
    set(handles.menu,'Enable','off')
end

function nSlices_Callback(hObject, eventdata, handles)

function nSlices_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function menu_Callback(hObject, eventdata, handles)

function menu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FDRclc_Callback(hObject, eventdata, handles)
