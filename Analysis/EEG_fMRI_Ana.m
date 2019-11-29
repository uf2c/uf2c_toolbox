function varargout = EEG_fMRI_Ana(varargin)
% UF²C - User Friendly Functional Connectivity
% Brunno Machado de Campos
% University of Campinas, 2017
%
% Copyright (c) 2017, Brunno Machado de Campos
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
                   'gui_OpeningFcn', @EEG_fMRI_Ana_OpeningFcn, ...
                   'gui_OutputFcn',  @EEG_fMRI_Ana_OutputFcn, ...
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

function EEG_fMRI_Ana_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

guidata(hObject, handles);

function varargout = EEG_fMRI_Ana_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function addfunc_Callback(hObject, eventdata, handles)
global Image4D Path
% clearvars -global Image4D Path

[Image4D,Path] = uigetfile({'*.nii','NIfTI file'},'Add all subject functional files','MultiSelect','on');

if ~iscell(Image4D)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
    Image4D = {Image4D};
end

if ~isequal(Image4D{:},0)
    set(handles.checkstatus,'foregroundColor',[1 0 0],'String','Pending!')
    set(handles.createB,'Enable','off')
    set(handles.addedF,'String','Added!')
    set(handles.addedF,'foregroundColor',[0.471 0.671 0.188])
else
    set(handles.addedF,'foregroundColor',[1 0 0])
    set(handles.addedF,'String','No Files')
    set(handles.checkstatus,'foregroundColor',[1 0 0],'String','Pending!')
    clear Image4D Path
    set(handles.createB,'Enable','off')
end

function addmarker_Callback(hObject, eventdata, handles)
global MarkerF MarkerPath
% clearvars -global MarkerF MarkerPath

[MarkerF,MarkerPath] = uigetfile({'*.*','MARKER file'},'Add all Marker files','MultiSelect','on');

if ~iscell(MarkerF)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
    MarkerF = {MarkerF};
end

if ~isequal(MarkerF{:},0)
    set(handles.checkstatus,'foregroundColor',[1 0 0],'String','Pending!')
    set(handles.createB,'Enable','off')
    set(handles.text10,'String','Added!')
    set(handles.text10,'foregroundColor',[0.471 0.671 0.188])
    set(handles.text18,'String','')
    set(handles.markers,'String','')
    set(handles.ProMark,'Enable','off')
else
    set(handles.text10,'foregroundColor',[1 0 0])
    set(handles.text10,'String','No Files')
    set(handles.checkstatus,'foregroundColor',[1 0 0],'String','Pending!')
    clear MarkerF MarkerPath
    set(handles.markers,'String','')
    set(handles.createB,'Enable','off')
    set(handles.ProMark,'Enable','off')
    set(handles.text18,'String','')

end

function multregre_Callback(hObject, eventdata, handles)
global rpF rpPath
% clearvars -global rpF rpPath


if get(handles.multregre,'Value')
    [rpF,rpPath] = uigetfile({'*.txt','TEXT file'},'Add all multiple regression files (e.g.: rp_***.txt)','MultiSelect','on');
    
    if ~iscell(rpF)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
        rpF = {rpF};
    end
    
    set(handles.checkstatus,'foregroundColor',[1 0 0],'String','Pending!')
    set(handles.createB,'Enable','off')
    set(handles.text9,'String','Added!')
    set(handles.text9,'foregroundColor',[0.471 0.671 0.188])
else
    set(handles.text9,'foregroundColor',[1 0 0])
    set(handles.text9,'String','No Files')
    set(handles.checkstatus,'foregroundColor',[1 0 0],'String','Pending!')
    clear rpF rpPath
end

function addstru_Callback(hObject, eventdata, handles)
global struF struPath
% clearvars -global struF struPath

if get(handles.addstru,'Value')
    [struF,struPath] = uigetfile({'*.nii','NIfTI file'},'Add subject structural image','MultiSelect','off');
    set(handles.text11,'foregroundColor',[0.471 0.671 0.188])
    set(handles.text11,'String','Added!')
else
    set(handles.text11,'foregroundColor',[1 0 0])
    set(handles.text11,'String','No Files')
    clear struF struPath
end

function pushbutton9_Callback(hObject, eventdata, handles)
checkcorres
set(handles.checkstatus,'foregroundColor',[0.471 0.671 0.188],'String','Done!')

if isequal(get(handles.text10,'String'),'Added!')
    set(handles.ProMark,'Enable','on')
else
    set(handles.ProMark,'Enable','off')
    set(handles.text18,'String','')
end

if  isequal(get(handles.text18,'String'),'Ok') && isequal(get(handles.addedF,'String'),'Added!')
    set(handles.createB,'Enable','on')
end

function ProMark_Callback(hObject, eventdata, handles)
global MarkerF MarkerPath

try
    handles = rmfield(handles.StrgsUni);
end

StrgsF = cell(0,1);
set(handles.markers,'String','')
set(handles.Conds,'String','')
set(handles.regres,'String','')
set(handles.text18,'String','')

for yy = 1:size(MarkerF,2)
    if ~isequal(MarkerF{yy},'Empty')
        EEGData = [MarkerPath filesep MarkerF{yy}];
        Strgs = cell(0,1);
        [Trigger{1:5}]   =  textread(EEGData,'%s %s %d %d %s',...
            'headerlines',2,'delimiter',',','endofline','\n');

        b = strfind(Trigger{1,1},'Stimulus');
        ix = cellfun(@isempty,b);
        ix = abs(ix-1);
        idxs = find(ix);

        for i = 1:size(idxs,1)
            Strgs{i,1} = Trigger{1,2}{idxs(i)};
        end
        StrgsF = [StrgsF;Strgs];
        clear idxs ix b Trigger
    end
end  

  [StrgsUnit,IA,IC] = unique(StrgsF);
  
for j = 1:size(IA,1)
    nOfinstan(j) = sum(IC(IA(j))==IC);
    StrgsUni{j} = [StrgsUnit{j} ': ' num2str(nOfinstan(j)) ' event(s)'];
end

set(handles.markers,'String',StrgsUni)

if ~isempty(StrgsUni)
    set(handles.text18,'String','Ok')
end

if ~isempty(StrgsUni) && isequal(get(handles.addedF,'String'),'Added!')
    set(handles.createB,'Enable','on')
end
handles.StrgsUni = StrgsUni;
guidata(hObject,handles);

function AddCond_Callback(hObject, eventdata, handles)
selections = get(handles.markers,'Value');
testStr = get(handles.markers,'String');
set(handles.Conds,'String',[get(handles.Conds,'String');testStr(selections)])
testStr(selections) = [];
set(handles.markers,'String',testStr)
set(handles.markers,'Value',1)

function AddRegre_Callback(hObject, eventdata, handles)
selections = get(handles.markers,'Value');
testStr = get(handles.markers,'String');
set(handles.regres,'String',[get(handles.regres,'String');testStr(selections)])
testStr(selections) = [];
set(handles.markers,'String',testStr)
set(handles.markers,'Value',1)

function Remcon_Callback(hObject, eventdata, handles)
selections = get(handles.Conds,'Value');
set(handles.Conds,'Value',1)
tmpStr = get(handles.Conds,'String');
set(handles.markers,'String',[get(handles.markers,'String');tmpStr(selections)])
tmpStr(selections) = [];
set(handles.Conds,'String',tmpStr)

function remreg_Callback(hObject, eventdata, handles)
selections = get(handles.regres,'Value');
set(handles.regres,'Value',1)
tmpStr = get(handles.regres,'String');
set(handles.markers,'String',[get(handles.markers,'String');tmpStr(selections)])
tmpStr(selections) = [];
set(handles.regres,'String',tmpStr)

function createB_Callback(hObject, eventdata, handles)
global struPath struF rpF rpPath MarkerF MarkerPath Image4D Path

condsStr = get(handles.Conds,'String');

if isempty(condsStr)
    warndlg('You should to select at least 1 condition to proced','Attention')
    return
end

set(handles.text17,'String','Running...')
drawnow

TR = str2num(get(handles.TRin,'String'));
SampInt = str2num(get(handles.edit8,'String'));

if get(handles.radiobutton4,'Value')
    mult_cor = 'none';
    strcor = 'Uncor';
else
    mult_cor = 'FWE';
    strcor = 'FWEcor';
end

pthresh = str2num(get(handles.edit2,'String'));
clust = str2num(get(handles.edit3,'String'));

if get(handles.ATD,'Value')
    StringDrift = get(handles.TDv,'String');
    if ~isequal(StringDrift(1),'[')
        StringDrift = ['[' StringDrift];
    end
    if ~isequal(StringDrift(end),']')
        StringDrift = [StringDrift ']'];
    end
    drifts = str2num(StringDrift);
else
    drifts = 0;
end

if get(handles.checkbox5,'Value') %Slice view

    SlicesNstr = get(handles.edit4,'String');

    if ~isequal(SlicesNstr(1),'[')
        SlicesNstr = ['[' SlicesNstr];
    end
    if ~isequal(SlicesNstr(end),']')
        SlicesNstr = [SlicesNstr ']'];
    end

    SlicesN = str2num(SlicesNstr);

    switch get(handles.popupmenu1,'Value')
        case 1
            Orie = 'axial';
        case 2
            Orie = 'sagittal';
        case 3
            Orie = 'sagittal';
    end

    if get(handles.addstru,'Value')
        file_T1 = [struPath struF];
    else
        spmDir = which('spm');
        spmDir = spmDir(1:end-5);
        file_T1 = [spmDir 'canonical' filesep 'avg152T1.nii'];
    end
end

spm('defaults','FMRI');

imgRR = getframe(EEG_fMRI_Ana);
imwrite(imgRR.cdata, [Path,filesep,'1-Your_Choices_Stats.png']);

for uy = 1:size(drifts,2)

    DirName = [Path num2str(uy) '-Stat_Drift_' num2str(drifts(uy)) 's_' 'p' num2str(pthresh) strcor '_c' num2str(clust)];
    mkdir(DirName);
    
    fid = fopen([DirName filesep '1-Experiment_definitions.txt'],'w+');
    fprintf(fid,'UF²C - EEG-fMRI Event Related Analysis\r\n\r\n');
    fprintf(fid,'Experiment Definitions\r\n\r\n');
    fprintf(fid,'Case directory: \t%s\r\n',DirName);
    fprintf(fid,'\r\n');
    
    fprintf(fid,'Number of Functional Images: \t%d\r\n',size(Image4D,2));
    for irt = 1:size(Image4D,2)
        fprintf(fid,'File %d: \t%s\r\n',irt,Image4D{irt});
    end
    fprintf(fid,'\r\n');
    
    fprintf(fid,'Functional Images Time to Repeat (TR):\t%d\tseconds\r\n',TR);
    fprintf(fid,'\r\n');
    
    fprintf(fid,'Number of Regressor Files: \t%d\r\n',size(rpF,2));
    for irt = 1:size(rpF,2)
        fprintf(fid,'File %d: \t%s\r\n',irt,rpF{irt});
    end
    fprintf(fid,'\r\n');
    
    fprintf(fid,'Number of EEG Markers Files:\t%d\r\n',size(MarkerF,2));
    for irt = 1:size(MarkerF,2)
        fprintf(fid,'File %d: \t%s\r\n',irt,MarkerF{irt});
    end
    fprintf(fid,'\r\n');
    
    fprintf(fid,'EEG Sampling Interval:\t%d\tmiliseconds\r\n',SampInt);
    fprintf(fid,'\r\n');
    
    if get(handles.addstru,'Value')
        fprintf(fid,'Structural file added?\tYES\r\n');
        fprintf(fid,'Structural file:\t%s\r\n',struF);
    else
        fprintf(fid,'Structural file added?\tNO\r\n');
    end
    fprintf(fid,'\r\n');
    
    clear matlabbatch
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = {DirName};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

    stopStr = strfind(condsStr,':');
    
    clear Final contra stopStr2 Trigger idxs b ix idxs2 b2 ix2 tmpOnset tmpDur OnsetSecs DurSecs findOver Diff
    Final.CondName = cell(size(Image4D,2),size(handles.StrgsUni,2));
    BinVetOnset = [];
    for yt = 1:size(Image4D,2)
        clear Funpre file1 EEGData
        
        
        Funpre = nifti([Path Image4D{yt}]);

        for tt = 1:Funpre.dat.dim(4)
            file1{tt,1} = [Path Image4D{yt},',',sprintf('%d',tt)];
        end
        
        matlabbatch{1}.spm.stats.fmri_spec.sess(yt).scans = file1;

        if ~isequal(MarkerF{yt},'Empty')
            EEGData = [MarkerPath filesep MarkerF{yt}];
            [Trigger{1:5}]   =  textread(EEGData,'%s %s %d %d %s',...
                'headerlines',2,'delimiter',',','endofline','\n');
        end

        condN = 0;

        if get(handles.allasone,'Value')
            if yt == 1
                Final.CondName = cell(size(Image4D,2),1);
            end
            condsStr = get(handles.Conds,'String');
            stopStr = strfind(condsStr,':');
            for tt2 = 1:size(condsStr,1)
                 if ~isequal(MarkerF{yt},'Empty')
                    b = strfind(Trigger{1,2},condsStr{tt2}(1:stopStr{tt2}-1));
                    ix = cellfun(@isempty,b);
                    ix = abs(ix-1);
                    idxs = find(ix);
                    Trigger{1,2}(idxs) = {'Spike'};
                 end
            end
            clear condsStr
            condsStr{1} = 'Spike:';
            stopStr = strfind(condsStr,':');
        end
        
        
        if ~isempty(condsStr)
            for tt2 = 1:size(condsStr,1)
                 if ~isequal(MarkerF{yt},'Empty')
                    b = strfind(Trigger{1,2},condsStr{tt2}(1:stopStr{tt2}-1));
                    ix = cellfun(@isempty,b);
                    ix = abs(ix-1);
                    idxs = find(ix);
                    
                    if ~isempty(idxs)
                        tmpOnset = Trigger{1,3}(idxs);
                        tmpDur = Trigger{1,4}(idxs);

                        OnsetSecs = ((tmpOnset-Trigger{1,3}(1))*SampInt)/1000;
                        OnsetSecs = OnsetSecs + drifts(uy);
                        
                        BinVetOnsetTMP = zeros(Funpre.dat.dim(4),1);
                        IdxBin = unique(round(OnsetSecs/TR));
                        IdxBin(IdxBin>Funpre.dat.dim(4)) = [];
                        IdxBin(IdxBin<=0) = [];
                        BinVetOnsetTMP(IdxBin) = 1;
                        BinVetOnset = [BinVetOnset;BinVetOnsetTMP];
                        
                        if any(OnsetSecs<1)
%                             tmpW = warning('query');
%                             fprintf('on')
                            fprintf('\r\n');
                            fprintf('++++++++++++++++++++++++||||||||||||||++++++++++++++++++++++++\n');
                            fprintf('UF²C: One or more onsets assumed values lower than 1 (first time point)\n');
                            fprintf('considering the one or more defined "Time Drifts".\n');
                            fprintf('UF²C: They will assume the value of 1 (first time point)\n');
                            fprintf('++++++++++++++++++++++++||||||||||||||++++++++++++++++++++++++\r\n');
%                             if isequal(tmpW.state,'off')
%                                 warning('off')
%                             end
                            OnsetSecs(OnsetSecs<1) = 1;
                        end    
                        if ~get(handles.MarkDur,'Value')
                            DurSecs = str2num(get(handles.edit1,'String'));
                        else
                            DurSecs = (tmpDur*SampInt)/1000;
                            DurSecs = round(DurSecs);
                        end
                        if any((OnsetSecs+DurSecs)>(Funpre.dat.dim(4)*TR))

%                             tmpW = warning('query');
%                             warning('on')
                            fprintf('\r\n');
                            fprintf('++++++++++++++++++++++++||||||||||||||++++++++++++++++++++++++\n');
                            fprintf('UF²C: One or more onsets assumed values higher than the\n');
                            fprintf('last data point considering one or more defined "Time Drifts".\n');
                            fprintf('UF²C: They will assume the closest possible value\n');
                            fprintf('++++++++++++++++++++++++||||||||||||||++++++++++++++++++++++++\r\n');
%                             if isequal(tmpW.state,'off')
%                                 warning('off')
%                             end

                            if any(OnsetSecs>=(Funpre.dat.dim(4)*TR))
                                findOver = find(OnsetSecs>(Funpre.dat.dim(4)*TR));
                                OnsetSecs(OnsetSecs>=(Funpre.dat.dim(4)*TR)) = Funpre.dat.dim(4)*TR-1;
                                DurSecs(findOver) = 0;
                            end
                            OnDur = OnsetSecs + DurSecs;
                            if any(OnDur>=(Funpre.dat.dim(4)*TR))
                                findOver = find(OnDur>(Funpre.dat.dim(4)*TR));
                                Diff = OnDur(findOver) - (Funpre.dat.dim(4)*TR);

                                DurSecs(findOver) = DurSecs(findOver) - Diff;
                            end
                        end
                    else
                        BinVetOnset = [BinVetOnset;zeros(Funpre.dat.dim(4),1)];
                        idxs = [];
                    end
                else
                     idxs = [];
                     BinVetOnset = [BinVetOnset;zeros(Funpre.dat.dim(4),1)];
                end
                
                if ~isempty(idxs) && ~isempty(OnsetSecs)
                    condN = condN + 1;
                    Final.CondName{yt,tt2} = condsStr{tt2}(1:stopStr{tt2}-1);

                    matlabbatch{1}.spm.stats.fmri_spec.sess(yt).cond(condN).name = condsStr{tt2}(1:stopStr{tt2}-1);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(yt).cond(condN).onset = OnsetSecs;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(yt).cond(condN).duration = DurSecs;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(yt).cond(condN).tmod = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(yt).cond(condN).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess(yt).cond(condN).orth = 1;
                else
                    Final.CondName{yt,tt2} = [];
                end
            end
        else
            matlabbatch{1}.spm.stats.fmri_spec.sess(yt).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        end

        matlabbatch{1}.spm.stats.fmri_spec.sess(yt).multi = {''};

        regreStr = get(handles.regres,'String');
        stopStr2 = strfind(regreStr,':');
        regrN = 0;
        if ~isempty(regreStr)
            for tt2 = 1:size(regreStr,1)

                b2 = strfind(Trigger{1,2},regreStr{tt2}(1:stopStr2{tt2}-1));
                ix2 = cellfun(@isempty,b2);
                ix2 = abs(ix2-1);
                idxs2 = find(ix2);

                tmpOnset2 = Trigger{1,3}(idxs2);
                tmpDur2 = Trigger{1,4}(idxs2);

                OnsetSecs2 = ((tmpOnset2-Trigger{1,3}(1))*SampInt)/1000;
                OnsetVols2 = ceil(OnsetSecs2/TR);
                
                if any(OnsetVols2>Funpre.dat.dim(4))
                    OnsetVols2(OnsetVols2>Funpre.dat.dim(4)) = [];
%                     tmpW = warning('query');
%                     warning('on')
                    fprintf('UF²C: One or more regressors assumed values higher than the');
                    fprintf('last data points.');
                    fprintf('UF²C: They will be excluded');

%                     if isequal(tmpW.state,'off')
%                         warning('off')
%                     end
                end
                
                if ~isempty(idxs2) && ~isempty(OnsetVols2)
                    regrN = regrN + 1;
                    regons = zeros(Funpre.dat.dim(4),1);
                    regons(OnsetVols2) = 1;

                    matlabbatch{1}.spm.stats.fmri_spec.sess(yt).regress(regrN).name = regreStr{tt2}(1:stopStr2{tt2}-1);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(yt).regress(regrN).val = regons;
                else
                    matlabbatch{1}.spm.stats.fmri_spec.sess(yt).regress = struct('name', {}, 'val', {});
                end

            end
        else
            matlabbatch{1}.spm.stats.fmri_spec.sess(yt).regress = struct('name', {}, 'val', {});
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(yt).multi_reg = {[rpPath rpF{yt}]};
        matlabbatch{1}.spm.stats.fmri_spec.sess(yt).hpf = 128;

        Final.Conds(yt,1) = condN;
        Final.Regres(yt,1) = regrN;
    end
    save([DirName filesep 'BinVetOnset'],'BinVetOnset')
    
    if get(handles.multregre,'Value')
        for gf = 1:yt
            regTest = textread([rpPath rpF{gf}]);
            SumReg = sum(regTest,1);
            Exreg = find(SumReg==0);
            regTest(:,Exreg) = [];
            Final.Regres(gf,1) = Final.Regres(gf,1) + size(regTest,2);
        end
    end

    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    
%%% No time Derivatives
    line = 0;
    for nc = 1:yt
        if Final.Conds(nc)>0
            CondNameTMP = Final.CondName(nc,:);
            CondNameTMP = CondNameTMP(~cellfun('isempty',CondNameTMP));
            for nc2 = 1:Final.Conds(nc)
                line = line + 1;
                contra = zeros(1,Final.Conds(nc));
                if Final.Conds(nc)>1
                    contra(1,(nc2)) = 1;
                else
                    contra(1,nc2) = 1;
                end
                contra = [contra,zeros(1,Final.Regres(nc))];
                Contrast{line,1} = contra;
                Contrast{line,2} = CondNameTMP{1,nc2};
            end
        end
    end
    
%%% With Time and Dispersion derivatives
%     line = 0;
%     for nc = 1:yt
%         if Final.Conds(nc)>0
%             CondNameTMP = Final.CondName(nc,:);
%             CondNameTMP = CondNameTMP(~cellfun('isempty',CondNameTMP));
%             for nc2 = 1:Final.Conds(nc)
%                 line = line + 1;
%                 contra = zeros(1,Final.Conds(nc)*3);
%                 if Final.Conds(nc)>1
%                     contra(1,(nc2*3)-2) = 1;
%                 else
%                     contra(1,nc2) = 1;
%                 end
%                 contra = [contra,zeros(1,Final.Regres(nc))];
%                 Contrast{line,1} = contra;
%                 Contrast{line,2} = CondNameTMP{1,nc2};
%             end
%         end
%     end
    
    condsName2 = unique(Contrast(:,2));
    cotr = 0;

    for rr = 1:size(condsName2,1)
        clear idxs2
        ctrF = [];
        b2 = strfind(Contrast(:,2),condsName2{rr});
        ix2 = cellfun(@isempty,b2);
        ix2 = abs(ix2-1);
        idxs2 = find(ix2);
        rr2F = 0;

        for rr2 = 1:yt

            rr2F = rr2F + 1;
            clear tster idxs
            tster = Final.CondName(rr2,:);
            tster = tster(~cellfun('isempty',tster));
%             prepCon = zeros(1,size(tster,2)*3);
            prepCon = zeros(1,size(tster,2));
            
            b = strfind(tster,condsName2{rr}); 
            ix = cellfun(@isempty,b);
            ix = abs(ix-1);
            idxs = find(ix);

            if isempty(idxs)
                rr2F = rr2F -1;
                ctrF = [ctrF,prepCon,zeros(1,Final.Regres(rr2))];
            else
                if isequal(size(idxs2,1),1)
                    ctrF = [ctrF,Contrast{idxs2(1),1}];
                else
                    ctrF = [ctrF,Contrast{idxs2(rr2F),1}];
                end
            end
        end

        cotr = cotr + 1;
        matlabbatch{3}.spm.stats.con.consess{cotr}.tcon.name = ['Positive Infere: ' condsName2{rr}];
        matlabbatch{3}.spm.stats.con.consess{cotr}.tcon.weights = ctrF;
        matlabbatch{3}.spm.stats.con.consess{cotr}.tcon.sessrep = 'none';

        cotr = cotr + 1;
        matlabbatch{3}.spm.stats.con.consess{cotr}.tcon.name = ['Negative Infere: ' condsName2{rr}];
        matlabbatch{3}.spm.stats.con.consess{cotr}.tcon.weights = ctrF*-1;
        matlabbatch{3}.spm.stats.con.consess{cotr}.tcon.sessrep = 'none';
    end
    matlabbatch{3}.spm.stats.con.delete = 0;

    save([DirName filesep 'BatchAna'],'matlabbatch')
    spm_jobman('run',matlabbatch)

    subloop = 0;

    for rr = 1:2:2*size(condsName2,1)
        clear matlabbatch
        matlabbatch{1}.spm.stats.results.spmmat(1) = {[DirName filesep 'SPM.mat']};
        matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
        matlabbatch{1}.spm.stats.results.conspec.contrasts = rr;
        matlabbatch{1}.spm.stats.results.conspec.threshdesc = mult_cor;
        matlabbatch{1}.spm.stats.results.conspec.thresh = pthresh;
        matlabbatch{1}.spm.stats.results.conspec.extent = clust;
        matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
        matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
        matlabbatch{1}.spm.stats.results.units = 1;
        matlabbatch{1}.spm.stats.results.print = 'tif';
        matlabbatch{1}.spm.stats.results.write.tspm.basename = [condsName2{rr-subloop} '_Posit_' num2str(pthresh) strcor '_c' num2str(clust) '_Drift_' num2str(drifts(uy))];
        [output_list, ~] = spm_jobman('run',matlabbatch);
        
        if rr>=10
            spmTprefix  = 'spmT_00';
        else
            if rr>=100
                spmTprefix  = 'spmT_0';
            else
                spmTprefix  = 'spmT_000';
            end
        end
        
        [~] = uf2c_AnatDescrip([DirName,filesep],...
            [spmTprefix num2str(rr) '_' condsName2{rr-subloop} '_Posit_' num2str(pthresh) strcor '_c' num2str(clust) '_Drift_' num2str(drifts(uy)),'_AnatDescrip.txt'],...
            {[DirName filesep spmTprefix num2str(rr) '_' condsName2{rr-subloop} '_Posit_' num2str(pthresh) strcor '_c' num2str(clust) '_Drift_' num2str(drifts(uy)),'.nii']},...
            'map');

        FDRcInfe = get(handles.FDRclust,'Value');
        
        if FDRcInfe
            clear matlabbatch
            FDRcV_Pos = round(output_list{1, 1}.TabDatvar.ftr{5,2}(4));
            matlabbatch{1}.spm.stats.results.spmmat(1) = {[DirName filesep 'SPM.mat']};
            matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
            matlabbatch{1}.spm.stats.results.conspec.contrasts = rr;
            matlabbatch{1}.spm.stats.results.conspec.threshdesc = mult_cor;
            matlabbatch{1}.spm.stats.results.conspec.thresh = pthresh;
            matlabbatch{1}.spm.stats.results.conspec.extent = FDRcV_Pos;
            matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
            matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
            matlabbatch{1}.spm.stats.results.units = 1;
            matlabbatch{1}.spm.stats.results.print = 'tif';
            matlabbatch{1}.spm.stats.results.write.tspm.basename = [condsName2{rr-subloop} '_Posit_iniThre' num2str(pthresh) '_FDRc_0.05_' num2str(FDRcV_Pos) 'vox_Drift_' num2str(drifts(uy))];
            spm_jobman('run',matlabbatch)
            
            [~] = uf2c_AnatDescrip([DirName,filesep],...
                [spmTprefix num2str(rr) '_' condsName2{rr-subloop} '_Posit_iniThre' num2str(pthresh) '_FDRc_0.05_' num2str(FDRcV_Pos) 'vox_Drift_' num2str(drifts(uy)),'_AnatDescrip.txt'],...
                {[DirName filesep spmTprefix num2str(rr) '_' condsName2{rr-subloop} '_Posit_iniThre' num2str(pthresh) '_FDRc_0.05_' num2str(FDRcV_Pos) 'vox_Drift_' num2str(drifts(uy)),'.nii']},...
                'map');
        end

        clear matlabbatch
        matlabbatch{1}.spm.stats.results.spmmat(1) = {[DirName filesep 'SPM.mat']};
        matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
        matlabbatch{1}.spm.stats.results.conspec.contrasts = rr+1;
        matlabbatch{1}.spm.stats.results.conspec.threshdesc = mult_cor;
        matlabbatch{1}.spm.stats.results.conspec.thresh = pthresh;
        matlabbatch{1}.spm.stats.results.conspec.extent = clust;
        matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
        matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
        matlabbatch{1}.spm.stats.results.units = 1;
        matlabbatch{1}.spm.stats.results.print = 'tif';
        matlabbatch{1}.spm.stats.results.write.tspm.basename = [condsName2{rr-subloop} '_Negat_' num2str(pthresh) strcor '_c' num2str(clust) '_Drift_' num2str(drifts(uy))];
        [output_list, ~] = spm_jobman('run',matlabbatch);
        
        if rr+1>=10
            spmTprefix  = 'spmT_00';
        else
            if rr+1>=100
                spmTprefix  = 'spmT_0';
            else
                spmTprefix  = 'spmT_000';
            end
        end
        
        [~] = uf2c_AnatDescrip([DirName,filesep],...
            [spmTprefix num2str(rr+1) '_' condsName2{rr-subloop} '_Negat_' num2str(pthresh) strcor '_c' num2str(clust) '_Drift_' num2str(drifts(uy)),'_AnatDescrip.txt'],...
            {[DirName filesep spmTprefix num2str(rr+1) '_' condsName2{rr-subloop} '_Negat_' num2str(pthresh) strcor '_c' num2str(clust) '_Drift_' num2str(drifts(uy)),'.nii']},...
            'map');

        if FDRcInfe
            FDRcV_Neg = round(output_list{1, 1}.TabDatvar.ftr{5,2}(4));
            matlabbatch{1}.spm.stats.results.spmmat(1) = {[DirName filesep 'SPM.mat']};
            matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
            matlabbatch{1}.spm.stats.results.conspec.contrasts = rr+1;
            matlabbatch{1}.spm.stats.results.conspec.threshdesc = mult_cor;
            matlabbatch{1}.spm.stats.results.conspec.thresh = pthresh;
            matlabbatch{1}.spm.stats.results.conspec.extent = FDRcV_Neg;
            matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
            matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
            matlabbatch{1}.spm.stats.results.units = 1;
            matlabbatch{1}.spm.stats.results.print = 'tif';
            matlabbatch{1}.spm.stats.results.write.tspm.basename = [condsName2{rr-subloop} '_Negat_iniThre' num2str(pthresh) '_FDRc_0.05_' num2str(FDRcV_Neg) 'vox_Drift_' num2str(drifts(uy))];
            spm_jobman('run',matlabbatch)
            
            [~] = uf2c_AnatDescrip([DirName,filesep],...
                [spmTprefix num2str(rr+1) '_' condsName2{rr-subloop} '_Negat_iniThre' num2str(pthresh) '_FDRc_0.05_' num2str(FDRcV_Neg) 'vox_Drift_' num2str(drifts(uy)),'_AnatDescrip.txt'],...
                {[DirName filesep spmTprefix num2str(rr+1) '_' condsName2{rr-subloop} '_Negat_iniThre' num2str(pthresh) '_FDRc_0.05_' num2str(FDRcV_Neg) 'vox_Drift_' num2str(drifts(uy)),'.nii']},...
                'map');
        end
        
        if get(handles.checkbox5,'Value') %Slice view
            if rr>=10
                spmTprefix  = 'spmT_00';
            else
                if rr>=100
                    spmTprefix  = 'spmT_0';
                else
                    spmTprefix  = 'spmT_000';
                end
            end
            
            namePos{rr-subloop,uy} = [DirName filesep,spmTprefix,num2str(rr),'_',...
                condsName2{rr-subloop},'_Posit_',num2str(pthresh),strcor,...
                '_c' num2str(clust),'_Drift_',num2str(drifts(uy)),'.nii'];
            
            SliceVPos{rr-subloop,uy} = [num2str(uy) '-SliceView_',spmTprefix,num2str(rr),'_',...
                condsName2{rr-subloop},'_Posit_',num2str(pthresh),strcor,...
                '_c',num2str(clust),'_Drift_',num2str(drifts(uy)),'.png'];
            
            SliceView_uf2c(file_T1,namePos{rr-subloop,uy},0.1,0.1,Orie,'hot',...
                SlicesN,['Time Drift ',num2str(drifts(uy))],...
                SliceVPos{rr-subloop,uy},[DirName filesep])
            
            if rr+1>=10
                spmTprefix  = 'spmT_00';
            else
                if rr+1>=100
                    spmTprefix  = 'spmT_0';
                else
                    spmTprefix  = 'spmT_000';
                end
            end
            
            nameNeg{rr-subloop,uy} = [DirName filesep,spmTprefix,num2str(rr+1),'_',...
                condsName2{rr-subloop},'_Negat_',num2str(pthresh),strcor,...
                '_c' num2str(clust),'_Drift_',num2str(drifts(uy)),'.nii'];
            
            SliceVNeg{rr-subloop,uy} = [num2str(uy) '-SliceView_',spmTprefix,num2str(rr+1),'_',...
                condsName2{rr-subloop},'_Negat_',num2str(pthresh),strcor,...
                '_c',num2str(clust),'_Drift_',num2str(drifts(uy)),'.png'];
            
            SliceView_uf2c(file_T1,nameNeg{rr-subloop,uy},0.1,0.1,Orie,'winter',...
                SlicesN,['Time Drift ',num2str(drifts(uy))],...
                SliceVNeg{rr-subloop,uy},[DirName filesep])
            
            if FDRcInfe
                    if rr>=10
                        spmTprefix  = 'spmT_00';
                    else
                        if rr>=100
                            spmTprefix  = 'spmT_0';
                        else
                            spmTprefix  = 'spmT_000';
                        end
                    end

                    namePos2{rr-subloop,uy} = [DirName filesep,spmTprefix,num2str(rr),'_',...
                        condsName2{rr-subloop},'_Posit_iniThre',num2str(pthresh),...
                        '_FDRc_0.05_',num2str(FDRcV_Pos),'vox_Drift_',num2str(drifts(uy)),'.nii'];

                    SliceVPos2{rr-subloop,uy} = [num2str(uy) '-SliceView_',spmTprefix,num2str(rr),'_',...
                        condsName2{rr-subloop},'_Posit_iniThre',num2str(pthresh),...
                        '_FDRc_0.05_',num2str(FDRcV_Pos),'vox_Drift_',num2str(drifts(uy)),'.png'];

                    SliceView_uf2c(file_T1,namePos2{rr-subloop,uy},0.1,0.1,Orie,'hot',...
                        SlicesN,['Time Drift ',num2str(drifts(uy))],...
                        SliceVPos2{rr-subloop,uy},[DirName filesep])
                    if rr+1>=10
                        spmTprefix  = 'spmT_00';
                    else
                        if rr+1>=100
                            spmTprefix  = 'spmT_0';
                        else
                            spmTprefix  = 'spmT_000';
                        end
                    end

                    nameNeg2{rr-subloop,uy} = [DirName filesep,spmTprefix,num2str(rr+1),'_',...
                        condsName2{rr-subloop},'_Negat_iniThre',num2str(pthresh),...
                        '_FDRc_0.05_',num2str(FDRcV_Neg),'vox_Drift_',num2str(drifts(uy)),'.nii'];

                    SliceVNeg2{rr-subloop,uy} = [num2str(uy) '-SliceView_',spmTprefix,num2str(rr+1),'_',...
                        condsName2{rr-subloop},'_Negat_iniThre',num2str(pthresh),...
                        '_FDRc_0.05_',num2str(FDRcV_Neg),'vox_Drift_',num2str(drifts(uy)),'.png'];

                    SliceView_uf2c(file_T1,nameNeg2{rr-subloop,uy},0.1,0.1,Orie,'winter',...
                        SlicesN,['Time Drift ',num2str(drifts(uy))],...
                        SliceVNeg2{rr-subloop,uy},[DirName filesep])
             end
        end
        subloop = subloop + 1;
    end
    
    fprintf(fid,'Number of Markers found:\t %d\r\n',size(handles.StrgsUni,2));
    for irt = 1:size(handles.StrgsUni,2)
        fprintf(fid,'Marker %d:\t %s\r\n',irt,handles.StrgsUni{irt});
    end
    fprintf(fid,'\r\n');
    
    fprintf(fid,'Number of Conditions Included:\t %d\r\n',size(condsStr,1));
    for irt = 1:size(condsStr,1)
        fprintf(fid,'Included Condition %d:\t %s\r\n',irt,condsStr{irt});
    end
    fprintf(fid,'\r\n');
    if ~get(handles.MarkDur,'Value')
        fprintf(fid,'Fixed Condition Duration:\t %s secs\r\n',get(handles.edit1,'String'));
        fprintf(fid,'\r\n');
    end
    
    fprintf(fid,'Number of Markers added as Regressors:\t %d\r\n',size(regreStr,1));
    for irt = 1:size(regreStr,1)
        fprintf(fid,'Included Marker''s Regressor  %d:\t %s\r\n',irt,regreStr{irt});
    end
    fprintf(fid,'\r\n');
    
    if get(handles.ATD,'VAlue')
        fprintf(fid,'Time drift to conditions onsets:\t YES\r\n');
        fprintf(fid,'Time drift to condition:\t %d\tseconds\r\n',drifts(uy));
    else
        fprintf(fid,'Time drift to condition onsets:\t NO\r\n');
    end
    fprintf(fid,'\r\n');
    
    fprintf(fid,'Conditions time and dispersion derivatives:\t NO\r\n');
    fprintf(fid,'Statistical correction:\t %s\r\n',strcor);
    fprintf(fid,'Statistical threshold:\t %s\r\n',num2str(pthresh));
    fprintf(fid,'Extend threshold (voxels):\t %s\r\n',num2str(clust));
    fclose(fid);
    DirNameFF{uy} = DirName;
end

if get(handles.ATD,'Value')
    if size(drifts,2)>1
        GenDdir = [Path '0-TotalDriftsResults'];
        mkdir(GenDdir)
        mkdir([GenDdir,filesep,'Cumulative_T-ScoreMaps'])
        mkdir([GenDdir,filesep,'Cumulative_ContrastMaps'])
        mkdir([GenDdir,filesep,'Average_ContrastMaps'])
        mkdir([GenDdir,filesep,'Average_T-ScoreMaps'])
        mkdir([GenDdir,filesep,'Merged_VoxWise_Thrs_T-score'])
        mkdir([GenDdir,filesep,'Merged_VoxWise_ConMaps'])
        mkdir([GenDdir,filesep,'Videos'])
        
        subloop = 0;
        subloop0 = 0;
        for rr = 1:2:2*size(condsName2,1)
            
            if rr>=10
                spmTprefix  = 'spmT_00';
                spmTprefix2  = 'con_00';
            else
                if rr>=100
                    spmTprefix  = 'spmT_0';
                    spmTprefix2  = 'con_0';
                else
                    spmTprefix  = 'spmT_000';
                    spmTprefix2  = 'con_000';
                end
            end

            % creating maps with the higher T-values among all Delays (mixing delays)
            % (voxel-wisely)
            for ff = 1:size(drifts,2)
                tmp1 = nifti(namePos{rr-subloop0,ff});
                mattmp1 = tmp1.dat(:,:,:);
                mattmp1 = reshape(mattmp1,prod(tmp1.dat.dim),1)';
                mattmp1(isnan(mattmp1)) = 0;

                tmp2 = nifti(namePos2{rr-subloop0,ff});
                mattmp2 = tmp2.dat(:,:,:);
                mattmp2 = reshape(mattmp2,prod(tmp2.dat.dim),1)';
                mattmp2(isnan(mattmp2)) = 0;

                tmp3 = nifti(nameNeg{rr-subloop0,ff});
                mattmp3 = tmp3.dat(:,:,:);
                mattmp3 = reshape(mattmp3,prod(tmp3.dat.dim),1)';
                mattmp3(isnan(mattmp3)) = 0;

                tmp4 = nifti(nameNeg2{rr-subloop0,ff});
                mattmp4 = tmp4.dat(:,:,:);
                mattmp4 = reshape(mattmp4,prod(tmp4.dat.dim),1)';
                mattmp4(isnan(mattmp4)) = 0;

                BVt_Pos1(ff,:) = mattmp1;
                BVt_Pos2(ff,:) = mattmp2;
                BVt_Neg1(ff,:) = mattmp3;
                BVt_Neg2(ff,:) = mattmp4;
                
                
                %Creating Merged voxel wisely contrast maps (mixing delays)
                tmp1con = nifti([DirNameFF{ff},filesep,spmTprefix2,num2str(rr),'.nii']);
                mattmp1con = tmp1con.dat(:,:,:);
                mattmp1con = reshape(mattmp1con,prod(tmp1con.dat.dim),1)';
                mattmp1con(isnan(mattmp1con)) = 0;

                tmp1conNeg = nifti([DirNameFF{ff},filesep,spmTprefix2,num2str(rr+1),'.nii']);
                mattmp1conNeg = tmp1conNeg.dat(:,:,:);
                mattmp1conNeg = reshape(mattmp1conNeg,prod(tmp1conNeg.dat.dim),1)';
                mattmp1conNeg(isnan(mattmp1conNeg)) = 0;

                BVcon_Pos1(ff,:) = mattmp1con;
                BVcon_Neg2(ff,:) = mattmp1conNeg;
            end
        
            BVt_PMax1 = max(BVt_Pos1,[],1);
                BVt_PMax1(BVt_PMax1==0) = NaN;
                mattmp1 = reshape(BVt_PMax1',tmp1.dat.dim);
                [~,baseName,~] = fileparts(namePos{rr-subloop0,ff});
                stoptmp = strfind(baseName,'_Drift');
                namePos3{1} = [GenDdir,filesep,'Merged_VoxWise_Thrs_T-score',filesep,'Merged_VoxWise_Thrs_T-score_',baseName(11:stoptmp-1),'.nii'];
                tmp1.dat.fname = namePos3{1};
                tmp1.dat(:,:,:) = mattmp1;
                create(tmp1)

            BVt_PMax2 = max(BVt_Pos2,[],1);
                BVt_PMax2(BVt_PMax2==0) = NaN;
                mattmp2 = reshape(BVt_PMax2',tmp2.dat.dim);
                [~,baseName,~] = fileparts(namePos2{rr-subloop0,ff});
                stoptmp = strfind(baseName,'_0.05');
                namePos3{2} = [GenDdir,filesep,'Merged_VoxWise_Thrs_T-score',filesep,'Merged_VoxWise_Thrs_T-score_',baseName(11:stoptmp+4),'.nii'];
                tmp2.dat.fname = namePos3{2};
                tmp2.dat(:,:,:) = mattmp2;
                create(tmp2)

            BVt_NMax1 = max(BVt_Neg1,[],1);
                BVt_NMax1(BVt_NMax1==0) = NaN;
                mattmp3 = reshape(BVt_NMax1',tmp3.dat.dim);
                [~,baseName,~] = fileparts(nameNeg{rr-subloop0,ff});
                stoptmp = strfind(baseName,'_Drift');
                nameNeg3{1} = [GenDdir,filesep,'Merged_VoxWise_Thrs_T-score',filesep,'Merged_VoxWise_Thrs_T-score_',baseName(11:stoptmp-1),'.nii'];
                tmp3.dat.fname = nameNeg3{1};
                tmp3.dat(:,:,:) = mattmp3;
                create(tmp3)

            BVt_NMax2 = max(BVt_Neg2,[],1);
                BVt_NMax2(BVt_NMax2==0) = NaN;
                mattmp4 = reshape(BVt_NMax2',tmp4.dat.dim);
                [~,baseName,~] = fileparts(nameNeg2{rr-subloop0,ff});
                stoptmp = strfind(baseName,'_0.05');
                nameNeg3{2} = [GenDdir,filesep,'Merged_VoxWise_Thrs_T-score',filesep,'Merged_VoxWise_Thrs_T-score_',baseName(11:stoptmp+4),'.nii'];
                tmp4.dat.fname = nameNeg3{2};
                tmp4.dat(:,:,:) = mattmp4;
                create(tmp4)
                
                
            % ConMap
                conMax1 = max(BVcon_Pos1,[],1);
                conMax1(conMax1==0) = NaN;
                mattmp1con = reshape(conMax1',tmp1con.dat.dim);
                [~,baseName,~] = fileparts(namePos{rr-subloop0,ff});
                stoptmp = strfind(baseName,'_Posit');
                namePosCon1 = [GenDdir,filesep,'Merged_VoxWise_ConMaps',filesep,'Merged_VoxWise_Posi_ConMap_',baseName(11:stoptmp-1),'.nii'];
                tmp1con.dat.fname = namePosCon1;
                tmp1con.dat(:,:,:) = mattmp1con;
                create(tmp1con)

                conMax2 = max(BVcon_Neg2,[],1);
                conMax2(conMax2==0) = NaN;
                mattmp1con2 = reshape(conMax2',tmp1conNeg.dat.dim);
                [~,baseName,~] = fileparts(nameNeg{rr-subloop0,ff});
                stoptmp = strfind(baseName,'_Negat');
                namePosCon2 = [GenDdir,filesep,'Merged_VoxWise_ConMaps',filesep,'Merged_VoxWise_Nega_ConMap_',baseName(11:stoptmp-1),'.nii'];
                tmp1conNeg.dat.fname = namePosCon2;
                tmp1conNeg.dat(:,:,:) = mattmp1con2;
                create(tmp1conNeg)
            % done!
                
               subloop0 = subloop0 + 1;
            
            if get(handles.checkbox5,'Value') % Slice View
                
                % Slice View of the maps with the higher T-values among all Delays
                % (voxel-wisely)
                [path100,tmp100,~] = fileparts(namePos3{1});
                SliceView_uf2c(file_T1,namePos3{1},0.1,0.1,Orie,'hot',...
                        SlicesN,'Merged_VoxWise_T-score',[tmp100,'.png'],path100)
                    
                [path100,tmp100,~] = fileparts(namePos3{2});
                SliceView_uf2c(file_T1,namePos3{2},0.1,0.1,Orie,'hot',...
                        SlicesN,'Merged_VoxWise_T-score',[tmp100,'.png'],path100)
                    
                [path100,tmp100,~] = fileparts(nameNeg3{1});
                SliceView_uf2c(file_T1,nameNeg3{1},0.1,0.1,Orie,'winter',...
                        SlicesN,'Merged_VoxWise_T-score',[tmp100,'.png'],path100)
                    
                [path100,tmp100,~] = fileparts(nameNeg3{2});
                SliceView_uf2c(file_T1,nameNeg3{2},0.1,0.1,Orie,'winter',...
                        SlicesN,'Merged_VoxWise_T-score',[tmp100,'.png'],path100)

                mkdir([GenDdir,filesep,'Temporary'])
                
%                 v1PC = VideoWriter([GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_ProgreCumu_Positive','.mp4'],'MPEG-4');
%                 v1PC.FrameRate = 0.25;
%                 v1PC.Quality = 100;
%                 open(v1PC)
                VFrame4D1 = imread([DirNameFF{1},filesep,SliceVPos{rr-subloop,1}]);
%                 writeVideo(v1PC,VFrame4D1)
                [VFrame4D1,cm] = rgb2ind(VFrame4D1,256);
                imwrite(VFrame4D1,cm,...
                   [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_ProgreCumu_Positive','.gif'],...
                    'gif','DelayTime',1.5,'Loopcount',inf);

                Stru1 = nifti(namePos{rr-subloop,1});
                MatPC1 = Stru1.dat(:,:,:);
                MatPC1(isnan(MatPC1)) = 0;
                for j = 2:size(drifts,2)
                    StruTMP = nifti(namePos{rr-subloop,j});
                    MatTMP = StruTMP.dat(:,:,:);
                    MatTMP(isnan(MatTMP)) = 0;
                    MatPC1 = MatPC1 + MatTMP;
                    StruTMP.dat.fname = [GenDdir,filesep,'Temporary',filesep,'TMP.nii'];
                    StruTMP.dat(:,:,:) = MatPC1;
                    create(StruTMP)
                    SliceView_uf2c(file_T1,[GenDdir,filesep,'Temporary',...
                        filesep,'TMP.nii'],0.1,0.1,Orie,'hot',...
                        SlicesN,['Time Drift ',num2str(drifts(j))],...
                        'TMP.png',[GenDdir,filesep,'Temporary',filesep])
                    VFrame4Dtmp = imread([GenDdir,filesep,'Temporary',filesep,'TMP.png']);
%                     writeVideo(v1PC,VFrame4Dtmp)
                    [VFrame4Dtmp,cm] = rgb2ind(VFrame4Dtmp,256);
                    imwrite(VFrame4Dtmp,cm,...
                       [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_ProgreCumu_Positive','.gif'],...
                        'gif','DelayTime',1.5,'WriteMode','append');
                end
%                 close(v1PC)
                clear v1PC VFrame4D1 VFrame4Dtmp StruTMP MatTMP
                
                if FDRcInfe
%                     v1PC = VideoWriter([GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_ProgreCumu_Positive_FDRcInfer','.mp4'],'MPEG-4');
%                     v1PC.FrameRate = 0.5;
%                     v1PC.Quality = 100;
%                     open(v1PC)
                    VFrame4D1 = imread([DirNameFF{1},filesep,SliceVPos2{rr-subloop,1}]);
%                     sfra = size(VFrame4D1);
%                     fnd2 = find(rem(sfra(1:2),2));
%                     if any(fnd2)
%                         if isequal(fnd2,1)
%                             VFrame4D1 = [VFrame4D1;zeros(1,sfra(2),3)];
%                         end
%                         if isequal(fnd2,2)
%                             VFrame4D1 = [VFrame4D1,zeros(sfra(1),1,3)];
%                         end
%                         if isequal(fnd2,[1 2])
%                             VFrame4D1 = [VFrame4D1;zeros(1,sfra(2),3)];
%                             VFrame4D1 = [VFrame4D1,zeros(sfra(1)+1,1,3)];
%                         end
%                     end
%                     writeVideo(v1PC,VFrame4D1)
                    [VFrame4D1,cm] = rgb2ind(VFrame4D1,256);
                    imwrite(VFrame4D1,cm,...
                        [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_ProgreCumu_Positive_FDRcInfer','.gif'],...
                        'gif','DelayTime',1.5,'Loopcount',inf);
                    
                    Stru1 = nifti(namePos2{rr-subloop,1});
                    MatPC1 = Stru1.dat(:,:,:);
                    MatPC1(isnan(MatPC1)) = 0;
                    
                    for j = 2:size(drifts,2)
                        StruTMP = nifti(namePos2{rr-subloop,j});
                        MatTMP = StruTMP.dat(:,:,:);
                        MatTMP(isnan(MatTMP)) = 0;
                        MatPC1 = MatPC1 + MatTMP;
                        StruTMP.dat.fname = [GenDdir,filesep,'Temporary',filesep,'TMP.nii'];
                        StruTMP.dat(:,:,:) = MatPC1;
                        create(StruTMP)
                        SliceView_uf2c(file_T1,[GenDdir,filesep,'Temporary',...
                            filesep,'TMP.nii'],0.1,0.1,Orie,'hot',...
                            SlicesN,['Time Drift ',num2str(drifts(j))],...
                            'TMP.png',[GenDdir,filesep,'Temporary',filesep]);
                        VFrame4Dtmp = imread([GenDdir,filesep,'Temporary',filesep,'TMP.png']);
%                         
%                         sfra = size(VFrame4Dtmp);
%                         fnd2 = find(rem(sfra(1:2),2));
%                         if any(fnd2)
%                             if isequal(fnd2,1)
%                                 VFrame4Dtmp = [VFrame4Dtmp;zeros(1,sfra(2),3)];
%                             end
%                             if isequal(fnd2,2)
%                                 VFrame4Dtmp = [VFrame4Dtmp,zeros(sfra(1),1,3)];
%                             end
%                             if isequal(fnd2,[1 2])
%                                 VFrame4Dtmp = [VFrame4Dtmp;zeros(1,sfra(2),3)];
%                                 VFrame4Dtmp = [VFrame4Dtmp,zeros(sfra(1)+1,1,3)];
%                             end
%                         end
%                         writeVideo(v1PC,VFrame4Dtmp)
                        [VFrame4Dtmp,cm] = rgb2ind(VFrame4Dtmp,256);
                        imwrite(VFrame4Dtmp,cm,...
                            [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_ProgreCumu_Positive_FDRcInfer','.gif'],...
                            'gif','DelayTime',1.5,'WriteMode','append');
                    end
%                     close(v1PC)
                end
            end
%             
%             v1 = VideoWriter([GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_Positive','.mp4'],'MPEG-4');
%             v1.FrameRate = 0.25;
%             v1.Quality = 100;
%             open(v1)
            for i = 1:size(drifts,2)
                if get(handles.checkbox5,'Value') % Slice View
                    VFrame4D = imread([DirNameFF{i},filesep,SliceVPos{rr-subloop,i}]);
%                     writeVideo(v1,VFrame4D)
                    [VFrame4D,cm] = rgb2ind(VFrame4D,256);
                    if i==1
                        imwrite(VFrame4D,cm,...
                            [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_Positive','.gif'],...
                            'gif','DelayTime',1.5,'Loopcount',inf);
                    else
                        imwrite(VFrame4D,cm,...
                           [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_Positive','.gif'],...
                            'gif','DelayTime',1.5,'WriteMode','append');
                    end
                end
                StruPoss = nifti([DirNameFF{i},filesep,spmTprefix,num2str(rr),'.nii']);
                Mat4DPoss(:,:,:,i) = StruPoss.dat(:,:,:);
                
                StruPossCon = nifti([DirNameFF{i},filesep,spmTprefix2,num2str(rr),'.nii']);
                Mat4DPossCon(:,:,:,i) = StruPossCon.dat(:,:,:);
            end
%             if get(handles.checkbox5,'Value') %Slice view
%                 close(v1)
%             end
            
            if FDRcInfe
%                 v1 = VideoWriter([GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_Positive_FDRcInfe','.mp4'],'MPEG-4');
%                 v1.FrameRate = 0.25;
%                 v1.Quality = 100;
%                 open(v1)
                for i = 1:size(drifts,2)
                    if get(handles.checkbox5,'Value') % Slice View
                        VFrame4D = imread([DirNameFF{i},filesep,SliceVPos2{rr-subloop,i}]);
%                         writeVideo(v1,VFrame4D)
                        [VFrame4D,cm] = rgb2ind(VFrame4D,256);
                        if i==1
                            imwrite(VFrame4D,cm,...
                                [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_Positive_FDRcInfe','.gif'],...
                                'gif','DelayTime',1.5,'Loopcount',inf);
                        else
                            imwrite(VFrame4D,cm,...
                               [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_Positive_FDRcInfe','.gif'],...
                                'gif','DelayTime',1.5,'WriteMode','append');
                        end
                    end
                end
%                 if get(handles.checkbox5,'Value') %Slice view
%                     close(v1)
%                 end
            end
            
            %Cumulative Positive
            Mat4DPoss(Mat4DPoss<0) = 0;
            SumMatPoss = sum(Mat4DPoss,4);
            SumMatPossStru = StruPoss;
            SumMatPossStru.dat.fname = [GenDdir,filesep,'Cumulative_T-ScoreMaps',filesep,'TScore_Cumulative_Positive_',condsName2{rr-subloop},'_map.nii'];
            SumMatPossStru.dat(:,:,:) = SumMatPoss;
            create(SumMatPossStru)
            
            Mat4DPossCon(Mat4DPossCon<0) = 0;
            SumMatPossCon = sum(Mat4DPossCon,4);
            SumMatPossStruCon = StruPossCon;
            SumMatPossStruCon.dat.fname = [GenDdir,filesep,'Cumulative_ContrastMaps',filesep,'ConMap_Cumulative_Positive_',condsName2{rr-subloop},'_map.nii'];
            SumMatPossStruCon.dat(:,:,:) = SumMatPossCon;
            create(SumMatPossStruCon)
            
            % Average Positive
            SumMatPoss2 = sum(Mat4DPoss,4)./size(Mat4DPoss,4);
            SumMatPossStru = StruPoss;
            SumMatPossStru.dat.fname = [GenDdir,filesep,'Average_T-ScoreMaps',filesep,'TScore_Average_Positive_',condsName2{rr-subloop},'_map.nii'];
            SumMatPossStru.dat(:,:,:) = SumMatPoss2;
            create(SumMatPossStru)
            
            SumMatPossCon2 = sum(Mat4DPossCon,4)./size(Mat4DPoss,4);
            SumMatPossStruCon = StruPossCon;
            SumMatPossStruCon.dat.fname = [GenDdir,filesep,'Average_ContrastMaps',filesep,'ConMap_Average_Positive_',condsName2{rr-subloop},'_map.nii'];
            SumMatPossStruCon.dat(:,:,:) = SumMatPossCon2;
            create(SumMatPossStruCon)

            
            clear Mat4DPoss VFrame4D v1 SumMatPossCon SumMatPossStruCon
            
            if rr+1>=10
                spmTprefix  = 'spmT_00';
                spmTprefix2  = 'con_00';
            else
                if rr+1>=100
                    spmTprefix  = 'spmT_0';
                    spmTprefix2  = 'con_0';
                else
                    spmTprefix  = 'spmT_000';
                    spmTprefix2  = 'con_000';
                end
            end
            
            if get(handles.checkbox5,'Value') %Slice view
%                 v2PC = VideoWriter([GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_ProgreCumul_Negative','.mp4'],'MPEG-4');
%                 v2PC.FrameRate = 0.25;
%                 v2PC.Quality = 100;
%                 open(v2PC)
                VFrame4D2 = imread([DirNameFF{1},filesep,SliceVNeg{rr-subloop,1}]);
%                 writeVideo(v2PC,VFrame4D2)
                [VFrame4D2,cm] = rgb2ind(VFrame4D2,256);
                imwrite(VFrame4D2,cm,...
                    [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_ProgreCumul_Negative','.gif'],...
                    'gif','DelayTime',1.5,'Loopcount',inf);
                
                Stru1 = nifti(nameNeg{rr-subloop,1});
                MatPC1 = Stru1.dat(:,:,:);
                MatPC1(isnan(MatPC1)) = 0;
                for j = 2:size(drifts,2)
                    StruTMP = nifti(nameNeg{rr-subloop,j});
                    MatTMP = StruTMP.dat(:,:,:);
                    MatTMP(isnan(MatTMP)) = 0;
                    MatPC1 = MatPC1 + MatTMP;
                    StruTMP.dat.fname = [GenDdir,filesep,'Temporary',filesep,'TMP.nii'];
                    StruTMP.dat(:,:,:) = MatPC1;
                    create(StruTMP)
                    SliceView_uf2c(file_T1,[GenDdir,filesep,'Temporary',...
                        filesep,'TMP.nii'],0.1,0.1,Orie,'winter',...
                        SlicesN,['Time Drift ',num2str(drifts(j))],...
                        'TMP.png',[GenDdir,filesep,'Temporary',filesep])
                    VFrame4Dtmp = imread([GenDdir,filesep,'Temporary',filesep,'TMP.png']);
                    [VFrame4Dtmp,cm] = rgb2ind(VFrame4Dtmp,256);
                    imwrite(VFrame4Dtmp,cm,...
                        [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_ProgreCumul_Negative','.gif'],...
                        'gif','DelayTime',1.5,'WriteMode','append');
%                     writeVideo(v2PC,VFrame4Dtmp)
                end
%                 close(v2PC)
                clear v2PC VFrame4D1 VFrame4Dtmp StruTMP MatTMP
                
                if FDRcInfe
%                     v2PC = VideoWriter([GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_ProgreCumul_Negative_FDRcInfer','.mp4'],'MPEG-4');
%                     v2PC.FrameRate = 0.25;
%                     v2PC.Quality = 100;
%                     open(v2PC)
                    VFrame4D2 = imread([DirNameFF{1},filesep,SliceVNeg2{rr-subloop,1}]);
%                     writeVideo(v2PC,VFrame4D2)
                    [VFrame4D2,cm] = rgb2ind(VFrame4D2,256);
                    imwrite(VFrame4D2,cm,...
                        [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_ProgreCumul_Negative_FDRcInfer','.gif'],...
                        'gif','DelayTime',1.5,'Loopcount',inf);
                    
                    Stru1 = nifti(nameNeg2{rr-subloop,1});
                    MatPC1 = Stru1.dat(:,:,:);
                    MatPC1(isnan(MatPC1)) = 0;
                    for j = 2:size(drifts,2)
                        StruTMP = nifti(nameNeg2{rr-subloop,j});
                        MatTMP = StruTMP.dat(:,:,:);
                        MatTMP(isnan(MatTMP)) = 0;
                        MatPC1 = MatPC1 + MatTMP;
                        StruTMP.dat.fname = [GenDdir,filesep,'Temporary',filesep,'TMP.nii'];
                        StruTMP.dat(:,:,:) = MatPC1;
                        create(StruTMP)
                        SliceView_uf2c(file_T1,[GenDdir,filesep,'Temporary',...
                            filesep,'TMP.nii'],0.1,0.1,Orie,'winter',...
                            SlicesN,['Time Drift ',num2str(drifts(j))],...
                            'TMP.png',[GenDdir,filesep,'Temporary',filesep])
                        VFrame4Dtmp = imread([GenDdir,filesep,'Temporary',filesep,'TMP.png']);
                        [VFrame4Dtmp,cm] = rgb2ind(VFrame4Dtmp,256);
                        imwrite(VFrame4Dtmp,cm,...
                            [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_ProgreCumul_Negative_FDRcInfer','.gif'],...
                            'gif','DelayTime',1.5,'WriteMode','append');
%                         writeVideo(v2PC,VFrame4Dtmp)
                    end
%                     close(v2PC)
                end
            end
            try
                rmdir([GenDdir,filesep,'Temporary'],'s')
            end
            
%             v2 = VideoWriter([GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_Negative','.mp4'],'MPEG-4');
%             v2.FrameRate = 0.25;
%             v2.Quality = 100;
%             open(v2)
            for i = 1:size(drifts,2)
                if get(handles.checkbox5,'Value') %Slice view
                    VFrame4D = imread([DirNameFF{i},filesep,SliceVNeg{rr-subloop,i}]);
                    [VFrame4D,cm] = rgb2ind(VFrame4D,256);
                    if i==1
                        imwrite(VFrame4D,cm,...
                            [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_Negative','.gif'],...
                            'gif','DelayTime',1.5,'Loopcount',inf);
                    else
                        imwrite(VFrame4D,cm,...
                           [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_Negative','.gif'],...
                            'gif','DelayTime',1.5,'WriteMode','append');
                    end
%                     writeVideo(v2,VFrame4D)
                end
                StruNegg = nifti([DirNameFF{i},filesep,spmTprefix,num2str(rr+1),'.nii']);
                Mat4DNegg(:,:,:,i) = StruNegg.dat(:,:,:);
                
                StruNeggCon = nifti([DirNameFF{i},filesep,spmTprefix2,num2str(rr+1),'.nii']);
                Mat4DNeggCon(:,:,:,i) = StruNeggCon.dat(:,:,:);
            end
%             if get(handles.checkbox5,'Value') %Slice view
%                 close(v2)
%             end
            
            if FDRcInfe
%                 v2 = VideoWriter([GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_Negative_FDRcInfer','.mp4'],'MPEG-4');
%                 v2.FrameRate = 0.25;
%                 v2.Quality = 100;
%                 open(v2)
                for i = 1:size(drifts,2)
                    if get(handles.checkbox5,'Value') %Slice view
                        VFrame4D = imread([DirNameFF{i},filesep,SliceVNeg2{rr-subloop,i}]);
                        [VFrame4D,cm] = rgb2ind(VFrame4D,256);
                        if i ==1
                            imwrite(VFrame4D,cm,...
                                [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_Negative_FDRcInfer','.gif'],...
                                'gif','DelayTime',1.5,'Loopcount',inf);
                        else
                            imwrite(VFrame4D,cm,...
                               [GenDdir,filesep,'Videos',filesep,condsName2{rr-subloop},'_Negative_FDRcInfer','.gif'],...
                                'gif','DelayTime',1.5,'WriteMode','append');
                        end
%                         writeVideo(v2,VFrame4D)
                    end
                end
%                 if get(handles.checkbox5,'Value') %Slice view
%                     close(v2)
%                 end
            end
            
            Mat4DNegg(Mat4DNegg<0) = 0;
            SumMatNegg = sum(Mat4DNegg,4);
            SumMatNeggStru = StruNegg;
            SumMatNeggStru.dat.fname = [GenDdir,filesep,'Cumulative_T-ScoreMaps',filesep,'TScore_Cumulative_Negative_',condsName2{rr-subloop},'_map.nii'];
            SumMatNeggStru.dat(:,:,:) = SumMatNegg;
            create(SumMatNeggStru)
            
            Mat4DNeggCon(Mat4DNeggCon<0) = 0;
            SumMatNeggCon = sum(Mat4DNeggCon,4);
            SumMatNeggStruCon = StruNeggCon;
            SumMatNeggStruCon.dat.fname = [GenDdir,filesep,'Cumulative_ContrastMaps',filesep,'ConMap_Cumulative_Negative_',condsName2{rr-subloop},'_map.nii'];
            SumMatNeggStruCon.dat(:,:,:) = SumMatNeggCon;
            create(SumMatNeggStruCon)

            SumMatNegg2 = sum(Mat4DNegg,4)./size(Mat4DNegg,4);
            SumMatNeggStru = StruNegg;
            SumMatNeggStru.dat.fname = [GenDdir,filesep,'Average_T-ScoreMaps',filesep,'TScore_Average_Negative_',condsName2{rr-subloop},'_map.nii'];
            SumMatNeggStru.dat(:,:,:) = SumMatNegg2;
            create(SumMatNeggStru)
            
            SumMatNeggCon2 = sum(Mat4DNeggCon,4)./size(Mat4DNeggCon,4);
            SumMatNeggStruCon = StruNeggCon;
            SumMatNeggStruCon.dat.fname = [GenDdir,filesep,'Average_ContrastMaps',filesep,'ConMap_Average_Negative_',condsName2{rr-subloop},'_map.nii'];
            SumMatNeggStruCon.dat(:,:,:) = SumMatNeggCon2;
            create(SumMatNeggStruCon)

            clear Mat4DNegg VFrame4D v2 SumMatNeggCon Mat4DNeggCon
            
            subloop = subloop +1;
        end
    end
end
fclose('all');
set(handles.text17,'String','Done!');
drawnow

function radiobutton4_Callback(hObject, eventdata, handles)
if get(handles.radiobutton4,'Value')
    set(handles.radiobutton3,'Value',0)
    set(handles.edit2,'String','0.001')
else
    set(handles.radiobutton3,'Value',1)
    set(handles.edit2,'String','0.05')
end
    
function radiobutton3_Callback(hObject, eventdata, handles)
if get(handles.radiobutton3,'Value')
    set(handles.radiobutton4,'Value',0)
    set(handles.edit2,'String','0.05')
else
    set(handles.radiobutton4,'Value',1)
    set(handles.edit2,'String','0.001')
end

function checkbox5_Callback(hObject, eventdata, handles)
if get(handles.radiobutton4,'Value')
    set(handles.text15,'Enable','on')
    set(handles.edit4,'Enable','on')
    set(handles.popupmenu1,'Enable','on')
else
    set(handles.text15,'Enable','off')
    set(handles.edit4,'Enable','off')
    set(handles.popupmenu1,'Enable','off')
end

function ATD_Callback(hObject, eventdata, handles)
if get(handles.ATD,'Value')
    set(handles.text19,'Enable','on')
    set(handles.TDv,'Enable','on')
else
    set(handles.text19,'Enable','off')
    set(handles.TDv,'Enable','off')
end

function markers_Callback(hObject, eventdata, handles)
function allasone_Callback(hObject, eventdata, handles)
function MarkDur_Callback(hObject, eventdata, handles)
if ~get(handles.MarkDur,'Value')
    set(handles.text2,'Enable','on')
    set(handles.edit1,'Enable','on')
else
    set(handles.text2,'Enable','off')
    set(handles.edit1,'Enable','off')
end

function edit1_Callback(hObject, eventdata, handles)
function markers_ButtonDownFcn(hObject, eventdata, handles)
function Conds_Callback(hObject, eventdata, handles)
function regres_Callback(hObject, eventdata, handles)
function popupmenu1_Callback(hObject, eventdata, handles)
function edit4_Callback(hObject, eventdata, handles)
function edit3_Callback(hObject, eventdata, handles)
function edit2_Callback(hObject, eventdata, handles)
function TDv_Callback(hObject, eventdata, handles)
function TRin_Callback(hObject, eventdata, handles)
function edit8_Callback(hObject, eventdata, handles)
function FDRclust_Callback(hObject, eventdata, handles)

function figure1_CloseRequestFcn(hObject, eventdata, handles)
clearvars -global Image4D MarkerF MarkerPath Path rpF rpPath struF struPath
delete(hObject);
function Conds_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function markers_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function regres_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function TRin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function TDv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



