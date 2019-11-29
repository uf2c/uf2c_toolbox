function varargout = uf2c_OutliersTool(varargin)
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
                   'gui_OpeningFcn', @uf2c_OutliersTool_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_OutliersTool_OutputFcn, ...
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

function uf2c_OutliersTool_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = uf2c_OutliersTool_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function addnumb_Callback(hObject, eventdata, handles)
global inputs

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='none';

inputs = inputdlg('Input numerical values','',5,{''},options);
inputs = str2num(cell2mat(inputs));
set(handles.text3,'String',sprintf('%d inputs added',numel(inputs)))

function runB2_Callback(hObject, eventdata, handles)
global inputs

set(handles.status,'String','Running...')
drawnow


switch get(handles.listbox1,'Value')
    case 1
        OutLierMode = 'severe';
    case 2 
        OutLierMode = 'minor';
    case 3
        OutLierMode = 'major';
end

switch OutLierMode
    case 'severe'
        whisker = 1;
    case 'minor'
        whisker = 1.5;
    case 'major'
        whisker = 3;
end

ar = uf2c_outlierdetec(inputs,OutLierMode);

if ~isempty(ar)
    disp('Outliers:')
    for i3 = 1:numel(ar)
        fprintf('input %d; value %d\n',ar(i3),inputs(ar(i3)))
    end
    fprintf('\n')
else
    disp('No outliers found')
    set(handles.status,'String','Done!')
end


figure;
boxplot(inputs,'Whisker',whisker);
ylabel('VALUES','FontSize',14);
mTextBox = uicontrol('style','text','HorizontalAlignment','left','position',round([1 1 450 30]));
set(mTextBox,'String',{'Click on the outliers to check the filenames. Right Click removes the label.';'Press ''Enter'' to close this tool'});
set(handles.status,'String','Done!')
drawnow
gname(inputs)



function addfileb_Callback(hObject, eventdata, handles)
global ReshMat1 matF filename ext Dx Dy Dz  Dv stru file2 typein a
 
file = {};
file = uipickfiles('Output','cell','Prompt',...
    'Add NIfTI files OR Matlab files with "double" precision variables',...
    'Type', {'*.nii','NIfTI';'*.mat','MAT-files'});

if isequal(file,0)
    return
end

set(handles.text3,'String','WAIT! Processing inputs...(can take a while)','ForegroundColor',[1 0 0])
drawnow

if ~iscell(file)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
    file = {file};
end

nOFsubjects = size(file,2);

if isequal(nOFsubjects,1)
    set(handles.text4,'String','Consider only values','ForegroundColor',[1 0 0])
    set(handles.threPop,'String',{'higher then:';'lower then:';'distinct of:'})
else
    set(handles.text4,'String','Turn it to ZERO values','ForegroundColor',[0 0 1])
    set(handles.threPop,'String','lower then:')
end

for yt = 1:nOFsubjects
    tmpDir = dir([file{yt}]);
    bitS(yt) = tmpDir.bytes;
end

sxs = unique(bitS);
nDiftypes = size(sxs,2);
if ~isequal(nDiftypes,1) && mean(diff(sxs))>5
    warndlg({sprintf('You added %d distinct types of files.',nDiftypes);...
    sprintf('Check the size of your files, it is the easiest way to identify the if they are equivalent')},'Attention!');
    set(handles.text3,'String','Aborted','ForegroundColor',[1 0 0])
    drawnow
    return
end
clear bitS

file2 = {};
file2 = file;

typein = '';

[pa,name,ext] = fileparts(file2{1});

filename = {};
ReshMat1 = [];
matF = [];
switch ext
    case '.nii'
        if isequal(nOFsubjects,1)
            set(handles.expmask,'Enable','off')
            set(handles.expmask,'Value',0)
            set(handles.woexpmask,'Enable','off')
            set(handles.woexpmask,'Value',0)
            set(handles.onlyshape,'Value',0)
            set(handles.onlyshape,'Enable','off')
            set(handles.checkbox7,'Enable','on')
            set(handles.checkbox7,'Value',1)
        else
            set(handles.expmask,'Enable','on')
            set(handles.expmask,'Value',1)
            set(handles.woexpmask,'Enable','on')
            set(handles.woexpmask,'Value',1)
            set(handles.onlyshape,'Value',1)
            set(handles.onlyshape,'Enable','on')
            set(handles.checkbox7,'Enable','on')
            set(handles.checkbox7,'Value',0)
        end
        set(handles.threshold,'Enable','on')
        set(handles.threshold,'String','')
        
        h = waitbar(0,'WAIT! Processing inputs...');
        
        stru = nifti(file2{1});
        if isequal(numel(stru.dat.dim),3)
            for j = 1:size(file2,2)
                stru = nifti(file2{j});
                mattmp = stru.dat(:,:,:);
                matF(:,:,:,j) = mattmp;
                [a,filename{j},c] = fileparts(file2{j});
                waitbar(j/size(file2,2),h)
            end
            typein = 'nii3d';
            set(handles.text3,'String',sprintf('Done! %d file(s) with a 3D NIfTI image was/were added.',size(file2,2)),'ForegroundColor',[0 0 0])
        end
        if isequal(numel(stru.dat.dim),4)
            disp('4D file(s): calculating and using the 3D median for intra or inter outliers computations.')
            for j = 1:size(file2,2)
                stru = nifti(file2{j});
                mattmp = stru.dat(:,:,:,:);
                matF(:,:,:,j) = median(mattmp,4);
                [a,filename{j},c] = fileparts(file2{j});
                waitbar(j/size(file2,2),h)
            end
            typein = 'nii4d';
            set(handles.text3,'String',sprintf('Done! %d file(s) with a 4D NIfTI image was/were added.',size(file2,2)),'ForegroundColor',[0 0 0])
        end
        Dx = size(matF,1);
        Dy = size(matF,2);
        Dz = size(matF,3);
        Dv = size(matF,4);
        
        ReshMat1 = reshape(matF,prod([Dx,Dy,Dz]),Dv)';
        close(h)
    
    case '.mat'
        set(handles.expmask,'Enable','off')
        set(handles.expmask,'Value',0)
        set(handles.woexpmask,'Enable','off')
        set(handles.woexpmask,'Value',0)
        set(handles.onlyshape,'Value',0)
        set(handles.onlyshape,'Enable','off')
        set(handles.checkbox7,'Value',1)
        set(handles.checkbox7,'Enable','on')
        
        a = load(file2{1});
        value = fields(a);
        value = value{1};
        eval(sprintf('%s = a.%s;',value,value));
        Trans = eval(value);
        
        if isvector(Trans)
            for j = 1:size(file2,2)
                a = load(file2{j});
                value = fields(a);
                value = value{1};
                eval(sprintf('%s = a.%s;',value,value));
                Trans = eval(value);
                matF = Trans;
                ReshMat1(j,:) = Trans;
                [a,filename{j},c] = fileparts(file2{j});
            end
            typein = 'mat1d';
            set(handles.text3,'String',sprintf('Done! %d file(s) with a double precision vector was/were added.',size(file2,2)),'ForegroundColor',[0 0 0])
        end
        
        if ~isvector(Trans) && isequal(numel(size(Trans)),2)
            for j = 1:size(file2,2)
                a = load(file2{j});
                value = fields(a);
                value = value{1};
                eval(sprintf('%s = a.%s;',value,value));
                Trans = eval(value);
%                     Trans(Trans<Theshold) = 0;
                matF(:,:,j) = Trans;
                [a,filename{j},c] = fileparts(file2{j});
            end
            typein = 'mat2d';
            set(handles.text3,'String',sprintf('Done! %d file(s) with a 2D double precision matrix was/were added.',size(file2,2)),'ForegroundColor',[0 0 0])

            Dx = size(matF,1);
            Dy = size(matF,2);
            Dv = size(matF,3);

            ReshMat1 = reshape(matF,prod([Dx,Dy]),Dv)';
        end
        if isequal(numel(size(Trans)),3)
            for j = 1:size(file2,2)
                a = load(file2{j});
                value = fields(a);
                value = value{1};
                eval(sprintf('%s = a.%s;',value,value));
                Trans = eval(value);
%                     Trans(Trans<Theshold) = 0;
                matF(:,:,:,j) = Trans;
                [a,filename{j},c] = fileparts(file2{j});
            end
            typein = 'mat3d';
            set(handles.text3,'String',sprintf('Done! %d file(s) with a 3D double precision matrix was/were added.',size(file2,2)),'ForegroundColor',[0 0 0])

            Dx = size(matF,1);
            Dy = size(matF,2);
            Dz = size(matF,3);
            Dv = size(matF,4);

            ReshMat1 = reshape(matF,prod([Dx,Dy,Dz]),Dv)';

        end
        if isequal(numel(size(Trans)),4)
            disp('4D file(s): calculating and using the 3D median for intra or inter outliers computations.')
            for j = 1:size(file2,2)
                a = load(file2{j});
                value = fields(a);
                value = value{1};
                eval(sprintf('%s = a.%s;',value,value));
                Trans = eval(value);
%                     Trans(Trans<Theshold) = 0;
                matF(:,:,:,j) = median(Trans,4);
                [a,filename{j},c] = fileparts(file2{j});
            end
            typein = 'mat4d';
            set(handles.text3,'String',sprintf('Done! %d file(s) with a 4D double precision matrix was/were added.',size(file2,2)),'ForegroundColor',[0 0 0])

            Dx = size(matF,1);
            Dy = size(matF,2);
            Dz = size(matF,3);
            Dv = size(matF,4);
            
            ReshMat1 = reshape(matF,prod([Dx,Dy,Dz]),Dv)';
        end
end

function runB_Callback(hObject, eventdata, handles)
global ReshMat1 filename ext Dx Dy Dz stru file2 typein matF a

set(handles.status,'String','Running...')
drawnow

ScreSize = get(0,'screensize');
ScreSize = ScreSize(3:end);


threshold = str2num(get(handles.threshold,'String'));

if get(handles.OutDetec,'Value')
    ReshMat = ReshMat1;

    switch get(handles.listbox1,'Value')
        case 1
            OutLierMode = 'severe';
        case 2 
            OutLierMode = 'minor';
        case 3
            OutLierMode = 'major';
    end

    switch OutLierMode
        case 'severe'
            whisker = 1;
        case 'minor'
            whisker = 1.5;
        case 'major'
            whisker = 3;
    end

    finalOut = [];
    MatOut = zeros(size(ReshMat,1),size(ReshMat,2));

    if isequal(size(ReshMat,1),1) % case user added only one file (intra-outliers will be search)
        
        fprintf('You input only one sample. Estimating intra-file outliers\n');
        threshold = str2num(get(handles.threshold,'String'));
        
        if ~isempty(threshold)
            switch get(handles.threPop,'Value')
                case 1
                    [yty,tyy] = find(ReshMat>threshold);
                case 2
                    [yty,tyy] = find(ReshMat<threshold);
                case 3
                    [yty,tyy] = find(ReshMat~=threshold);
            end
            
            ReshMattmp = ReshMat(tyy);
            
            ReshMattmp = ReshMattmp./max(ReshMattmp);
            
            ar = uf2c_outlierdetec(ReshMattmp,OutLierMode);
            artmp = tyy(ar);
            s2 = ReshMattmp;
            OutInputs = artmp;
        else
            ReshMat = ReshMat./max(ReshMat);
            ar = uf2c_outlierdetec(ReshMat,OutLierMode);
            s2 = ReshMat(1,:);
            OutInputs = ar;
        end
        
        for gg = 1:size(ReshMat,2)
            switch typein
                case {'mat1d','mat2d'}
                    [Xi,Yi] = ind2sub(size(matF),gg);
                    filename2{gg} = ['[' num2str(Xi) ',' num2str(Yi) ']'];
                case {'nii3d','nii4d','mat3d','mat4d'}
                    [Xi,Yi,Zi] = ind2sub(size(matF),gg);
                    filename2{gg} = ['[' num2str(Xi) ',' num2str(Yi) ',' num2str(Zi) ']'];
            end
        end

        if ~isempty(OutInputs)
            if ~isequal(ext,'.nii')
                disp('Outliers:')
                for i3 = 1:numel(OutInputs)
                    fprintf('index: %s\n',filename2{OutInputs(i3)})
                end
                fprintf('\n')
            end
        else
            disp('No intra-file outliers found')
            set(handles.status,'String','Done!')
        end
        
    else % case user added more than one file
        
        if ~isempty(threshold)
            ReshMat(ReshMat<threshold) = 0;
        end
        
        tmpVet = 0;
        if get(handles.onlyshape,'Value')
            for u2 = 1:size(ReshMat,2)
                tmpVet = ReshMat(:,u2);
                tmpVet(tmpVet~=0) = 1;
                ReshMat(:,u2) = tmpVet;
            end
        else
            for u = 1:size(ReshMat,1)
                ReshMat(u,:) = ReshMat(u,:)./max(ReshMat(u,:));
            end
        end
        
        for j = 1:size(ReshMat,2)
            if ~isequal(sum(ReshMat(:,j)),0)
                a2 = uf2c_outlierdetec(ReshMat(:,j),OutLierMode,'twosides');
                MatOut(1:numel(a2),j) = a2;
                clear a2
            end
        end
        
        volsoccur = 0;
        volsoccur = unique(MatOut);
        volsoccur = volsoccur(2:end);

        s2 = [];
        inds = [];
        if ~isempty(volsoccur)
            s = [];
            for k = 1:size(volsoccur,1)
                s(k) = sum(sum(volsoccur(k)==MatOut));
            end

            [s2,inds] = sort(s,'descend');
            volsoccur = volsoccur(inds);
            filename2 = filename(volsoccur);

            OutInputs = 0;
            OutInputs = uf2c_outlierdetec(s2,OutLierMode,'upperside');
            OutInputs = volsoccur(OutInputs);
        else
            OutInputs = [];
        end

        if ~isempty(OutInputs)
            disp('Outliers:')
            for i3 = 1:numel(OutInputs)
                fprintf('%s: due %d outlier(s) point(s) or voxel(s)\n',filename{OutInputs(i3)},s2(i3))
            end
            fprintf('\n')
        else
            disp('No inter-files outliers found')
            set(handles.status,'String','Done!')
        end
    end

    if ~isempty(OutInputs)
        switch typein
            case {'nii3d','nii4d'}
                mkdir([a filesep '1-Outliers']);
                disp('Creating:')
                if isequal(size(ReshMat,1),1)
                    bin = zeros(1,size(ReshMat,2));
                    bin(OutInputs) = 1;
                    BinMat = reshape(bin',[Dx,Dy,Dz]);
                    [dp,dn,de] = fileparts(file2{1});
                    
                    strutmp = stru;
                    strutmp.dat.fname = [dp filesep '1-Outliers' filesep filename{1} '__OutilierVoxels.nii'];
                    fprintf('%s\n',[filename{1} '__OutilierVoxels.nii'])
                    strutmp.dat.dim = [Dx,Dy,Dz];
                    strutmp.dat.dtype = 'FLOAT32-LE';
                    strutmp.dat(:,:,:) = BinMat;
                    create(strutmp)
                else
                    for l = 1:numel(OutInputs)
                        bin = MatOut==OutInputs(l);
                        bin = sum(bin,1);
                        BinMat = reshape(bin',[Dx,Dy,Dz]);
                        
                        [dp,dn,de] = fileparts(file2{OutInputs(l)});
                        
                        strutmp = stru;
                        strutmp.dat.fname = [a filesep '1-Outliers' filesep filename{OutInputs(l)} '__OutilierVoxels.nii'];
                        fprintf('%s\n',[filename{OutInputs(l)} '__OutilierVoxels.nii'])
                        strutmp.dat.dim = [Dx,Dy,Dz];
                        strutmp.dat.dtype = 'FLOAT32-LE';
                        strutmp.dat(:,:,:) = BinMat;
                        create(strutmp)
                    end
                end
            case 'mat1d'
                if isequal(size(ReshMat,1),1)
                    s2 = ReshMat;
                else

                end
            case 'mat2d'
                if isequal(size(ReshMat,1),1)
                    s2 = ReshMat;
                else
                    mkdir([a filesep '1-Outliers']);
                    for l = 1:numel(OutInputs)
                        bin = MatOut==OutInputs(l);
                        bin = sum(bin,1);
                        BinMat = reshape(bin',[Dx,Dy]);
                        [dp,dn,de] = fileparts(file2{OutInputs(l)});
                        save([dp filesep '1-Outliers' filesep filename{OutInputs(l)} '__OutilierVoxels.mat'],'BinMat')
                    end
                end
            case 'mat3d'
                if isequal(size(ReshMat,1),1)
                    s2 = ReshMat;
                else
                    mkdir([a filesep '1-Outliers']);
                    for l = 1:numel(OutInputs)
                        bin = MatOut==OutInputs(l);
                        bin = sum(bin,1);
                        BinMat = reshape(bin',[Dx,Dy,Dz]);
                        [dp,dn,de] = fileparts(file2{OutInputs(l)});
                        save([dp filesep '1-Outliers' filesep filename{OutInputs(l)} '__OutilierVoxels.mat'],'BinMat')
                    end
                end
            case 'mat4d'
                if isequal(size(ReshMat,1),1)
                    s2 = ReshMat;
                else
                    mkdir([a filesep '1-Outliers']);
                    for l = 1:numel(OutInputs)
                        bin = MatOut==OutInputs(l);
                        bin = sum(bin,1);
                        BinMat = reshape(bin',[Dx,Dy,Dz]);
                        [dp,dn,de] = fileparts(file2{OutInputs(l)});
                        save([dp filesep '1-Outliers' filesep filename{OutInputs(l)} '__OutilierVoxels.mat'],'BinMat')
                    end
                end
        end
    end
end

if get(handles.expmask,'Value')

    spm_jobman('initcfg')
    spm('defaults','fmri');
    
    mkdir([a filesep '1-Masks']);

    if isequal(typein,'nii4d') || isequal(typein,'nii3d')
        matF2 = matF;
        struExp = stru;
        if isequal(typein,'nii4d')
            struExp.dat.fname = [a filesep '1-Masks' filesep '4D_Medians.nii'];
        else
            struExp.dat.fname = [a filesep '1-Masks' filesep '4D_Files.nii'];
        end
        
        if ~isempty(threshold)
            matF2(matF2<threshold) = 0;
        end
        struExp.dat.dim = size(matF2);
        struExp.dat(:,:,:,:) = matF2;
        create(struExp)

        file11 = cell(struExp.dat.dim(4),1);
        for tt = 1:struExp.dat.dim(4)
            if isequal(typein,'nii4d')
                file11{tt,1} = [a filesep '1-Masks' filesep '4D_Medians.nii,',sprintf('%d',tt)];
            else
                file11{tt,1} = [a filesep '1-Masks' filesep '4D_files.nii,',sprintf('%d',tt)];
            end
        end
    else
        file11 = file2;
    end

    [~,~] = uf2c_avgimg(file2,threshold,'PREPOS',[a filesep '1-Masks' filesep],'TotalAverage.nii',1);
    
    matlabbatch{1}.spm.stats.factorial_design.dir = {[a filesep '1-Masks']};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = file11;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tma.athresh = 0.01;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'OneSampleTTest';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 0;

    matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
    matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
    matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
    matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
    matlabbatch{4}.spm.stats.results.conspec.extent = 0;
    matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
    matlabbatch{4}.spm.stats.results.units = 1;
    matlabbatch{4}.spm.stats.results.print = 'tif';
    matlabbatch{4}.spm.stats.results.write.tspm.basename = 'Mask_Original';

    matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
    matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
    matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'FWE';
    matlabbatch{5}.spm.stats.results.conspec.thresh = 0.05;
    matlabbatch{5}.spm.stats.results.conspec.extent = 0;
    matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
    matlabbatch{5}.spm.stats.results.units = 1;
    matlabbatch{5}.spm.stats.results.print = 'tif';
    matlabbatch{5}.spm.stats.results.write.binary.basename = 'BinaryMask_Original';
    try
        spm_jobman('run',matlabbatch)
        delete([a filesep '1-Masks' filesep 'beta_0001.nii']);
        delete([a filesep '1-Masks' filesep 'con_0001.nii']);
        delete([a filesep '1-Masks' filesep 'mask.nii']);
        delete([a filesep '1-Masks' filesep 'ResMS.nii']);
        delete([a filesep '1-Masks' filesep 'RPV.nii']);
        delete([a filesep '1-Masks' filesep 'spmT_0001.nii']);
        delete([a filesep '1-Masks' filesep 'SPM.mat']);
        movefile([a filesep '1-Masks' filesep 'spmT_0001_Mask_Original.nii'],[a filesep '1-Masks' filesep 'Mask_Original.nii'])
        movefile([a filesep '1-Masks' filesep 'spmT_0001_BinaryMask_Original.nii'],[a filesep '1-Masks' filesep 'BinaryMask_Original.nii'])

        tmpdirn =  dir('*.tif');
        movefile([a filesep '1-Masks' filesep tmpdirn(1).name],[a filesep '1-Masks' filesep 'Mask_Original.tif'])
        delete([a filesep '1-Masks' filesep tmpdirn(2).name])
    end
end
    
if get(handles.woexpmask,'Value')
    if ~isempty(OutInputs)
        
        spm_jobman('initcfg')
        spm('defaults','fmri');
        mkdir([a filesep '1-Masks']);
        
        if isequal(typein,'nii4d')  || isequal(typein,'nii3d')
            matF3 = matF;
            matF3(:,:,:,OutInputs) = [];
            struExp2 = stru;
            if isequal(typein,'nii4d')
                struExp2.dat.fname = [a filesep '1-Masks' filesep '4D_Medians_NoOutliers.nii'];
            else
                struExp2.dat.fname = [a filesep '1-Masks' filesep '4D_files_NoOutliers.nii'];
            end
            if ~isempty(threshold)
                matF3(matF3<threshold) = 0;
            end
            struExp2.dat.dim = size(matF3);
            struExp2.dat(:,:,:,:) = matF3;
            create(struExp2)

            file12 = cell(struExp2.dat.dim(4),1);
            for tt = 1:struExp2.dat.dim(4)
                if isequal(typein,'nii4d')
                    file12{tt,1} = [a filesep '1-Masks' filesep '4D_Medians_NoOutliers.nii,',sprintf('%d',tt)];
                else
                    file12{tt,1} = [a filesep '1-Masks' filesep '4D_Files_NoOutliers.nii,',sprintf('%d',tt)];
                end
            end
        else
            file12 = file2;
        end

        matlabbatch{1}.spm.stats.factorial_design.dir = {[a filesep '1-Masks']};
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = file12;
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tma.athresh = 0.01;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'OneSampleTTest';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.delete = 0;

        matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
        matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
        matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
        matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
        matlabbatch{4}.spm.stats.results.conspec.extent = 0;
        matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
        matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
        matlabbatch{4}.spm.stats.results.units = 1;
        matlabbatch{4}.spm.stats.results.print = 'tif';
        matlabbatch{4}.spm.stats.results.write.tspm.basename = 'Mask_OutlierRemmoved';

        matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
        matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
        matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'FWE';
        matlabbatch{5}.spm.stats.results.conspec.thresh = 0.05;
        matlabbatch{5}.spm.stats.results.conspec.extent = 0;
        matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
        matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
        matlabbatch{5}.spm.stats.results.units = 1;
        matlabbatch{5}.spm.stats.results.print = 'tif';
        matlabbatch{5}.spm.stats.results.write.binary.basename = 'BinaryMask_OutlierRemmoved';
        spm_jobman('run',matlabbatch)
        
        delete([a filesep '1-Masks' filesep 'beta_0001.nii']);
        delete([a filesep '1-Masks' filesep 'con_0001.nii']);
        delete([a filesep '1-Masks' filesep 'mask.nii']);
        delete([a filesep '1-Masks' filesep 'ResMS.nii']);
        delete([a filesep '1-Masks' filesep 'RPV.nii']);
        delete([a filesep '1-Masks' filesep 'spmT_0001.nii']);
        delete([a filesep '1-Masks' filesep 'SPM.mat']);
        movefile([a filesep '1-Masks' filesep 'spmT_0001_Mask_OutlierRemmoved.nii'],[a filesep '1-Masks' filesep 'Mask_OutlierRemmoved.nii'])
        movefile([a filesep '1-Masks' filesep 'spmT_0001_BinaryMask_OutlierRemmoved.nii'],[a filesep '1-Masks' filesep 'BinaryMask_OutlierRemmoved.nii'])
        
        tmpdirn =  dir('*.tif');
        movefile([a filesep '1-Masks' filesep tmpdirn(2).name],[a filesep '1-Masks' filesep 'Mask_OutlierRemmoved.tif'])
        delete([a filesep '1-Masks' filesep tmpdirn(3).name])
    else
        disp('No outliers to compute the Explicit Mask without outliers')
    end
end

    
set(handles.status,'String','Done!')

if get(handles.OutDetec,'Value')
    if isequal(size(ReshMat,1),1)
        if ~isequal(ext,'.nii')
            boxfig = figure;
            set(boxfig,'Name','Data Box Plot',...
                'Position', round([ScreSize(2)*.1 ScreSize(2)*.1 ScreSize(2)*.4 ScreSize(2)*.7]),...
                'Color',[1 0.94 0.86]);
            
            boxplot(s2,'Whisker',whisker);
            ylabel('HIGH <--- QUALITY ---> POOR','FontSize',14);
            mTextBox = uicontrol('style','text','HorizontalAlignment','left','position',round([1 1 450 30]));
            set(mTextBox,'String',{'Click on the outliers to check the filenames. Right Click removes the label.';'Press ''Enter'' to close this tool'});
            gname(filename2)
        end
    else
        if ~isempty(volsoccur)
            boxfig = figure;
            set(boxfig,'Name','Data Box Plot',...
                'Position', round([ScreSize(2)*.1 ScreSize(2)*.1 ScreSize(2)*.4 ScreSize(2)*.7]),...
                'Color',[1 0.94 0.86]);
            
            boxplot(s2,'Whisker',whisker);
            ylabel('HIGH <--- QUALITY ---> POOR','FontSize',14);
            mTextBox = uicontrol('style','text','HorizontalAlignment','left','position',round([1 1 450 30]));
            set(mTextBox,'String',{'Click on the outliers to check the filenames. Right Click removes the label.';'Press ''Enter'' to close this tool'});
            gname(filename2)
        end
    end
end

function checkbox7_Callback(hObject, eventdata, handles)
function woexpmask_Callback(hObject, eventdata, handles)
function onlyshape_Callback(hObject, eventdata, handles)
function threshold_Callback(hObject, eventdata, handles)
function expmask_Callback(hObject, eventdata, handles)
function threPop_Callback(hObject, eventdata, handles)
function listbox1_Callback(hObject, eventdata, handles)

function cbnumb_Callback(hObject, eventdata, handles)

if get(handles.cbnumb,'Value')
    set(handles.addnumb,'Enable','on')
    set(handles.addfileb,'Enable','off')
    set(handles.cbfiles,'Value',0)
    set(handles.expmask,'Value',0)
    set(handles.expmask,'Enable','off')
    set(handles.woexpmask,'Value',0)
    set(handles.woexpmask,'Enable','off')
    set(handles.onlyshape,'Enable','off')
    set(handles.onlyshape,'Value',0)
    set(handles.checkbox7,'Enable','off')
    set(handles.checkbox7,'Value',1)
    set(handles.threshold,'Enable','off')
    set(handles.threshold,'String','')
    set(handles.threPop,'Enable','off')
    set(handles.runB,'Visible','off')
    set(handles.runB2,'Visible','on')
    set(handles.OutDetec,'Value',1)
    set(handles.listbox1,'Enable','on')
else
    set(handles.addnumb,'Enable','off')
    set(handles.addfileb,'Enable','on')
    set(handles.cbfiles,'Value',1)
    set(handles.expmask,'Value',1)
    set(handles.expmask,'Enable','on')
    set(handles.woexpmask,'Value',1)
    set(handles.woexpmask,'Enable','on')
    set(handles.threshold,'Enable','on')
    set(handles.threshold,'String','')
    set(handles.onlyshape,'Enable','on')
    set(handles.onlyshape,'Value',1)
    set(handles.checkbox7,'Enable','on')
    set(handles.checkbox7,'Value',0)
    set(handles.threPop,'Enable','on')
    set(handles.runB,'Visible','on')
    set(handles.runB2,'Visible','off')

end

function cbfiles_Callback(hObject, eventdata, handles)
if get(handles.cbfiles,'Value')
    set(handles.expmask,'Enable','on')
    set(handles.expmask,'Value',1)
    set(handles.runB,'Visible','on')
    set(handles.runB2,'Visible','off')
    set(handles.woexpmask,'Enable','on')
    set(handles.woexpmask,'Value',1)
    set(handles.text4,'Enable','on')
    set(handles.onlyshape,'Value',1)
    set(handles.onlyshape,'Enable','on')
    set(handles.checkbox7,'Enable','on')
    set(handles.checkbox7,'Value',0)
    set(handles.threshold,'Enable','on')
    set(handles.threshold,'String','')
    set(handles.cbnumb,'Value',0)
    set(handles.addnumb,'Enable','off')
    set(handles.threPop,'Enable','on')
else
    set(handles.expmask,'Enable','off')
    set(handles.expmask,'Value',0)
    set(handles.runB,'Visible','off')
    set(handles.runB2,'Visible','on')
    set(handles.woexpmask,'Enable','off')
    set(handles.woexpmask,'Value',0)
    set(handles.text4,'Enable','off')
    set(handles.onlyshape,'Value',0)
    set(handles.onlyshape,'Enable','off')
    set(handles.checkbox7,'Enable','off')
    set(handles.checkbox7,'Value',1)
    set(handles.threshold,'Enable','off')
    set(handles.threshold,'String','')
    set(handles.cbnumb,'Value',1)
    set(handles.addnumb,'Enable','on')
    set(handles.threPop,'Enable','off')
    set(handles.OutDetec,'Value',1)
    set(handles.listbox1,'Enable','on')
end

function OutDetec_Callback(hObject, eventdata, handles)

if get(handles.OutDetec,'Value')
    if get(handles.cbnumb,'Value')
        set(handles.OutDetec,'Value',1)
        set(handles.listbox1,'Enable','on')
    else
        set(handles.woexpmask,'Value',1)
        set(handles.woexpmask,'Enable','on')
        set(handles.listbox1,'Enable','on')
        set(handles.onlyshape,'Enable','on')
        set(handles.onlyshape,'Value',1)
        set(handles.checkbox7,'Enable','on')
        set(handles.checkbox7,'Value',0)
    end
else
    if get(handles.cbnumb,'Value')
        set(handles.OutDetec,'Value',1)
        set(handles.listbox1,'Enable','on')
    else
        set(handles.OutDetec,'Value',0)
        set(handles.woexpmask,'Value',0)
        set(handles.woexpmask,'Enable','off')
        set(handles.listbox1,'Enable','off')
        set(handles.onlyshape,'Enable','off')
        set(handles.onlyshape,'Value',0)
        set(handles.checkbox7,'Enable','off')
        set(handles.checkbox7,'Value',1)
    end
end

function threshold_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function threPop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
