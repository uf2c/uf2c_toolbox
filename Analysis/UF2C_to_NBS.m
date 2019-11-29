function varargout = uf2c_to_NBS(varargin)
% UF²C M-file for uf2c_to_NBS.fig
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
                   'gui_OpeningFcn', @uf2c_to_NBS_OpeningFcn, ...
                   'gui_OutputFcn',  @uf2c_to_NBS_OutputFcn, ...
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

function uf2c_to_NBS_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = uf2c_to_NBS_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function push1_Callback(hObject, eventdata, handles)
global g1_VAR g2OK g1OK dirOK nodeOK

[fileG1,pathG1] = uigetfile({'*.mat','Matlab file'},...
        'Select the result file (All_Subjs-VAR.mat) for group 1','MultiSelect','off');
if ~isequal(fileG1,0)
    load([pathG1 fileG1])
    try
        g1_VAR = AllSubj3D;
    catch
        g1_VAR = TotalPairWise;
    end
    set(handles.text4,'String',sprintf('Data from %d subjects were added',size(g1_VAR,3)))
    g1OK = 1;
    if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
        set(handles.runB,'Enable','on')
    end
end

function push2_Callback(hObject, eventdata, handles)
global g2_VAR g2OK g1OK dirOK nodeOK

[fileG2,pathG2] = uigetfile({'*.mat','Matlab file'},...
        'Select the result file (All_Subjs-VAR.mat) for group 2','MultiSelect','off');
if ~isequal(fileG2,0)
    load([pathG2 fileG2])
    try
        g2_VAR = AllSubj3D;
    catch
        g2_VAR = TotalPairWise;
    end
    set(handles.text5,'String',sprintf('Data from %d subjects were added',size(g2_VAR,3)))
    g2OK = 1;
    if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
        set(handles.runB,'Enable','on')
    end
end

function COV_Callback(hObject, eventdata, handles)
global COVDATA nOfLines nOfCovs g2OK g1OK dirOK nodeOK g1_VAR g2_VAR

if get(handles.COV,'Value')
    [COVlist,COVpath] = uigetfile({'*.txt','Text files';'*.mat','3D mat file'},'Select a text file or a 3D mat (ROI-to-ROI cov) file with all covariates','MultiSelect','off');
    if ~isequal(COVlist,0)
        COVDATA = importdata([COVpath COVlist]);
        nOfLines = size(COVDATA,1);
        nOfCovs = size(COVDATA,2);
        if isequal(g1OK,1) && isequal(g2OK,1)
            if isequal((size(g1_VAR,3)+size(g2_VAR,3)),nOfLines)
                set(handles.text9,'String',sprintf('%d covariate(s) added!',nOfCovs))
            else
                warndlg('The number of inputs in each covariate column should to match with the sum of Group 1 and 2 subjects','Ops!')
                set(hanldes.COV,'Value',0)
                set(hanldes.text9,'String','')
            end
        else
            set(handles.text9,'String',sprintf('%d covariate(s) added!',nOfCovs))
        end
        if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
            set(handles.runB,'Enable','on')
        end
    else
        set(hanldes.COV,'Value',0)
        set(hanldes.text9,'String','')
    end
else
    set(hanldes.text9,'String','')
end

function OutDir_Callback(hObject, eventdata, handles)
global Dir g2OK g1OK  dirOK nodeOK

Dir = uigetdir('','Define the output directory');
if ~isequal(Dir,0)
    set(handles.text10,'String',Dir)
    drawnow
    Dir = ([Dir filesep 'NBS_Inputs' filesep]);
    dirOK = 1;
    if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
        set(handles.runB,'Enable','on')
    end
end

function addNODES_Callback(hObject, eventdata, handles)
global CoordMS seedList nodeOK g1OK g2OK dirOK nEle nOfNets fileROI EROI EList

if get(handles.CoodList,'Value')
    set(handles.text12,'String','Wait.......')
    drawnow
    EROI = 0;
    EList = 1;

    [fileCor,pathCor] = uigetfile({'*.txt','Text file'},...
            'Select the cordinates list file','MultiSelect','off');
        if ~isequal(fileCor,0)
            seedListT = importdata([pathCor fileCor]);
            nOfNets = str2num(seedListT{end,1}(1:2));
            xs=1;
            nE = 1;
            while xs <= size(seedListT,1)
                xStr = seedListT{xs,1}(1:3);
                net1 = strncmp(xStr,seedListT,3);
                nEle(nE) = sum(net1);
                xs = xs + sum(net1);
                nE = nE+1;
            end
            for i = 1:size(seedListT,1)
                seedList(i,:) = str2num(seedListT{i}(4:end));
            end
            RootDir = which('uf2c');
            RootDir = ([RootDir(1:end-6) 'Image_Examples' filesep]);
            for jh = 1:size(seedList,1)
                Vhdr = spm_vol([RootDir,'FiltRegrSW_Example.nii']);
                TcM = Vhdr.mat;
                EPIcoord = [seedList(jh,1) seedList(jh,2) seedList(jh,3) 1]*(inv(TcM))';
                EPIcoord = round(EPIcoord);
                CoordMS(jh,:) = EPIcoord(1,:);
            end
            set(handles.text12,'String',sprintf('%d coordinates of %d netorworks were added!',size(CoordMS,1),nOfNets))
            nodeOK = 1;
            if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
                set(handles.runB,'Enable','on')
            end
        end
else
    EROI = 1;
    EList = 0;
    set(handles.text12,'String','Wait.......')
    drawnow
    [fileROI,pathROI] = uigetfile({'*.nii','NIfTI files'},'Select 3D ROI masks',...
        'MultiSelect','on');
    if ~isequal(fileROI,0)
        if ~iscell(fileROI)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
            fileROI = {fileROI};
        end
        fileROI = sort(fileROI);
        xs=1;
        nE = 1;
        loopSize = size(fileROI,2);
        while xs <= loopSize
           xStr = fileROI{1,xs}(1:3);
           net1 = strncmp(xStr,fileROI,3);
           nEle(nE) = sum(net1);
           xs = xs + sum(net1);
           nE = nE+1;
        end
        nOfNets = size(nEle,2);

        for jh = 1:loopSize
            struRoi = nifti([pathROI fileROI{jh}]);
            ROIMAt = struRoi.dat(:,:,:);
            ROIMAt(ROIMAt>0)=1;
            roi1 = ROIMAt;
%             [iCo,jCo,kCo]= ind2sub(size(roi1), find(roi1>0));
%             coords = [iCo,jCo,kCo];
            ncentr = 1;
            for yt = 1:size(roi1,3)
                roiFind = roi1(:,:,yt);
                roiFind = fliplr(roiFind);
                roiFind = permute(roiFind,[2,1]);
                roiFind = fliplr(roiFind);
                roiFind = flipud(roiFind);
                if ~isequal(sum(sum(roiFind)),0)
                    roiFind2 = logical(roiFind);
                    s = regionprops(roiFind2,'centroid');
                    if size(s,1)>1
                        for ytr = 1:size(s,1)
                            centNs(ytr,:) = s(ytr).Centroid;
                        end
                        sFnl(ncentr,1:2) = mean(centNs,1);
                        sFnl(ncentr,3) = yt;  
                    else
                        sFnl(ncentr,1:2) = s.Centroid(1,:);
                        sFnl(ncentr,3) = yt;  
                    end
                    ncentr = ncentr+1;
                    clear s
                end
            end
            FinalCentrCoord(jh,:) = round(median(sFnl,1));
            clear sFnl
        end
        nodeOK = 1;
        CoordMS = FinalCentrCoord;
        
        DimROI = struRoi.dat.dim;
        DimROIDiff = [91 109 91]./DimROI;
        
        for hg = 1:size(CoordMS,1)
            CoordMS(hg,:) = CoordMS(hg,:).*DimROIDiff;
        end
        CoordMS = round(CoordMS);
        
        seedList = transpose(fileROI);
        set(handles.text12,'String',sprintf('%d coordinates of %d netorworks were added!',size(CoordMS,1),nOfNets))
        if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
            set(handles.runB,'Enable','on')
        end
    end
end
    
function runB_Callback(hObject, eventdata, handles)
global CoordMS Dir COVDATA g1_VAR g2_VAR nOfLines EROI fileROI
set(handles.text13,'String','Running.....')
drawnow
mkdir(Dir)

%%%%%% Checking...
if ~isequal(size(CoordMS,1),size(g1_VAR,1))
    warndlg(sprintf('The number of nodes added (%d) should to match with the used on first level (%d)',size(CoordMS,1),size(g1_VAR,1)),'Ops! Process Aborted')
    return
end
if get(handles.COV,'Value')
    if ~isequal((size(g1_VAR,3)+size(g2_VAR,3)),nOfLines)
        warndlg(sprintf('The number of data on each covariate (%d) should to match with the sum of the number of subjects on Group 1 and 2',nOfLines,(size(g1_VAR,3)+size(g2_VAR,3))),'Ops! Process Aborted')
        return
    end
end

Mat = cat(3,g1_VAR,g2_VAR);

save([Dir,'matrices'],'Mat')

loopSize = size(CoordMS,1);

if get(handles.COV,'Value')
    design = zeros(size(Mat,3),3+size(COVDATA,2));
    design(:,1) = 1;
    design(1:size(g1_VAR,3),2) = 1;
    design(size(g1_VAR,3)+1:end,3) = 1;
    design(:,4:4+size(COVDATA,2)-1) = COVDATA;
    ContrastVetG1_G2 = [0 1 -1 zeros(1,size(COVDATA,2))];
    ContrastVetG2_G1 = [0 -1 1 zeros(1,size(COVDATA,2))];
else
    design = zeros(size(Mat,3),2);
    design(1:size(g1_VAR,3),1) = 1;
    design(size(g1_VAR,3)+1:end,2) = 1;
    ContrastVetG1_G2 = [1 -1];
    ContrastVetG2_G1 = [-1 1];
end

save([Dir,'designMatrix'],'design')

fideG = fopen([Dir,'nodeLabels.txt'],'w+'); % CHARGE OUTPUT LOG FILE

if EROI
    for jh = 1:loopSize
        fprintf(fideG,'%s\r\n',fileROI{jh});
    end
else
    for jh = 1:loopSize
        fprintf(fideG,'%s\r\n',jh);
    end
end

T = [1.9999   -0.0112   -0.0104  -91.3707
    0.0113    1.9999    0.0120 -131.9940
    0.0103   -0.0120    1.9999  -69.0593
         0         0         0    1.0000];

cor = round(CoordMS);
mni = T*[cor(:,1) cor(:,2) cor(:,3) ones(size(cor,1),1)]';
mni = mni';
mni(:,4) = [];
mni = round(mni);

save([Dir,'NodeCoordinates'],'mni')
fclose('all');
set(handles.text13,'String','Done!')
disp('Possible Contrasts:')
fprintf('For G1>G2: [%s]\r\n',num2str(ContrastVetG1_G2))
fprintf('For G2>G1: [%s]\r\n',num2str(ContrastVetG2_G1))
disp('------------')

function ROIfiles_Callback(hObject, eventdata, handles)
if get(handles.ROIfiles,'Value')
    set(handles.CoodList,'Value',0)
    set(handles.addNODES,'String','Add ROIs images')
else
    set(handles.CoodList,'Value',1)
    set(handles.addNODES,'String','Add Coordinate list')
end

function CoodList_Callback(hObject, eventdata, handles)
if get(handles.CoodList,'Value')
    set(handles.ROIfiles,'Value',0)
    set(handles.addNODES,'String','Add Coordinate list')
else
    set(handles.ROIfiles,'Value',1)
    set(handles.addNODES,'String','Add ROIs images')
end

function figure1_CreateFcn(hObject, eventdata, handles)
global g2OK g1OK dirOK nodeOK
g2OK=0; g1OK=0; dirOK=0; nodeOK=0;
