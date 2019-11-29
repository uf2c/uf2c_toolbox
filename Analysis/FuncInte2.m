function varargout = FuncInte2(varargin)
% UF²C M-file for FuncInte2.fig
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

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FuncInte2_OpeningFcn, ...
                   'gui_OutputFcn',  @FuncInte2_OutputFcn, ...
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

function FuncInte2_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = FuncInte2_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function AddFunc_Callback(hObject, eventdata, handles)
global fileFunc pathFunc numofVoluF nOFsubjects matrixF VoxSizeF

if ~get(handles.FazPrep,'Value') 
    set(handles.checkReg,'Enable','on')
end

try
    clear('regVets');
end
try
    delete([tmpDIR 'additional_Reg.mat']);
end

[fileFunc,pathFunc] = uigetfile({'*.nii','NIfTI files'},'Select all the functional images','MultiSelect','on');

if ~isequal(fileFunc,0)
    
    set(handles.txtFunc,'String','Analyzing files. Wait....')
    drawnow

    if ~iscell(fileFunc)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
        fileFunc = {fileFunc};
    end
    fileFunc = sort(fileFunc); % SORT FILES IN THE ALPHABETIC ORDER
    nOFsubjects = size(fileFunc,2);

    for yt = 1:nOFsubjects
        tmpDir = dir([pathFunc fileFunc{yt}]);
        bitS(yt) = tmpDir.bytes;
    end
    sxs = unique(bitS);
    nDiftypes = size(sxs,2);
    if ~isequal(nDiftypes,1)
        warndlg({sprintf('There are %d distinct protocol types among the functional images that you added.',nDiftypes);...
          sprintf('Check the size of your files, it is the easiest way to identify the protocol homogeneity');'We can continue, but some errors and/or interpolations with distinct weights can occur.'},'Attention!');
    end
    clear bitS

    set(handles.txtFunc,'String',sprintf('%d functional image(s) added',nOFsubjects))
    drawnow

    previewF = nifti([pathFunc fileFunc{1}]);
    matrixF = previewF.dat.dim(1:3);
    numofVoluF = previewF.dat.dim(4);
    VoxSizeF = previewF.hdr.pixdim(2:4);
    VoxS_Res = [num2str(VoxSizeF(1,1)),'x',num2str(VoxSizeF(1,2)),'x',num2str(VoxSizeF(1,3))];
    MatS_Res = [num2str(matrixF(1,1)),'x',num2str(matrixF(1,2)),'x',num2str(matrixF(1,3))];

    set(handles.numDyna,'String',numofVoluF)
    set(handles.wsIn,'Enable','on')
    set(handles.wsIn,'String',numofVoluF)
    set(handles.SizeVoxFunc,'String',VoxS_Res)
    set(handles.SizeMatFunc,'String',MatS_Res)
    set(handles.txtregress,'String','0 added')
    set(handles.checkReg,'Value',0)
    set(handles.pushbutton6,'Enable','on')
    
    if get(handles.FazPrep,'Value') && get(handles.AddROIm,'Value')
        set(handles.Run,'Enable','on')
    end
    if get(handles.FazPrep,'Value') && get(handles.Addlist,'Value')
        set(handles.Run,'Enable','on')
    end
    clear previewF
end

function AddStru_Callback(hObject, eventdata, handles)
global fileStru pathStru

[fileStru,pathStru] = uigetfile({'*.nii','NIfTI files'},'Select all the Structural images','MultiSelect','on');

if ~isequal(fileStru,0)
    if ~iscell(fileStru)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
        fileStru = {fileStru};
    end

    fileStru = sort(fileStru);

    for yt = 1:size(fileStru,2)
        tmpDir = dir([pathStru fileStru{yt}]);
        bitS(yt) = tmpDir.bytes;
    end
    sxs = unique(bitS);
    nDiftypes = size(sxs,2);
    if ~isequal(nDiftypes,1)
        warndlg({sprintf('There are %d distinct protocol types among the structural images that you added.',nDiftypes);...
          sprintf('Check the size of your files, this is the easiest way to check the protocol homogeneity');'We can continue, but some errors and/or interpolations with distinct weights can occur.'},'Attention!');
    end
    clear bitS

    set(handles.txtstru,'String',sprintf('%d Structural image(s) added',size(fileStru,2)))

    previewS = nifti([pathStru fileStru{1}]);
    matrixS = previewS.dat.dim;
    VoxSizeS = previewS.hdr.pixdim(2:4);

    VoxS_Res = [num2str(VoxSizeS(1,1)),'x',num2str(VoxSizeS(1,2)),'x',num2str(VoxSizeS(1,3))];
    MatS_Res = [num2str(matrixS(1,1)),'x',num2str(matrixS(1,2)),'x',num2str(matrixS(1,3))];

    set(handles.edit9,'String',num2str(VoxS_Res))
    set(handles.edit10,'String',num2str(MatS_Res))
    set(handles.pushbutton7,'Enable','on')
    clear previewS
end

function AddROIm_Callback(hObject, eventdata, handles)
global VoxSizeF pathROI  matrixF

if get(handles.AddROIm,'Value')
    set(handles.Addlist,'Enable','off')
    set(handles.Addlist,'Value',0)
    set(handles.text24,'String','')
    [fileROI,pathROI] = uigetfile({'*.nii','NIfTI files'},'Select 3D ROI masks','MultiSelect','on');
    
    if ~isequal(fileROI,0) && size(fileROI,2)>1
        if ~iscell(fileROI)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
            fileROI = {fileROI};
        end
        fileROI = sort(fileROI);

        previewF2 = nifti([pathROI fileROI{1}]);
        ROImaskMAt = previewF2.dat(:,:,:);
        matrixF2 = previewF2.dat.dim(1:3);
        VoxSizeF2 = previewF2.hdr.pixdim(2:4);
        VoxSizeF2 = num2str(VoxSizeF2);
        VoxSizeF2 = double(str2num(VoxSizeF2));
        
        VoxSizeFtmp = num2str(VoxSizeF);
        VoxSizeFtmp = double(str2num(VoxSizeFtmp));
        
        if get(handles.FazPrep,'Value')
            if isequal(matrixF2,matrixF) && isequal(VoxSizeF2,VoxSizeFtmp)
               set(handles.text7,'String',sprintf('%d ROI mask(s) added!',size(fileROI,2)))
               set(handles.Run,'Enable','on')
               xs=1;
               nE = 1;
               while xs <= size(fileROI,2)
                   xStr = fileROI{1,xs}(1:3);
                   net1 = strncmp(xStr,fileROI,3);
                   nEle(nE) = sum(net1);
                   xs = xs + sum(net1);
                   nE = nE+1;
               end
               nOfNets = size(nEle,2);
               set(handles.textNETS,'String',sprintf('%d network(s) identified!',nOfNets))
            else
                warndlg('The ROI mask added do not match to the preprocessed functional images that you previously added!', 'Warning');
                set(handles.AddROIm,'Value',0)
                set(handles.Addlist,'Enable','on')
                set(handles.Run,'Enable','off')
           end
        else
            if get(handles.SPMbb,'Value')
               if isequal(matrixF2,[53 63 46]) && isequal(VoxSizeF2,VoxSizeF)
                   set(handles.text7,'String',sprintf('%d ROI mask(s) added!',size(fileROI,2)))
                   set(handles.Run,'Enable','on')
                   xs=1;
                   nE = 1;
                   while xs <= size(fileROI,2)
                       xStr = fileROI{1,xs}(1:3);
                       net1 = strncmp(xStr,fileROI,3);
                       nEle(nE) = sum(net1);
                       xs = xs + sum(net1);
                       nE = nE+1;
                   end
                   nOfNets = size(nEle,2);
                   set(handles.textNETS,'String',sprintf('%d network(s) identified!',nOfNets))
               else
                   warndlg('The ROI mask added do not match to the functional images post-processing parameters (SPM default)!', 'My Warn Dialog');
                   set(handles.AddROIm,'Value',0)
                   set(handles.Addlist,'Enable','on')
                   set(handles.Run,'Enable','off')
               end
            else
                tmpStr = num2str(VoxSizeF2);
                VoxSizeF2 = str2num(tmpStr);
                if isequal(matrixF2,[91 109 91]) && isequal(double(VoxSizeF2),[2 2 2])
                   set(handles.text7,'String',sprintf('%d ROI mask(s) added!',size(fileROI,2)))
                   set(handles.Run,'Enable','on')
                   xs=1;
                   nE = 1;
                   while xs <= size(fileROI,2)
                       xStr = fileROI{1,xs}(1:3);
                       net1 = strncmp(xStr,fileROI,3);
                       nEle(nE) = sum(net1);
                       xs = xs + sum(net1);
                       nE = nE+1;
                   end
                   nOfNets = size(nEle,2);
                   set(handles.textNETS,'String',sprintf('%d network(s) identified!',nOfNets))
                else
                    warndlg('The ROI mask added do not match to the functional images that you previously added!', 'Warning');
                    set(handles.AddROIm,'Value',0)
                    set(handles.Addlist,'Enable','on')
                    set(handles.Run,'Enable','off')
               end
            end
        end
        
        try
            delete('ddROIs.mat')
        catch
            disp('Seeds Added!')
        end
        ddROIs.prevF2 = previewF2;
        ddROIs.mtxF2 = matrixF2;
        ddROIs.VoxSF2 = VoxSizeF2;
        ddROIs.ROImMat = ROImaskMAt;
        ddROIs.fROI = fileROI;
        ddROIs.PathRoi = pathROI;
        ddROIs.NoNet = nOfNets;
        ddROIs.NoEle = nEle;
        uf2cdir = fileparts(which('uf2c'));
        save([pathROI 'ddROIs.mat'],'ddROIs')
    else
        if size(fileROI,2) == 1
            warndlg('You should to add at least 2 ROIs.', 'Warning');
        end
        set(handles.AddROIm,'Value',0)
        set(handles.Addlist,'Enable','on')
        set(handles.Addlist,'Value',0)
        set(handles.text24,'String','')
    end
else
    set(handles.Addlist,'Enable','on')
    set(handles.text7,'String','')
    set(handles.textNETS,'String','')
end

function Addlist_Callback(hObject, eventdata, handles)
%global seedList nOfNets nEle

if isequal(get(handles.Addlist,'Value'),1)
    set(handles.AddROIm,'Enable','off')
    set(handles.AddROIm,'Value',0)
    set(handles.text7,'String','')
    [ROIlist,pathlist] = uigetfile({'*.txt','Text files'},'Select text file with all coordinates','MultiSelect','off');
    if ~isequal(ROIlist,0)
        seedListT = importdata([pathlist ROIlist]);
        if iscell(seedListT)
            if size(seedListT,1) > 1
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
                set(handles.text24,'String',sprintf('%d coord(s) of %d nets were added!',size(seedListT,1),nOfNets))
                set(handles.Run,'Enable','on')
                for i = 1:size(seedListT,1)
                    seedList(i,:) = str2num(seedListT{i}(4:end));
                end
                try
                    delete('Global_SL.mat')
                catch
                    disp('Seeds Added!')
                end
                Global_SL.SL = seedList;
                Global_SL.NoN = nOfNets;
                Global_SL.NoEl = nEle;
                save('Global_SL.mat','Global_SL')
            else
                warndlg('You should to add at least 2 coordinates.','Warning');
                set(handles.Addlist,'Value',0)
                set(handles.AddROIm,'Enable','on')
                set(handles.AddROIm,'Value',0)
                set(handles.text7,'String','')
            end
        else
            if size(seedListT.data,1) > 1
                nOfNets = str2num(seedListT.textdata{end,1}(1:2));
                xs=1;
                nE = 1;
                while xs <= size(seedListT.data,1)
                   xStr = seedListT.textdata{xs,1}(1:3);
                   net1 = strncmp(xStr,seedListT.textdata,3);
                   nEle(nE) = sum(net1);
                   xs = xs + sum(net1);
                   nE = nE+1;
                end
                set(handles.text24,'String',sprintf('%d coord(s) of %d nets were added!',size(seedListT.data,1),nOfNets))
                set(handles.Run,'Enable','on')
                seedList(:,:) = seedListT.data;
                try
                    delete('Global_SL.mat')
                catch
                    disp('Seeds Added!')
                end
                Global_SL.SL = seedList;
                Global_SL.NoN = nOfNets;
                Global_SL.NoEl = nEle;
                save('Global_SL.mat','Global_SL')
            else
                warndlg('You should to add at least 2 coordinates.','Warning');
                set(handles.Addlist,'Value',0)
                set(handles.AddROIm,'Enable','on')
                set(handles.AddROIm,'Value',0)
                set(handles.text7,'String','')
            end
        end
    else
        set(handles.Addlist,'Value',0)
        set(handles.AddROIm,'Enable','on')
        set(handles.AddROIm,'Value',0)
        set(handles.text7,'String','')
    end
else
    set(handles.AddROIm,'Enable','on')
    set(handles.text24,'String','')
    set(handles.Run,'Enable','off')
end
    
function checkReg_Callback(hObject, eventdata, handles)
if isequal(get(handles.checkReg,'Value'),1)
    set(handles.addreg,'Enable','on')
    set(handles.refreshb,'Enable','on')
else
    set(handles.addreg,'Enable','off')
    set(handles.refreshb,'Enable','off')
    set(handles.txtregress,'String','0 added')
end

function addreg_Callback(hObject, eventdata, handles)
global regVets
try
    clear(regVets)
end
Addreg

function refreshb_Callback(hObject, eventdata, handles)
global nrg 
try
    if nrg < 0
        nrg = 0;
    end
    
    nrT = num2str(nrg);
    set(handles.txtregress,'String',[nrT ' added'])
end

if nrg == 0
     set(handles.checkReg,'Value',0)
     set(handles.addreg,'Enable','off')
     set(handles.refreshb,'Enable','off')
end


function Run_Callback(hObject, eventdata, handles)
global fileStru pathStru fileFunc pathFunc numofVoluF nrg tmpDIR pathROI

SPMdir2 = which('spm');
SPMdir2 = SPMdir2(1:end-5);

if isequal(get(handles.FazPrep,'Value'),0)
    if ~isequal(size(fileStru,2),size(fileFunc,2))
        warndlg({sprintf('The numbers of functional (%d) and structural (%d) images are different.',size(fileFunc,2),size(fileStru,2));...
          'You need to add one functional image for each strutural and vice versa.';...
          'If you have more than one functional section from the same subject, make copies of the structural image and add all!'},...
          'Attention: process aborted!');
        return
    end
end

set(handles.status,'String','Running....')
drawnow

fprintf('\r\n')
fprintf('UF²C =============== %s\n',spm('time'));
fprintf('Process started\n');
fprintf('==========================================\r\n');

if get(handles.Addlist,'Value')
    load('Global_SL.mat')
    %delete('Global_SL.mat')
    seedList = Global_SL.SL;
    nOfNets = Global_SL.NoN;
    nEle = Global_SL.NoEl;
else
    load([pathROI 'ddROIs.mat'])
    % delete('ddROIs.mat')
    previewF2 = ddROIs.prevF2;
    matrixF2 = ddROIs.mtxF2;
    VoxSizeF2 = ddROIs.VoxSF2;
    % VoxSizeF = ddROIs.VoxSF;
    ROImaskMAt = ddROIs.ROImMat;
    fileROI = ddROIs.fROI;
    pathROI = ddROIs.PathRoi;
    nOfNets = ddROIs.NoNet;
    nEle = ddROIs.NoEle;
end
    
seedSizeMat = [str2num(get(handles.SSInp,'String'))];
WS = str2num(get(handles.wsIn,'String')); %Window Size (SIZE OF THE MOVING AVERAGE TO CALC DE CORRELATIONS)
TR = str2num(get(handles.TRInp,'String'));    % repetition time (s)

x = seedSizeMat(1); % ROI size in voxels (HALF VALUE FOR EACH SIDE)
y = seedSizeMat(2);
z = seedSizeMat(3);

smWin = str2double(get(handles.wsIn, 'String'));
ratioVol = numofVoluF/smWin;

multiWaitbar('Total Progress', 0, 'Color', 'g' );  %CHARGE WAIT BAR
if isequal(get(handles.FazPrep,'Value'),0)
    multiWaitbar('Filtering & Regression', 0, 'Color', 'y' ); %CHARGE WAIT BAR
end
multiWaitbar('Extracting and processing ROIs time-series', 0, 'Color', 'b' ); %CHARGE WAIT BAR

dateNow = clock;
foldername = sprintf('1-Total_Log_%d_%d_%d--%d_%d', dateNow(3),dateNow(2),dateNow(1),dateNow(4),dateNow(5));

mkdir(pathFunc,foldername)

fideSJ = fopen([pathFunc,filesep,foldername,filesep,'1-Subjects.txt'],'w+'); 
for yTy = 1:size(fileFunc,2)
    fprintf(fideSJ,'%d: \t%s\r\n',yTy,fileFunc{1,yTy});
end
fclose(fideSJ);

imgRR = getframe(FuncInte2);
imwrite(imgRR.cdata, [pathFunc,filesep,foldername,filesep,'1-Your_Choices.png']);

fideG = fopen([pathFunc,filesep,foldername,filesep,'1-Seeds_Order.txt'],'w+'); % CHARGE OUTPUT LOG FILE
pathFLOG =  [pathFunc,foldername,filesep];
fprintf(fideG,'                                UF²C - Cross_Correlation ROI Analysis\r\n\r\n');
fprintf(fideG,'Number of subjects included: %d \r\n\r\n', size(fileFunc,2));
if get(handles.Addlist,'Value')
    fprintf(fideG,'Number of seeds included: \t%d \r\n\r\n', size(seedList,1));
    for jh = 1:size(seedList,1)
        fprintf(fideG,'%d: \t%d \t%d \t%d\r\n',jh,seedList(jh,:));
    end
    fclose(fideG);
else
    fprintf(fideG,'Number of seeds included: \t%d \r\n\r\n', size(fileROI,2));
    for jh = 1:size(fileROI,2)
        fprintf(fideG,'%d: \t%s\r\n',jh,fileROI{jh});
    end
    fclose(fideG);
end

fideGx = fopen([pathFunc,filesep,foldername,filesep,'1-Seeds_Description.txt'],'w+'); % CHARGE OUTPUT LOG FILE

if isequal(get(handles.AddROIm,'Value'),1)
    fprintf('\r\n')
    fprintf('UF²C =============== %s\n',spm('time'));
    fprintf('Processing ROIs Files\n');
        fprintf('\tPerforming ROIs anatomical labeling\n');

    fprintf(fideGx,'UF²C: ROIs Anatomical Labelling and Quantifications\r\n\r\n');
    fprintf(fideGx,'These data are better visualized if imported to a MS Excel type software\r\n\r\n');
    fprintf(fideGx,'Credits:\r\n');
    fprintf(fideGx,'Maximum probability tissue labels derived from:\r\n');
    fprintf(fideGx,'''MICCAI 2012 Grand Challenge and Workshop on Multi-Atlas Labeling''\r\n');
    fprintf(fideGx,'These data were released under the Creative Commons Attribution-NonCommercial (CC BY-NC) with no end date.\r\n');
    fprintf(fideGx,'Labeled data provided by Neuromorphometrics, Inc. (http://Neuromorphometrics.com/) under academic subscription".\r\n\r\n\r\n');

    try
        loopSize = size(fileROI,2);

        stIndx = 1;

        for hg = 1:size(nEle,2)
            str_nEle(hg) = stIndx;
            stIndx = stIndx+nEle(hg);
        end

        for jh = 1:loopSize
            fprintf(fideGx,'ROI %d:\t%s\r\n', jh,fileROI{jh});
            tester = jh<=(str_nEle+nEle)-1;
            networkNumb = find(tester==1,1);

            Net = networkNumb;

            if Net>1
                ROIn = jh - sum(nEle(1:Net-1));
            else
                ROIn = jh;
            end

            fprintf(fideGx,'\tROI(r)/network(n) short name:\tr%dn%d\r\n',ROIn,Net);

            struRoi = nifti([pathROI fileROI{jh}]);
            roi11 = double(struRoi.dat(:,:,:));
            ROIMAt(:,:,:,jh) = roi11;
            roi11(roi11>0)=1;
            nofVoxTOT = sum(sum(sum(roi11)));
            RoiVol = nofVoxTOT.*prod(struRoi.hdr.pixdim(2:4));
            fprintf(fideGx,'\tNumber of Voxels:\t%d\r\n',nofVoxTOT);
            fprintf(fideGx,'\tTotal Volume:\t%d\tmm³\r\n',RoiVol);

            %%%%%%%% Find centroi of all ROIs
            ncentr = 1;
            for yt = 1:size(roi11,3)
                roiFind = roi11(:,:,yt);
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

            if ~get(handles.SPMbb,'Value')
                %%%%%%%% Neuromorphometric Labeling
                DirUF2C = which('uf2c');
                try
                    LabelImgS = nifti([DirUF2C(1:end-6) 'Utilities' filesep 'labels_Neuromorphometrics_91-109-91.nii']);
                    load([DirUF2C(1:end-6) 'Utilities' filesep 'NeuromorphometricsLabels.mat']); 
                catch
                    LabelImgS = nifti([DirUF2C(1:end-13) 'Utilities' filesep 'labels_Neuromorphometrics_91-109-91.nii']);
                    load([DirUF2C(1:end-13) 'Utilities' filesep 'NeuromorphometricsLabels.mat']); 
                end

                LabelImgMAT = LabelImgS.dat(:,:,:);
                MultROILab = LabelImgMAT.*roi11;
                UniQ = unique(MultROILab);
                UniQ = UniQ(2:end);
                nofregions = numel(UniQ);

                IDXcentroid = LabelImgMAT((size(roi11,1)-FinalCentrCoord(jh,1)),FinalCentrCoord(jh,2),FinalCentrCoord(jh,3));
                if IDXcentroid == 0
                    CentroidRegions{jh,1} = 'Extra Brain Region';
                else
                    Posi = find(NeuromorphometricsLabels.Index == IDXcentroid);
                    CentroidRegions{jh,1} = NeuromorphometricsLabels.Label{Posi};
                end
                fprintf(fideGx,'\tCentroid Coordinate (voxel space):\t%dx%dx%d \r\n',FinalCentrCoord(jh,1),FinalCentrCoord(jh,2),FinalCentrCoord(jh,3));
                fprintf(fideGx,'\tCentroid Exact Regions:\t%s \r\n\r\n',CentroidRegions{jh,1});
                fprintf(fideGx,'\tRegions included on the ROI mask\r\n');
                fprintf(fideGx,'\tRegion Name\tNumber of Voxels\tPercentage of the Total\r\n');

                for nr = 1:nofregions
                    Bint = double(MultROILab==UniQ(nr));
                    nOfvox(nr,1) = sum(sum(sum(Bint)));
                    Posi = find(NeuromorphometricsLabels.Index == UniQ(nr));
                    Regions{nr,1} = NeuromorphometricsLabels.Label{Posi};
                end
                [nOfvox,www] = sort(nOfvox,'descend');
                Regions = Regions(www);

                for nr = 1:nofregions
                    fprintf(fideGx,'\t%s\t%d\t%.2f\r\n',Regions{nr,1},nOfvox(nr,1),(100*(nOfvox(nr,1)/nofVoxTOT)));
                end
                fprintf(fideGx,'\r\n\r\n');

                clear sFnl nOfvox Regions Posi
            end
        end
        CoordMSori = FinalCentrCoord;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        seedList = transpose(fileROI);
        fclose(fideGx);
        fprintf('\tAnatomical labeling file saved\n');
        fprintf('Done! ============== %s\r\n',spm('time'));
    catch
        fprintf(fideGx,'FAILURE TO PERFORM THE ANATOMICAL LABELING.\r\n');
        fprintf(fideGx,'Probbably, the images you addded do not match our anatomical temlplate parameters.\r\n');
        fprintf('\tFail to perform natomical labeling\n');
        fprintf('Continuing! ============== %s\r\n',spm('time'));
        fclose(fideGx);
    end
end    


if isequal(get(handles.Addlist,'Value'),1)
    fprintf('\r\n')
    fprintf('UF²C =============== %s\n',spm('time'));
    fprintf('Processing Coordinate List\r\n');
    fprintf('\tPerforming ROIs anatomical labeling\n');

    fprintf(fideGx,'UF²C: ROIs Anatomical Labelling and Quantifications\r\n\r\n');
    fprintf(fideGx,'These data are better visualized if imported to a MS Excel type software\r\n\r\n');
    fprintf(fideGx,'Credits:\r\n');
    fprintf(fideGx,'Maximum probability tissue labels derived from:\r\n');
    fprintf(fideGx,'''MICCAI 2012 Grand Challenge and Workshop on Multi-Atlas Labeling''\r\n');
    fprintf(fideGx,'These data were released under the Creative Commons Attribution-NonCommercial (CC BY-NC) with no end date.\r\n');
    fprintf(fideGx,'Labeled data provided by Neuromorphometrics, Inc. (http://Neuromorphometrics.com/) under academic subscription".\r\n\r\n\r\n');

    try
        UF2Cdir = which('uf2c');
        tmpDIR = [UF2Cdir(1:end-6) 'Utilities' filesep 'labels_Neuromorphometrics_91-109-91.nii'];

        Vhdr = nifti(tmpDIR);

        Dx = 91;
        Dy = 109;
        Dz = 91;

        loopSize = size(seedList,1);

        stIndx = 1;
        for hg = 1:size(nEle,2)
            str_nEle(hg) = stIndx;
            stIndx = stIndx+nEle(hg);
        end

        for jh = 1:loopSize
            fprintf(fideGx,'ROI:\t%d\r\n', jh);

            tester = jh<=(str_nEle+nEle)-1;
            networkNumb = find(tester==1,1);

            Net = networkNumb;

            if Net>1
                ROIn = jh - sum(nEle(1:Net-1));
            else
                ROIn = jh;
            end

            fprintf(fideGx,'\tROI(r)/network(n) short name:\tr%dn%d\r\n',ROIn,Net);

            TcM = Vhdr.mat;
            EPIcoord = [seedList(jh,1) seedList(jh,2) seedList(jh,3) 1]*(inv(TcM))';
            EPIcoord = round(EPIcoord);

            CoordMS(jh,:) = EPIcoord(1,:);

            Sxx = x/2;
            Syy = y/2;
            Szz = z/2;
            if Sxx<1
                Sxx=0;
            end
            if Syy<1
                Syy=0;
            end
            if Szz<1
                Szz=0;
            end
            nofVoxTOT = (x+1)*(y+1)*(z+1);

            RoiVol = nofVoxTOT.*prod(Vhdr.hdr.pixdim(2:4));

            fprintf(fideGx,'\tNumber of Voxels:\t%d\tTotal Volume:\t%d\tmm³\r\n',nofVoxTOT,RoiVol);

            roi1 = zeros(Dx,Dy,Dz);%Prelocation
            %%%%%%%%%%%%%%%%%%%% Create A cubic ROI centered in the input coordinate %%%%%%%%%%%%%%%%%%%%%%
            roi1(EPIcoord(1):EPIcoord(1)+Sxx, EPIcoord(2):EPIcoord(2)+Syy, EPIcoord(3):EPIcoord(3)+Szz) = 1;
            roi1(EPIcoord(1)-Sxx:EPIcoord(1), EPIcoord(2):EPIcoord(2)+Syy, EPIcoord(3):EPIcoord(3)+Szz) = 1;

            roi1(EPIcoord(1):EPIcoord(1)+Sxx, EPIcoord(2)-Syy:EPIcoord(2), EPIcoord(3):EPIcoord(3)+Szz) = 1;
            roi1(EPIcoord(1):EPIcoord(1)+Sxx, EPIcoord(2):EPIcoord(2)+Syy, EPIcoord(3)-Szz:EPIcoord(3)) = 1;

            roi1(EPIcoord(1)-Sxx:EPIcoord(1), EPIcoord(2):EPIcoord(2)+Syy, EPIcoord(3)-Szz:EPIcoord(3)) = 1;
            roi1(EPIcoord(1):EPIcoord(1)+Sxx, EPIcoord(2)-Syy:EPIcoord(2), EPIcoord(3)-Szz:EPIcoord(3)) = 1;

            roi1(EPIcoord(1)-Sxx:EPIcoord(1), EPIcoord(2)-Syy:EPIcoord(2), EPIcoord(3)-Szz:EPIcoord(3)) = 1;
            roi1(EPIcoord(1)-Sxx:EPIcoord(1), EPIcoord(2)-Syy:EPIcoord(2), EPIcoord(3):EPIcoord(3)+Szz) = 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            roi1Fall(:,:,:,jh) = roi1;

            DirUF2C = which('uf2c');
            try
                LabelImgS = nifti([DirUF2C(1:end-6) 'Utilities' filesep 'labels_Neuromorphometrics_91-109-91.nii']);
                load([DirUF2C(1:end-6) 'Utilities' filesep 'NeuromorphometricsLabels.mat']); 
            catch
                LabelImgS = nifti([DirUF2C(1:end-13) 'Utilities' filesep 'labels_Neuromorphometrics_91-109-91.nii']);
                load([DirUF2C(1:end-13) 'Utilities' filesep 'NeuromorphometricsLabels.mat']); 
            end

            LabelImgMAT = LabelImgS.dat(:,:,:);
            MultROILab = LabelImgMAT.*roi1;
            UniQ = unique(MultROILab);
            UniQ = UniQ(2:end);
            nofregions = numel(UniQ);

            IDXcentroid = LabelImgMAT(EPIcoord(1),EPIcoord(2),EPIcoord(3));
            if IDXcentroid == 0
                CentroidRegions{jh,1} = 'Extra Brain Region';
            else
                Posi = find(NeuromorphometricsLabels.Index == IDXcentroid);
                CentroidRegions{jh,1} = NeuromorphometricsLabels.Label{Posi};
            end
            fprintf(fideGx,'\tCentroid Coordinate (Vox Space):\t%dx%dx%d \r\n',(size(roi1,1)-EPIcoord(1)),EPIcoord(2),EPIcoord(3));
            fprintf(fideGx,'\tCentroid Exact Regions:\t%s \r\n\r\n',CentroidRegions{jh,1});
            fprintf(fideGx,'\tRegions included on the ROI mask\r\n');
            fprintf(fideGx,'\tRegion Name\tNumber of Voxels\tPercentage of the Total\r\n');

            for nr = 1:nofregions
                Bint = MultROILab==UniQ(nr);
                nOfvox(nr,1) = sum(sum(sum(Bint)));
                Posi = find(NeuromorphometricsLabels.Index == UniQ(nr));
                Regions{nr,1} = NeuromorphometricsLabels.Label{Posi};
            end

            [nOfvox,www] = sort(nOfvox,'descend');
            Regions = Regions(www);

            for nr = 1:nofregions
                fprintf(fideGx,'\t%s\t%d\t%.2f\r\n',Regions{nr,1},nOfvox(nr,1),(100*(nOfvox(nr,1)/nofVoxTOT)));
            end

            fprintf(fideGx,'\r\n\r\n');
            clear sFnl
        end
        fclose(fideGx);
        fprintf('\tAnatomical labeling file saved\n');
        fprintf('Done! ============== %s\r\n',spm('time'));
    catch
        fprintf(fideGx,'FAILURE TO PERFORM THE ANATOMICAL LABELING.\r\n');
        fprintf(fideGx,'Probbably, the images you addded do not match our anatomical temlplate parameters.\r\n');
        fprintf('\tFail to perform natomical labeling\n');
        fprintf('Continuing! ============== %s\r\n',spm('time'));
        fclose(fideGx);
    end
end

if isequal(get(handles.FazPrep,'Value'),0) % just if you need preprocessing
    
    Totmovtxt = fopen([pathFunc,filesep,foldername,filesep,'Total_Movement.txt'],'w+'); % CHARGE OUTPUT LOG FILE
    fprintf(Totmovtxt,'Subject \t Maximum Displacemente (mm) \t Maximum Rotation (degree) \t Avg Framiwise Displacemente (FD) (mm) \t numb of FD Censored scans \t Avg DVARs (%%) \t numb of DVAR Censored scans \t Total Censored scans \r\n');
    if get(handles.FazFilt,'Value') % just if you want pass-band filtering
        fprintf('\r\n')
        fprintf('UF²C =============== %s\n',spm('time'));
        fprintf('Designing filters\n');
        %%%%%%%%%%%%%%%%%%%%%% FILTER DESIGN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Fs = 1/TR;     % sampling frequency (Hz)

        FpassLOW = str2double(get(handles.LPF,'String'));   % passband frequency   (Hz) (default: 0.1)
        FstopLOW = str2double(get(handles.LPF,'String')) + (Fs/numofVoluF); % stopband frequency   (Hz) (default: 0.15)
        ApassLOW = 1;                                       % passband ripple      (dB) (default: 1)
        AstopLOW = str2double(get(handles.stba,'String'));  % stopband attenuation (dB) (default: 40)

        FpassHIGH = str2double(get(handles.HPF,'String'));  % passband frequency   (Hz) (default: 0.008)
        FstopHIGH = str2double(get(handles.HPF,'String')) - (Fs/numofVoluF); % stopband frequency   (Hz) (default: 0.005)
                        
        if FstopHIGH<0
            FstopHIGH = str2double(get(handles.HPF,'String'))-(str2double(get(handles.HPF,'String')).*0.9);
        end
        
        AstopHIGH = str2double(get(handles.stba,'String')); % stopband attenuation (dB) (default: 40)
        ApassHIGH = 1;
        
        hLOW = fdesign.lowpass('Fp,Fst,Ap,Ast',FpassLOW,FstopLOW,ApassLOW,AstopLOW,Fs);
        HdLOW = design(hLOW,'equiripple');

        hHIGH  = fdesign.highpass('Fst,Fp,Ast,Ap',FstopHIGH,FpassHIGH,AstopHIGH,ApassHIGH,Fs);
        HdHIGH = design(hHIGH, 'equiripple');
        fprintf('Done! ============== %s\r\n',spm('time'));
        
        aFil = 1;
        tmp1 = 3*(max(length(HdLOW.Numerator),length(aFil))-1);
        tmp2 = 3*(max(length(HdHIGH.Numerator),length(aFil))-1);
    else
        FpassLOW = [];  % passband frequency   (Hz) (default: 0.1)
        FstopLOW = []; % stopband frequency   (Hz) (default: 0.15)
        ApassLOW = [];                                     % passband ripple      (dB) (default: 1)
        AstopLOW = [];  % stopband attenuation (dB) (default: 40)
        FpassHIGH = [];  % passband frequency   (Hz) (default: 0.008)
        FstopHIGH = []; % stopband frequency   (Hz) (default: 0.005)
        AstopHIGH = []; % stopband attenuation (dB) (default: 40)
        ApassHIGH = [];
        hLOW = [];
        HdLOW = [];
        hHIGH  = [];
        HdHIGH = [];
        aFil = [];
        tmp1 = [];
        tmp2 = [];
    end 
end

ScreSize = get(0,'screensize');
ScreSize = ScreSize(3:end);

WinAvSave = struct;
ExclSub = 0;

FunpreT = nifti([pathFunc,fileFunc{1}]);

if isequal(WS,FunpreT.dat.dim(4))
    dynStu= 0;
else
    dynStu= 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% START OF SUBJECT LOOP %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% START OF SUBJECT LOOP %%%%%%%%%%%%%%%%%%%
for yy = 1:size(fileFunc,2)
    
    multiWaitbar('Total Progress', (yy/size(fileFunc,2))*0.95 );
    multiWaitbar('Extracting and processing ROIs time-series', 'Reset' );
    
    fprintf('\r\n')
    fprintf('UF²C ===================== %s\n',spm('time'));
    if isequal(get(handles.FazPrep,'Value'),0)
        fprintf('Starting subj. %d  (%s...)\n',yy,fileFunc{yy}(1:round(size(fileFunc{yy},2)/3)));
    else
        fprintf('Starting subj. %d  (%s...)\n',yy,fileFunc{yy}(12:12+round(size(fileFunc{yy},2)/3)));
    end
    fprintf('================================================\r\n\r\n');
    
    try
        close(fig)
    end
    
    Funpre = nifti([pathFunc,fileFunc{yy}]);
    
    if dynStu
         if mod(Funpre.dat.dim(4),WS)
             ExclSub = ExclSub+1;
                 if isequal(get(handles.FazPrep,'Value'),0)
                     fprintf('subj. %d  (%s...) process aborted!\n',yy,fileFunc{yy}(1:round(size(fileFunc{yy},2)/3)));
                 else
                     fprintf('subj. %d  (%s...) process aborted!\n',yy,fileFunc{yy}(12:12+round(size(fileFunc{yy},2)/3)));
                 end
             fprintf('The number of dynamics and the time blocks are not multiple\n');
             continue
         end
    else
        WS = Funpre.dat.dim(4);
    end
             
    fvs = double(Funpre.hdr.pixdim(2:4)); %GET PIXEL SIZE (FUNCTIONAL)
    
    if isequal(get(handles.FazPrep,'Value'),0)
        multiWaitbar('Filtering & Regression', 'Reset' );
        dirname = Funpre.dat.fname;
        [ax,dirname,cx] = fileparts(dirname);
        
        if numel(dirname)>30
            fprintf('Attention!\n');
            fprintf('The filename is long. A shorter version will be used.\r\n')
            dirname = dirname(1:30);
        end

        nfidx = 1;
        while isequal(exist([pathFunc,dirname],'dir'),7)
            dirname = [dirname '_' num2str(nfidx)];
            nfidx = nfidx + 1;
        end

        mkdir(pathFunc,dirname)
        copyfile([pathStru,fileStru{yy}],[pathFunc,dirname,filesep,fileStru{yy}])
        file2 = fullfile(pathFunc,dirname,filesep,fileStru{yy});
        Strpre = nifti(file2);
        Svs = double(Strpre.hdr.pixdim(2:4)); % GET PIXEL SIZE (STRUCTURAL)
        copyfile([pathFunc,fileFunc{yy}],[pathFunc,dirname,filesep,fileFunc{yy}])
        
        file1 = cell(Funpre.dat.dim(4),1);
        for tt = 1:Funpre.dat.dim(4)
            file1{tt,1} = [pathFunc,dirname,filesep,fileFunc{yy},',',sprintf('%d',tt)];
        end
        
        fprintf(Totmovtxt,'%s \t',dirname);

    else
        dirname = Funpre.dat.fname;
        [ax,dirname,cx] = fileparts(dirname);
        dirname = dirname(12:end);
        
        if numel(dirname)>30
            fprintf('Attention!\n');
            fprintf('The filename is long. A shorter version will be used.\r\n')
            dirname = dirname(1:30);
        end
        
        nfidx = 1;
        while isequal(exist([pathFunc,dirname],'dir'),7)
            dirname = [dirname '_' num2str(nfidx)];
            nfidx = nfidx + 1;
        end

        mkdir(pathFunc,dirname)
    end
    

    mkdir([pathFunc,dirname,filesep,'ROI_mask'])
    mkdir([pathFunc,dirname],'Results_Log')
    mkdir([pathFunc,dirname],'Correlation_Map')

    fideSubP = fopen([pathFunc,dirname,filesep,'Results_Log',filesep,'Positive_Corr_Values.txt'],'w+'); % CHARGE OUTPUT LOG FILE
    fprintf(fideSubP,'                           UF²C - Cross_Correlation ROI Analysis\r\n\r\n');
    fprintf(fideSubP,'Case name: %s \r\n', fileFunc{yy});
    fprintf(fideSubP,'Number of temporal blocks: \t%d \r\n',ratioVol);
    fprintf(fideSubP,'Number of networks identified: \t%d \r\n',nOfNets);
    
    fideSubN = fopen([pathFunc,dirname,filesep,'Results_Log',filesep,'Negative_Corr_Values.txt'],'w+'); % CHARGE OUTPUT LOG FILE
    fprintf(fideSubN,'                           UF²C - Cross_Correlation ROI Analysis\r\n\r\n');
    fprintf(fideSubN,'Case name: %s \r\n', fileFunc{yy});
    fprintf(fideSubN,'Number of temporal blocks: \t%d \r\n',ratioVol);
    fprintf(fideSubN,'Number of networks identified: \t%d \r\n',nOfNets);

    
    mkdir([pathFunc,dirname],'Correlation_map')
    
    if isequal(get(handles.FazPrep,'Value'),0)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Preprocessing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
            STC = get(handles.STCo,'Value');
            BounB = get(handles.SPMbb,'Value');
            fileFuncX = fileFunc{yy};
            fileStruX = fileStru{yy};
            Preproc_uf2c(file1,file2,BounB,SPMdir2,Svs,fvs,pathFunc,dirname,STC,Funpre,TR,fileFuncX,fileStruX)
        catch
             set(handles.status,'String','Error...')
             warndlg(sprintf('An error occured during subject %d preprocessing. This is a SPM error. Check your data',yy), 'Process Aborted')
             return
        end
        try
            close(fig);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot and save the motion series
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mkdir([pathFunc,dirname,filesep,'Motion_Controls_Params'])
        Motfid = fopen([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'1-Motion_Quantifications.txt'],'w+');
        fprintf(Motfid,'Motion Quantifications\r\n\r\n');
        fprintf(Motfid,'Subject Name: \t%s\r\n',dirname);

        [Mot_fig,HDvV,HTvV] = uf2c_plot_motion([pathFunc,dirname,filesep,'rp_',fileFunc{yy}(1:end-3),'txt'],'on');
        imgRR = getframe(Mot_fig);
        imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'Realignment_Parameters_Plot.png']);
        saveas(Mot_fig,[pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'Realignment_Parameters_Plot'],'fig')

        fprintf(Motfid,'Maximum Displacemente (mm):\t %s\r\n',HDvV);
        fprintf(Motfid,'Maximum Rotation (degree):\t %s\r\n',HTvV);
        fprintf('Done! ============== %s\r\n',spm('time'));   
        
        fprintf(Totmovtxt,'%s \t %s \t',HDvV,HTvV);
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BRAMILA's Framewise Displacemente (FD)
        % Code from: Brain and Mind Lab at Aalto University
        % Power et al. (2012) doi:10.1016/j.neuroimage.2011.10.018 and also 
        % Power et al. (2014) doi:10.1016/j.neuroimage.2013.08.048
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\r\n')
        fprintf('UF²C =============== %s\n',spm('time'));
        fprintf('Quantifying Framewise Displacemente (FD)\r\n');

        if get(handles.FDcheck,'Value')
            cfg.motionparam = [pathFunc,dirname,filesep,'rp_',fileFunc{yy}(1:end-3),'txt'];
            cfg.prepro_suite = 'spm';
            cfg.radius = 50;
            
            [FDts,rms] = bramila_framewiseDisplacement(cfg);
            
            save([pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'FramewiseDisplacement.mat'],'FDts')

            figuFD = figure('Visible','off');
            set(figuFD,'Name','Framewise Displacemente (FD) TS',...
                'Position', round([ScreSize(1)*.15 ScreSize(2)*.15 ScreSize(1)*.4 ScreSize(1)*.2]),...
                'Color',[1 0.94 0.86]);
            plot(FDts);
            hold on
            plot(ones(1,numel(FDts)).*str2num(get(handles.FDthres,'String')));
            ylabel('FD (mm)')
            xlabel('Time Points')
            title('Average Framewise Displacement TS','FontSize', 14);
            drawnow
            
            imgRR = getframe(figuFD);
            imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'FramewiseDisplacement_avgTS.png']);
            saveas(figuFD,[pathFunc,dirname,filesep,'Motion_Controls_Params',filesep,'FramewiseDisplacement_avgTS'],'fig')
            close(figuFD)
            fprintf(Motfid,'BRAMILA''s average Framiwise Displacemente (FD):\t %.3f\r\n',mean(FDts));
            fprintf('Done! ============== %s\r\n',spm('time'));
            
            fprintf(Totmovtxt,'%.5f \t',mean(FDts));
        else
            fprintf(Totmovtxt,'N/A \t');
        end
        
        if get(handles.TM_FD,'Value')
            fprintf('\r\n')
            fprintf('UF²C =============== %s\n',spm('time'));
            fprintf('Creating FD Temporal Mask\r\n');

            FD_TM = FDts>str2num(get(handles.FDthres,'String'));
            fprintf('Done! ============== %s\r\n',spm('time'));
            
            fprintf(Totmovtxt,'%d \t',sum(FD_TM));
        else
            FD_TM = [];
            fprintf(Totmovtxt,'N/A \t');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Extracting Globals and Masking
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
            FI = 1;
            [sizeFunc,finalEPI,MeanWM,MeanCSF,Func] = mask_globals_uf2c(fvs,pathFunc,dirname,fileFunc{yy},fileStru{yy},get(handles.SPMbb,'Value'),FI);
        catch
            set(handles.status,'String','Error...')
            warndlg(sprintf('An error occured during subject %d globals and mask extraction. Check your data.',yy), 'Process Aborted')
            return
        end
        try
            close(Mot_fig)
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Reshapeing and Thresholding
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
            Dx = size(finalEPI,1);
            Dy = size(finalEPI,2);
            Dz = size(finalEPI,3);
            Dt = size(finalEPI,4);

            [finalEPI,histCut] = imgthres_uf2c(finalEPI,Dx,Dy,Dz,ScreSize,pathFunc,dirname);
            imgRR = getframe(histCut);
            imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Avg_Img_Histogram.tif']);
            
            newEPI1 = Func;
            newEPI1.dat.fname = [pathFunc,dirname,filesep,'sw',fileFunc{yy}];
            newEPI1.descrip = 'UF²C Thresholded Masked swEPI';
            newEPI1.dat.dtype = 'INT16-LE';
            newEPI1.dat(:,:,:,:) = finalEPI;
            create(newEPI1)
        catch
            set(handles.status,'String','Error...')
            warndlg(sprintf('An error occured during subject %d thresholding. Check your data.',yy), 'Process Aborted')
            return
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BRAMILA DVARS - Computes Derivative VARiance
        % Code from: Brain and Mind Lab at Aalto University
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        extGMvar = get(handles.extGM,'Value');
        extWMvar = get(handles.extWM,'Value');
        extCSFvar = get(handles.extCSF,'Value');
        DVARthresV = str2num(get(handles.DVARthres,'String'));

        if extGMvar || extWMvar || extCSFvar
            [dvarsC1,dvarsC2,dvarsC3,MeanDvars] = bramila_dvars_uf2c(pathFunc,dirname,finalEPI,extGMvar,extWMvar,extCSFvar,DVARthresV);
            fprintf(Motfid,'BRAMILA''s Average Derivative VARiance (DVAR) GM masked: %.3f \r\n',mean(dvarsC1));
            fprintf(Motfid,'BRAMILA''s Average Derivative VARiance (DVAR) WM masked: %.3f \r\n',mean(dvarsC2));
            fprintf(Motfid,'BRAMILA''s Average Derivative VARiance (DVAR) CSF masked: %.3f \r\n',mean(dvarsC3));
            fprintf(Totmovtxt,'%.4f \t',mean(MeanDvars));
        else
            MeanDvars = [];
            fprintf(Totmovtxt,'N/A \t');
        end
        
        if get(handles.TM_DVARs,'Value')
            fprintf('\r\n')
            fprintf('UF²C =============== %s\n',spm('time'));
            fprintf('Creating DVARs Temporal Mask\r\n');

            dvars_TM = zeros(size(finalEPI,4),1);
            dvars_TM = MeanDvars>str2num(get(handles.DVARthres,'String'));
            fprintf('Done! ============== %s\r\n',spm('time'));
            fprintf(Totmovtxt,'%d \t',sum(dvars_TM));
        else
            dvars_TM = [];
            fprintf(Totmovtxt,'N/A \t');
        end

        fclose(Motfid);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Creating the regression Matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~get(handles.checkReg,'Value')
             nrg = 0;
        end

        if get(handles.checkMOV,'Value') || get(handles.checkOsc,'Value') || get(handles.checkReg,'Value')
            try
                [Frgval,Regstr,figu,thicksL] = regrprep_uf2c(yy,get(handles.checkMOV,'Value'),pathFunc,dirname,fileFunc{yy},get(handles.checkOsc,'Value'),get(handles.checkReg,'Value'),ScreSize,MeanWM,MeanCSF,size(finalEPI,4),nrg);
                if ~isempty(figu)
                    if isequal(get(handles.TM_FD,'Value'),0) && isequal(get(handles.TM_DVARs,'Value'),0)
                        imgRR = getframe(figu);
                        imwrite(imgRR.cdata, [pathFunc,dirname,filesep,'Regression Matrix.tif']);
                    end
                end
            catch
                set(handles.status,'String','Error...')
                warndlg(sprintf('An error occured during subject %d regressors creation. Check your data.',yy), 'Process  Aborted')
                return
            end
        else
            Regstr = '';
            Frgval = [];
            thicksL = {};
        end            

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Adding the FD and DVAR regressors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TM_DVARsV = get(handles.TM_DVARs,'Value');
        TM_FDV = get(handles.TM_FD,'Value');
        DVAR_tsV = get(handles.DVAR_TS,'Value');

        if TM_DVARsV || TM_FDV || DVAR_tsV
           [Frgval,figu,Censu] = Add_FD_DVAR_uf2c(figu,thicksL,Frgval,pathFunc,dirname,TM_DVARsV,TM_FDV,DVAR_tsV,MeanDvars,FD_TM,dvars_TM,size(finalEPI,4));
            if size(Frgval,2)>1
                Regstr = [Regstr,'_tmpMask'];
            end
            if TM_DVARsV || TM_FDV
                nOfSensu = size(Censu,2);
                fprintf(Totmovtxt,'%d \r\n',nOfSensu);
            end
        else
            fprintf(Totmovtxt,'N/A \r\n');
        end

        try
            close(figu)
        end
        try
            close(histCut)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Filtering and regressions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  get(handles.FazFilt,'Value') || ~isequal(Regstr,'') || TM_DVARsV || TM_FDV
            fprintf('\r\n')
            fprintf('UF²C =============== %s\n',spm('time'));
            fprintf('Performing regressions and/or filtering\r\n');

            BinFazFilt = get(handles.FazFilt,'Value');
            finalEPI = filtregr_uf2c(finalEPI,BinFazFilt,Dx,Dy,Dz,Dt,Frgval,HdLOW,HdHIGH,aFil,tmp1,tmp2,Regstr);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Creating the FiltRegre Image
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
            if  get(handles.FazFilt,'Value') || ~isequal(Regstr,'')
                if  get(handles.FazFilt,'Value') && ~isequal(Regstr,'')
                    newEPI = Func;
                    newEPI.dat.fname = [pathFunc,dirname,filesep,'FiltRegrSW_',fileFunc{yy}];
                    newEPI.descrip = 'UF²C-Filtered_Regressed_sw_EPI';
                    newEPI.dat.dtype = 'INT16-LE';
                    newEPI.dat(:,:,:,:) = finalEPI;
                    create(newEPI)
                    fprintf('Done! ============== %s\n\r',spm('time'));
                else if get(handles.FazFilt,'Value') && isequal(Regstr,'')
                        newEPI = Func;
                        newEPI.dat.fname = [pathFunc,dirname,filesep,'FiltSW_',fileFunc{yy}];
                        newEPI.descrip = 'UF²C-Filtered_sw_EPI';
                        newEPI.dat.dtype = 'INT16-LE';
                        newEPI.dat(:,:,:,:) = finalEPI;
                        create(newEPI)
                        fprintf('Done! ============== %s\n\r',spm('time'));
                    else
                        newEPI = Func;
                        newEPI.dat.fname = [pathFunc,dirname,filesep,'RegrSW_',fileFunc{yy}];
                        newEPI.descrip = 'UF²C-Regressed_sw_EPI';
                        newEPI.dat.dtype = 'INT16-LE';
                        newEPI.dat(:,:,:,:) = finalEPI;
                        create(newEPI)
                        fprintf('Done! ============== %s\n\r',spm('time'));
                    end
                end
            else
                newEPI = Func;
                newEPI.dat.fname = [pathFunc,dirname,filesep,'ThreSW_',fileFunc{yy}];
                newEPI.descrip = 'UF²C-Masked and Thresholded_sw_EPI';
                newEPI.dat.dtype = 'INT16-LE';
                newEPI.dat(:,:,:,:) = finalEPI;
                create(newEPI)
                fprintf('Done! ============== %s\n\r',spm('time'));
            end
        catch
            set(handles.status,'String','Error...')
            warndlg(sprintf('An error occured during subject %d image creation. Check the avaiable disk space.',yy), 'Process Aborted')
            return
        end
        
    else   % case preprocessing was not performed
        try
            Func = nifti([pathFunc,fileFunc{yy}]);
            finalEPI = Func.dat(:,:,:,:);
            Dx = size(finalEPI,1);
            Dy = size(finalEPI,2);
            Dz = size(finalEPI,3);
            Dt = size(finalEPI,4);
        catch
            set(handles.status,'String','Error...')
            warndlg(sprintf('An error occured trying to read the subject %d processed image. Check your data.',yy), 'Process Aborted')
            return
        end
    end
        
    if isequal(get(handles.AddROIm,'Value'),1)
        
        fprintf('\r\n')
        fprintf('UF²C =============== %s\n',spm('time'));
        fprintf('Extracting and processing ROIs time-series\r\n');
        
        MatrixSeeds = zeros(Dt,size(fileROI,2));
        roiF = zeros(Dx,Dy,Dz); %Prelocation
        Final4Dmap = zeros(Dx,Dy,Dz,loopSize);
        finalEPItmp = reshape(finalEPI,prod([Dx,Dy,Dz]),Dt)';

        for jh = 1:loopSize
            multiWaitbar('Extracting and processing ROIs time-series', jh/loopSize);
            roi1 = squeeze(ROIMAt(:,:,:,jh));
            %%%% Loop option
            [iCo,jCo,kCo]= ind2sub(size(roi1), find(roi1>0));
            
            IniNofVox = numel(iCo);

            coords = [iCo,jCo,kCo];
            for vx = 1:Dt % get the ROI voxels time series. This step exclude from the ROI eventuals with matter voxels.
                final3D(:,:,:) = roi1.*finalEPI(:,:,:,vx);
                for hj = 1:size(coords,1)
                    vetsCor(vx,hj) = final3D(coords(hj,1),coords(hj,2),coords(hj,3));
                end
            end
            %%%%
            
            %%%%% Tester of corregistering - Uncomment to create image
            roiF =  (roi1.*finalEPI(:,:,:,1))+roiF; % add matrix value to the NIfTI object
            RoiNii = Func;   % Creates a NIfTI object
            RoiNii.dat.dim = [Dx Dy Dz]; % Apply the correct matrix size
            RoiNii.dat.fname = [pathFunc,dirname,filesep,'ROI_mask',filesep,dirname,'_Roi.nii']; % Change file name
            RoiNii.dat.dtype = 'FLOAT32-LE'; % verify matrix datatype
            RoiNii.descrip = dirname; % Apply Name information in the Nii struct (description field)
            RoiNii.dat(:,:,:) = roiF; % add matriz value to the NIfTI object
            create(RoiNii) % creates the NIfTI file from the Object

            vetsCor(:,~any(vetsCor,1) ) = []; %removes time series of zeros (outside util area)
            verifCor = mean(vetsCor,2);
            xxX = corr(vetsCor(:,:),verifCor);
            xxX(isnan(xxX)) = 0;
            
            % Option 1: Outliers (lower than 1 IQR, for "severe") would be removed         
            xxXxx = outlierdetec_uf2c(xxX,'severe','lowerside');
            vetsCorF = vetsCor;
            vetsCorF(:,xxXxx) = [];
            
            finalNofVox = size(vetsCorF,2);
            
            mean4D = (mean(vetsCorF,2))';
            
            fprintf('\r\n')
            fprintf('UF²C =============== %s\n',spm('time'));
            fprintf('Performing Cross-Correlations\n');

            map = zeros(Dt/WS,prod([Dx,Dy,Dz]));%Prelocation
            Vn = 1;
            for MA = 1:WS:Dt
                map(Vn,:) = corr(finalEPItmp(MA:MA+(WS-1),:),mean4D(MA:MA+(WS-1))');
                Vn = Vn+1;
            end
            map(isnan(map)) = 0;
            map = reshape(map',Dx,Dy,Dz,Dt/WS);
            
            Final4Dmap(:,:,:,jh) = squeeze(map);
            
            fprintf('Done! ============== %s\r\n',spm('time'));
            
            Rois_N_of_Voxels__Init_vs_Final(jh,1) = IniNofVox;
            Rois_N_of_Voxels__Init_vs_Final(jh,2) = finalNofVox;

            clear coords iCo jCo kCo vetsCor verifCor xxX xxX_Z stdxxX meanxxX lowerCut CutFinal mean4D
        end
        fprintf('Done! ============== %s\r\n',spm('time'));
    else
        fprintf('\r\n')
        fprintf('UF²C =============== %s\n',spm('time'));
        fprintf('Extracting and processing ROIs Time-series\r\n');

         loopSize = size(seedList,1);
         MatrixSeeds = zeros(Dt,size(seedList,1));
         roiF = zeros(Dx,Dy,Dz); %Prelocation
        Final4Dmap = zeros(Dx,Dy,Dz,loopSize);
        finalEPItmp = reshape(finalEPI,prod([Dx,Dy,Dz]),Dt)';
         for jh = 1:loopSize
            [iCo,jCo,kCo]= ind2sub(size(roi1Fall(:,:,:,jh)), find(roi1Fall(:,:,:,jh)>0));
            coords = [iCo,jCo,kCo];
            IniNofVox = numel(iCo);
            
            for vx = 1:Dt % get the ROI voxels time series. This step exclude from the ROI eventuals with matter voxels.
                final3D(:,:,:) = roi1Fall(:,:,:,jh).*finalEPI(:,:,:,vx);
                for hj = 1:size(coords,1)
                    vetsCor(vx,hj) = final3D(coords(hj,1),coords(hj,2),coords(hj,3));
                end
            end

            %%%%% Tester of corregistering
            roiF =  (roi1Fall(:,:,:,jh).*finalEPI(:,:,:,1))+roiF; % add matriz value to the NIfTI object
            %%%%% Tester of corregistering
            
            vetsCor(:,~any(vetsCor,1)) = [];
            verifCor = mean(vetsCor,2);
            xxX = corr(vetsCor(:,:),verifCor);
            xxX(isnan(xxX)) = 0;
            
            % Option 1: Outliers (lower than 1 IQR, for "severe") would be removed         
            xxXxx = outlierdetec_uf2c(xxX,'severe','lowerside');
            vetsCorF = vetsCor;
            vetsCorF(:,xxXxx) = [];
            
            finalNofVox = size(vetsCorF,2);

            mean4D = (mean(vetsCorF,2))';
            
            fprintf('\r\n')
            fprintf('UF²C =============== %s\n',spm('time'));
            fprintf('Performing Cross-Correlations\n');

            map = zeros(Dt/WS,prod([Dx,Dy,Dz]));%Prelocation
            Vn = 1;
            for MA = 1:WS:Dt
                map(Vn,:) = corr(finalEPItmp(MA:MA+(WS-1),:),mean4D(MA:MA+(WS-1))');
                Vn = Vn+1;
            end
            map(isnan(map)) = 0;
            map = reshape(map',Dx,Dy,Dz,Dt/WS);
            
            Final4Dmap(:,:,:,jh) = squeeze(map);
            
            fprintf('Done! ============== %s\r\n',spm('time'));
            
            Rois_N_of_Voxels__Init_vs_Final(jh,1) = IniNofVox;
            Rois_N_of_Voxels__Init_vs_Final(jh,2) = finalNofVox;

            clear coords iCo jCo kCo vetsCor verifCor xxX stdxxX meanxxX lowerCut CutFinal mean4D xxXxx
        end
        
        RoiNii = Func;   % Creates a NIfTI object
        RoiNii.dat.dim = [Dx Dy Dz]; % Apply the correct matrix size
        RoiNii.dat.fname = [pathFunc,dirname,filesep,'ROI_mask',filesep,dirname,'_Roi.nii']; % Change file name
        RoiNii.dat.dtype = 'FLOAT32-LE'; % verify matrix datatype
        RoiNii.descrip = dirname; % Apply Name information in the Nii struct (description field)
        RoiNii.dat(:,:,:) = roiF; % add matriz value to the NIfTI object
        create(RoiNii) % creates the NIfTI file from the Object
        fprintf('Done! ============== %s\r\n',spm('time'));
        
        
        FM4DNii = Func;   % Creates a NIfTI object
        FM4DNii.dat.dim = [Dx,Dy,Dz,loopSize]; % Apply the correct matrix size
        FM4DNii.dat.fname = [pathFunc,dirname,filesep,'Correlation_map',filesep,dirname,'_AllSeeds4D.nii']; % Change file name
        FM4DNii.dat.dtype = 'FLOAT32-LE'; % verify matrix datatype
        FM4DNii.descrip = dirname; % Apply Name information in the Nii struct (description field)
        FM4DNii.dat(:,:,:,:) = Final4Dmap; % add matriz value to the NIfTI object
        create(FM4DNii) % creates the NIfTI file from the Object
        fprintf('Done! ============== %s\r\n',spm('time'));

    end
    
    save([pathFunc,dirname,filesep,'Rois_NofVoxels_Init_vs_Final.mat'],'Rois_N_of_Voxels__Init_vs_Final')
    save([pathFunc,dirname,filesep,'Seeds_TimeSeries.mat'],'MatrixSeeds')
    Vn = 1;
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Done! ============== %s\r\n',spm('time'));
    clear MatrixSeeds mean4D vetsCorF

    if isequal(get(handles.OutputFil, 'Value'),1)
        if isequal(get(handles.FazPrep,'Value'),0)
            delete([pathFunc dirname filesep fileStru{yy}])
            delete([pathFunc dirname filesep fileFunc{yy}])
            delete([pathFunc dirname filesep fileFunc{yy}(1:end-4) '.mat'])
            delete([pathFunc dirname filesep fileStru{yy}(1:end-4) '_seg_inv_sn.mat'])
            delete([pathFunc dirname filesep fileStru{yy}(1:end-4) '_seg_sn.mat'])
            delete([pathFunc dirname filesep fileStru{yy}(1:end-4) '_seg8.mat'])
            delete([pathFunc dirname filesep 'mean' fileFunc{yy}])
            delete([pathFunc dirname filesep 'sw' fileFunc{yy}])
            delete([pathFunc dirname filesep 'y_' fileFunc{yy}])
            delete([pathFunc dirname filesep 'c1' fileStru{yy}])
            delete([pathFunc dirname filesep 'c2' fileStru{yy}])
            delete([pathFunc dirname filesep 'c3' fileStru{yy}])
            delete([pathFunc dirname filesep 'm' fileStru{yy}])
            delete([pathFunc dirname filesep 'mwc1' fileStru{yy}])
            delete([pathFunc dirname filesep 'mwc2' fileStru{yy}])
            delete([pathFunc dirname filesep 'mwc3' fileStru{yy}])
            delete([pathFunc dirname filesep 'wc1' fileStru{yy}])
            delete([pathFunc dirname filesep 'wc2' fileStru{yy}])
            delete([pathFunc dirname filesep 'wc3' fileStru{yy}])
            delete([pathFunc dirname filesep 'w' fileFunc{yy}])
            if STC
                delete([pathFunc dirname filesep 'a' fileFunc{yy}])
                delete([pathFunc dirname filesep 'wa' fileFunc{yy}])
            end
            delete([pathFunc dirname filesep 'w' fileFunc{yy}(1:end-4) '.mat'])
            delete([pathFunc dirname filesep 'sw' fileFunc{yy}(1:end-4) '.mat'])
            delete([pathFunc dirname filesep 'c1interp.nii'])
            delete([pathFunc dirname filesep 'sc1interp.nii'])
            delete([pathFunc dirname filesep 'c1interpF.nii'])
            delete([pathFunc dirname filesep 'c2interp.nii'])
            delete([pathFunc dirname filesep 'c3interp.nii'])
            delete([pathFunc dirname filesep 'Template.nii'])
            delete([pathFunc dirname filesep 'Template.mat'])
        end
    end
    
    fprintf('\r\n')
    fprintf('UF²C =============== %s\n',spm('time'));
    fprintf('Subject %d Done!\r\n',yy);

    try
        close(fig)
    end
    try
        close(figi)
    end
    try
        close(figxT)
    end
    try
        delete('Global_SL.mat')
    end
    clear('ex1_1','ex2_1','ex1_2','ex2_2','ex3_1','ex3_2','ex4_2','ex4_1')
    fprintf('\r\n')

end

%%%%%%%%%%%% END OF THE SUBJECT LOOP
        %%%%%%% 
        %%%%%%% 
        %%%%%%% 
        %%%%%%% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yytmp = yy;        
yy = yy-ExclSub; %to exclude possible volunteers that were excluded from the number of subjects

if yy == 0
    set(handles.status,'String','Aborted');
    fprintf('Process aborted because there is nothing to do!\n')
    fprintf('All subjects analyses were aborted due ROIs spatial mismatch.\n')
    fprintf('Number of subjects included: %d. Number of aborted subjects: %d.\n',yytmp,ExclSub)
    multiWaitbar('CloseAll'); %close the progressbar window
    fclose('all');
    fprintf('==========================================\r\n');
    return
end

fprintf('\r\n')
fprintf('UF²C =============== %s\n',spm('time'));
fprintf('Generating overall results\r\n');

multiWaitbar('CloseAll'); %close the progressbar window
fclose('all');
save([pathFunc,foldername,filesep,'Total_Variable'],'WinAvSave')
set(handles.status,'String','Done!!!')
fprintf('\r\n')
fprintf('All Done! ========== %s\n',spm('time'));
fprintf('==========================================\r\n');

function FazPrep_Callback(hObject, eventdata, handles)
if get(handles.FazPrep,'Value')
    set(handles.AddStru,'Enable','off')
    set(handles.FazFilt,'Enable','off')
    set(handles.LPF,'Enable','off')
    set(handles.HPF,'Enable','off')
    set(handles.stba,'Enable','off')
    set(handles.checkMOV,'Enable','off')
    set(handles.checkOsc,'Enable','off')
    set(handles.checkReg,'Enable','off')
    set(handles.text8,'Enable','off')
    set(handles.text9,'Enable','off')
    set(handles.text20,'Enable','off')
    set(handles.text19,'Enable','off')
    set(handles.text21,'Enable','off')
    set(handles.STCo,'Enable','off')
    set(handles.AddFunc,'String','Add FiltRegrSW Files')
    set(handles.FDcheck,'Enable','off')
    set(handles.extGM,'Enable','off')
    set(handles.extWM,'Enable','off')
    set(handles.extCSF,'Enable','off')
    set(handles.TM_FD,'Enable','off')
    set(handles.TM_DVARs,'Enable','off')
    set(handles.DVAR_TS,'Enable','off')
else
    set(handles.AddStru,'Enable','on')
    set(handles.FazFilt,'Enable','on')
    set(handles.LPF,'Enable','on')
    set(handles.HPF,'Enable','on')
    set(handles.stba,'Enable','on')
    set(handles.checkMOV,'Enable','on')
    set(handles.checkOsc,'Enable','on')
    set(handles.checkReg,'Enable','on')
    set(handles.text8,'Enable','on')
    set(handles.text9,'Enable','on')
    set(handles.text20,'Enable','on')
    set(handles.text19,'Enable','on')
    set(handles.text21,'Enable','on')
    set(handles.STCo,'Enable','on')
    set(handles.AddFunc,'String','Add Functional Files')
    set(handles.FDcheck,'Enable','on')
    set(handles.extGM,'Enable','on')
    set(handles.extWM,'Enable','on')
    set(handles.extCSF,'Enable','on')
    set(handles.TM_FD,'Enable','on')
    set(handles.TM_DVARs,'Enable','on')
    set(handles.DVAR_TS,'Enable','on')

end

function wsIn_Callback(hObject, eventdata, handles)
global numofVoluF ratioVol
smWin = str2num(get(handles.wsIn,'String'));
ratioVol = numofVoluF/smWin;

if ~isequal(ratioVol,round(ratioVol)) || smWin==0
    warndlg('The "Windows size" value can not be zero and should to be multiple of the "Number of Dynamics"','Attention')
    set(handles.wsIn,'String',num2str(numofVoluF))
    set(handles.fancyRes,'Enable','on')
    set(handles.fancyRes,'Value',1)
else if ~isequal(numofVoluF,smWin)
        set(handles.fancyRes,'Enable','off')
        set(handles.fancyRes,'Value',0)
    else
        set(handles.fancyRes,'Enable','on')
        set(handles.fancyRes,'Value',1)
    end
end

function SPMbb_Callback(hObject, eventdata, handles)
if isequal(get(handles.SPMbb,'Value'),1)
    set(handles.SSInp,'String','2 2 2')
else
    set(handles.SSInp,'String','4 4 4')
end


function pushbutton6_Callback(hObject, eventdata, handles)
global fileFunc

UF2Cdir = which('uf2c');
tmpDIR = [UF2Cdir(1:end-6) 'Analysis' filesep 'FC_tmp' filesep];

try
    delete([tmpDIR 'Subject_List_Func.txt']);
end

subL = transpose(fileFunc);
sList = size(subL,1);
fideSL = fopen([tmpDIR 'Subject_List_Func.txt'],'w+');

for lt = 1:sList
    fprintf(fideSL,'%d - %s\r\n',lt,fileFunc{lt});
end

fclose(fideSL)
open([tmpDIR 'Subject_List_Func.txt'])

function pushbutton7_Callback(hObject, eventdata, handles)
global fileStru

UF2Cdir = which('uf2c');
tmpDIR = [UF2Cdir(1:end-6) 'Analysis' filesep 'FC_tmp' filesep];

try
    delete([tmpDIR 'Subject_List_Stru.txt']);
end

subL = transpose(fileStru);
sList = size(subL,1);
fideSL = fopen([tmpDIR 'Subject_List_Stru.txt'],'w+');

for lt = 1:sList
    fprintf(fideSL,'%d - %s\r\n',lt,fileStru{lt});
end

fclose(fideSL)
open([tmpDIR 'Subject_List_Stru.txt'])

function HPF_Callback(hObject, eventdata, handles)
function FazFilt_Callback(hObject, eventdata, handles)
function TRInp_Callback(hObject, eventdata, handles)
function OutputFil_Callback(hObject, eventdata, handles)
function checkMOV_Callback(hObject, eventdata, handles)
function checkOsc_Callback(hObject, eventdata, handles)
function LPF_Callback(hObject, eventdata, handles)
function stba_Callback(hObject, eventdata, handles)
function SSInp_Callback(hObject, eventdata, handles)
function fancyRes_Callback(hObject, eventdata, handles)
function STCo_Callback(hObject, eventdata, handles)
function DVARthres_Callback(hObject, eventdata, handles)
function FDthres_Callback(hObject, eventdata, handles)

function FDcheck_Callback(hObject, eventdata, handles)
if get(handles.FDcheck,'Value')
    set(handles.TM_FD,'Value',0)
    set(handles.TM_FD,'Enable','on')
    set(handles.text28,'Enable','on')
    set(handles.FDthres,'Enable','on')
else
    set(handles.TM_FD,'Value',0)
    set(handles.TM_FD,'Enable','off')
    set(handles.text28,'Enable','off')
    set(handles.FDthres,'Enable','off')
end

function extGM_Callback(hObject, eventdata, handles)
if get(handles.extGM,'Value')
    set(handles.TM_DVARs,'Enable','on')
    set(handles.DVAR_TS,'Enable','on')
    set(handles.DVAR_TS,'Value',1)
else
    if get(handles.extWM,'Value') || get(handles.extCSF,'Value')
        set(handles.TM_DVARs,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)
    else
        set(handles.TM_DVARs,'Enable','off')
        set(handles.TM_DVARs,'Value',0)
        set(handles.text27,'Enable','off')
        set(handles.DVARthres,'Enable','off')
        set(handles.DVAR_TS,'Enable','off')
        set(handles.DVAR_TS,'Value',0)
    end
end
    
function extWM_Callback(hObject, eventdata, handles)
if get(handles.extWM,'Value')
    set(handles.TM_DVARs,'Enable','on')
    set(handles.DVAR_TS,'Enable','on')
    set(handles.DVAR_TS,'Value',1)
else
    if get(handles.extGM,'Value') || get(handles.extCSF,'Value')
        set(handles.TM_DVARs,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
        set(handles.DVAR_TS,'Value',1)
    else
        set(handles.TM_DVARs,'Enable','off')
        set(handles.TM_DVARs,'Value',0)
        set(handles.text27,'Enable','off')
        set(handles.DVARthres,'Enable','off')
        set(handles.DVAR_TS,'Enable','off')
        set(handles.DVAR_TS,'Value',0)
    end
end

function extCSF_Callback(hObject, eventdata, handles)
if get(handles.extCSF,'Value')
    set(handles.TM_DVARs,'Enable','on')
    set(handles.TM_DVARs,'Enable','on')
    set(handles.DVAR_TS,'Enable','on')
else
    if get(handles.extGM,'Value') || get(handles.extWM,'Value')
        set(handles.TM_DVARs,'Enable','on')
        set(handles.TM_DVARs,'Enable','on')
        set(handles.DVAR_TS,'Enable','on')
    else
        set(handles.TM_DVARs,'Enable','off')
        set(handles.TM_DVARs,'Value',0)
        set(handles.text27,'Enable','off')
        set(handles.DVARthres,'Enable','off')
        set(handles.DVAR_TS,'Enable','off')
        set(handles.DVAR_TS,'Value',0)
    end
end

function TM_DVARs_Callback(hObject, eventdata, handles)
if get(handles.TM_DVARs,'Value')
    set(handles.text27,'Enable','on')
    set(handles.DVARthres,'Enable','on')
else
    set(handles.text27,'Enable','off')
    set(handles.DVARthres,'Enable','off')
end

function TM_FD_Callback(hObject, eventdata, handles)
if get(handles.TM_FD,'Value')
    set(handles.text28,'Enable','on')
    set(handles.FDthres,'Enable','on')
else
    set(handles.text28,'Enable','off')
    set(handles.FDthres,'Enable','off')
end

function SSInp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function TRInp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function HPF_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function LPF_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function stba_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function wsIn_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function FDthres_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function DVARthres_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DVAR_TS_Callback(hObject, eventdata, handles)
