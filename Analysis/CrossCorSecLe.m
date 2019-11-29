function varargout = CrossCorSecLe(varargin)
% UF²C M-file for CrossCorSecLe.fig
% UF²C - User Friendly Functional Connectivity
% Brunno Machado de Campos
% University of Campinas, 2018
%
% Copyright (c) 2018, Brunno Machado de Campos
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
                   'gui_OpeningFcn', @CrossCorSecLe_OpeningFcn, ...
                   'gui_OutputFcn',  @CrossCorSecLe_OutputFcn, ...
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

function CrossCorSecLe_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = CrossCorSecLe_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function push1_Callback(hObject, eventdata, handles)
global g1_VAR g2OK g1OK dirOK nodeOK

[fileG1,pathG1] = uigetfile({'*.mat','Matlab file'},...
        'Select the result file (All_Subjs-VAR.mat) for G1nameT','MultiSelect','off');
if ~isequal(fileG1,0)
    load([pathG1 fileG1])
    try
        g1_VAR = AllSubj3D;
    catch
        g1_VAR = TotalPairWise;
    end
    set(handles.text4,'String',sprintf('Data from %d seeds and %d subjects were added',size(g1_VAR,1),size(g1_VAR,3)))
    drawnow
    g1OK = 1;
    if g1OK && g2OK && dirOK && nodeOK
        set(handles.runB,'Enable','on')
    else
        if g1OK && g2OK && dirOK && get(handles.checkIfNet,'Value')
            set(handles.runB,'Enable','on')
        end
    end
    if g1OK && get(handles.checkIfNet,'Value')
        if size(g1_VAR,1)>1
            for nt = 1:size(g1_VAR,1)
                strNetList{nt,:} = sprintf('Network %d',nt);
            end
            set(handles.TarNetPop,'String',strNetList)
            if get(handles.tarNet,'Value')
                set(handles.TarNetPop,'Enable','on')
            end
        end
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
    set(handles.text5,'String',sprintf('Data from %d seeds and %d subjects were added',size(g2_VAR,1),size(g2_VAR,3)))
    drawnow
    g2OK = 1;
    if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
        set(handles.runB,'Enable','on')
    else
        if g1OK && g2OK && dirOK && get(handles.checkIfNet,'Value')
            set(handles.runB,'Enable','on')
        end
    end
    if g1OK && get(handles.checkIfNet,'Value')
        if size(g1_VAR,1)>1
            for nt = 1:size(g1_VAR,1)
                strNetList{nt,:} = sprintf('Network %d',nt);
            end
            set(handles.TarNetPop,'String',strNetList)
            drawnow
            if get(handles.tarNet,'Value')
                set(handles.TarNetPop,'Enable','on')
            end
        end
    end
    drawnow
end

function COV_Callback(hObject, eventdata, handles)
global COVDATA nOfLines nOfCovs g2OK g1OK dirOK nodeOK g1_VAR g2_VAR

if get(handles.COV,'Value')
    set(handles.TSTTbtm,'String','ANCOVA')
    [COVlist,COVpath] = uigetfile({'*.txt','Text files';'*.mat','3D mat file'},'Select a text file or a 3D mat (ROI-to-ROI cov) file with all covariates','MultiSelect','off');
    if ~isequal(COVlist,0)
        COVDATA = importdata([COVpath COVlist]);
        nOfLines = size(COVDATA,1);
        nOfCovs = size(COVDATA,2);
        if isequal(g1OK,1) && isequal(g2OK,1)
            if isequal((size(g1_VAR,3)+size(g2_VAR,3)),nOfLines)
                set(handles.text9,'String',sprintf('%d covariate(s) added!',nOfCovs))
                drawnow
            else
                warndlg('The number of inputs in each covariate column should to match with the sum of Group 1 and 2 subjects','Ops!')
                set(handles.COV,'Value',0)
                set(handles.text9,'String','')
                drawnow
            end
        else
            set(handles.text9,'String',sprintf('%d covariate(s) added!',nOfCovs))
            drawnow
        end
        if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
            set(handles.runB,'Enable','on')
        end
    else
        set(handles.COV,'Value',0)
        set(handles.text9,'String','')
        set(handles.TSTTbtm,'String','Two Sample t-test')
        drawnow
    end
else
    set(handles.TSTTbtm,'String','Two Sample t-test')
    set(handles.text9,'String','')
    drawnow
end

function OutDir_Callback(hObject, eventdata, handles)
global Dir g2OK g1OK  dirOK nodeOK DirTmp

DirTmp = uigetdir('','Define the output directory');
if ~isequal(DirTmp,0)
    Dir = [DirTmp filesep 'UF2C_SL_' get(handles.G1name,'String') 'x' get(handles.G2name,'String')];
    set(handles.text10,'String',Dir)
    dirOK = 1;
    set(handles.text10,'Enable','on')
    
    if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
        set(handles.runB,'Enable','on')
    else
        if g1OK && g2OK && dirOK && get(handles.checkIfNet,'Value')
            set(handles.runB,'Enable','on')
        end
    end
    drawnow
end

function addNODES_Callback(hObject, eventdata, handles)
global CoordMS seedList nodeOK g1OK g2OK dirOK nEle nOfNets REtype

if get(handles.CoodList,'Value')
    set(handles.text12,'String','Wait.......')
    drawnow
    [fileCor,pathCor] = uigetfile({'*.txt','Text file'},...
            'Select the cordinates list file','MultiSelect','off');
        if ~isequal(fileCor,0)
            seedListT = importdata([pathCor fileCor]);
            if iscell(seedListT)
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

            else
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
                
                seedListTMP = seedListT.data;
                clear seedListT
                seedList = seedListTMP;
            end
            
            TcM = [-2 0 0 92;0 2 0 -128;0 0 2 -74;0 0 0 1];
            
            for jh = 1:size(seedList,1)
                EPIcoord = [seedList(jh,1) seedList(jh,2) seedList(jh,3) 1]*(inv(TcM))';
                EPIcoord = round(EPIcoord);
                CoordMS(jh,:) = EPIcoord(1,:);
            end
            
            CoordMS(:,1) = abs(CoordMS(:,1)-91); % because Matlab plots images in radiological orientation (inverting L-R)
            
            set(handles.text12,'String',sprintf('%d coordinates of %d networks were added!',size(CoordMS,1),nOfNets))
            drawnow
            if nOfNets>1
                if get(handles.NetLab,'Value')
                    for nt = 1:nOfNets
                        strNetList{nt,:} = NetLabels{nt};
                    end
                else
                    for nt = 1:nOfNets
                        strNetList{nt,:} = sprintf('Network %d',nt);
                    end
                end
                set(handles.TarNetPop,'String',strNetList)
                set(handles.NetLab,'Enable','on')
            end
            nodeOK = 1;
            if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
                set(handles.runB,'Enable','on')
            end
            
            RNq = questdlg('How do you want to enumerate rois ROIs?',...
                'ROIs numbering',...
                'Ascending order from the first to last ROI',...
                'Ascending order inside each network',...
                'Ascending order from the first to last ROI');
            switch RNq
                case 'Ascending order from the first to last ROI'
                    REtype = 1;
                case 'Ascending order inside each network'
                    REtype = 2;
            end
            
        end
        drawnow
else
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
        CoordMS = FinalCentrCoord;
        
        DimROI = struRoi.dat.dim;
        DimROIDiff = [91 109 91]./DimROI;
        
        for hg = 1:size(CoordMS,1)
            CoordMS(hg,:) = CoordMS(hg,:).*DimROIDiff;
        end
        CoordMS = round(CoordMS);
        
        seedList = transpose(fileROI);
        set(handles.text12,'String',sprintf('%d coordinates of %d networks were added!',size(CoordMS,1),nOfNets))
        drawnow
        
        if nOfNets>1
            if get(handles.NetLab,'Value')
                for nt = 1:nOfNets
                    strNetList{nt,:} = NetLabels{nt};
                end
            else
                for nt = 1:nOfNets
                    strNetList{nt,:} = sprintf('Network %d',nt);
                end
            end
            set(handles.TarNetPop,'String',strNetList)
            set(handles.NetLab,'Enable','on')
        end
        nodeOK = 1;
        
        if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
            set(handles.runB,'Enable','on')
        end
        RNq = questdlg('How do you want to enumerate rois ROIs?',...
            'ROIs numbering',...
            'Ascending order from the first to last ROI',...
            'Ascending order inside each network',...
            'Ascending order from the first to last ROI');
        
        switch RNq
            case 'Ascending order from the first to last ROI'
                REtype = 1;
            case 'Ascending order inside each network'
                REtype = 2;
        end
        
    end
    drawnow
end
    
function runB_Callback(hObject, eventdata, handles)
global CoordMS seedList Dir COVDATA g1_VAR g2_VAR nOfLines nEle nOfNets REtype NetLabels

set(handles.text13,'String','Running.....')
drawnow
%%%%%% Checking...
if get(handles.COV,'Value')
    if ~isequal((size(g1_VAR,3)+size(g2_VAR,3)),nOfLines)
        warndlg(sprintf('The number of data on each covariate (%d) should to match with the sum of the number of subjects on Group 1 and 2',nOfLines,(size(g1_VAR,3)+size(g2_VAR,3))),'Ops! Process Aborted')
        return
    end
end

if ~isequal(size(g1_VAR,1),size(g2_VAR,1))
    warndlg(sprintf('The number of ROIs included on Group 1 first level analysis (%d) should to match with the number of Group 2 (%d)',size(g1_VAR,1),size(g2_VAR,1)),'Ops! Process Aborted')
    return
end
if ~isequal(size(CoordMS,1),size(g1_VAR,1))
    if logical(get(handles.Grapho2D,'Value')) || logical(get(handles.Graphos3D,'Value'))
        warndlg(sprintf('The number of nodes added (%d) should to match with the used on first level (%d)',size(CoordMS,1),size(g1_VAR,1)),'Ops! Process Aborted')
    return
    end
end

if get(handles.checkIfNet,'Value')
    AnaTypeStr = '_Netwise';
    nEle = ones(1,size(g1_VAR,1));
    nOfNets = size(g1_VAR,1);
else
    AnaTypeStr = '_Seedwise';
end

if ~get(handles.NetLab,'Value')
    NetLabels = strcat('n',strsplit(num2str(1:nOfNets)));
end

if get(handles.tarNet,'Value')
    TarNetC = get(handles.TarNetPop,'Value');
    TarNetCnEle = nEle(TarNetC);
    TarNetCInit = sum(nEle(1:TarNetC-1))+1;
    g1_VAR2 = g1_VAR(TarNetCInit:TarNetCInit+nEle(TarNetC)-1,:,:);
    g2_VAR2 = g2_VAR(TarNetCInit:TarNetCInit+nEle(TarNetC)-1,:,:);
    fact = TarNetCInit;
    if ~get(handles.NetLab,'Value')
        TarNetStr = ['_TargetNet_' NetLabels{TarNetC}(2:end)];
    else
        TarNetStr = ['_TargetNet_' NetLabels{TarNetC}];
    end
else
    g1_VAR2 = g1_VAR;
    g2_VAR2 = g2_VAR;
    fact = 1;
    TarNetStr = [''];
    TarNetC = [];
end

DirF = [Dir AnaTypeStr TarNetStr];

mkdir(DirF)

imgRR = getframe(CrossCorSecLe);
imwrite(imgRR.cdata, [DirF filesep '1-Your_Choices.png']);


% % Find Critial r-score value
% dg = 180;
% rvet = 0:0.001:1;
% yind = 0;
% TY = 0;
% 
% while TY == 0
%     yind = yind+1;
%     t = rvet(yind)*sqrt((dg-2)/(1-rvet(yind)^2));
%     p = 2*(1-tcdf(t,dg-2));
%     if p < 0.05
%         TY = 1;
%     end
% end
% 
% Crit_Z = 0.5.*(log(1+rvet(yind-1)) - log(1-rvet(yind-1))) %Z transf
% 
% g1_VAR(g1_VAR<Crit_Z & g1_VAR>-Crit_Z) = 0;
% g2_VAR(g2_VAR<Crit_Z & g2_VAR>-Crit_Z) = 0;

if get(handles.TSTTbtm,'Value') %twoSample Ttest
    for i = 1:size(g1_VAR2,1)
        for j = 1:size(g1_VAR2,2)
            if i == j - (fact - 1)
                Pmap(i,j) = 1;
            else if sum(squeeze(g1_VAR2(i,j,:)))==0 || sum(squeeze(g2_VAR2(i,j,:)))==0
                   Pmap(i,j) = 1;
                else
                    if get(handles.COV,'Value')
                        [T,P] = mancovan([squeeze(g1_VAR2(i,j,:));squeeze(g2_VAR2(i,j,:))],...
                            [ones(size(g1_VAR2,3),1);zeros(size(g2_VAR2,3),1)],COVDATA,'group-group');
                         Pmap(i,j) = P(1);
                    else
                        [T,P] = mancovan([squeeze(g1_VAR2(i,j,:));squeeze(g2_VAR2(i,j,:))],...
                            [ones(size(g1_VAR2,3),1); zeros(size(g2_VAR2,3),1)],[],'group-group');
                        Pmap(i,j) = P;
                    end
                 end
            end
        end
    end
end

if get(handles.PTTbtm,'Value') %Paired ttest
    for i = 1:size(g1_VAR2,1)
        for j = 1:size(g1_VAR2,2)
            if i==j
                Pmap(i,j) = 1;
            else
                if sum(squeeze(g1_VAR2(i,j,:)))==0 || sum(squeeze(g2_VAR2(i,j,:)))==0
                     Pmap(i,j) = 1;
                else
                     [H,P,CI,tvalue] = ttest(squeeze(g1_VAR2(i,j,:)),squeeze(g2_VAR2(i,j,:)));
                     Pmap(i,j) = P;
                 end
            end
        end
    end
end

g1_mean = mean(g1_VAR2,3);
g2_mean = mean(g2_VAR2,3);

MultMean = g1_mean.*g2_mean;
MultMean(MultMean > 0) = 1;
MultMean(MultMean < 0) = -1;

ScreSize = get(0,'screensize');
ScreSize = ScreSize(3:end);

G1nameT = get(handles.G1name,'String');
G2nameT = get(handles.G2name,'String');

% Quantifying ROIs for each networks
stIndx = 1;
for hg = 1:size(nEle,2)
    str_nEle(hg) = stIndx;
    stIndx = stIndx+nEle(hg);
end

% Uncorrected P Map
if get(handles.chec1,'Value')
    
    unCorPmap = Pmap;
    unCorPmap(unCorPmap>str2num(get(handles.edit1,'String'))) = 1;
    unCorPmap = -1.*unCorPmap;
    if get(handles.tarNet,'Value')
        
        unCorPmapTrim = unCorPmap(:,str_nEle(TarNetC):str_nEle(TarNetC)+nEle(TarNetC)-1);
        unCorPmapTrim = tril(unCorPmapTrim,-1);
        unCorPmapTrim(unCorPmapTrim==0) = -1;
        unCorPmap(:,str_nEle(TarNetC):str_nEle(TarNetC)+nEle(TarNetC)-1) = unCorPmapTrim;
    end
    
        
    if ~isequal(numel(unique(unCorPmap)),1)
        if get(handles.CorMats,'Value')   %%%%% Correlation Matrix
            fig1 = figure;
            imagesc(unCorPmap)

            set(fig1,'Name','Uncorrected Correlation Matrix (p-value)','Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.45 ScreSize(1)*.4]),...
                    'Color',[1 .949 .867]);
            title('Uncorrected Correlation Matrix','FontSize',14);
            xlabel('Seeds'  ,'FontSize',14);
            ylabel('Seeds','FontSize',14);
            set(gca, 'color', [0 0 0])
            
            try  set(gca, 'gridcolor', [1 1 1]); end
            try set(gca, 'gridalpha', 1); end
            
            daspect([1 1 1])
            colorbar
            caxis([-1.*(str2num(get(handles.edit1,'String'))) (max(max(unCorPmap))-(0.2*max(max(unCorPmap))))]) %reescala a volorbar para fugir dos extremos
            
            drawnow
            imgRR = getframe(fig1);
            imwrite(imgRR.cdata, [DirF filesep 'Uncorrected_P_Map' '.png']);
            clear imgRR
            saveas(fig1,[DirF filesep 'Uncorrected_P_Map'],'fig')
            save([DirF filesep 'unCorRmap-VAR'],'unCorPmap')
            close(fig1)
        end
    
        matrixR = unCorPmap;
        matrixR(matrixR==-1) = 0;
        matrixR = abs(matrixR);

        % Creating Output text file
        matrix_bin = matrixR>0;

        fidRes = fopen([DirF,filesep,'Uncorr_Results_Description.txt'],'w+');
        
        for r = 1:size(g1_VAR2,1)
            Ridx = find(matrix_bin(r,:));
            if sum(Ridx)>0
                if get(handles.tarNet,'Value')
                    fprintf(fidRes,'ROI %d: r%d%s: \r\n',r+fact-1,r,NetLabels{TarNetC});
                else
                    tester1 = r<=(str_nEle+nEle)-1;
                    networkNumb1 = find(tester1==1,1);
                    try
                        RnUmb1 = r - sum(nEle(1:networkNumb1-1));
                    catch
                        RnUmb1 = r;
                    end
                    if isequal(REtype,2)
                        fprintf(fidRes,'ROI %d: r%d%s: \r\n',r,RnUmb1,NetLabels{networkNumb1});
                    else
                        fprintf(fidRes,'ROI %d: r%d%s: \r\n',r,r,NetLabels{networkNumb1});
                    end
                end
                for ixR = 1:numel(Ridx)
                    if MultMean(r,Ridx(ixR))>0
                         if abs(g1_mean(r,Ridx(ixR)))>abs(g2_mean(r,Ridx(ixR)))
                             strDir = [G1nameT ' > ' G2nameT];
                         else
                             strDir = [G2nameT ' > ' G1nameT];
                         end
                    else
                        strDir = [G2nameT ' ' '<->' ' -(' G1nameT ')'];
                    end
                    tester2 = Ridx(ixR)<=(str_nEle+nEle)-1;
                    networkNumb2 = find(tester2==1,1);
                    try
                        RnUmb2 = Ridx(ixR) - sum(nEle(1:networkNumb2-1));
                    catch
                        RnUmb2 = Ridx(ixR);
                    end
                    if isequal(REtype,2)
                        fprintf(fidRes,'\t ROI %d: r%d%s \t %s\r\n',Ridx(ixR),RnUmb2,NetLabels{networkNumb2},strDir);
                    else
                        fprintf(fidRes,'\t ROI %d: r%d%s \t %s\r\n',Ridx(ixR),Ridx(ixR),NetLabels{networkNumb2},strDir);
                    end
                end
            end
        end
        fclose(fidRes);
        
        limiteG = (size(g1_VAR2,2)/2).*1.1;

        RAIO = limiteG/2.2;
        ang=0:0.01:2*pi; 
        xp=RAIO*cos(ang);
        yp=RAIO*sin(ang);
        xpt=(limiteG/2.15)*cos(ang);
        ypt=(limiteG/2.15)*sin(ang);

        AnoRa = size(xp,2)/size(matrixR,2);
        xp2 = xp(1:AnoRa:end);
        yp2 = yp(1:AnoRa:end);

        xp2t = xpt(1:AnoRa:end);
        yp2t = ypt(1:AnoRa:end);
        xp2t = xp2t+(limiteG/2);
        yp2t = yp2t+(limiteG/2);
        xp2t = xp2t(end:-1:1);
        yp2t = yp2t(end:-1:1);

        xp2 = xp2+(limiteG/2);
        yp2 = yp2+(limiteG/2);
        yp2 = yp2(end:-1:1);
        xp2 = xp2(end:-1:1);
        
        minH = (pdist([xp2(1),yp2(1);xp2(2),yp2(2)],'euclidean'));
        
        stIndx = 1;

         for hg = 1:size(nEle,2)
             str_nEle(hg) = stIndx;
             stIndx = stIndx+nEle(hg);
         end

         cc=hsv(nOfNets);
         ccV1 = cc(1:2:end,:);
         ccV2 = cc(2:2:end,:);
         cc = [ccV1;ccV2];

        if get(handles.tarNet,'Value')
            LoopAte = fact+TarNetCnEle-1;
            TranspMin = 0.4;
        else
            LoopAte = size(matrixR,1);
            TranspMin = 0;
        end
        Diam = max(xp)-min(xp);
        
        if get(handles.CircConnB,'Value') 

            circ = figure('Visible','on');
            set(circ,'Name','Connectome2D',...
                'Position',round([ScreSize(1)*.05 ScreSize(2)*.05 ScreSize(1)*.45 ScreSize(1)*.45]),...
                      'Color',[0 0 0]);
            set(gca,'xticklabels',[])
            set(gca,'yticklabels',[])
            axis off
            hold on

            xlim([-1 limiteG+1]);
            ylim([-1 limiteG+1]);

            daspect([1 1 1])
            cp1 = plot((limiteG/2)+xp,(limiteG/2)+yp,'LineWidth',1.5,'Color',[1 1 1]);

             mTextBox1 = uicontrol('style','text',...
                 'position',round([ScreSize(1)*.001 ScreSize(1)*.4375 ScreSize(1)*.2  ScreSize(1)*.0125]),...
                 'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
                 'ForegroundColor',[1 .949 .867],'FontWeight','bold');
             set(mTextBox1,'String','Circular Connectome - DRAWING........');

             mTextBox2 = uicontrol('style','text',...
                 'position',round([ScreSize(1)*.001 ScreSize(1)*.422 ScreSize(1)*.2  ScreSize(1)*.0125]),...
                 'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
                 'ForegroundColor',[1 .949 .867],'FontWeight','bold');
             set(mTextBox2,'String',['Uncorrected (alpha ' get(handles.edit1,'String') ')']);
            drawnow

            for i = fact:LoopAte
                drawnow
                tester = i<=(str_nEle+nEle)-1;
                networkNumb = find(tester==1,1);
                colorV = cc(networkNumb,:);
                if i<10
                    str = ['0' num2str(i)];
                else
                    str = num2str(i);
                end
                vet = matrixR(i-(fact-1),:);
                vet = find(vet);

                if ~isempty(vet)
                    bt = text(xp2t(i),yp2t(i),[str '-' NetLabels{networkNumb}],'FontSize',9.5,'Color',colorV,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
                    Angle = atan2(yp2t(i)-(limiteG/2),xp2t(i)-(limiteG/2));
                    if abs(Angle)>pi/2
                        Angle = Angle - pi;
                        set(bt,'HorizontalAlignment','right')
                    end
                    set(bt,'Rotation',rad2deg(Angle))
                    for j = 1:size(vet,2)

                        P1 = [xp2(i),yp2(i)];
                        P2 = [xp2(vet(j)),yp2(vet(j))];

                        ARCalpha = pdist([P1(1),P1(2);P2(1),P2(2)],'euclidean');

                        h = minH + (((Diam/ARCalpha)*(minH/Diam)));
                        R = ((ARCalpha^2) + (4*(h^2)))/(8*h);
        %                     r = (sqrt(4*(R^2) - (alpha^2)))/2;
        %                     p = [xp2(i),xp2(vet(j));yp2(i),yp2(vet(j))];
        %                     Teta = (2*acos(r/R));
                        hold on
                        NewC = circ_cent(P1,P2,R);
    %                     [C1,C2] = TanCircCenters(P1,P2,R,0);
    %                     NewC = [C1;C2];

                        Marg = (limiteG-2*RAIO)/2;
                        DisC1 = pdist([NewC(1,1),NewC(1,2);RAIO+Marg,RAIO+Marg],'euclidean');
                        DisC2 = pdist([NewC(2,1),NewC(2,2);RAIO+Marg,RAIO+Marg],'euclidean');

                        if DisC1>DisC2
                            centerIdx = 1;
                        else
                            centerIdx = 2;
                        end

                        AngRad(1) = atan2(P1(2)-NewC(centerIdx,2),P1(1)-NewC(centerIdx,1));
                        AngRad(2) = atan2(P2(2)-NewC(centerIdx,2),P2(1)-NewC(centerIdx,1));

                        if sum(abs(AngRad))>pi && any(AngRad<0)
                            idxCic = find(AngRad<0);
                            AngRad(idxCic) = pi +  pi-abs(AngRad(idxCic));
                        end

                        [AngRad,ksord] = sort(AngRad);

                        if isequal(ksord,[1 2])
                            TranOrd = 'descend';
                        else
                            TranOrd = 'ascend';
                        end

                        ang2 = AngRad(1):0.03:AngRad(2);
                        ang2 = [ang2,AngRad(2)];

                        xunit = R * cos(ang2) + NewC(centerIdx,1);
                        yunit = R * sin(ang2) + NewC(centerIdx,2);

                        TrnaspSca = TranspMin:(1-TranspMin)/size(xunit,2):1;
                        TrnaspSca = sort(TrnaspSca,TranOrd);
                        figure(circ)
                        for hg = 1:size(xunit,2)-1
                            h = plot(xunit(1,hg:hg+1),yunit(1,hg:hg+1),'Color',colorV,'LineWidth',1.5);
                            h.Color(4) = TrnaspSca(hg);
                        end

                        if get(handles.tarNet,'Value')
                            testerj = vet(j)<=(str_nEle+nEle)-1;
                            networkNumbj = find(testerj==1,1);
                            colorVj = cc(networkNumbj,:);
                            if networkNumbj ~= TarNetC
                                if vet(j)<10
                                    strf = ['0' num2str(vet(j))];
                                else
                                    strf = num2str(vet(j));
                                end
                                bt = text(xp2t(vet(j)),yp2t(vet(j)),[strf '-' NetLabels{networkNumbj}],'FontSize',9.5,'Color',...
                                    colorVj,'FontWeight','bold','HorizontalAlignment','left',...
                                    'VerticalAlignment','middle');

                                Angle = atan2(yp2t(vet(j))-(limiteG/2),xp2t(vet(j))-(limiteG/2));

                                if abs(Angle)>pi/2
                                    Angle = Angle - pi;
                                    set(bt,'HorizontalAlignment','right')
                                end
                                set(bt,'Rotation',rad2deg(Angle))
                            end
                        end
                    end
                end
            end
            set(mTextBox1,'String','Circular Connectome');
            drawnow

            set(circ, 'InvertHardCopy', 'off');
        %     saveas(circ,[DirF filesep 'Uncorrected_Circ_Map'],'fig')
            imgRR = getframe(circ);
            imwrite(imgRR.cdata, [DirF filesep 'Uncorrected_Circ_Map' '.png']);
            clear imgRR
            close(circ)
            
            circTS = figure('Visible','on');
            set(circTS,'Name','Connectome2D',...
                'Position',round([ScreSize(1)*.05 ScreSize(2)*.05 ScreSize(1)*.45 ScreSize(1)*.45]),...
                      'Color',[0 0 0],'Units','Normalized');
            set(gca,'xticklabels',[])
            set(gca,'yticklabels',[])
            axis off
            hold on
            
            xlim([-1 limiteG+1]);
            ylim([-1 limiteG+1]);
            
            pos1 = [0.15 0.15 0.7 0.7];
            subplot('Position',pos1,'Color',[0 0 0]);
            daspect([1 1 1])
            cp1 = plot((limiteG/2)+xp,(limiteG/2)+yp,'LineWidth',1.5,'Color',[1 1 1]);
            axis off
            
             mTextBox1 = uicontrol('style','text',...
                 'position',round([ScreSize(1)*.001 ScreSize(1)*.4375 ScreSize(1)*.2  ScreSize(1)*.0125]),...
                 'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
                 'ForegroundColor',[1 .949 .867],'FontWeight','bold');
             set(mTextBox1,'String','Circular Connectome - DRAWING........');
             
            TvaCrit = tinv(1-(str2num(get(handles.edit1,'String'))*.5),(size(g1_VAR2,3)+size(g2_VAR2,3)-2));
            
             mTextBox2 = uicontrol('style','text',...
                 'position',round([ScreSize(1)*.001 ScreSize(1)*.422 ScreSize(1)*.3  ScreSize(1)*.0125]),...
                 'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
                 'ForegroundColor',[1 .949 .867],'FontWeight','bold');
             set(mTextBox2,'String',['Uncorrected (alpha ' get(handles.edit1,'String') ', |Tscore|>' num2str(TvaCrit) ')']);
            drawnow

            ccJet = jet(64);
            ccGre = gray(64);
            for xcv = 1:size(g1_VAR2,1)
                for ycv = 1:size(g1_VAR2,2)
                    Tva(xcv,ycv) = tinv(1-(matrixR(xcv,ycv)*.5),(size(g1_VAR2,3)+size(g2_VAR2,3)-2));
                end
            end
            Tva(Tva==inf) = 0;
            TvaMax = max(max(Tva));
            
            cbm1 = 0;
            cbm2 = 0;
            cbm3 = 0;
            
            for i = fact:LoopAte
                drawnow
                tester = i<=(str_nEle+nEle)-1;
                networkNumb = find(tester==1,1);
                colorV = cc(networkNumb,:);
                
                if i<10
                    str = ['0' num2str(i)];
                else
                    str = num2str(i);
                end
                
                vet = matrixR(i-(fact-1),:);
                vet = find(vet);

                if ~isempty(vet)
                    bt = text(xp2t(i),yp2t(i),[str '-' NetLabels{networkNumb}],'FontSize',9.5,'Color',colorV,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
                    Angle = atan2(yp2t(i)-(limiteG/2),xp2t(i)-(limiteG/2));
                    if abs(Angle)>pi/2
                        Angle = Angle - pi;
                        set(bt,'HorizontalAlignment','right')
                    end
                    set(bt,'Rotation',rad2deg(Angle))
                    for j = 1:size(vet,2)
                        
                        if (g1_mean(i-(fact-1),vet(j))*g2_mean(i-(fact-1),vet(j)))>0
                            if abs(g1_mean(i-(fact-1),vet(j)))>abs(g2_mean(i-(fact-1),vet(j)))
                                colorVlidx = (round((Tva(i-(fact-1),vet(j))/TvaMax)*32))+32;
                                colorVl = ccJet(colorVlidx,:);
                            else
                                colorVlidx = round(((-Tva(i-(fact-1),vet(j))/TvaMax)*32))+33;
                                colorVl = ccJet(colorVlidx,:);
                            end
                        else
                            colorVlidx = round((Tva(i-(fact-1),vet(j))/TvaMax)*64);
                            colorVl = ccGre(colorVlidx,:);
                        end

                        P1 = [xp2(i),yp2(i)];
                        P2 = [xp2(vet(j)),yp2(vet(j))];

                        ARCalpha = pdist([P1(1),P1(2);P2(1),P2(2)],'euclidean');

                        h = minH + (((Diam/ARCalpha)*(minH/Diam)));
                        R = ((ARCalpha^2) + (4*(h^2)))/(8*h);
                        
                        hold on
                        NewC = circ_cent(P1,P2,R);
                        Marg = (limiteG-2*RAIO)/2;
                        DisC1 = pdist([NewC(1,1),NewC(1,2);RAIO+Marg,RAIO+Marg],'euclidean');
                        DisC2 = pdist([NewC(2,1),NewC(2,2);RAIO+Marg,RAIO+Marg],'euclidean');

                        if DisC1>DisC2
                            centerIdx = 1;
                        else
                            centerIdx = 2;
                        end

                        AngRad(1) = atan2(P1(2)-NewC(centerIdx,2),P1(1)-NewC(centerIdx,1));
                        AngRad(2) = atan2(P2(2)-NewC(centerIdx,2),P2(1)-NewC(centerIdx,1));

                        if sum(abs(AngRad))>pi && any(AngRad<0)
                            idxCic = find(AngRad<0);
                            AngRad(idxCic) = pi +  pi-abs(AngRad(idxCic));
                        end

                        [AngRad,ksord] = sort(AngRad);

                        if isequal(ksord,[1 2])
                            TranOrd = 'descend';
                        else
                            TranOrd = 'ascend';
                        end

                        ang2 = AngRad(1):0.03:AngRad(2);
                        ang2 = [ang2,AngRad(2)];

                        xunit = R * cos(ang2) + NewC(centerIdx,1);
                        yunit = R * sin(ang2) + NewC(centerIdx,2);
                        figure(circTS)
                        h1 = plot(xunit,yunit,'Color',colorVl,'LineWidth',3);
                        h1.Color(4) = 0.4;

                        if get(handles.tarNet,'Value')
                            testerj = vet(j)<=(str_nEle+nEle)-1;
                            networkNumbj = find(testerj==1,1);
                            colorVj = cc(networkNumbj,:);
                            if networkNumbj ~= TarNetC
                                if vet(j)<10
                                    strf = ['0' num2str(vet(j))];
                                else
                                    strf = num2str(vet(j));
                                end
                                bt = text(xp2t(vet(j)),yp2t(vet(j)),[strf '-' NetLabels{networkNumbj}],'FontSize',9.5,'Color',...
                                    colorVj,'FontWeight','bold','HorizontalAlignment','left',...
                                    'VerticalAlignment','middle');

                                Angle = atan2(yp2t(vet(j))-(limiteG/2),xp2t(vet(j))-(limiteG/2));

                                if abs(Angle)>pi/2
                                    Angle = Angle - pi;
                                    set(bt,'HorizontalAlignment','right')
                                end
                                set(bt,'Rotation',rad2deg(Angle))
                            end
                        end
                    end
                end
            end
            
            clear A_cell Labx1
            subplot('Position',[0.04 0.01 0.9 0.015])
            imagesc(-TvaMax:(TvaMax/32):TvaMax);
            colormap(gca,'jet');
            set(gca, 'XAxisLocation', 'top','YColor',[0 0 0],'XColor',[1 1 1])
            xticks([1:64/16:64]);
            Labx1 = [(round(10*(-TvaMax:(TvaMax/8):TvaMax)))/10 TvaMax];
            A_cell = split(cellstr(num2str(Labx1)));
            xticklabels(A_cell)
            yticklabels('')
            title(gca,'T-Score G1>G2 (hot colors) G2>G1 (cool colors)','Color','w')
            
            set(mTextBox1,'String','Circular Connectome');
            drawnow

            set(circTS, 'InvertHardCopy', 'off');
            imgRR = getframe(circTS);
            imwrite(imgRR.cdata, [DirF filesep 'Uncorrected_T-Score_Circ_Map' '.png']);
            clear imgRR
            close(circTS)
        end

        if get(handles.OrgCircConnB,'Value')

            circ2 = figure('Visible','on');
            set(circ2,'Name','Connectome2D',...
                'Position',round([ScreSize(1)*.05 ScreSize(2)*.05 ScreSize(1)*.45 ScreSize(1)*.45]),...
                      'Color',[0 0 0]);
            set(gca,'xticklabels',[])
            set(gca,'yticklabels',[])
            axis off
            hold on

            xlim([-1 limiteG+1]);
            ylim([-1 limiteG+1]);

            daspect([1 1 1])
            cp1 = plot((limiteG/2)+xp,(limiteG/2)+yp,'LineWidth',1.5,'Color',[1 1 1]);

             mTextBox1 = uicontrol('style','text',...
                 'position',round([ScreSize(1)*.001 ScreSize(1)*.4375 ScreSize(1)*.3  ScreSize(1)*.0125]),...
                 'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
                 'ForegroundColor',[1 .949 .867],'FontWeight','bold');
             set(mTextBox1,'String','Organic Circular Connectome - DRAWING........');

             mTextBox2 = uicontrol('style','text',...
                 'position',round([ScreSize(1)*.001 ScreSize(1)*.422 ScreSize(1)*.2  ScreSize(1)*.0125]),...
                 'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
                 'ForegroundColor',[1 .949 .867],'FontWeight','bold');
             set(mTextBox2,'String',['Uncorrected (alpha ' get(handles.edit1,'String') ')']);
            drawnow

            daspect([1 1 1])
            RAIO2 = RAIO*0.15;        
            xp3=RAIO2*cos(ang);
            yp3=RAIO2*sin(ang);

            xp3 = xp3(1:AnoRa:end);
            yp3 = yp3(1:AnoRa:end);
            xp3 = flip(xp3+(limiteG/2));
            yp3 = flip(yp3+(limiteG/2));


            RAIO3 = RAIO*0.85;        
            xp4=RAIO3*cos(ang);
            yp4=RAIO3*sin(ang);


            idxNR = ceil((str_nEle+(nEle./2))-1);
            xp4 = flip(xp4(1:AnoRa:end));
            yp4 = flip(yp4(1:AnoRa:end));
            xp4 = xp4(idxNR);
            yp4 = yp4(idxNR);

            xp4 = xp4+(limiteG/2);
            yp4 = yp4+(limiteG/2);


            RAIO4 = RAIO*0.95;        
            xp5=RAIO4*cos(ang);
            yp5=RAIO4*sin(ang);

            xp5 = xp5((AnoRa/2):AnoRa:end);
            yp5 = yp5((AnoRa/2):AnoRa:end);
            xp5 = flip(xp5+(limiteG/2));
            yp5 = flip(yp5+(limiteG/2));

            for i = fact:LoopAte
                drawnow
                tester = i<=(str_nEle+nEle)-1;
                networkNumb = find(tester==1,1);
                colorV = cc(networkNumb,:);
                if i<10
                    str = ['0' num2str(i)];
                else
                    str = num2str(i);
                end
                vet = matrixR(i-(fact-1),:);
                vet = find(vet);

                if ~isempty(vet)
                    bt = text(xp2t(i),yp2t(i),[str '-' NetLabels{networkNumb}],'FontSize',9.5,'Color',colorV,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
                    Angle = atan2(yp2t(i)-(limiteG/2),xp2t(i)-(limiteG/2));
                    if abs(Angle)>pi/2
                        Angle = Angle - pi;
                        set(bt,'HorizontalAlignment','right')
                    end
                    set(bt,'Rotation',rad2deg(Angle))
                    for j = 1:size(vet,2)

                        P1 = [xp2(i),yp2(i)];
                        P2 = [xp2(vet(j)),yp2(vet(j))];

                        testerj = vet(j)<=(str_nEle+nEle)-1;
                        networkNumbj = find(testerj==1,1);
                        colorVj = cc(networkNumbj,:);
                        xxF = [];
                        yyF = [];
                        if networkNumb~= networkNumbj
                            opt1 =  round(i-(i-vet(j))/2);
                            opt2 = round(opt1+(LoopAte/2));
                            if opt2<=LoopAte
                                dopt1 = pdist([xp3(opt1),yp3(opt1);P1(1),P1(2)]);
                                dopt2 = pdist([xp3(opt2),yp3(opt2);P1(1),P1(2)]);
                                if dopt1>dopt2
                                    dopt = opt2;
                                else
                                    dopt = opt1;
                                end
                            else
                                dopt = opt1;
                            end

                            [xx,yy]=fillline([P1(1),P1(2)],[xp4(networkNumb),yp4(networkNumb)],15);
                            [xx2,yy2]=fillline([xp4(networkNumb),yp4(networkNumb)],[xp3(dopt),yp3(dopt)],20);
                            [xx3,yy3]=fillline([xp3(dopt),yp3(dopt)],[xp4(networkNumbj),yp4(networkNumbj)],20);
                            [xx4,yy4]=fillline([xp4(networkNumbj),yp4(networkNumbj)],[P2(1),P2(2)],15);
                            
                            xxF = [xx,xx2,xx3,xx4];
                            yyF = [yy,yy2,yy3,yy4];
                            xxF = [xxF(1);smoothLNI(xxF(2:end-1),30);xxF(end)]';
                            yyF = [yyF(1);smoothLNI(yyF(2:end-1),30);yyF(end)]';
                            TranOrd = 'descend';
                            
                            TrnaspSca = [zeros(1,size(xxF,2)*0.3)+TranspMin,TranspMin:(1-TranspMin)/(size(xxF,2)*0.7):1];
                            TrnaspSca = sort(TrnaspSca,TranOrd);
                            figure(circ2)
                            if networkNumb~=networkNumbj
                                if get(handles.tarNet,'Value')
                                    for hg = 1:(size(xxF,2)-1)
                                        h = plot(xxF(1,hg:hg+1),yyF(1,hg:hg+1),'Color',colorV,'LineWidth',1.5);
                                        h.Color(4) = TrnaspSca(hg);
                                    end
                                else
                                    for hg = 1:(size(xxF,2)-1)*0.7
                                        h = plot(xxF(1,hg:hg+1),yyF(1,hg:hg+1),'Color',colorV,'LineWidth',1.5);
                                        h.Color(4) = TrnaspSca(hg);
                                    end
                                end
                            else
                                if i<vet(j)
                                    h = plot(xxF,yyF,'Color',colorV,'LineWidth',1);
                                end
                            end

                        else

                            ARCalpha = pdist([P1(1),P1(2);P2(1),P2(2)],'euclidean');

                            h = minH + (((Diam/ARCalpha)*(minH/Diam)));
                            R = ((ARCalpha^2) + (4*(h^2)))/(8*h);
                            NewC = circ_cent(P1,P2,R);
                            Marg = (limiteG-2*RAIO)/2;
                            DisC1 = pdist([NewC(1,1),NewC(1,2);RAIO+Marg,RAIO+Marg],'euclidean');
                            DisC2 = pdist([NewC(2,1),NewC(2,2);RAIO+Marg,RAIO+Marg],'euclidean');

                            if DisC1>DisC2
                                centerIdx = 1;
                            else
                                centerIdx = 2;
                            end

                            AngRad(1) = atan2(P1(2)-NewC(centerIdx,2),P1(1)-NewC(centerIdx,1));
                            AngRad(2) = atan2(P2(2)-NewC(centerIdx,2),P2(1)-NewC(centerIdx,1));

                            if sum(abs(AngRad))>pi && any(AngRad<0)
                                idxCic = find(AngRad<0);
                                AngRad(idxCic) = pi +  pi-abs(AngRad(idxCic));
                            end

                            ang2 = AngRad(1):0.03:AngRad(2);
                            ang2 = [ang2,AngRad(2)];

                            xunit = R * cos(ang2) + NewC(centerIdx,1);
                            yunit = R * sin(ang2) + NewC(centerIdx,2);

                            h = plot(xunit,yunit,'Color',colorV.*0.5,'LineWidth',2);
                        end

                        if vet(j)<10
                            strf = ['0' num2str(vet(j))];
                        else
                            strf = num2str(vet(j));
                        end
                        bt = text(xp2t(vet(j)),yp2t(vet(j)),[strf '-' NetLabels{networkNumbj}],'FontSize',9.5,'Color',...
                            colorVj,'FontWeight','bold','HorizontalAlignment','left',...
                            'VerticalAlignment','middle');

                        Angle = atan2(yp2t(vet(j))-(limiteG/2),xp2t(vet(j))-(limiteG/2));

                        if abs(Angle)>pi/2
                            Angle = Angle - pi;
                            set(bt,'HorizontalAlignment','right')
                        end
                        set(bt,'Rotation',rad2deg(Angle))
                    end
                end
            end


                set(mTextBox1,'String','Organic Circular Connectome');
                drawnow

                set(circ2, 'InvertHardCopy', 'off');
            %     saveas(circ,[DirF filesep 'Uncorrected_Circ_Map'],'fig')
                imgRR = getframe(circ2);
                imwrite(imgRR.cdata, [DirF filesep 'Uncorrected_OrganicCirc_Map' '.png']);
                clear imgRR
                close(circ2)
        end
        
        if get(handles.Graphos3D,'Value')  %%%%% 3D Graphos
            GraphoFig = figure;

            set(GraphoFig,'Name','3DConnectome',...
            'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.6 ScreSize(1)*.4]),...
                        'Color',[1 .949 .867]);

            title(['Uncorrected (alpha ' get(handles.edit1,'String') ') 3D Connectome'],'FontSize',14);
            xlabel('X axis'  ,'FontSize',14);
            ylabel('Y axis','FontSize',14);
            zlabel('Z axis','FontSize',14);
            set(gca, 'color', [1 1 1])
            try set(gca, 'gridcolor', [0.8 0.8 0.8]); end

            try set(gca, 'gridalpha', 1); end
            axis off
            xlim([0 91]);
            ylim([0 109]);
            zlim([0 91]);
            view([-1 -1 1])
            hold on

            cc=hsv(nOfNets);
            ccV1 = cc(1:2:end,:);
            ccV2 = cc(2:2:end,:);
            cc = [ccV1;ccV2];
            [rws,res,rdw] = sphere(10);

            if  get(handles.tarNet,'Value')
                for huy = fact:size(g1_VAR2,1)+(fact-1)
                    CoordQ = CoordMS(huy,:);
                    CoordList = CoordMS(1:end,:);
                    for hgf = 1:size(CoordList,1)
                        isF = [CoordQ(1,1);CoordList(hgf,1)];
                        jsF = [CoordQ(1,2);CoordList(hgf,2)];
                        ksF = [CoordQ(1,3);CoordList(hgf,3)];
                       if unCorPmap(huy-(fact-1),hgf)>-1.*(str2num(get(handles.edit1,'String')))
                            [xre,yre,zre] = tubeplot([isF(1) isF(2);jsF(1) jsF(2);ksF(1) ksF(2)],0.4,8);
                            if (g1_mean(huy-(fact-1),hgf)*g2_mean(huy-(fact-1),hgf))>0
                                if abs(g1_mean(huy-(fact-1),hgf))>abs(g2_mean(huy-(fact-1),hgf))
                                    ySurface1 = surf(xre,yre,zre);
                                    set(ySurface1,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                        'gouraud','FaceColor',[0.1 0.5 1],'FaceLighting','gouraud','LineSmoothing','on');
                                    ex1_1 = 1;
                                else
                                    ySurface2 = surf(xre,yre,zre);
                                    set(ySurface2,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                        'gouraud','FaceColor',[1 0.1 0.1],'FaceLighting','gouraud','LineSmoothing','on');
                                    ex1_2 = 1;
                                end
                            else
                                ySurface3 = surf(xre,yre,zre);
                                set(ySurface3,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                        'gouraud','FaceColor',[1 0.8 0],'FaceLighting','gouraud','LineSmoothing','on');
                                ex1_3 = 1;
                            end
                       end   
                    end
                end
            else
                for huy = 1:size(CoordMS,1)
                    CoordQ = CoordMS(huy,:);
                    CoordList = CoordMS(huy+1:end,:);
                    for hgf = 1:size(CoordList,1)
                        isF = [CoordQ(1,1);CoordList(hgf,1)];
                        jsF = [CoordQ(1,2);CoordList(hgf,2)];
                        ksF = [CoordQ(1,3);CoordList(hgf,3)];
                       if unCorPmap(hgf+huy,huy)>-1.*(str2num(get(handles.edit1,'String')))
                            [xre,yre,zre] = tubeplot([isF(1) isF(2);jsF(1) jsF(2);ksF(1) ksF(2)],0.4,8);
                            if (g1_mean(hgf+huy,huy)*g2_mean(hgf+huy,huy))>0
                                if abs(g1_mean(hgf+huy,huy))>abs(g2_mean(hgf+huy,huy))
                                    ySurface1 = surf(xre,yre,zre);
                                    set(ySurface1,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                        'gouraud','FaceColor',[0.1 0.5 1],'FaceLighting','gouraud','LineSmoothing','on');
                                    ex1_1 = 1;
                                else
                                    ySurface2 = surf(xre,yre,zre);
                                    set(ySurface2,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                        'gouraud','FaceColor',[1 0.1 0.1],'FaceLighting','gouraud','LineSmoothing','on');
                                    ex1_2 = 1;
                                end
                            else
                                ySurface3 = surf(xre,yre,zre);
                                set(ySurface3,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                        'gouraud','FaceColor',[1 0.8 0],'FaceLighting','gouraud','LineSmoothing','on');
                                ex1_3 = 1;
                            end
                       end   
                    end
               end
           end
            %%%%%%%%%%%%%%%%%%%%%%%  Include only spheres with alterations
            
            teste = unCorPmap;
            teste(teste==-1)=0;
            teste(teste~=0)=1;
            vetTTT = sum(teste,1);
            VetFF = vetTTT;
            VetFF(VetFF>1) = 1;
            idxROIs = find(VetFF);

            vetTT2 = sum(teste,2);
            VetFF2 = vetTT2;
            VetFF2(VetFF2>1) = 1;
            idxROIs2 = find(VetFF2);
            idxROIs2 = idxROIs2+(fact-1);
            vetTT2(vetTT2==0) = [];
            VetFF(idxROIs2) = 1;
            vetTTT(idxROIs2) = vetTT2;

            idxROIs = [idxROIs,idxROIs2'];
            idxROIs = unique(idxROIs);

            for i = 1:numel(idxROIs)
                tester = idxROIs(i)<=(str_nEle+nEle)-1;
                networkNumb = find(tester==1,1);
                colorV = cc(networkNumb,:);
                CoordQ = CoordMS(idxROIs(i),:);
                CoordV4 = (vetTTT(idxROIs(i))/max(vetTTT))+1;
                rNumb = find(VetFF~=0,i);
                rNumb = rNumb(end);

                if isequal(networkNumb,1)
                    if isequal(REtype,2)
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(3*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(rNumb) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    else
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(3*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(i) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    end
                else
                    if isequal(REtype,2)
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(5*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(rNumb-sum(nEle(1:networkNumb-1))) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    else
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(5*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(idxROIs(i)) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    end
                end
                hSurface = surf(CoordV4.*rws+CoordQ(1,1),CoordV4.*res+CoordQ(1,2),CoordV4.*rdw+CoordQ(1,3));
                set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             if exist('ex1_1') || exist('ex1_2')
                if exist('ex1_1') && exist('ex1_2')
                    if exist('ex1_3')
                        legend([ySurface1,ySurface2,ySurface3],[G1nameT ' > ' G2nameT],[G2nameT ' > ' G1nameT],[G2nameT ' ' char(8596) ' -(' G1nameT ')'])
                    else
                        legend([ySurface1,ySurface2],[G1nameT ' > ' G2nameT],[G2nameT ' > ' G1nameT])
                    end
                else if exist('ex1_1')
                        if exist('ex1_3')
                            legend([ySurface1,ySurface3],[G1nameT ' > ' G2nameT],[G2nameT ' ' char(8596) ' -(' G1nameT ')'])
                        else
                            legend(ySurface1,[G1nameT ' > ' G2nameT])
                        end
                    else
                         if exist('ex1_3')
                             legend([ySurface2,ySurface3],[G2nameT ' > ' G1nameT],[G2nameT ' ' char(8596) ' -(' G1nameT ')'])
                         else
                            legend(ySurface2,[G2nameT ' > ' G1nameT])
                         end
                     end
                end
            else if exist('ex1_3')
                     legend(ySurface3,[G2nameT ' ' char(8596) ' -(' G1nameT ')']) % Simbolos apenas copiando e colando...
                 end
            end

            mTextBox = uicontrol('style','text','position',round([ScreSize(1)*.01 ScreSize(2)*.01 ScreSize(1)*.15 ScreSize(1)*.01]));
            set(mTextBox,'String','r = Roi number and n = network number');

            DirUF2C = which('uf2c');
            try
                aRT = nifti([DirUF2C(1:end-6) 'Utilities' filesep 'WFUteste2.nii']);
            catch
                aRT = nifti([DirUF2C(1:end-13) 'Utilities' filesep 'WFUteste2.nii']);
            end

            matI = aRT.dat(:,:,:);
            matI = permute(matI,[1,3,2]);
            matI = permute(matI,[2,1,3]);
            matI = permute(matI,[3,2,1]);
            matI = matI.*100;
            matI = int8(matI);
            Ds = smooth3(matI);
            Ds = int8(Ds);
            Xiso = isosurface(Ds,10);
            hiso = patch(Xiso,'FaceColor',[0.6,0.6,0.6],'EdgeColor','none','FaceAlpha',0.3);%,'FaceVertexAlphaData',0.8,'AlphaDataMapping','Direct','CDataMapping','Direct');
            daspect([1,0.83,0.83])
            view(113,14)
            lighting gouraud
            light('Style','infinite'); %'Position',[1 0 1]);
            
            try
                legTmp = legend;
                IndexLegC = strfind(legTmp.String,'data1');
                IndexLeg = find(not(cellfun('isempty',IndexLegC)));
                legTmp.String(IndexLeg) = []; 
            end
            
            saveas(GraphoFig,[DirF filesep 'Uncorrected_Graph_3D'],'fig')
            drawnow
            imgRR = getframe(GraphoFig);
            imwrite(imgRR.cdata, [DirF filesep 'Uncorrected_Graph_3D' '.png']);
            clear imgRR
            
            view(90,0)
            axis off
            drawnow
            set(gcf, 'InvertHardCopy', 'off');
            imgRR = getframe(GraphoFig);
            imwrite(imgRR.cdata, [DirF filesep 'Uncorrected_Graph_3D_Sag' '.png']);
            clear imgRR
            
            mTextBoxOri1 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.05 ScreSize(2)*.101 ScreSize(1)*.03 ScreSize(1)*.05]),...
             'FontSize',60,'BackgroundColor',[1 .949 .867],...
             'ForegroundColor',[0 0 0]);
            set(mTextBoxOri1,'String','L');
            
            mTextBoxOri2 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.4 ScreSize(2)*.101 ScreSize(1)*.03 ScreSize(1)*.05]),...
             'FontSize',60,'BackgroundColor',[1 .949 .867],...
             'ForegroundColor',[0 0 0]);
            set(mTextBoxOri2,'String','R');
            
            view(0,90)
            drawnow
            set(gcf, 'InvertHardCopy', 'off');
            imgRR = getframe(GraphoFig);
            imwrite(imgRR.cdata, [DirF filesep 'Uncorrected_Graph_3D_Ax' '.png']);
            clear imgRR
            
            view(0,0)
            axis off
            drawnow
            set(gcf, 'InvertHardCopy', 'off');
            imgRR = getframe(GraphoFig);
            imwrite(imgRR.cdata, [DirF filesep 'Uncorrected_Graph_3D_Cor' '.png']);
            clear imgRR

            close(GraphoFig)

            CoLeg = figure('Visible','off');
            set(CoLeg,'Name','Color Legend','Position', [400 400 250 400],...
                            'Color',[1 .949 .867]);
            title('Color Legend','FontSize',14);
            xlim([0 0.6]);
            ylim([0 nOfNets]);

            for uiC = 1:nOfNets
                rectangle('Position',[0,uiC-1,0.6,1],'LineWidth',2,'FaceColor',cc(uiC,:))
                string = sprintf('Network %d color',uiC);
                text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
            end

            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])

    %             saveas(CoLeg,[DirF filesep 'Color_Legend'],'png')
            drawnow
            imgRR = getframe(CoLeg);
            imwrite(imgRR.cdata, [DirF filesep 'Color_Legend' '.png']);

            saveas(CoLeg,[DirF filesep 'Color_Legend'],'fig')
            close(CoLeg)

        end
        
        if get(handles.Grapho2D,'Value')  %%%%% 2D Graphos
            CoordMS2 = 3.*CoordMS;
            GraphoFig2 = figure;

            set(GraphoFig2,'Name','Connectome2D',...
                'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.45 ScreSize(1)*.45]),...
                            'Color',[1 .949 .867]);

            title(['Uncorrected (alpha ' get(handles.edit1,'String') ') 2D Connectome'],'FontSize',14);
            xlabel('X axis'  ,'FontSize',14);
            ylabel('Y axis','FontSize',14);
            zlabel('Z axis','FontSize',14);

            set(gca, 'color', [1 1 1])
            try set(gca, 'gridcolor', [0 0 0]); end
            try set(gca, 'gridalpha', 1); end
            hold on
            xlim([0 3*91]);
            ylim([0 3*109]);
            zlim([0 3*91]);
            view([0 0 1])

             stIndx = 1;
             for hg = 1:size(nEle,2)
                 str_nEle(hg) = stIndx;
                 stIndx = stIndx+nEle(hg);
             end

            cc=hsv(nOfNets);
            ccV1 = cc(1:2:end,:);
            ccV2 = cc(2:2:end,:);
            cc = [ccV1;ccV2];
            [rws,res,rdw] = sphere(15);

             if get(handles.tarNet,'Value')
                for huy = fact:size(g1_VAR2,1)+(fact-1)
                    CoordQ2 = CoordMS2(huy,:);
                    CoordList = CoordMS2;
                    for hgf = 1:size(CoordList,1)
                        isF = [CoordQ2(1,1);CoordList(hgf,1)];
                        jsF = [CoordQ2(1,2);CoordList(hgf,2)];
                        ksF = [CoordQ2(1,3);CoordList(hgf,3)];
                       if unCorPmap(huy-(fact-1),hgf)>-1.*(str2num(get(handles.edit1,'String')))
                            if (g1_mean(huy-(fact-1),hgf)*g2_mean(huy-(fact-1),hgf))>0
                                if abs(g1_mean(huy-(fact-1),hgf))>abs(g2_mean(huy-(fact-1),hgf))
                                   pLine1 =  line(isF,jsF,ksF,'Color',[0.2 0.5 1],'LineStyle','-','LineWidth',0.8);
                                   ex2_1 = 1;
                                else
                                   pLine2 = line(isF,jsF,ksF,'Color',[1 0.3 0.3],'LineStyle','-','LineWidth',0.8);
                                   ex2_2 = 1;
                                end
                            else
                                pLine3 =  line(isF,jsF,ksF,'Color',[0 0 0],'LineStyle','--','LineWidth',0.7);
                                ex2_3 = 1;
                            end
                       end   
                    end
                end
            else
                for huy = 1:size(CoordMS2,1)
                    CoordQ2 = CoordMS2(huy,:);
                    CoordList = CoordMS2(huy+1:end,:);
                    for hgf = 1:size(CoordList,1)
                        isF = [CoordQ2(1,1);CoordList(hgf,1)];
                        jsF = [CoordQ2(1,2);CoordList(hgf,2)];
                        ksF = [CoordQ2(1,3);CoordList(hgf,3)];
                        if unCorPmap(hgf+huy,huy)>-1.*(str2num(get(handles.edit1,'String')))
                            if (g1_mean(hgf+huy,huy)*g2_mean(hgf+huy,huy))>0
                                if abs(g1_mean(hgf+huy,huy))>abs(g2_mean(hgf+huy,huy))
                                   pLine1 =  line(isF,jsF,ksF,'Color',[0.2 0.5 1],'LineStyle','-','LineWidth',0.8);
                                   ex2_1 = 1;
                                else
                                   pLine2 = line(isF,jsF,ksF,'Color',[1 0.3 0.3],'LineStyle','-','LineWidth',0.8);
                                   ex2_2 = 1;
                                end
                            else
                                pLine3 =  line(isF,jsF,ksF,'Color',[0 0 0],'LineStyle','--','LineWidth',0.7);
                                ex2_3 = 1;
                            end
                        end   
                    end
                end
             end
            %%%%%%%%%%%%%%%%%%%%%%%  include only  spheres with alterations
            
            teste = unCorPmap;
            teste(teste==-1)=0;
            teste(teste~=0)=1;
            vetTTT = sum(teste,1);
            VetFF = vetTTT;
            VetFF(VetFF>1) = 1;
            idxROIs = find(VetFF);

            vetTT2 = sum(teste,2);
            VetFF2 = vetTT2;
            VetFF2(VetFF2>1) = 1;
            idxROIs2 = find(VetFF2);
            idxROIs2 = idxROIs2+(fact-1);
            vetTT2(vetTT2==0) = [];
            VetFF(idxROIs2) = 1;
            vetTTT(idxROIs2) = vetTT2;

            idxROIs = [idxROIs,idxROIs2'];
            idxROIs = unique(idxROIs);


            for i = 1:sum(VetFF~=0)
                tester = idxROIs(i)<=(str_nEle+nEle)-1;
                networkNumb = find(tester==1,1);
                colorV = cc(networkNumb,:);
                CoordQ = CoordMS2(idxROIs(i),:);
                CoordV4 = (vetTTT(idxROIs(i))/max(vetTTT))+1;
                rNumb = find(VetFF~=0,i);
                rNumb = rNumb(end);
                
                if isequal(networkNumb,1)
                    if isequal(REtype,2)
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(3*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(rNumb) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    else
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(3*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(i) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    end
                else
                    if isequal(REtype,2)
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(5*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(rNumb-sum(nEle(1:networkNumb-1))) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    else
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(5*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(idxROIs(i)) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    end
                end
                hSurface = surf(CoordV4*3.*rws+CoordQ(1,1),CoordV4*3.*res+CoordQ(1,2),CoordV4*3.*rdw+CoordQ(1,3));
                set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if exist('ex2_1') || exist('ex2_2')
                if exist('ex2_1') && exist('ex2_2')
                    if exist('ex2_3')
                        legend([pLine1,pLine2,pLine3],[G1nameT ' > ' G2nameT],[G2nameT ' > ' G1nameT],[G2nameT ' ' char(8596) ' -(' G1nameT ')'])
                    else
                        legend([pLine1,pLine2],[G1nameT ' > ' G2nameT],[G2nameT ' > ' G1nameT])
                    end
                else if exist('ex2_1')
                        if exist('ex2_3')
                            legend([pLine1,pLine3],[G1nameT ' > ' G2nameT],[G2nameT ' ' char(8596) ' -(' G1nameT ')'])
                        else
                            legend(pLine1,[G1nameT ' > ' G2nameT])
                        end
                     else
                        if exist('ex2_3')
                            legend([pLine2,pLine3],[G2nameT ' > ' G1nameT],[G2nameT ' ' char(8596) ' -(' G1nameT ')'])
                        else
                            legend(pLine2,[G2nameT  ' > ' G1nameT])
                        end
                     end
                end
            else if exist('ex2_3')
                     legend(pLine3,[G2nameT ' ' char(8596) ' -(' G1nameT ')'])
                 end
            end

            mTextBox = uicontrol('style','text','position',round([ScreSize(1)*.01 ScreSize(2)*.01 ScreSize(1)*.15 ScreSize(1)*.01]));
            set(mTextBox,'String','r = Roi number and n = network number');

            DirUF2C = which('uf2c');

            try
                load([DirUF2C(1:end-6) 'Utilities' filesep 'Brain_Contour.mat'])
            catch
                load([DirUF2C(1:end-13) 'Utilities' filesep 'Brain_Contour.mat'])
            end

            daspect([1,0.95,1])
            ater2 = slice(mat3d,[],[],1);
            set(ater2, 'FaceAlpha',0.7);
            colormap gray
            set(ater2, 'LineStyle','none');
            matConto = mat3d(:,:,1);%round(size(final3D,3)/2));
            matConto = squeeze(matConto);
            matConto(matConto<1)=0.0001;
            set(ater2(1),'alphadata',matConto,'facealpha','flat')
            grid off
            axis off
            
            try
                legTmp = legend;
                IndexLegC = strfind(legTmp.String,'data1');
                IndexLeg = find(not(cellfun('isempty',IndexLegC)));
                legTmp.String(IndexLeg) = []; 
            end

            light('Style','infinite');
                        
            mTextBoxOri12 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.05 ScreSize(2)*.101 ScreSize(1)*.03 ScreSize(1)*.05]),...
             'FontSize',60,'BackgroundColor',[1 .949 .867],...
             'ForegroundColor',[0 0 0]);
            set(mTextBoxOri12,'String','L');
            
            mTextBoxOri22 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.4 ScreSize(2)*.101 ScreSize(1)*.03 ScreSize(1)*.05]),...
             'FontSize',60,'BackgroundColor',[1 .949 .867],...
             'ForegroundColor',[0 0 0]);
            set(mTextBoxOri22,'String','R');

            
            set(gcf, 'InvertHardCopy', 'off');
            saveas(GraphoFig2,[DirF filesep 'Uncorrected_Graph_2D'],'fig')
            drawnow
            imgRR = getframe(GraphoFig2);
            imwrite(imgRR.cdata, [DirF filesep 'Uncorrected_Graph_2D' '.png']);
            clear imgRR
            
            close(GraphoFig2)

            CoLeg = figure('Visible','off');
            set(CoLeg,'Name','Color Legend','Position', [400 400 250 400],...
                            'Color',[1 .949 .867]);
            title('Color Legend','FontSize',14);
            xlim([0 0.6]);
            ylim([0 nOfNets]);

            for uiC = 1:nOfNets
                rectangle('Position',[0,uiC-1,0.6,1],'LineWidth',2,'FaceColor',cc(uiC,:))
                string = sprintf('Network %d color',uiC);
                text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
            end
            set(gcf, 'InvertHardCopy', 'off');
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])

            drawnow
            imgRR = getframe(CoLeg);
            imwrite(imgRR.cdata, [DirF filesep 'Color_Legend' '.png']);
            clear imgRR

            saveas(CoLeg,[DirF filesep 'Color_Legend'],'fig')
            close(CoLeg)

        end
    else
        warndlg('No significant results with the ''uncorrected'' defined threshold ','Ops!')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%  FDR
if get(handles.chec2,'Value')  %%%%% FDR corrected Matrix
    FDRcor = Pmap;
    
    if get(handles.AnaLe,'Value')
        FigNT = 'Analysis Level - FDR';
        FigNTpx = 'AL';
        
        if get(handles.tarNet,'Value')
            FDRcor(:,TarNetCInit:(TarNetCInit+TarNetCnEle-1)) = tril(FDRcor(:,TarNetCInit:(TarNetCInit+TarNetCnEle-1)),-1);
            reshFDRcor = reshape(FDRcor,1,prod(size(FDRcor)));
            ZeroFind = find(reshFDRcor<=0);
            tPmap = FDRcor(FDRcor>0);
            
            [~,crit_p,adj_p] = fdr_bh(tPmap,str2num(get(handles.edit2,'String')),'pdep','yes');
            if size(adj_p,1)<size(adj_p,2)
                adj_p = adj_p';
            end
            
            for i = 1:numel(ZeroFind)
                if isequal(ZeroFind(i),1)
                    adj_p = [0;adj_p(ZeroFind(i):end)];
                else
                    adj_p = [adj_p(1:ZeroFind(i)-1);0;adj_p(ZeroFind(i):end)];
                end
            end
            resF = reshape(adj_p,size(FDRcor));
        else
            tPmap = tril(FDRcor,-1);
            tPmap = tPmap(tPmap>0);
            [~,crit_p,adj_p] = fdr_bh(tPmap,str2num(get(handles.edit2,'String')),'pdep','yes');
            res = zeros(size(g1_VAR2,1),size(g1_VAR2,2));
            stt = 1;

            for i = 1:size(g1_VAR2,2)
                res(1:i,i) = zeros(i,1);
                res(i+1:end,i) = adj_p(stt:stt+size(g1_VAR2,1)-i-1);
                stt = stt+(size(g1_VAR2,1)-i);
            end
            resF = res+transpose(res);
        end
    else
        FigNT = 'Node Level - FDR';
        FigNTpx = 'NL';

        res = zeros(size(g1_VAR2,1),size(g1_VAR2,1));
        res(1,1) = 1;
        for rl = 1:size(g1_VAR2,1)
            tPmapV = FDRcor(:,rl);
            tPmapV(rl) = [];
            [~,crit_p,adj_p] = fdr_bh(tPmapV,str2num(get(handles.edit2,'String')),'pdep','no');
            res(1:rl-1,rl) = adj_p(1:rl-1);
            res(rl,rl) = 1;
            res(rl+1:end,rl) = adj_p(rl:end);
        end
        
        res1 = tril(res,-1);
        res1(res1>(str2num(get(handles.edit2,'String')))) = 1;
        res2 = tril(res',-1);
        res2(res2>(str2num(get(handles.edit2,'String')))) = 1;
        for ggg = 1:size(g1_VAR,1)
            for nnn = 1:size(g1_VAR,1)
                if isequal(res1(ggg,nnn),res2(ggg,nnn))
                    resF(ggg,nnn) = res1(ggg,nnn);
                else
                    if isequal(res1(ggg,nnn),1)
                        resF(ggg,nnn) = res2(ggg,nnn);
                    else if isequal(res2(ggg,nnn),1)
                            resF(ggg,nnn) = res1(ggg,nnn);
                        else
                            resF(ggg,nnn) = (res1(ggg,nnn)+res2(ggg,nnn))/2;
                        end
                    end
                end
            end
        end
        resF = resF+resF';
    end
    
    FDRcor = resF;

%     FDRcor(FDRcor>(str2num(get(handles.edit2,'String')))) = 1;
    FDRcor(FDRcor>0.05) = 1; %0.05 because we are using the adjust p-values, regarding the input initial threshold.
    FDRcor(FDRcor==0) = 1;
    FDRcor = -1.*FDRcor;
    matrix_bin2 = FDRcor>-1;
    
    clear Ridx tester1 RnUmb1 networkNumb1 ixR r strDir RnUmb2 networkNumb2 tester2 ixR
    
    fidRes2 = fopen([DirF,filesep,FigNTpx,'_FDR_Results_Description.txt'],'w+');
    
    for r = 1:size(g1_VAR2,1)
        Ridx = find(matrix_bin2(r,:));
        if sum(Ridx)>0
            if get(handles.tarNet,'Value')
                fprintf(fidRes2,'ROI %d: r%d%s: \r\n',r+fact-1,r,NetLabels{TarNetC});
            else
                tester1 = r<=(str_nEle+nEle)-1;
                networkNumb1 = find(tester1==1,1);
                try
                    RnUmb1 = r - sum(nEle(1:networkNumb1-1));
                catch
                    RnUmb1 = r;
                end
                if isequal(REtype,2)
                    fprintf(fidRes2,'ROI %d: r%d%s: \r\n',r,RnUmb1,NetLabels{networkNumb1});
                else
                    fprintf(fidRes2,'ROI %d: r%d%s: \r\n',r,r,NetLabels{networkNumb1});
                end
            end
            for ixR = 1:numel(Ridx)
                if MultMean(r,Ridx(ixR))>0
                     if abs(g1_mean(r,Ridx(ixR)))>abs(g2_mean(r,Ridx(ixR)))
                         strDir = [G1nameT ' > ' G2nameT];
                     else
                         strDir = [G2nameT ' > ' G1nameT];
                     end
                else
                    strDir = [G2nameT ' ' '<->' ' -(' G1nameT ')'];
                end
                tester2 = Ridx(ixR)<=(str_nEle+nEle)-1;
                networkNumb2 = find(tester2==1,1);
                try
                    RnUmb2 = Ridx(ixR) - sum(nEle(1:networkNumb2-1));
                catch
                    RnUmb2 = Ridx(ixR);
                end
                if isequal(REtype,2)
                    fprintf(fidRes2,'\t ROI %d: r%d%s \t %s\r\n',Ridx(ixR),RnUmb2,NetLabels{networkNumb2},strDir);
                else
                    fprintf(fidRes2,'\t ROI %d: r%d%s \t %s\r\n',Ridx(ixR),Ridx(ixR),NetLabels{networkNumb2},strDir);
                end
            end
        end
    end
    fclose(fidRes2);
    
    if ~isequal(numel(unique(FDRcor)),1)
        if get(handles.CorMats,'Value')
            fig2 = figure;
            imagesc(FDRcor)
            
            set(fig2,'Name',[FigNT ' FDR-Corrected-Correlation-Matrix'],'Position',...
                round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.45 ScreSize(1)*.4]),...
                    'Color',[1 .949 .867]);
                
            title([FigNT ' (alpha ' get(handles.edit2,'String') ') Corrected Correlation Matrix'],'FontSize',14);
            
            xlabel('Seeds'  ,'FontSize',14);
            ylabel('Seeds','FontSize',14);
            set(gca, 'color', [0 0 0])
            try set(gca, 'gridcolor', [1 1 1]); end
            try set(gca, 'gridalpha', 1); end
            colorbar
            daspect([1 1 1])
            if get(handles.AnaLe,'Value')
                caxis([-1.*(str2num(get(handles.edit2,'String'))) (-1*(min(adj_p)-(0.2*min(adj_p))))])
            else
                tmpFDRcor = FDRcor;
                tmpFDRcor(tmpFDRcor==-1) = NaN;
                caxis([-1.*(str2num(get(handles.edit2,'String'))) ((max(max(tmpFDRcor))))])
            end
            set(gcf, 'InvertHardCopy', 'off');
            
            drawnow
            imgRR = getframe(fig2);
            imwrite(imgRR.cdata, [DirF filesep FigNTpx '_FDR_corrected_P_Map' '.png']);
            clear imgRR

            saveas(fig2,[DirF filesep filesep FigNTpx '_FDR_corrected_P_Map'],'fig')
            save([DirF filesep filesep FigNTpx '_FDRcorPmap-VAR'],'FDRcor')
            close(fig2)
        end
        
        
        if get(handles.CircConnB,'Value') || get(handles.OrgCircConnB,'Value')
            
            matrixR = FDRcor;
            matrixR(matrixR==-1) = 0;
            matrixR = abs(matrixR);
            
            limiteG = (size(g1_VAR2,2)/2).*1.1;
            
            RAIO = limiteG/2.2;
            ang=0:0.01:2*pi; 
            xp=RAIO*cos(ang);
            yp=RAIO*sin(ang);
            xpt=(limiteG/2.15)*cos(ang);
            ypt=(limiteG/2.15)*sin(ang);

            AnoRa = size(xp,2)/size(matrixR,2);
            xp2 = xp(1:AnoRa:end);
            yp2 = yp(1:AnoRa:end);

            xp2t = xpt(1:AnoRa:end);
            yp2t = ypt(1:AnoRa:end);
            xp2t = xp2t+(limiteG/2);
            yp2t = yp2t+(limiteG/2);
            xp2t = xp2t(end:-1:1);
            yp2t = yp2t(end:-1:1);

            xp2 = xp2+(limiteG/2);
            yp2 = yp2+(limiteG/2);
            yp2 = yp2(end:-1:1);
            xp2 = xp2(end:-1:1);

            minH = (pdist([xp2(1),yp2(1);xp2(2),yp2(2)],'euclidean'));
            Diam = max(xp) - min(xp);

            stIndx = 1;

            for hg = 1:size(nEle,2)
                str_nEle(hg) = stIndx;
                stIndx = stIndx+nEle(hg);
            end

            cc=hsv(nOfNets);
            ccV1 = cc(1:2:end,:);
            ccV2 = cc(2:2:end,:);
            cc = [ccV1;ccV2];

            if get(handles.tarNet,'Value')
                LoopAte = fact+TarNetCnEle-1;
                TranspMin = 0.5;
            else
                LoopAte = size(matrixR,1);
                TranspMin = 0;
            end
        end
        
        if get(handles.CircConnB,'Value')
            
            CircConnFuncFDR(matrixR,limiteG,ScreSize,cc,...
                get(handles.edit2,'String'),nEle,LoopAte,fact,...
                str_nEle,get(handles.tarNet,'Value'),TarNetC,...
                DirF,FigNTpx,NetLabels,TranspMin);
            
            if get(handles.AnaLe,'Value')
                % T-scored Circular
                matrixR2 = Pmap;
                matrixR2(matrixR2>crit_p) = 0;
                
                CircConnTscoreFDR(matrixR2,limiteG,ScreSize,cc,g1_VAR2,...
                    g2_VAR2,nEle,LoopAte,fact,str_nEle,get(handles.edit2,'String'),...
                    get(handles.tarNet,'Value'),TarNetC,...
                    DirF,NetLabels,crit_p,g1_mean,g2_mean);
            end
        end

        if get(handles.OrgCircConnB,'Value')

            circ2 = figure('Visible','on');
            set(circ2,'Name','Connectome2D',...
                'Position',round([ScreSize(1)*.05 ScreSize(2)*.05 ScreSize(1)*.45 ScreSize(1)*.45]),...
                      'Color',[0 0 0]);
            set(gca,'xticklabels',[])
            set(gca,'yticklabels',[])
            axis off
            hold on

            xlim([-1 limiteG+1]);
            ylim([-1 limiteG+1]);

            daspect([1 1 1])
            cp1 = plot((limiteG/2)+xp,(limiteG/2)+yp,'LineWidth',1.5,'Color',[1 1 1]);

            daspect([1 1 1])
            RAIO2 = RAIO*0.15;        
            xp3=RAIO2*cos(ang);
            yp3=RAIO2*sin(ang);

            xp3 = xp3(1:AnoRa:end);
            yp3 = yp3(1:AnoRa:end);
            xp3 = flip(xp3+(limiteG/2));
            yp3 = flip(yp3+(limiteG/2));

            RAIO3 = RAIO*0.85;        
            xp4=RAIO3*cos(ang);
            yp4=RAIO3*sin(ang);

            idxNR = ceil((str_nEle+(nEle./2))-1);
            xp4 = flip(xp4(1:AnoRa:end));
            yp4 = flip(yp4(1:AnoRa:end));
            xp4 = xp4(idxNR);
            yp4 = yp4(idxNR);

            xp4 = xp4+(limiteG/2);
            yp4 = yp4+(limiteG/2);

            RAIO4 = RAIO*0.95;        
            xp5=RAIO4*cos(ang);
            yp5=RAIO4*sin(ang);

            xp5 = xp5((AnoRa/2):AnoRa:end);
            yp5 = yp5((AnoRa/2):AnoRa:end);
            xp5 = flip(xp5+(limiteG/2));
            yp5 = flip(yp5+(limiteG/2));
            
             mTextBo21 = uicontrol('style','text',...
                 'position',round([ScreSize(1)*.001 ScreSize(1)*.4375 ScreSize(1)*.3  ScreSize(1)*.0125]),...
                 'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
                 'ForegroundColor',[1 .949 .867],'FontWeight','bold');
             set(mTextBo21,'String','Organic Circular Connectome - DRAWING........');

             mTextBo22 = uicontrol('style','text',...
                 'position',round([ScreSize(1)*.001 ScreSize(1)*.422 ScreSize(1)*.2  ScreSize(1)*.0125]),...
                 'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
                 'ForegroundColor',[1 .949 .867],'FontWeight','bold');
            set(mTextBo22,'String',['FDR-Corrected (alpha ' get(handles.edit2,'String') ')']);
            drawnow
            for i = fact:LoopAte
                drawnow
                tester = i<=(str_nEle+nEle)-1;
                networkNumb = find(tester==1,1);
                colorV = cc(networkNumb,:);
                if i<10
                    str = ['0' num2str(i)];
                else
                    str = num2str(i);
                end
                vet = matrixR(i-(fact-1),:);
                vet = find(vet);

                if ~isempty(vet)
                    bt = text(xp2t(i),yp2t(i),[str '-' NetLabels{networkNumb}],'FontSize',9.5,'Color',colorV,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
                    Angle = atan2(yp2t(i)-(limiteG/2),xp2t(i)-(limiteG/2));
                    if abs(Angle)>pi/2
                        Angle = Angle - pi;
                        set(bt,'HorizontalAlignment','right')
                    end
                    set(bt,'Rotation',rad2deg(Angle))
                    for j = 1:size(vet,2)

                        P1 = [xp2(i),yp2(i)];
                        P2 = [xp2(vet(j)),yp2(vet(j))];

                        testerj = vet(j)<=(str_nEle+nEle)-1;
                        networkNumbj = find(testerj==1,1);
                        colorVj = cc(networkNumbj,:);
                        xxF = [];
                        yyF = [];
                        if networkNumb ~= networkNumbj
                            opt1 =  round(i-(i-vet(j))/2);
                            opt2 = round(opt1+(LoopAte/2));
                            if opt2<=LoopAte
                                dopt1 = pdist([xp3(opt1),yp3(opt1);P1(1),P1(2)]);
                                dopt2 = pdist([xp3(opt2),yp3(opt2);P1(1),P1(2)]);
                                if dopt1>dopt2
                                    dopt = opt2;
                                else
                                    dopt = opt1;
                                end
                            else
                                dopt = opt1;
                            end

                            [xx,yy]=fillline([P1(1),P1(2)],[xp4(networkNumb),yp4(networkNumb)],15);
                            [xx2,yy2]=fillline([xp4(networkNumb),yp4(networkNumb)],[xp3(dopt),yp3(dopt)],20);
                            [xx3,yy3]=fillline([xp3(dopt),yp3(dopt)],[xp4(networkNumbj),yp4(networkNumbj)],20);
                            [xx4,yy4]=fillline([xp4(networkNumbj),yp4(networkNumbj)],[P2(1),P2(2)],15);
                            xxF = [xx,xx2,xx3,xx4];
                            yyF = [yy,yy2,yy3,yy4];
                            xxF = [xxF(1);smoothLNI(xxF(2:end-1),30);xxF(end)]';
                            yyF = [yyF(1);smoothLNI(yyF(2:end-1),30);yyF(end)]';
                            TranOrd = 'descend';
                            TrnaspSca = [zeros(1,size(xxF,2)*0.3)+TranspMin,TranspMin:(1-TranspMin)/(size(xxF,2)*0.7):1];
                            TrnaspSca = sort(TrnaspSca,TranOrd);
                            figure(circ2)
                            if networkNumb~=networkNumbj
                                if get(handles.tarNet,'Value')
                                    for hg = 1:(size(xxF,2)-1)
                                        h = plot(xxF(1,hg:hg+1),yyF(1,hg:hg+1),'Color',colorV,'LineWidth',1.5);
                                        h.Color(4) = TrnaspSca(hg);
                                    end
                                else
                                    for hg = 1:(size(xxF,2)-1)*0.7
                                        h = plot(xxF(1,hg:hg+1),yyF(1,hg:hg+1),'Color',colorV,'LineWidth',1.5);
                                        h.Color(4) = TrnaspSca(hg);
                                    end
                                end
                            else
                                if i<vet(j)
                                    h = plot(xxF,yyF,'Color',colorV,'LineWidth',1);
                                end
                            end

                        else

                            ARCalpha = pdist([P1(1),P1(2);P2(1),P2(2)],'euclidean');

                            h = minH + (((Diam/ARCalpha)*(minH/Diam)));
                            R = ((ARCalpha^2) + (4*(h^2)))/(8*h);
                            NewC = circ_cent(P1,P2,R);
                            Marg = (limiteG-2*RAIO)/2;
                            DisC1 = pdist([NewC(1,1),NewC(1,2);RAIO+Marg,RAIO+Marg],'euclidean');
                            DisC2 = pdist([NewC(2,1),NewC(2,2);RAIO+Marg,RAIO+Marg],'euclidean');

                            if DisC1>DisC2
                                centerIdx = 1;
                            else
                                centerIdx = 2;
                            end

                            AngRad(1) = atan2(P1(2)-NewC(centerIdx,2),P1(1)-NewC(centerIdx,1));
                            AngRad(2) = atan2(P2(2)-NewC(centerIdx,2),P2(1)-NewC(centerIdx,1));
                            
                            if sum(abs(AngRad))>pi && any(AngRad<0)
                                idxCic = find(AngRad<0);
                                AngRad(idxCic) = pi +  pi-abs(AngRad(idxCic));
                            end

                            [AngRad,ksord] = sort(AngRad);

                            if sum(abs(AngRad))>pi && any(AngRad<0)
                                idxCic = find(AngRad<0);
                                AngRad(idxCic) = pi +  pi-abs(AngRad(idxCic));
                            end

                            ang2 = AngRad(1):0.03:AngRad(2);
                            ang2 = [ang2,AngRad(2)];

                            xunit = R * cos(ang2) + NewC(centerIdx,1);
                            yunit = R * sin(ang2) + NewC(centerIdx,2);
                            figure(circ2)
                            h = plot(xunit,yunit,'Color',colorV.*0.5,'LineWidth',2);
                        end

                        if vet(j)<10
                            strf = ['0' num2str(vet(j))];
                        else
                            strf = num2str(vet(j));
                        end
                        bt = text(xp2t(vet(j)),yp2t(vet(j)),[strf '-' NetLabels{networkNumbj}],'FontSize',9.5,'Color',...
                            colorVj,'FontWeight','bold','HorizontalAlignment','left',...
                            'VerticalAlignment','middle');

                        Angle = atan2(yp2t(vet(j))-(limiteG/2),xp2t(vet(j))-(limiteG/2));

                        if abs(Angle)>pi/2
                            Angle = Angle - pi;
                            set(bt,'HorizontalAlignment','right')
                        end
                        set(bt,'Rotation',rad2deg(Angle))
                    end
                end
            end
            
            set(mTextBo21,'String','Organic Circular Connectome');
            drawnow

            set(gcf,'InvertHardCopy', 'off');
            drawnow
            imgRR = getframe(circ2);
            imwrite(imgRR.cdata, [DirF filesep FigNTpx '_FDR_corrected_OrganicCirc_Map' '.png']);
            clear imgRR
%             saveas(circ,[DirF filesep FigNTpx '_FDR_corrected_Circ_Map'],'fig')
            close(circ2)
        end
        
        if get(handles.Graphos3D,'Value') %%%%% FDR corrected 3D Grapho
            
            GraphoFig3 = figure;
            set(GraphoFig3,'Name','Connectome3D',...
                'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.6 ScreSize(1)*.4]),...
                            'Color',[1 .949 .867]);
            
            title([FigNT ' (alpha ' get(handles.edit2,'String') ') Corrected 3D Connectome'],'FontSize',14);
            
            xlabel('X axis'  ,'FontSize',14);
            ylabel('Y axis','FontSize',14);
            zlabel('Z axis','FontSize',14);
            set(gca, 'color', [1 1 1])
            try set(gca, 'gridcolor', [0.8 0.8 0.8]); end
            try set(gca, 'gridalpha', 1); end
            xlim([0 91]);
            ylim([0 109]);
            zlim([0 91]);
            hold on
            axis off
            view([-1 -1 1])

             stIndx = 1;
             for hg = 1:size(nEle,2)
                 str_nEle(hg) = stIndx;
                 stIndx = stIndx+nEle(hg);
             end

            cc=hsv(nOfNets);
            ccV1 = cc(1:2:end,:);
            ccV2 = cc(2:2:end,:);
            cc = [ccV1;ccV2];
            
             [rws,res,rdw] = sphere(15);
             
             if  get(handles.tarNet,'Value')
                for huy = fact:size(g1_VAR2,1)+(fact-1)
                    CoordQ = CoordMS(huy,:);
                    CoordList = CoordMS(1:end,:);
                    for hgf = 1:size(CoordList,1)
                        isF = [CoordQ(1,1);CoordList(hgf,1)];
                        jsF = [CoordQ(1,2);CoordList(hgf,2)];
                        ksF = [CoordQ(1,3);CoordList(hgf,3)];
                       if FDRcor(huy-(fact-1),hgf)>-1.*(str2num(get(handles.edit2,'String')))
                            [xre,yre,zre] = tubeplot([isF(1) isF(2);jsF(1) jsF(2);ksF(1) ksF(2)],0.4,8);
                            if (g1_mean(huy-(fact-1),hgf)*g2_mean(huy-(fact-1),hgf))>0
                                if abs(g1_mean(huy-(fact-1),hgf))>abs(g2_mean(huy-(fact-1),hgf))
                                    ySurface1 = surf(xre,yre,zre);
                                    set(ySurface1,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                        'gouraud','FaceColor',[0.1 0.5 1],'FaceLighting','gouraud','LineSmoothing','on');
                                    ex3_1 = 1;
                                else
                                    ySurface2 = surf(xre,yre,zre);
                                    set(ySurface2,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                        'gouraud','FaceColor',[1 0.1 0.1],'FaceLighting','gouraud','LineSmoothing','on');
                                    ex3_2 = 1;
                                end
                            else
                                ySurface3 = surf(xre,yre,zre);
                                set(ySurface3,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                        'gouraud','FaceColor',[1 0.8 0],'FaceLighting','gouraud','LineSmoothing','on');
                                ex3_3 = 1;
                            end
                       end   
                    end
                end
            else
                for huy = 1:size(CoordMS,1)
                    CoordQ = CoordMS(huy,:);
                    CoordList = CoordMS(huy+1:end,:);
                    for hgf = 1:size(CoordList,1)
                        isF = [CoordQ(1,1);CoordList(hgf,1)];
                        jsF = [CoordQ(1,2);CoordList(hgf,2)];
                        ksF = [CoordQ(1,3);CoordList(hgf,3)];
                       if FDRcor(hgf+huy,huy)>-1.*(str2num(get(handles.edit2,'String')))
                            [xre,yre,zre] = tubeplot([isF(1) isF(2);jsF(1) jsF(2);ksF(1) ksF(2)],0.4,8);
                            if (g1_mean(hgf+huy,huy)*g2_mean(hgf+huy,huy))>0
                                if abs(g1_mean(hgf+huy,huy))>abs(g2_mean(hgf+huy,huy))
                                    ySurface1 = surf(xre,yre,zre);
                                    set(ySurface1,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                        'gouraud','FaceColor',[0.1 0.5 1],'FaceLighting','gouraud','LineSmoothing','on');
                                    ex3_1 = 1;
                                else
                                    ySurface2 = surf(xre,yre,zre);
                                    set(ySurface2,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                        'gouraud','FaceColor',[1 0.1 0.1],'FaceLighting','gouraud','LineSmoothing','on');
                                    ex3_2 = 1;
                                end
                            else
                                ySurface3 = surf(xre,yre,zre);
                                set(ySurface3,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                        'gouraud','FaceColor',[1 0.8 0],'FaceLighting','gouraud','LineSmoothing','on');
                                ex3_3 = 1;
                            end
                       end   
                    end
               end
           end
            %%%%%%%%%%%%%%%%%%%%%%%  Include only spheres with alterations
            teste = FDRcor;
            teste(teste==-1)=0;
            teste(teste~=0)=1;
            vetTTT = sum(teste,1);
            VetFF = vetTTT;
            VetFF(VetFF>1) = 1;
            idxROIs = find(VetFF);

            vetTT2 = sum(teste,2);
            VetFF2 = vetTT2;
            VetFF2(VetFF2>1) = 1;
            idxROIs2 = find(VetFF2);
            idxROIs2 = idxROIs2+(fact-1);
            vetTT2(vetTT2==0) = [];
            VetFF(idxROIs2) = 1;
            vetTTT(idxROIs2) = vetTT2;

            idxROIs = [idxROIs,idxROIs2'];
            idxROIs = unique(idxROIs);

            for i = 1:numel(idxROIs)
                tester = idxROIs(i)<=(str_nEle+nEle)-1;
                networkNumb = find(tester==1,1);
                colorV = cc(networkNumb,:);
                CoordQ = CoordMS(idxROIs(i),:);
                CoordV4 = (vetTTT(idxROIs(i))/max(vetTTT))+1;
                rNumb = find(VetFF~=0,i);
                rNumb = rNumb(end);

                if isequal(networkNumb,1)
                    if isequal(REtype,2)
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(3*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(rNumb) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    else
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(3*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(i) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    end
                else
                    if isequal(REtype,2)
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(5*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(rNumb-sum(nEle(1:networkNumb-1))) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    else
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(5*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(idxROIs(i)) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    end
                end
                hSurface = surf(CoordV4.*rws+CoordQ(1,1),CoordV4.*res+CoordQ(1,2),CoordV4.*rdw+CoordQ(1,3));
                set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if exist('ex3_1') || exist('ex3_2')
                if exist('ex3_1') && exist('ex3_2')
                    if exist('ex3_3')
                        legend([ySurface1,ySurface2,ySurface3],[G1nameT ' > ' G2nameT],[G2nameT  ' > ' G1nameT],[G2nameT ' ' char(8596) ' -(' G1nameT ')'])
                    else
                        legend([ySurface1,ySurface2],[G1nameT ' > ' G2nameT],[G2nameT  ' > ' G1nameT])
                    end
                else if exist('ex3_1')
                        if exist('ex3_3')
                            legend([ySurface1,ySurface3],[G1nameT ' > ' G2nameT],[G2nameT ' ' char(8596) ' -(' G1nameT ')'])
                        else
                            legend(ySurface1,[G1nameT ' > ' G2nameT])
                        end
                    else
                         if exist('ex3_3')
                             legend([ySurface2,ySurface3],[G2nameT ' > ' G1nameT],[G2nameT ' ' char(8596) ' -(' G1nameT ')'])
                         else
                            legend(ySurface2,[G2nameT ' > ' G1nameT])
                         end
                     end
                end
            else
                if exist('ex3_3')
                    legend(ySurface3,[G2nameT ' ' char(8596) ' -(' G1nameT ')']) % Simbolos apenas copiando e colando...
                end
            end
            
            mTextBox = uicontrol('style','text','position',round([ScreSize(1)*.01 ScreSize(2)*.01 ScreSize(1)*.15 ScreSize(1)*.01]));
            set(mTextBox,'String','r = Roi number and n = network number');
            
            DirUF2C = which('uf2c');
            try
                aRT = nifti([DirUF2C(1:end-6) 'Utilities' filesep 'WFUteste2.nii']);
            catch
                aRT = nifti([DirUF2C(1:end-13) 'Utilities' filesep 'WFUteste2.nii']);
            end

            matI = aRT.dat(:,:,:);
            matI = permute(matI,[1,3,2]);
            matI = permute(matI,[2,1,3]);
            matI = permute(matI,[3,2,1]);
            matI = matI.*100;
            matI = int8(matI);
            Ds = smooth3(matI);
            Ds = int8(Ds);
            Xiso = isosurface(Ds,10);
            hiso = patch(Xiso,'FaceColor',[0.6,0.6,0.6],'EdgeColor','none','FaceAlpha',0.3);%,'FaceVertexAlphaData',0.8,'AlphaDataMapping','Direct','CDataMapping','Direct');
            daspect([1,0.83,0.83])
            view(113,14)
            lighting gouraud
            set(gca,'xticklabels',[])
            set(gca,'yticklabels',[])

            light('Style','infinite');%'Position',[1 0 1]);
            
            try
                legTmp = legend;
                IndexLegC = strfind(legTmp.String,'data1');
                IndexLeg = find(not(cellfun('isempty',IndexLegC)));
                legTmp.String(IndexLeg) = []; 
            end
            
            saveas(GraphoFig3,[DirF filesep FigNTpx '_FDR_Corrected_Graph_3D'],'fig')
            set(gcf, 'InvertHardCopy', 'off');
            
            drawnow
            imgRR = getframe(GraphoFig3);
            imwrite(imgRR.cdata, [DirF filesep FigNTpx '_FDR_Corrected_Graph_3D' '.png']);
            clear imgRR
            
            view(90,0)
            drawnow
            set(gcf, 'InvertHardCopy', 'off');
            imgRR = getframe(GraphoFig3);
            imwrite(imgRR.cdata, [DirF filesep FigNTpx '_FDR_Corrected_Graph_3D_Sag' '.png']);
            clear imgRR
            
            mTextBoxOri13 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.05 ScreSize(2)*.101 ScreSize(1)*.03 ScreSize(1)*.05]),...
             'FontSize',60,'BackgroundColor',[1 .949 .867],...
             'ForegroundColor',[0 0 0]);
            set(mTextBoxOri13,'String','L');
            
            mTextBoxOri23 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.4 ScreSize(2)*.101 ScreSize(1)*.03 ScreSize(1)*.05]),...
             'FontSize',60,'BackgroundColor',[1 .949 .867],...
             'ForegroundColor',[0 0 0]);
            set(mTextBoxOri23,'String','R');
            
            view(0,90)
            drawnow
            set(gcf, 'InvertHardCopy', 'off');
            imgRR = getframe(GraphoFig3);
            imwrite(imgRR.cdata, [DirF filesep FigNTpx '_FDR_Corrected_Graph_3D_Ax' '.png']);
            clear imgRR
            
            view(0,0)
            axis off
            drawnow
            set(gcf, 'InvertHardCopy', 'off');
            imgRR = getframe(GraphoFig3);
            imwrite(imgRR.cdata, [DirF filesep FigNTpx '_FDR_Corrected_Graph_3D_Cor' '.png']);
            clear imgRR

            close(GraphoFig3)
            
            CoLeg = figure('Visible','off');
            set(CoLeg,'Name','Color Legend','Position', [400 400 250 400],...
                            'Color',[1 .949 .867]);
            title('Color Legend','FontSize',14);
            xlim([0 0.6]);
            ylim([0 nOfNets]);

            for uiC = 1:nOfNets
                rectangle('Position',[0,uiC-1,0.6,1],'LineWidth',2,'FaceColor',cc(uiC,:))
                string = sprintf('Network %d color',uiC);
                text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
            end
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])

            drawnow
            imgRR = getframe(CoLeg);
            imwrite(imgRR.cdata, [DirF filesep 'Color_Legend' '.png']);
            clear imgRR

            saveas(CoLeg,[DirF filesep 'Color_Legend'],'fig')
            close(CoLeg)
            
        end
        
        if get(handles.Grapho2D,'Value')  %%%%% FDR corrected 2D Grapho
            CoordMS2 = 3.*CoordMS;
            GraphoFig4 = figure;
            set(GraphoFig4,'Name','Connectome2D',...
                'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.45 ScreSize(1)*.45]),...
                            'Color',[1 .949 .867]);
                        
            title([FigNT ' (alpha ' get(handles.edit2,'String') ') Corrected 2D Connectome'],'FontSize',14);
            
            xlabel('X axis'  ,'FontSize',14);
            ylabel('Y axis','FontSize',14);
            zlabel('Z axis','FontSize',14);
            set(gca, 'color', [1 1 1])
            try set(gca, 'gridcolor', [0 0 0]); end
            try set(gca, 'gridalpha', 1); end
            hold on
            xlim([0 3*91]);
            ylim([0 3*109]);
            zlim([0 3*91]);
            view([0 0 1])
            
             stIndx = 1;
             for hg = 1:size(nEle,2)
                 str_nEle(hg) = stIndx;
                 stIndx = stIndx+nEle(hg);
             end

            cc=hsv(nOfNets);
            ccV1 = cc(1:2:end,:);
            ccV2 = cc(2:2:end,:);
            cc = [ccV1;ccV2];
            [rws,res,rdw] = sphere(15);

             if get(handles.tarNet,'Value')
                for huy = fact:size(g1_VAR2,1)+(fact-1)
                    CoordQ2 = CoordMS2(huy,:);
                    CoordList = CoordMS2;
                    for hgf = 1:size(CoordList,1)
                        isF = [CoordQ2(1,1);CoordList(hgf,1)];
                        jsF = [CoordQ2(1,2);CoordList(hgf,2)];
                        ksF = [CoordQ2(1,3);CoordList(hgf,3)];
                       if FDRcor(huy-(fact-1),hgf)>-1.*(str2num(get(handles.edit2,'String')))
                            if (g1_mean(huy-(fact-1),hgf)*g2_mean(huy-(fact-1),hgf))>0
                                if abs(g1_mean(huy-(fact-1),hgf))>abs(g2_mean(huy-(fact-1),hgf))
                                   pLine1 =  line(isF,jsF,ksF,'Color',[0.2 0.5 1],'LineStyle','-','LineWidth',0.8);
                                   ex4_1 = 1;
                                else
                                   pLine2 = line(isF,jsF,ksF,'Color',[1 0.3 0.3],'LineStyle','-','LineWidth',0.8);
                                   ex4_2 = 1;
                                end
                            else
                                pLine3 =  line(isF,jsF,ksF,'Color',[0 0 0],'LineStyle','--','LineWidth',0.7);
                                ex4_3 = 1;
                            end
                       end   
                    end
                end
            else
                for huy = 1:size(CoordMS2,1)
                    CoordQ2 = CoordMS2(huy,:);
                    CoordList = CoordMS2(huy+1:end,:);
                    for hgf = 1:size(CoordList,1)
                        isF = [CoordQ2(1,1);CoordList(hgf,1)];
                        jsF = [CoordQ2(1,2);CoordList(hgf,2)];
                        ksF = [CoordQ2(1,3);CoordList(hgf,3)];
                        if FDRcor(hgf+huy,huy)>-1.*(str2num(get(handles.edit2,'String')))
                            if (g1_mean(hgf+huy,huy)*g2_mean(hgf+huy,huy))>0
                                if abs(g1_mean(hgf+huy,huy))>abs(g2_mean(hgf+huy,huy))
                                   pLine1 =  line(isF,jsF,ksF,'Color',[0.2 0.5 1],'LineStyle','-','LineWidth',0.8);
                                   ex4_1 = 1;
                                else
                                   pLine2 = line(isF,jsF,ksF,'Color',[1 0.3 0.3],'LineStyle','-','LineWidth',0.8);
                                   ex4_2 = 1;
                                end
                            else
                                pLine3 =  line(isF,jsF,ksF,'Color',[0 0 0],'LineStyle','--','LineWidth',0.7);
                                ex4_3 = 1;
                            end
                        end   
                    end
                end
             end
            %%%%%%%%%%%%%%%%%%%%%%%  include only  spheres with alterations
            teste = FDRcor;
            teste(teste==-1)=0;
            teste(teste~=0)=1;
            vetTTT = sum(teste,1);
            VetFF = vetTTT;
            VetFF(VetFF>1) = 1;
            idxROIs = find(VetFF);

            vetTT2 = sum(teste,2);
            VetFF2 = vetTT2;
            VetFF2(VetFF2>1) = 1;
            idxROIs2 = find(VetFF2);
            idxROIs2 = idxROIs2+(fact-1);
            vetTT2(vetTT2==0) = [];
            VetFF(idxROIs2) = 1;
            vetTTT(idxROIs2) = vetTT2;

            idxROIs = [idxROIs,idxROIs2'];
            idxROIs = unique(idxROIs);


            for i = 1:sum(VetFF~=0)
                tester = idxROIs(i)<=(str_nEle+nEle)-1;
                networkNumb = find(tester==1,1);
                colorV = cc(networkNumb,:);
                CoordQ = CoordMS2(idxROIs(i),:);
                CoordV4 = (vetTTT(idxROIs(i))/max(vetTTT))+1;
                rNumb = find(VetFF~=0,i);
                rNumb = rNumb(end);

                if isequal(networkNumb,1)
                    if isequal(REtype,2)
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(3*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(rNumb) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    else
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(3*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(i) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    end
                else
                    if isequal(REtype,2)
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(5*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(rNumb-sum(nEle(1:networkNumb-1))) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    else
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(5*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),...
                            ['r' num2str(idxROIs(i)) NetLabels{networkNumb}],'HorizontalAlignment','Right','FontSize',9);
                    end
                end
                hSurface = surf(CoordV4*3.*rws+CoordQ(1,1),CoordV4*3.*res+CoordQ(1,2),CoordV4*3.*rdw+CoordQ(1,3));
                set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            if exist('ex4_1') || exist('ex4_2')
                if exist('ex4_1') && exist('ex4_2')
                    if exist('ex4_3')
                        legend([pLine1,pLine2,pLine3],[G1nameT ' > ' G2nameT],[G2nameT ' > ' G1nameT],[G2nameT ' ' char(8596) ' -(' G1nameT ')'])
                    else
                        legend([pLine1,pLine2],[G1nameT ' > ' G2nameT],[G2nameT ' > ' G1nameT])
                    end
                else if exist('ex4_1')
                        if exist('ex4_3')
                            legend([pLine1,pLine3],[G1nameT ' > ' G2nameT],[G2nameT ' ',char(8596),' -(' G1nameT ')'])
                        else
                            legend(pLine1,[G1nameT ' > ' G2nameT])
                        end
                     else
                        if exist('ex4_3')
                            legend([pLine2,pLine3],[G2nameT ' > ' G1nameT],[G2nameT ' ' char(8596) ' -(' G1nameT ')'])
                        else
                            legend(pLine2,[G2nameT ' > ' G1nameT])
                        end
                    end
                 end
            else if exist('ex4_3')
                     legend(pLine3,[G2nameT ' ' char(8596) ' -(' G1nameT ')'])
                 end
            end

            mTextBox = uicontrol('style','text','position',round([ScreSize(1)*.01 ScreSize(2)*.01 ScreSize(1)*.15 ScreSize(1)*.01]));
            set(mTextBox,'String','r = Roi number and n = network number');
            
            DirUF2C = which('uf2c');

            try
                load([DirUF2C(1:end-6) 'Utilities' filesep 'Brain_Contour.mat'])
            catch
                load([DirUF2C(1:end-13) 'Utilities' filesep 'Brain_Contour.mat'])
            end

            daspect([1,0.95,1])
            ater2 = slice(mat3d,[],[],1);
            set(ater2, 'FaceAlpha',0.7);
            colormap gray
            set(ater2, 'LineStyle','none');
            matConto = mat3d(:,:,1);%round(size(final3D,3)/2));
            matConto = squeeze(matConto);
            matConto(matConto<1)=0.0001;
            set(ater2(1),'alphadata',matConto,'facealpha','flat')
            grid off
            axis off
            set(gca,'xticklabels',[])
            set(gca,'yticklabels',[])

            light('Style','infinite');
            
            try
                legTmp = legend;
                IndexLegC = strfind(legTmp.String,'data1');
                IndexLeg = find(not(cellfun('isempty',IndexLegC)));
                legTmp.String(IndexLeg) = []; 
            end
            
            mTextBoxOri14 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.05 ScreSize(2)*.101 ScreSize(1)*.03 ScreSize(1)*.05]),...
             'FontSize',60,'BackgroundColor',[1 .949 .867],...
             'ForegroundColor',[0 0 0]);
            set(mTextBoxOri14,'String','L');
            
            mTextBoxOri24 = uicontrol('style','text',...
             'position',round([ScreSize(1)*.4 ScreSize(2)*.101 ScreSize(1)*.03 ScreSize(1)*.05]),...
             'FontSize',60,'BackgroundColor',[1 .949 .867],...
             'ForegroundColor',[0 0 0]);
            set(mTextBoxOri24,'String','R');
            
            set(gcf, 'InvertHardCopy', 'off');
            saveas(GraphoFig4,[DirF filesep FigNTpx '_FDR_Corrected_Graph_2D'],'fig')
            drawnow
            imgRR = getframe(GraphoFig4);
            imwrite(imgRR.cdata, [DirF filesep FigNTpx '_FDR_Corrected_Graph_2D' '.png']);
            clear imgRR

            close(GraphoFig4)
            
            CoLeg = figure('Visible','off');
            set(CoLeg,'Name','Color Legend','Position', [400 400 250 400],...
                            'Color',[1 .949 .867]);
            title('Color Legend','FontSize',14);
            xlim([0 0.6]);
            ylim([0 nOfNets]);

            for uiC = 1:nOfNets
                rectangle('Position',[0,uiC-1,0.6,1],'LineWidth',2,'FaceColor',cc(uiC,:))
                string = sprintf('Network %d color',uiC);
                text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
            end
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            
            drawnow
            imgRR = getframe(CoLeg);
            imwrite(imgRR.cdata, [DirF filesep 'Color_Legend' '.png']);
            clear imgRR

            saveas(CoLeg,[DirF filesep 'Color_Legend'],'fig')
            close(CoLeg)
            
        end
    else
        warndlg('No significant results with FDR correction','Ops!')
    end
end
fclose('all');

set(handles.text13,'String','Done!')
disp('------------')
disp('You can open the *.fig files on the Matlab to interact (rotate, zoom, edit...) with the figures!')
disp('Done!')

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

function Grapho2D_Callback(hObject, eventdata, handles)
global g2OK g1OK  dirOK nodeOK

if isequal(get(handles.Grapho2D,'Value'),0)
    if logical(g2OK) && logical(g1OK) && logical(dirOK)
        if isequal(get(handles.Graphos3D,'Value'),0)
            set(handles.runB,'Enable','on')
        end
    end
else
    if logical(g2OK) && logical(g1OK) && logical(dirOK) && logical(nodeOK)
        set(handles.runB,'Enable','on')
    else
        set(handles.runB,'Enable','off')
    end
end

function Graphos3D_Callback(hObject, eventdata, handles)
global g2OK g1OK dirOK nodeOK

if isequal(get(handles.Graphos3D,'Value'),0)
    if logical(g2OK) && logical(g1OK) && logical(dirOK)
        if isequal(get(handles.Grapho2D,'Value'),0)
            set(handles.runB,'Enable','on')
        end
    end
else
    if logical(g2OK) && logical(g1OK) && logical(dirOK) && logical(nodeOK)
        set(handles.runB,'Enable','on')
    else
        set(handles.runB,'Enable','off')
    end
end

function TSTTbtm_Callback(hObject, eventdata, handles)
if get(handles.TSTTbtm,'Value')
    set(handles.COV,'Enable','on')
    set(handles.PTTbtm,'Value',0)
else
    set(handles.COV,'Enable','off')
    set(handles.COV,'Value',0)
    set(handles.TSTTbtm,'String','Two Sample t-test')
    set(handles.PTTbtm,'Value',1)
end


function PTTbtm_Callback(hObject, eventdata, handles)
if get(handles.PTTbtm,'Value')
    set(handles.COV,'Enable','off')
    set(handles.COV,'Value',0)
    set(handles.TSTTbtm,'Value',0)
    set(handles.TSTTbtm,'String','Two Sample t-test')
else
    set(handles.COV,'Enable','on')
    set(handles.TSTTbtm,'Value',1)
end

function figure1_CreateFcn(hObject, eventdata, handles)
global g2OK g1OK dirOK nodeOK
g2OK=0; g1OK=0; dirOK=0; nodeOK=0;

function AnaLe_Callback(hObject, eventdata, handles)
if get(handles.AnaLe,'Value')
    set(handles.RoiLe,'Value',0)
else
    set(handles.RoiLe,'Value',1)
end

function RoiLe_Callback(hObject, eventdata, handles)
if get(handles.RoiLe,'Value')
    set(handles.AnaLe,'Value',0)
else
    set(handles.AnaLe,'Value',1)
end

function chec1_Callback(hObject, eventdata, handles)
if get(handles.chec1,'Value')
    set(handles.edit1,'Enable','on')
    set(handles.text7,'Enable','on')
else
    set(handles.edit1,'Enable','off')
    set(handles.text7,'Enable','off')
end

function chec2_Callback(hObject, eventdata, handles)
if get(handles.chec2,'Value')
    set(handles.AnaLe,'Value',1)
    set(handles.RoiLe,'Value',0)
    set(handles.AnaLe,'Enable','on')
    set(handles.RoiLe,'Enable','on')
    set(handles.edit2,'Enable','on')
    set(handles.text8,'Enable','on')

else
    set(handles.AnaLe,'Enable','off')
    set(handles.RoiLe,'Enable','off')
    set(handles.AnaLe,'Value',1)
    set(handles.RoiLe,'Value',0)
    set(handles.edit2,'Enable','off')
    set(handles.text8,'Enable','off')
end



function CorMats_Callback(hObject, eventdata, handles)
function edit1_Callback(hObject, eventdata, handles)
function edit2_Callback(hObject, eventdata, handles)

function G1name_Callback(hObject, eventdata, handles)
global Dir dirOK DirTmp
if exist('dirOK','var')
    Dir = [DirTmp filesep 'UF2C_SL_' get(handles.G1name,'String') 'x' get(handles.G2name,'String')];
    set(handles.text10,'String',Dir)
    drawnow
end

function G2name_Callback(hObject, eventdata, handles)
global Dir dirOK DirTmp
if exist('dirOK','var')
    Dir = [DirTmp filesep 'UF2C_SL_' get(handles.G1name,'String') 'x' get(handles.G2name,'String')];
    set(handles.text10,'String',Dir)
    drawnow
end

function tarNet_Callback(hObject, eventdata, handles)
global g2OK g1OK nodeOK

if get(handles.tarNet,'Value')
    if get(handles.checkIfNet,'Value')
        if g1OK || g2OK
            set(handles.TarNetPop,'Enable','on')
        end
    else
        if nodeOK
            set(handles.TarNetPop,'Enable','on')
        end
    end
    set(handles.AnaLe,'Value',1)
    set(handles.RoiLe,'Value',0)
    set(handles.AnaLe,'Enable','off')
    set(handles.RoiLe,'Enable','off')
else
    set(handles.TarNetPop,'Enable','off')
    set(handles.AnaLe,'Value',1)
    set(handles.RoiLe,'Value',0)
    set(handles.AnaLe,'Enable','on')
    set(handles.RoiLe,'Enable','on')
end

function checkIfNet_Callback(hObject, eventdata, handles)
global g2OK dirOK g1OK nodeOK g1_VAR g2_VAR nOfNets NetLabels

if get(handles.checkIfNet,'Value')
    
    set(handles.ROIfiles,'Enable','off')
    set(handles.CoodList,'Enable','off')
    set(handles.addNODES,'Enable','off')
    set(handles.Grapho2D,'Enable','off')
    set(handles.Graphos3D,'Enable','off')
    set(handles.Grapho2D,'Value',0)
    set(handles.Graphos3D,'Value',0)
    if g2OK && dirOK && g1OK
        set(handles.runB,'Enable','on')
    end
    if g1OK
        if size(g1_VAR,1)>1
            nOfNets = size(g1_VAR,1);
            if get(handles.NetLab,'Value')
                for nt = 1:nOfNets
                    strNetList{nt,:} = NetLabels{nt};
                end
            else
                for nt = 1:nOfNets
                    strNetList{nt,:} = sprintf('Network %d',nt);
                end
            end
            set(handles.TarNetPop,'String',strNetList)
            if get(handles.tarNet,'Value')
                set(handles.TarNetPop,'Enable','on')
            end
            set(handles.NetLab,'Enable','on')
            
        end
    end
    if g2OK
        if size(g2_VAR,1)>1
            nOfNets = size(g1_VAR,1);
            if get(handles.NetLab,'Value')
                for nt = 1:nOfNets
                    strNetList{nt,:} = NetLabels{nt};
                end
            else
                for nt = 1:nOfNets
                    strNetList{nt,:} = sprintf('Network %d',nt);
                end
            end
            set(handles.TarNetPop,'String',strNetList)
            if get(handles.tarNet,'Value')
                set(handles.TarNetPop,'Enable','on')
            end
            set(handles.NetLab,'Enable','on')
        end
    end
else
    set(handles.NetLab,'Enable','off')
    set(handles.ROIfiles,'Enable','on')
    set(handles.CoodList,'Enable','on')
    set(handles.addNODES,'Enable','on')
    set(handles.Grapho2D,'Enable','on')
    set(handles.Graphos3D,'Enable','on')
    set(handles.Grapho2D,'Value',1)
    set(handles.Graphos3D,'Value',1)
    if g2OK && dirOK && g1OK && nodeOK
        set(handles.runB,'Enable','on')
    else
        set(handles.runB,'Enable','off')
    end
    if nodeOK && get(handles.tarNet,'Value')
        set(handles.TarNetPop,'Enable','on')
    else
         set(handles.TarNetPop,'Enable','off')
    end
end

function NetLab_Callback(hObject, eventdata, handles)
global nOfNets NetLabels
if get(handles.NetLab,'Value')
    [fileNetLab,pathNetLab] = uigetfile({'*.txt','Text file'},...
            'Add a text file witha Network Label on each line (avoid spaces)','MultiSelect','off');
    if ~isequal(fileNetLab,0)
        NetLabels = importdata([pathNetLab fileNetLab]);
        if isequal(nOfNets,numel(NetLabels))
            set(handles.NetLab,'ForeGroundColor',[0 0.8 0]);
            if nOfNets>1
                if get(handles.NetLab,'Value')
                    for nt = 1:nOfNets
                        strNetList{nt,:} = NetLabels{nt};
                    end
                else
                    for nt = 1:nOfNets
                        strNetList{nt,:} = sprintf('Network %d',nt);
                    end
                end
                set(handles.TarNetPop,'String',strNetList)
                set(handles.NetLab,'Enable','on')
            end

        else
            set(handles.NetLab,'ForeGroundColor',[1 0 0]);
            set(handles.NetLab,'Value',0);
            warndlg(sprintf('The number of labels added (%d) does not matches the numbers of Networks identfied (%d)',numel(NetLabels),nOfNets),'Attention');
        end
    else
        set(handles.NetLab,'Value',0);
    end
else
    set(handles.NetLab,'ForeGroundColor',[0 0 0]);
    set(handles.NetLab,'Value',0);
    for nt = 1:nOfNets
        strNetList{nt,:} = sprintf('Network %d',nt);
    end
    set(handles.TarNetPop,'String',strNetList)
end

function TarNetPop_Callback(hObject, eventdata, handles)

function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function G1name_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function G2name_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TarNetPop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CircConnB_Callback(hObject, eventdata, handles)
function OrgCircConnB_Callback(hObject, eventdata, handles)
