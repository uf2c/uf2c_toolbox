function varargout = CrossCorSecLe3(varargin)
% UF²C M-file for CrossCorSecLe3.fig
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
                   'gui_OpeningFcn', @CrossCorSecLe3_OpeningFcn, ...
                   'gui_OutputFcn',  @CrossCorSecLe3_OutputFcn, ...
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

function CrossCorSecLe3_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = CrossCorSecLe3_OutputFcn(hObject, eventdata, handles) 
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
    Dir = [Dir filesep 'UF2C_SecondLevel'];
    set(handles.text10,'String',Dir)
    dirOK = 1;
    if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
        set(handles.runB,'Enable','on')
    end
end

function addNODES_Callback(hObject, eventdata, handles)
global CoordMS seedList nodeOK g1OK g2OK dirOK nEle nOfNets

if get(handles.CoodList,'Value')
    set(handles.text12,'String','Wait.......')
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
    set(handles.text12,'String','Wait.......')
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
        seedList = transpose(fileROI);
        set(handles.text12,'String',sprintf('%d coordinates of %d netorworks were added!',size(CoordMS,1),nOfNets))
        if isequal(g1OK,1) && isequal(g2OK,1) && isequal(dirOK,1) && isequal(nodeOK,1)
            set(handles.runB,'Enable','on')
        end
    end
end
    
function runB_Callback(hObject, eventdata, handles)
global CoordMS seedList Dir COVDATA g1_VAR g2_VAR nOfLines nEle nOfNets G1covv G2covv
set(handles.text13,'String','Running.....')
%%%%%% Checking...
if ~isequal(size(CoordMS,1),size(g1_VAR,1))
    if logical(get(handles.Grapho2D,'Value')) || logical(get(handles.Graphos3D,'Value'))
        warndlg(sprintf('The number of nodes added (%d) should to match with the used on first level (%d)',size(CoordMS,1),size(g1_VAR,1)),'Ops! Process Aborted')
    return
    end
end
if get(handles.COV,'Value')
    if ~isequal((size(g1_VAR,3)+size(g2_VAR,3)),nOfLines)
        warndlg(sprintf('The number of data of each covariate (%d) should to match with the sum of Group 1 and 2 subjects',nOfLines,(size(g1_VAR,3)+size(g2_VAR,3))),'Ops! Process Aborted')
        return
    end
end
%%%%%%%%%%%%%%%%
mkdir(Dir)
%  g1_VAR(g1_VAR<0)=0;
%  g2_VAR(g2_VAR<0)=0;

COVnumber = (get(handles.COV,'Value') + get(handles.addROIwCov,'Value'));

if get(handles.addROIwCov,'Value')
    ROIwCov = cat(3,G1covv,G2covv);
%     ROIwCov = 1 - ROIwCov;
end

switch COVnumber
    case 0 % no covariates
        for i = 1:size(g1_VAR,1)
            for j = 1:size(g1_VAR,1)
                if i==j
                    Pmap(i,j) = 1;
                else
                    %[H,P,CI,tvalue] = ttest2(squeeze(g1_VAR(i,j,:)),squeeze(g2_VAR(i,j,:)));
                    [T,P] = mancovan([squeeze(g1_VAR(i,j,:)); squeeze(g2_VAR(i,j,:))],...
                    [ones(size(g1_VAR,3),1); zeros(size(g2_VAR,3),1)],[],'group-group');
                    Pmap(i,j) = P;
                end
            end
        end
     case 1 %one type of covariate
         for i = 1:size(g1_VAR,1)
             for j = 1:size(g1_VAR,1)
                if i==j
                    Pmap(i,j) = 1;
                else
                    if get(handles.COV,'Value')
                        [T,P] = mancovan([squeeze(g1_VAR(i,j,:)); squeeze(g2_VAR(i,j,:))],...
                        [ones(size(g1_VAR,3),1); zeros(size(g2_VAR,3),1)],COVDATA);
                        Pmap(i,j) = P(1);
                    else
                        [T,P] = mancovan([squeeze(g1_VAR(i,j,:)); squeeze(g2_VAR(i,j,:))],...
                        [ones(size(g1_VAR,3),1); zeros(size(g2_VAR,3),1)],squeeze(ROIwCov(i,j,:)));
                        Pmap(i,j) = P(1);
                    end
                end
             end
         end
      case 2
            for i = 1:size(g1_VAR,1)
                for j = 1:size(g1_VAR,1)
                    if i==j
                        Pmap(i,j) = 1;
                    else
                        [T,P] = mancovan([squeeze(g1_VAR(i,j,:)); squeeze(g2_VAR(i,j,:))],...
                        [ones(size(g1_VAR,3),1); zeros(size(g2_VAR,3),1)],[squeeze(ROIwCov(i,j,:)),COVDATA]);
                        Pmap(i,j) = P(1);
                    end
                end
            end
end

g1_mean = mean(g1_VAR,3);
g2_mean = mean(g2_VAR,3);
% Uncorrected P Map
if get(handles.chec1,'Value')
    
    unCorPmap = Pmap;
    unCorPmap(unCorPmap>str2num(get(handles.edit1,'String'))) = 1;
    unCorPmap = -1.*unCorPmap;

    if ~isequal(numel(unique(unCorPmap)),1)
        if get(handles.CorMats,'Value')   %%%%% Correlation Matrix
            fig1 = figure;
            imagesc(unCorPmap)
            title('Uncorrected P Map')
            colorbar
            caxis([-1.*(str2num(get(handles.edit1,'String'))) (max(max(unCorPmap))-(0.2*max(max(unCorPmap))))]) %reescala a volorbar para fugir dos extremos

            saveas(fig1,[Dir filesep 'Uncorrected_P_Map'],'tif')
            saveas(fig1,[Dir filesep 'Uncorrected_P_Map'],'fig')
            save([Dir filesep 'unCorPmap-VAR'],'unCorPmap')
            close(fig1)
        end
        if get(handles.Graphos3D,'Value')  %%%%% 3D Graphos
            GraphoFig = figure;
            set(GraphoFig,'Name','Connecticome3D','Position', [130 130 1350 800],...
                            'Color',[1 .949 .867]);
            title('Connecticome3D','FontSize',14);
            xlabel('X axis'  ,'FontSize',14);
            ylabel('Y axis','FontSize',14);
            zlabel('Z axis','FontSize',14);
            scatter3(CoordMS(:,1),CoordMS(:,2),CoordMS(:,3),'Marker','square','MarkerFaceColor','red','MarkerEdgeColor','k');
            set(gca, 'color', [1 1 1])
            set(gca, 'gridcolor', [0.8 0.8 0.8])
            set(gca, 'gridalpha', 1)
            xlim([0 91]);
            ylim([0 109]);
            zlim([0 91]);
            hold on

             stIndx = 1;
             for hg = 1:size(nEle,2)
                 str_nEle(hg) = stIndx;
                 stIndx = stIndx+nEle(hg);
             end

             cc=hsv(nOfNets);
             [rws,res,rdw] = sphere(15);

            for huy = 1:size(CoordMS,1)
                CoordQ = CoordMS(huy,:);
                CoordV4 = (sum(abs(unCorPmap(:,huy)))-1)/(size(CoordMS,1)-1);
                CoordList = CoordMS(huy+1:end,:);
                CoordV4 = sum(unCorPmap(:,huy)~=-1);
                CoordV4 = 1+ (2*CoordV4/size(CoordMS,1));
                tester = huy <=(str_nEle+nEle)-1;
                networkNumb = find(tester==1,1);
                colorV = cc(networkNumb,:);
                if isequal(networkNumb,1)
                    text(CoordQ(1)+(0.8*(2*CoordV4)),CoordQ(2)+(0.8*(2*CoordV4)),CoordQ(3)+(0.8*(2*CoordV4)),['r' num2str(huy) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                else
                    text(CoordQ(1)+(0.8*(2*CoordV4)),CoordQ(2)+(0.8*(2*CoordV4)),CoordQ(3)+(0.8*(2*CoordV4)),['r' num2str(huy-sum(nEle(1:networkNumb-1))) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                end
                for hgf = 1:size(CoordList,1)
                    isF = [CoordQ(1,1);CoordList(hgf,1)];
                    jsF = [CoordQ(1,2);CoordList(hgf,2)];
                    ksF = [CoordQ(1,3);CoordList(hgf,3)];
                    if unCorPmap(hgf+huy,huy)>-1.*(str2num(get(handles.edit1,'String')))
                        hSurface = surf(CoordV4.*2.*rws+CoordQ(1,1),CoordV4.*2.*res+CoordQ(1,2),CoordV4.*2.*rdw+CoordQ(1,3));
                        set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0,'LineSmoothing','on');
                        [xre,yre,zre] = tubeplot([isF(1) isF(2);jsF(1) jsF(2);ksF(1) ksF(2)],0.4,10);
                        if g1_mean(hgf+huy,huy)>g2_mean(hgf+huy,huy)
                            ySurface1 = surf(xre,yre,zre);
                            set(ySurface1,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                'gouraud','FaceColor',[0.2 0.5 1],'FaceLighting','gouraud','LineSmoothing','on');
                             ex1_1 = 1;

                        else
                            ySurface2 = surf(xre,yre,zre);
                            set(ySurface2,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                'gouraud','FaceColor',[1 0.3 0.3],'FaceLighting','gouraud','LineSmoothing','on');
                             ex1_2 = 1;
                        end
                    else
                         hSurface = surf(CoordV4.*2.*rws+CoordQ(1,1),CoordV4.*2.*res+CoordQ(1,2),CoordV4.*2.*rdw+CoordQ(1,3));
                         set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
                    end   
                 end
            end

            hSurface = surf(CoordV4.*2.*rws+CoordMS(end,1),CoordV4.*2.*res+CoordMS(end,2),CoordV4.*2.*rdw+CoordMS(end,3));
            set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);


            if exist('ex1_1') && exist('ex1_2')
                    legend([ySurface1,ySurface2],'Group 1 > Group 2', 'Group 2 > Group 1')
               else if exist('ex1_1')
                        legend(ySurface1,'Group 1 > Group 2')
                    else
                        legend(ySurface2,'Group 2 > Group 1')
                    end
            end

            DirUF2C = which('uf2c');
            try
                aRT = nifti([DirUF2C(1:end-6) 'Utilities' filesep 'rControls_Template.nii']);
            catch
                aRT = nifti([DirUF2C(1:end-13) 'Utilities' filesep 'rControls_Template.nii']);
            end
            mTextBox = uicontrol('style','text','position',[20    20    200    20]);
            set(mTextBox,'String','r = Roi number and n = network number');

            matI = aRT.dat(:,:,:);
            matI(isnan(matI))=0;
            mat2 = permute(matI,[1,3,2]);
            mat2 = permute(mat2,[2,1,3]);
            mat2 = permute(mat2,[3,2,1]);

            ater = slice(mat2,46,[],20);

            colormap gray
            set(ater, 'LineStyle','none');
            set(ater, 'FaceAlpha',1);

            aphamap3d = smooth3(mat2,'box',[11 11 11]);    
            aphamap1 = aphamap3d(:,46,:);
            aphamap1 = squeeze(aphamap1);
            aphamap1(aphamap1<300)=0.0001;
            aphamap1(aphamap1>=300)=0.7;
            set(ater(1),'alphadata',aphamap1,'facealpha','flat')

            aphamap3 = aphamap3d(:,:,20);
            aphamap3 = squeeze(aphamap3);
            aphamap3(aphamap3<300)=0.0001;
            aphamap3(aphamap3>=300)=0.7;
            set(ater(2),'alphadata',aphamap3,'facealpha','flat')
            light('Style','infinite');%'Position',[1 0 1]);
            
            saveas(GraphoFig,[Dir filesep 'Uncorrected_Graph_3D'],'fig')
            tmpFig = myaa('publish');
            saveas(tmpFig,[Dir filesep 'Uncorrected_Graph_3D'],'tif')
            close(tmpFig)

            CoLeg = figure('Visible','off');
            set(CoLeg,'Name','Color Legend','Position', [400 400 250 400],...
                            'Color',[1 .949 .867]);
            title('Color Legend','FontSize',14);
            xlim([0 0.6]);
            ylim([0 nOfNets]);

            for uiC = 1:nOfNets
                rectangle('Position',[0,uiC-1,0.6,1],'LineWidth',2,'FaceColor',cc(uiC,:))
                string = sprintf('Netowork %d color',uiC);
                text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
            end
            saveas(CoLeg,[Dir filesep 'Color_Legend'],'tif')
            saveas(CoLeg,[Dir filesep 'Color_Legend'],'fig')
            close(CoLeg)

        end
        if get(handles.Grapho2D,'Value')  %%%%% 2D Graphos
            GraphoFig2 = figure;
            set(GraphoFig2,'Name','Connecticome2D','Position', [130 130 800 800],...
                            'Color',[1 .949 .867]);
            title('Connecticome2D','FontSize',14);
            xlabel('X axis'  ,'FontSize',14);
            ylabel('Y axis','FontSize',14);
            zlabel('Z axis','FontSize',14);
            scatter3(CoordMS(:,1),CoordMS(:,2),CoordMS(:,3),'Marker','square','MarkerFaceColor','red','MarkerEdgeColor','k');
            set(gca, 'color', [1 1 1])
            set(gca, 'gridcolor', [0 0 0])
            set(gca, 'gridalpha', 1)
            hold on
            xlim([0 91]);
            ylim([0 109]);
            zlim([0 91]);
            view([0 0 1])

             stIndx = 1;
             for hg = 1:size(nEle,2)
                 str_nEle(hg) = stIndx;
                 stIndx = stIndx+nEle(hg);
             end

             cc=hsv(nOfNets);
             [rws,res,rdw] = sphere(15);

            for huy = 1:size(CoordMS,1)
                CoordQ2 = CoordMS(huy,:);
                CoordList = CoordMS(huy+1:end,:);
                CoordV4 = sum(unCorPmap(:,huy)~=-1);
                CoordV4 = 1+ (2*CoordV4/size(CoordMS,1));
                tester2 = huy<=(str_nEle+nEle)-1;
                networkNumb2 = find(tester2==1,1);
                colorV2 = cc(networkNumb2,:);
                if isequal(networkNumb2,1)
                    text(CoordQ2(1)+(0.8*(2*CoordV4)),CoordQ2(2)+(0.8*(2*CoordV4)),CoordQ2(3)+(0.8*(2*CoordV4)),['r' num2str(huy) 'n' num2str(networkNumb2)],'HorizontalAlignment','Right','FontSize',9);
                else
                    text(CoordQ2(1)+(0.8*(2*CoordV4)),CoordQ2(2)+(0.8*(2*CoordV4)),CoordQ2(3)+(0.8*(2*CoordV4)),['r' num2str(huy-sum(nEle(1:networkNumb2-1))) 'n' num2str(networkNumb2)],'HorizontalAlignment','Right','FontSize',9);
                end
                for hgf = 1:size(CoordList,1)
                    isF = [CoordQ2(1,1);CoordList(hgf,1)];
                    jsF = [CoordQ2(1,2);CoordList(hgf,2)];
                    ksF = [CoordQ2(1,3);CoordList(hgf,3)];
                    if unCorPmap(hgf+huy,huy)>-1.*(str2num(get(handles.edit1,'String')))
                        hSurface = surf(CoordV4.*rws+CoordQ2(1,1),CoordV4.*res+CoordQ2(1,2),CoordV4.*rdw+CoordQ2(1,3));
                        set(hSurface, 'FaceColor',colorV2, 'FaceAlpha',1, 'EdgeAlpha',0);
                        if g1_mean(hgf+huy,huy)>g2_mean(hgf+huy,huy)
                            pLine1 = line(isF,jsF,ksF,'Color',[0.2 0.5 1],'LineStyle','-','LineWidth',0.5);
                            ex2_1 = 1;
                        else
                            pLine2 = line(isF,jsF,ksF,'Color',[1 0.3 0.3],'LineStyle','-','LineWidth',0.5);
                            ex2_2 = 1;
                        end
                    else
                        hSurface = surf(CoordV4.*rws+CoordQ2(1,1),CoordV4.*res+CoordQ2(1,2),CoordV4.*rdw+CoordQ2(1,3));
                        set(hSurface, 'FaceColor',colorV2, 'FaceAlpha',1, 'EdgeAlpha', 0);
                    end   
                end
            end
            hSurface = surf(CoordV4.*rws+CoordQ2(end,1),CoordV4.*res+CoordQ2(end,2),CoordV4.*rdw+CoordQ2(end,3));
            set(hSurface, 'FaceColor',colorV2, 'FaceAlpha',1, 'EdgeAlpha', 0);

            if exist('ex2_1') && exist('ex2_2')
                    legend([pLine1,pLine2],'Group 1 > Group 2', 'Group 2 > Group 1')
               else if exist('ex2_1')
                        legend(pLine1,'Group 1 > Group 2')
                    else
                        legend(pLine2,'Group 2 > Group 1')
                    end
            end

            mTextBox = uicontrol('style','text','position',[20    20    200    20]);
            set(mTextBox,'String','r = Roi number and n = network number');
            
            DirUF2C = which('uf2c');

            try
                load([DirUF2C(1:end-6) 'Utilities' filesep 'Brain_Countour.mat'])
            catch
                load([DirUF2C(1:end-13) 'Utilities' filesep 'Brain_Countour.mat'])
            end

            aRR = Brain_Contour(:,:,1);
            aRR = double(aRR);
            aRR(aRR<250)=1;
            aRR(aRR>1)=0;
            aRR = flipud(aRR);

            for i = 1:91
                mat3d(:,:,i) = aRR;
            end

            aRR = permute(aRR,[1,3,2]);
            ater2 = slice(mat3d,[],[],1);
            set(ater2, 'FaceAlpha',0.7);
            colormap gray
            set(ater2, 'LineStyle','none');
            matConto = mat3d(:,:,46);
            matConto = squeeze(matConto);
            matConto(matConto<1)=0.0001;
            set(ater2(1),'alphadata',matConto,'facealpha','flat')
            grid off
            light('Style','infinite');
            
            saveas(GraphoFig2,[Dir filesep 'Uncorrected_Graph_2D'],'fig')

            tmpFig = myaa('publish');
            saveas(tmpFig,[Dir filesep 'Uncorrected_Graph_2D'],'tif')
            close(tmpFig)

            CoLeg = figure('Visible','off');
            set(CoLeg,'Name','Color Legend','Position', [400 400 250 400],...
                            'Color',[1 .949 .867]);
            title('Color Legend','FontSize',14);
            xlim([0 0.6]);
            ylim([0 nOfNets]);

            for uiC = 1:nOfNets
                rectangle('Position',[0,uiC-1,0.6,1],'LineWidth',2,'FaceColor',cc(uiC,:))
                string = sprintf('Netowork %d color',uiC);
                text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
            end
            saveas(CoLeg,[Dir filesep 'Color_Legend'],'tif')
            saveas(CoLeg,[Dir filesep 'Color_Legend'],'fig')
            close(CoLeg)

        end
    else
        warndlg('No significant results with the ''uncorrected'' defined threshold ','Ops!')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%  FDR
%%%%%%%%%%%%%%%%%%%%%%%%%%%  FDR

if get(handles.chec2,'Value')  %%%%% FDR corrected Matrix
    FDRcor = Pmap;
    tPmap = tril(FDRcor,-1);
    tPmap = tPmap(tPmap>0);
    [h, crit_p, adj_p] = fdr_bh(tPmap,str2num(get(handles.edit2,'String')),'pdep','yes');
    res = zeros(size(g1_VAR,1),size(g1_VAR,1));
    stt = 1;

    for i = 1:size(g1_VAR,1)
        res(1:i,i) = zeros(i,1);
        res(i+1:end,i) = adj_p(stt:stt+size(g1_VAR,1)-i-1);
        stt = stt+(size(g1_VAR,1)-i);
    end
    resF = res+transpose(res);
    FDRcor = resF;

    FDRcor(FDRcor>(str2num(get(handles.edit2,'String')))) = 1;
    FDRcor(FDRcor==0) = 1;
    FDRcor = -1.*FDRcor;
    
    if ~isequal(numel(unique(FDRcor)),1)
        if get(handles.CorMats,'Value')
            fig2 = figure;
            imagesc(FDRcor)
            title('FDR corrected P Map')
            colorbar
            caxis([-1.*(str2num(get(handles.edit2,'String'))) (-1*(min(adj_p)-(0.2*min(adj_p))))])

            saveas(fig2,[Dir filesep filesep 'FDR_corrected_P_Map'],'tif')
            saveas(fig2,[Dir filesep filesep 'FDR_corrected_P_Map'],'fig')
            save([Dir filesep filesep 'unCorPmap-VAR'],'FDRcor')
            close(fig2)
        end
        if get(handles.Graphos3D,'Value') %%%%% FDR corrected 3D Grapho
            GraphoFig3 = figure;
            set(GraphoFig3,'Name','Connecticome3D','Position', [130 130 1350 800],...
                            'Color',[1 .949 .867]);
            title('Connecticome3D','FontSize',14);
            xlabel('X axis'  ,'FontSize',14);
            ylabel('Y axis','FontSize',14);
            zlabel('Z axis','FontSize',14);
            scatter3(CoordMS(:,1),CoordMS(:,2),CoordMS(:,3),'Marker','square','MarkerFaceColor','red','MarkerEdgeColor','k');
            set(gca, 'color', [1 1 1])
            set(gca, 'gridcolor', [0.8 0.8 0.8])
            set(gca, 'gridalpha', 1)
            xlim([0 91]);
            ylim([0 109]);
            zlim([0 91]);
            hold on

             stIndx = 1;
             for hg = 1:size(nEle,2)
                 str_nEle(hg) = stIndx;
                 stIndx = stIndx+nEle(hg);
             end

             cc=hsv(nOfNets);
             [rws,res,rdw] = sphere(15);

             for huy = 1:size(CoordMS,1)
                CoordQ = CoordMS(huy,:);
                CoordList = CoordMS(huy+1:end,:);
                CoordV4 = sum(FDRcor(:,huy)~=-1);
                CoordV4 = 1+ (2*CoordV4/size(CoordMS,1));
                tester = huy<=(str_nEle+nEle)-1;
                networkNumb = find(tester==1,1);
                colorV = cc(networkNumb,:);
                if isequal(networkNumb,1)
                    text(CoordQ(1)+(0.8*(2*CoordV4)),CoordQ(2)+(0.8*(2*CoordV4)),CoordQ(3)+(0.8*(2*CoordV4)),['r' num2str(huy) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                else
                    text(CoordQ(1)+(0.8*(2*CoordV4)),CoordQ(2)+(0.8*(2*CoordV4)),CoordQ(3)+(0.8*(2*CoordV4)),['r' num2str(huy-sum(nEle(1:networkNumb-1))) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                end
                for hgf = 1:size(CoordList,1)
                    isF = [CoordQ(1,1);CoordList(hgf,1)];
                    jsF = [CoordQ(1,2);CoordList(hgf,2)];
                    ksF = [CoordQ(1,3);CoordList(hgf,3)];
                    if FDRcor(hgf+huy,huy)>-1.*(str2num(get(handles.edit2,'String')))
                        hSurface = surf(CoordV4.*rws+CoordQ(1,1),CoordV4.*res+CoordQ(1,2),CoordV4.*rdw+CoordQ(1,3));
                        set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
                        [xre,yre,zre] = tubeplot([isF(1) isF(2);jsF(1) jsF(2);ksF(1) ksF(2)],0.4,10);
                        if g1_mean(hgf+huy,huy)>g2_mean(hgf+huy,huy)
                            ySurface1 = surf(xre,yre,zre);
                            set(ySurface1,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                'gouraud','FaceColor',[0.2 0.5 1],'FaceLighting','gouraud','LineSmoothing','on');
                            ex3_1=1;
                        else
                            ySurface2 = surf(xre,yre,zre);
                            set(ySurface2,'EdgeColor','none','LineStyle','none','FaceLighting',...
                                'gouraud','FaceColor',[1 0.3 0.3],'FaceLighting','gouraud','LineSmoothing','on');
                            ex3_2=1;
                        end
                    else
                         hSurface = surf(CoordV4.*rws+CoordQ(1,1),CoordV4.*res+CoordQ(1,2),CoordV4.*rdw+CoordQ(1,3));
                         set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
                    end   
                end
             end

            hSurface = surf(CoordV4.*rws+CoordMS(end,1),CoordV4.*res+CoordMS(end,2),CoordV4.*rdw+CoordMS(end,3));
            set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
            
            if exist('ex3_1') && exist('ex3_2')
                legend([ySurface1,ySurface2],'Group 1 > Group 2', 'Group 2 > Group 1')
            else if exist('ex3_1')
                    legend(ySurface1,'Group 1 > Group 2')
                else
                    legend(ySurface2,'Group 2 > Group 1')
                end
            end
            
            mTextBox = uicontrol('style','text','position',[20    20    200    20]);
            set(mTextBox,'String','r = Roi number and n = network number');
            
            DirUF2C = which('uf2c');
            try
                aRT = nifti([DirUF2C(1:end-6) 'Utilities' filesep 'rControls_Template.nii']);
            catch
                aRT = nifti([DirUF2C(1:end-13) 'Utilities' filesep 'rControls_Template.nii']);
            end

            matI = aRT.dat(:,:,:);
            matI(isnan(matI))=0;
            mat2 = permute(matI,[1,3,2]);
            mat2 = permute(mat2,[2,1,3]);
            mat2 = permute(mat2,[3,2,1]);

            ater = slice(mat2,46,[],20);

            colormap gray
            set(ater, 'LineStyle','none');
            set(ater, 'FaceAlpha',1);

            aphamap3d = smooth3(mat2,'box',[11 11 11]);    
            aphamap1 = aphamap3d(:,46,:);
            aphamap1 = squeeze(aphamap1);
            aphamap1(aphamap1<300)=0.0001;
            aphamap1(aphamap1>=300)=0.6;
            set(ater(1),'alphadata',aphamap1,'facealpha','flat')

            aphamap3 = aphamap3d(:,:,20);
            aphamap3 = squeeze(aphamap3);
            aphamap3(aphamap3<300)=0.0001;
            aphamap3(aphamap3>=300)=0.6;
            set(ater(2),'alphadata',aphamap3,'facealpha','flat')
            light('Style','infinite');
            saveas(GraphoFig3,[Dir filesep 'FDR_Corrected_Graph_3D'],'fig')
            
            tmpFig = myaa('publish');
            saveas(tmpFig,[Dir filesep 'FDR_Corrected_Graph_3D'],'tiff')
            close(tmpFig)
            
            CoLeg = figure('Visible','off');
            set(CoLeg,'Name','Color Legend','Position', [400 400 250 400],...
                            'Color',[1 .949 .867]);
            title('Color Legend','FontSize',14);
            xlim([0 0.6]);
            ylim([0 nOfNets]);

            for uiC = 1:nOfNets
                rectangle('Position',[0,uiC-1,0.6,1],'LineWidth',2,'FaceColor',cc(uiC,:))
                string = sprintf('Netowork %d color',uiC);
                text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
            end
            saveas(CoLeg,[Dir filesep 'Color_Legend'],'tif')
            saveas(CoLeg,[Dir filesep 'Color_Legend'],'fig')
            close(CoLeg)
            
        end
        
        if get(handles.Grapho2D,'Value')  %%%%% FDR corrected 2D Grapho
            GraphoFig4 = figure;
            set(GraphoFig4,'Name','Connecticome2D','Position', [130 130 800 800],...
                            'Color',[1 .949 .867]);
            title('Connecticome2D','FontSize',14);
            xlabel('X axis'  ,'FontSize',14);
            ylabel('Y axis','FontSize',14);
            zlabel('Z axis','FontSize',14);
            scatter3(CoordMS(:,1),CoordMS(:,2),CoordMS(:,3),'Marker','square','MarkerFaceColor','red','MarkerEdgeColor','k');
            set(gca, 'color', [1 1 1])
            set(gca, 'gridcolor', [0 0 0])
            set(gca, 'gridalpha', 1)
            hold on
            xlim([0 91]);
            ylim([0 109]);
            zlim([0 91]);
            view([0 0 1])
            
             stIndx = 1;
             for hg = 1:size(nEle,2)
                 str_nEle(hg) = stIndx;
                 stIndx = stIndx+nEle(hg);
             end

            cc=hsv(nOfNets);
            [rws,res,rdw] = sphere(15);
            
            for huy = 1:size(CoordMS,1)
                CoordQ2 = CoordMS(huy,:);
                CoordList = CoordMS(huy+1:end,:);
                CoordV4 = sum(FDRcor(:,huy)~=-1);
                CoordV4 = 1+ (2*CoordV4/size(CoordMS,1));
                tester2 = huy<=(str_nEle+nEle)-1;
                networkNumb2 = find(tester2==1,1);
                colorV2 = cc(networkNumb2,:);
                if isequal(networkNumb2,1)
                    text(CoordQ2(1)+(0.8*(2*CoordV4)),CoordQ2(2)+(0.8*(2*CoordV4)),CoordQ2(3)+(0.8*(2*CoordV4)),['r' num2str(huy) 'n' num2str(networkNumb2)],'HorizontalAlignment','Right','FontSize',9);
                else
                    text(CoordQ2(1)+(0.8*(2*CoordV4)),CoordQ2(2)+(0.8*(2*CoordV4)),CoordQ2(3)+(0.8*(2*CoordV4)),['r' num2str(huy-sum(nEle(1:networkNumb2-1))) 'n' num2str(networkNumb2)],'HorizontalAlignment','Right','FontSize',9);
                end
                for hgf = 1:size(CoordList,1)
                    isF = [CoordQ2(1,1);CoordList(hgf,1)];
                    jsF = [CoordQ2(1,2);CoordList(hgf,2)];
                    ksF = [CoordQ2(1,3);CoordList(hgf,3)];
                    if FDRcor(hgf+huy,huy)>-1.*(str2num(get(handles.edit2,'String')))
                        hSurface = surf(CoordV4.*rws+CoordQ2(1,1),CoordV4.*res+CoordQ2(1,2),CoordV4.*rdw+CoordQ2(1,3));
                        set(hSurface, 'FaceColor',colorV2, 'FaceAlpha',1, 'EdgeAlpha',0);
                        if g1_mean(hgf+huy,huy)>g2_mean(hgf+huy,huy)
                           pLine1 =  line(isF,jsF,ksF,'Color',[0.2 0.5 1],'LineStyle','-','LineWidth',0.5);
                           ex4_1 = 1;
                        else
                           pLine2 = line(isF,jsF,ksF,'Color',[1 0.3 0.3],'LineStyle','-','LineWidth',0.5);
                           ex4_2 = 1;
                        end
                    else
                        hSurface = surf(CoordV4.*rws+CoordQ2(1,1),CoordV4.*res+CoordQ2(1,2),CoordV4.*rdw+CoordQ2(1,3));
                        set(hSurface, 'FaceColor',colorV2, 'FaceAlpha',1, 'EdgeAlpha', 0);
                    end   
                end
            end
            hSurface = surf(CoordV4.*rws+CoordQ2(end,1),CoordV4.*res+CoordQ2(end,2),CoordV4.*rdw+CoordQ2(end,3));
            set(hSurface, 'FaceColor',colorV2, 'FaceAlpha',1, 'EdgeAlpha', 0);
            
            if exist('ex4_1') && exist('ex4_2')
                    legend([pLine1,pLine2],'Group 1 > Group 2', 'Group 2 > Group 1')
               else if exist('ex4_1')
                        legend(pLine1,'Group 1 > Group 2')
                    else
                        legend(pLine2,'Group 2 > Group 1')
                    end
            end

            mTextBox = uicontrol('style','text','position',[20    20    200    20]);
            set(mTextBox,'String','r = Roi number and n = network number');
            
            DirUF2C = which('uf2c');
            try
                load([DirUF2C(1:end-6) 'Utilities' filesep 'Brain_Countour.mat'])
            catch
                load([DirUF2C(1:end-13) 'Utilities' filesep 'Brain_Countour.mat'])
            end

            aRR = Brain_Contour(:,:,1);
            aRR = double(aRR);
            aRR(aRR<250)=1;
            aRR(aRR>1)=0;
            aRR = flipud(aRR);

            for i = 1:91
                mat3d(:,:,i) = aRR;
            end

            aRR = permute(aRR,[1,3,2]);
            ater2 = slice(mat3d,[],[],1);
            set(ater2, 'FaceAlpha',0.7);
            colormap gray
            set(ater2, 'LineStyle','none');
            matConto = mat3d(:,:,46);
            matConto = squeeze(matConto);
            matConto(matConto<1)=0.0001;
            set(ater2(1),'alphadata',matConto,'facealpha','flat')
            grid off
            light('Style','infinite');
            
            saveas(GraphoFig4,[Dir filesep 'FDR_Corrected_Graph_2D'],'fig')
            tmpFig = myaa('publish');
            saveas(tmpFig,[Dir filesep 'FDR_Corrected_Graph_2D'],'tif')
            close(tmpFig)
            
            CoLeg = figure('Visible','off');
            set(CoLeg,'Name','Color Legend','Position', [400 400 250 400],...
                            'Color',[1 .949 .867]);
            title('Color Legend','FontSize',14);
            xlim([0 0.6]);
            ylim([0 nOfNets]);

            for uiC = 1:nOfNets
                rectangle('Position',[0,uiC-1,0.6,1],'LineWidth',2,'FaceColor',cc(uiC,:))
                string = sprintf('Netowork %d color',uiC);
                text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
            end
            saveas(CoLeg,[Dir filesep 'Color_Legend'],'tif')
            saveas(CoLeg,[Dir filesep 'Color_Legend'],'fig')
            close(CoLeg)

            
        end
    else
        warndlg('No significant results with FDR correction','Ops!')
    end
end
set(handles.text13,'String','Done!')

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

function chec1_Callback(hObject, eventdata, handles)
function chec2_Callback(hObject, eventdata, handles)
function CorMats_Callback(hObject, eventdata, handles)
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

function edit1_Callback(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function figure1_CreateFcn(hObject, eventdata, handles)
global g2OK g1OK dirOK nodeOK
g2OK=0; g1OK=0; dirOK=0; nodeOK=0;


function addROIwCov_Callback(hObject, eventdata, handles)
if isequal(get(handles.addROIwCov,'Value'),1)
    set(handles.G1covarBt,'Enable','on')
    set(handles.G2covarBt,'Enable','on')
else
    set(handles.G1covarBt,'Enable','off')
    set(handles.G2covarBt,'Enable','off')
end


function G1covarBt_Callback(hObject, eventdata, handles)
global g1_VAR g2OK g1OK dirOK nodeOK G1covv

[G1cov,pathG1cov] = uigetfile({'*.mat','Matlab file'},...
        'Select a .mat file for covariate','MultiSelect','off');
load([pathG1cov G1cov]);

G1covv = Total_DifGM_Matrix;
clear Total_DifGM_Matrix
if isequal(size(G1covv),size(g1_VAR))
    set(handles.G1covtxt,'String','Cov. added!')
else
    warndlg('The sizes of the covariate column should to match with the Group 1 file sizes','Ops!')
    set(handles.G1covtxt,'String','Empty')
end



function G2covarBt_Callback(hObject, eventdata, handles)
global g2_VAR g2OK g1OK dirOK nodeOK G2covv

[G2cov,pathG2cov] = uigetfile({'*.mat','Matlab file'},...
        'Select a .mat file for covariate','MultiSelect','off');
load([pathG2cov G2cov]);

G2covv = Total_DifGM_Matrix;
clear Total_DifGM_Matrix
if isequal(size(G2covv),size(g2_VAR))
    set(handles.G2covtxt,'String','Cov. added!')
else
    warndlg('The sizes of the covariate column should to match with the Group 2 file sizes','Ops!')
    set(handles.G2covtxt,'String','Empty')
    clear G2cov
end

function checkbox12_Callback(hObject, eventdata, handles)

function checkbox13_Callback(hObject, eventdata, handles)

function edit8_Callback(hObject, eventdata, handles)

function edit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit9_Callback(hObject, eventdata, handles)

function edit9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
