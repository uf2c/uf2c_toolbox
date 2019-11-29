function varargout = CrossCorSecLe_Ztest(varargin)
% UF²C M-file for CrossCorSecLe_Ztest.fig
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
                   'gui_OpeningFcn', @CrossCorSecLe_Ztest_OpeningFcn, ...
                   'gui_OutputFcn',  @CrossCorSecLe_Ztest_OutputFcn, ...
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

function CrossCorSecLe_Ztest_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function varargout = CrossCorSecLe_Ztest_OutputFcn(hObject, eventdata, handles) 
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
    drawnow
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
  a = 1;
  
function runB_Callback(hObject, eventdata, handles)
global CoordMS seedList Dir COVDATA g1_VAR g2_VAR nOfLines nEle nOfNets
set(handles.text13,'String','Running.....')
drawnow
%%%%%% Checking...
if ~isequal(size(CoordMS,1),size(g1_VAR,1))
    if logical(get(handles.Grapho2D,'Value')) || logical(get(handles.Graphos3D,'Value'))
        warndlg(sprintf('The number of nodes added (%d) should to match with the used on first level (%d)',size(CoordMS,1),size(g1_VAR,1)),'Ops! Process Aborted')
    return
    end
end
%%%%%%%%%%%%%%%%
mkdir(Dir)
for Sjo = 1:size(g2_VAR,3)
    Pmap = zeros(size(g1_VAR,1),size(g1_VAR,1));
    tic
    for i = 1:size(g1_VAR,1)
        for j = 1:size(g1_VAR,1)
            if i==j
                Pmap(i,j) = 1;
            else if sum(squeeze(g1_VAR(i,j,:)))==0 || g2_VAR(i,j,Sjo)==0
                   Pmap(i,j) = 1;
                else
                    T = (g2_VAR(i,j,Sjo)-mean(squeeze(g1_VAR(i,j,:))))/(std(squeeze(g1_VAR(i,j,:)))/sqrt(numel(squeeze(g1_VAR(i,j,:)))));
                    if T>25 || T<-25% At this value the P explode to zero
                        Pmap(i,j) = 5.5511e-17;
                    else
                        Pmap(i,j) = tcdf(-1*abs(T/2),size(g1_VAR,3));
                    end
                 end
            end
        end
    end
    
    g1_mean = mean(g1_VAR,3);
    ScreSize = get(0,'screensize');
    ScreSize = ScreSize(3:end);
    
%     Pmap = Pmap.*435;

    % Uncorrected P Map
    if get(handles.chec1,'Value')
       
        unCorPmap = Pmap;
        unCorPmap(unCorPmap>str2num(get(handles.edit1,'String'))) = 1;
        unCorPmap = -1.*unCorPmap;

        if ~isequal(numel(unique(unCorPmap)),1)
            if get(handles.CorMats,'Value')   %%%%% Correlation Matrix
                fig1 = figure;
                imagesc(unCorPmap)

                set(fig1,'Name','Uncorrected Correlation Matrix','Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.45 ScreSize(1)*.4]),...
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
                imgRR = getframe(gcf);
                imwrite(imgRR.cdata, [Dir filesep 'Uncorrected_R_Map_Subj_' num2str(Sjo)  '.png']);
                clear imgRR
                saveas(fig1,[Dir filesep 'Uncorrected_R_Map_Subj_' num2str(Sjo)],'fig')
                save([Dir filesep 'unCorRmap-VAR_Subj_' num2str(Sjo)],'unCorPmap')
                close(fig1)
            end

        matrixR = unCorPmap;
        matrixR(matrixR==-1) = 0;
        matrixR = abs(matrixR);

        circ = figure('Visible','on');

        set(circ,'Name','Connectome2D',...
            'Position',round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.45 ScreSize(1)*.4]),...
                  'Color',[1 .949 .867]);

        title('Uncorrected Circular Connectome','FontSize',18);
        xlabel('X axis'  ,'FontSize',14);
        ylabel('Y axis','FontSize',14);
        set(gca, 'color', [0 0 0])
        try set(gca, 'gridcolor', [1 1 1]); end
        try set(gca, 'gridalpha', 1); end
        hold on
        limiteG = (size(g1_VAR,1)/2).*1.1;

        xlim([0 limiteG]);
        ylim([0 limiteG]);

        daspect([1 1 1])

        ang=0:0.01:2*pi; 
        xp=(limiteG/2.2)*cos(ang);
        yp=(limiteG/2.2)*sin(ang);
        xpt=(limiteG/2.1)*cos(ang);
        ypt=(limiteG/2.1)*sin(ang);

        plot((limiteG/2)+xp,(limiteG/2)+yp,'LineWidth',1);

        AnoRa = size(xp,2)/size(matrixR,1);
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


        stIndx = 1;

         for hg = 1:size(nEle,2)
             str_nEle(hg) = stIndx;
             stIndx = stIndx+nEle(hg);
         end

         cc=hsv(nOfNets);

        for i = 1:size(matrixR,1)

            tester = i<=(str_nEle+nEle)-1;
            networkNumb = find(tester==1,1);
            colorV = cc(networkNumb,:);

            if i<10
                str = ['0' num2str(i)];
            else
                str = num2str(i);
            end

            text(xp2t(i),yp2t(i),str,'FontSize',9,'Color',colorV,'HorizontalAlignment','center','VerticalAlignment','middle');

            vet = matrixR(:,i);
            vet = find(vet);
            for j = 1:size(vet)
                line([xp2(i);xp2(vet(j))],[yp2(i);yp2(vet(j))],'Color',[0.2 0.5 1],'LineStyle','-','LineWidth',1);
            end
        end

        set(gcf, 'InvertHardCopy', 'off');
    %     saveas(circ,[Dir filesep 'Uncorrected_Circ_Map'],'png')
        saveas(circ,[Dir filesep 'Uncorrected_Circ_Map_Subj_' num2str(Sjo)],'fig')
        drawnow
        imgRR = getframe(gcf);
        imwrite(imgRR.cdata, [Dir filesep 'Uncorrected_Circ_Map_Subj_' num2str(Sjo) '.png']);
        clear imgRR
        close(circ)

            if get(handles.Graphos3D,'Value')  %%%%% 3D Graphos
                GraphoFig = figure;

                set(GraphoFig,'Name','3DConnectome',...
                    'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.6 ScreSize(1)*.4]),...
                                'Color',[1 .949 .867]);

                title('Uncorrected 3D Connectome','FontSize',14);
                xlabel('X axis'  ,'FontSize',14);
                ylabel('Y axis','FontSize',14);
                zlabel('Z axis','FontSize',14);
    %             scatter3(CoordMS(:,1),CoordMS(:,2),CoordMS(:,3),'Marker','square','MarkerFaceColor','red','MarkerEdgeColor','k');
                set(gca, 'color', [1 1 1])
                try set(gca, 'gridcolor', [0.8 0.8 0.8]); end

                try set(gca, 'gridalpha', 1); end

                xlim([0 91]);
                ylim([0 109]);
                zlim([0 91]);
                view([-1 -1 1])
                hold on

               stIndx = 1;
               for hg = 1:size(nEle,2)
                     str_nEle(hg) = stIndx;
                     stIndx = stIndx+nEle(hg);
               end

               cc=hsv(nOfNets);
               [rws,res,rdw] = sphere(10);
               
               UnCorrLat_Right(Sjo) = 0;
               UnCorrLat_Left(Sjo) = 0;
               
               for huy = 1:size(CoordMS,1)
                    CoordQ = CoordMS(huy,:);
                    CoordV4 = (sum(abs(unCorPmap(:,huy)))-1)/(size(CoordMS,1)-1);
                    CoordList = CoordMS(huy+1:end,:);
                    for hgf = 1:size(CoordList,1)
                        isF = [CoordQ(1,1);CoordList(hgf,1)];
                        jsF = [CoordQ(1,2);CoordList(hgf,2)];
                        ksF = [CoordQ(1,3);CoordList(hgf,3)];
                        if unCorPmap(hgf+huy,huy)>-1.*(str2num(get(handles.edit1,'String')))
                            [xre,yre,zre] = tubeplot([isF(1) isF(2);jsF(1) jsF(2);ksF(1) ksF(2)],0.4,8);
                            if CoordQ(1)>=47
                                if CoordList(hgf,1)>=47
                                    UnCorrLat_Right(Sjo) = UnCorrLat_Right(Sjo) + 1;
                                end
                            end
                            if CoordQ(1)<=45
                                if CoordList(hgf,1)<=45
                                    UnCorrLat_Left(Sjo) = UnCorrLat_Left(Sjo) + 1;
                                end
                            end

                            if (g1_mean(hgf+huy,huy)*g2_VAR(hgf+huy,huy,Sjo))>0;
                                if abs(g1_mean(hgf+huy,huy))>abs(g2_VAR(hgf+huy,huy,Sjo))
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

                %%%%%%%%%%%%%%%%%%%%%%%  added to just include spheres with
                %%%%%%%%%%%%%%%%%%%%%%%  alterations
                teste = unCorPmap;
                teste(teste==-1)=0;
                teste(teste~=0)=1;
                vetTTT = sum(teste,1);
                VetFF = vetTTT;
                VetFF(VetFF>1) = 1;
                binVet = (VetFF~=0);
                idxROIs = find(binVet);

                for i = 1:sum(VetFF~=0)
                    tester = idxROIs(i)<=(str_nEle+nEle)-1;
                    networkNumb = find(tester==1,1);
                    colorV = cc(networkNumb,:);
                    CoordQ = CoordMS(idxROIs(i),:);
                    CoordV4 = (vetTTT(idxROIs(i))/max(vetTTT))+1;
                    rNumb = find(VetFF~=0,i);
                    rNumb = rNumb(end);
                    
                    if isequal(networkNumb,1)
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(3*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),['r' num2str(rNumb) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                    else
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(5*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),['r' num2str(rNumb-sum(nEle(1:networkNumb-1))) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                    end
                    hSurface = surf(CoordV4.*rws+CoordQ(1,1),CoordV4.*res+CoordQ(1,2),CoordV4.*rdw+CoordQ(1,3));
                    set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                 if exist('ex1_1') || exist('ex1_2')
                    if exist('ex1_1') && exist('ex1_2')
                        if exist('ex1_3')
                            legend([ySurface1,ySurface2,ySurface3],'Group 1 > Group 2','Group 2 > Group 1',['Group2 ',char(8596),' -(Group1)'])
                        else
                            legend([ySurface1,ySurface2],'Group 1 > Group 2','Group 2 > Group 1')
                        end
                    else if exist('ex1_1')
                            if exist('ex1_3')
                                legend([ySurface1,ySurface3],'Group 1 > Group 2',['Group2 ',char(8596),' -(Group1)'])
                            else
                                legend(ySurface1,'Group 1 > Group 2')
                            end
                        else
                             if exist('ex1_3')
                                 legend([ySurface2,ySurface3],'Group 2 > Group 1',['Group2 ',char(8596),' -(Group1)'])
                             else
                                legend(ySurface2,'Group 2 > Group 1')
                             end
                         end
                    end
                else if exist('ex1_3')
                         legend(ySurface3,['Group2 ',char(8596),' -(Group1)']) % Simbolos apenas copiando e colando...
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
                light('Style','infinite');%'Position',[1 0 1]);

                saveas(GraphoFig,[Dir filesep 'Uncorrected_Graph_3D_Subj_' num2str(Sjo)],'fig')
    %             saveas(GraphoFig,[Dir filesep 'Uncorrected_Graph_3D'],'png')
                drawnow
                imgRR = getframe(gcf);
                imwrite(imgRR.cdata, [Dir filesep 'Uncorrected_Graph_3D_Subj_' num2str(Sjo) '.png']);
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
                    string = sprintf('Netowork %d color',uiC);
                    text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
                end

                set(gca,'xtick',[])
                set(gca,'xticklabel',[])
                set(gca,'ytick',[])
                set(gca,'yticklabel',[])

    %             saveas(CoLeg,[Dir filesep 'Color_Legend'],'png')
                drawnow
                imgRR = getframe(gcf);
                imwrite(imgRR.cdata, [Dir filesep 'Color_Legend' '.png']);

                saveas(CoLeg,[Dir filesep 'Color_Legend'],'fig')
                close(CoLeg)

            end
            if get(handles.Grapho2D,'Value')  %%%%% 2D Graphos
                CoordMS2 = 3.*CoordMS;
                GraphoFig2 = figure;

                set(GraphoFig2,'Name','Connectome2D',...
                    'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.45 ScreSize(1)*.45]),...
                                'Color',[1 .949 .867]);

                title('Uncorrected 2D Connectome','FontSize',14);
                xlabel('X axis'  ,'FontSize',14);
                ylabel('Y axis','FontSize',14);
                zlabel('Z axis','FontSize',14);
    %             scatter3(CoordMS(:,1),CoordMS(:,2),CoordMS(:,3),'Marker','square','MarkerFaceColor','red','MarkerEdgeColor','k');
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
                 [rws,res,rdw] = sphere(15);


                for huy = 1:size(CoordMS2,1)
                    CoordQ2 = CoordMS2(huy,:);
                    CoordList = CoordMS2(huy+1:end,:);
                    tester2 = huy<=(str_nEle+nEle)-1;
                    networkNumb2 = find(tester2==1,1);
                    colorV2 = cc(networkNumb2,:);
                    for hgf = 1:size(CoordList,1)
                        isF = [CoordQ2(1,1);CoordList(hgf,1)];
                        jsF = [CoordQ2(1,2);CoordList(hgf,2)];
                        ksF = [CoordQ2(1,3);CoordList(hgf,3)];
                        if unCorPmap(hgf+huy,huy)>-1.*(str2num(get(handles.edit1,'String')))
                            if (g1_mean(hgf+huy,huy)*g2_VAR(hgf+huy,huy,Sjo))>0;
                                if abs(g1_mean(hgf+huy,huy))>abs(g2_VAR(hgf+huy,huy,Sjo))
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

                %%%%%%%%%%%%%%%%%%%%%%%  added to just include spheres with
                %%%%%%%%%%%%%%%%%%%%%%%  alterations
                teste = unCorPmap;
                teste(teste==-1)=0;
                teste(teste~=0)=1;
                vetTTT = sum(teste,1);
                VetFF = vetTTT;
                VetFF(VetFF>1) = 1;
                binVet = (VetFF~=0);
                idxROIs = find(binVet);

                for i = 1:sum(VetFF~=0)
                    tester = idxROIs(i)<=(str_nEle+nEle)-1;
                    networkNumb = find(tester==1,1);
                    colorV = cc(networkNumb,:);
                    CoordQ = CoordMS2(idxROIs(i),:);
                    CoordV4 = (vetTTT(idxROIs(i))/max(vetTTT))+1;
                    rNumb = find(VetFF~=0,i);
                    rNumb = rNumb(end);

                    if isequal(networkNumb,1)
                        text(CoordQ(1)+(0.8*(5*CoordV4)),CoordQ(2)+(0.8*(5*CoordV4)),CoordQ(3)+(0.8*(5*CoordV4)),['r' num2str(rNumb) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                    else
                        text(CoordQ(1)+(0.8*(5*CoordV4)),CoordQ(2)+(0.8*(5*CoordV4)),CoordQ(3)+(0.8*(5*CoordV4)),['r' num2str(rNumb-sum(nEle(1:networkNumb-1))) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                    end
                    hSurface = surf(CoordV4*3.*rws+CoordQ(1,1),CoordV4*3.*res+CoordQ(1,2),CoordV4*3.*rdw+CoordQ(1,3));
                    set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if exist('ex2_1') || exist('ex2_2')
                    if exist('ex2_1') && exist('ex2_2')
                        if exist('ex2_3')
                            legend([pLine1,pLine2,pLine3],'Group 1 > Group 2','Group 2 > Group 1',['Group2 ',char(8596),' -(Group1)'])
                        else
                            legend([pLine1,pLine2],'Group 1 > Group 2','Group 2 > Group 1')
                        end
                    else if exist('ex2_1')
                            if exist('ex2_3')
                                legend([pLine1,pLine3],'Group 1 > Group 2',['Group2 ',char(8596),' -(Group1)'])
                            else
                                legend(pLine1,'Group 1 > Group 2')
                            end
                         else
                            if exist('ex2_3')
                                legend([pLine2,pLine3],'Group 1 > Group 2',['Group2 ',char(8596),' -(Group1)'])
                            else
                                legend(pLine2,'Group 1 > Group 2')
                            end
                         end
                    end
                else if exist('ex2_3')
                         legend(pLine3,'Group2 = -(Group1)')
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
                light('Style','infinite');

                set(gcf, 'InvertHardCopy', 'off');
                saveas(GraphoFig2,[Dir filesep 'Uncorrected_Graph_2D_Subj_' num2str(Sjo)],'fig')
    %             saveas(GraphoFig2,[Dir filesep 'Uncorrected_Graph_2D'],'png')
                drawnow
                imgRR = getframe(gcf);
                imwrite(imgRR.cdata, [Dir filesep 'Uncorrected_Graph_2D_Subj_' num2str(Sjo) '.png']);
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
                    string = sprintf('Netowork %d color',uiC);
                    text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
                end
                set(gcf, 'InvertHardCopy', 'off');
                set(gca,'xtick',[])
                set(gca,'xticklabel',[])
                set(gca,'ytick',[])
                set(gca,'yticklabel',[])


    %             saveas(CoLeg,[Dir filesep 'Color_Legend'],'png')
                drawnow
                imgRR = getframe(gcf);
                imwrite(imgRR.cdata, [Dir filesep 'Color_Legend' '.png']);
                clear imgRR

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

                set(fig2,'Name','FDR-Corrected-Correlation-Matrix','Position',...
                    round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.45 ScreSize(1)*.4]),...
                        'Color',[1 .949 .867]);

                title('FDR Corrected Correlation Matrix','FontSize',14);
                xlabel('Seeds'  ,'FontSize',14);
                ylabel('Seeds','FontSize',14);
                set(gca, 'color', [0 0 0])
                try set(gca, 'gridcolor', [1 1 1]); end
                try set(gca, 'gridalpha', 1); end
                colorbar
                daspect([1 1 1])

                caxis([-1.*(str2num(get(handles.edit2,'String'))) (-1*(min(adj_p)-(0.2*min(adj_p))))])
                set(gcf, 'InvertHardCopy', 'off');
    %             saveas(fig2,[Dir filesep filesep 'FDR_corrected_P_Map'],'png')
                drawnow
                imgRR = getframe(gcf);
                imwrite(imgRR.cdata, [Dir filesep 'FDR_corrected_R_Map_Subj_' num2str(Sjo) '.png']);
                clear imgRR

                saveas(fig2,[Dir filesep filesep 'FDR_corrected_R_Map_Subj_' num2str(Sjo)],'fig')
                save([Dir filesep filesep 'FDRcorRmap-VAR_Subj_' num2str(Sjo)],'FDRcor')
                close(fig2)

                matrixR = FDRcor;
                matrixR(matrixR==-1) = 0;
                matrixR = abs(matrixR);

                circ = figure('Visible','on');
                set(circ,'Name','Connectome2D',...
                    'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.45 ScreSize(1)*.4]),...
                        'Color',[1 .949 .867]);

                title('FDR Corrected Circular Connectome','FontSize',14);
                xlabel('X axis'  ,'FontSize',14);
                ylabel('Y axis','FontSize',14);
                set(gca, 'color', [0 0 0])
                try set(gca, 'gridcolor', [1 1 1]); end
                try set(gca, 'gridalpha', 1); end
                hold on

                limiteG = (size(g1_VAR,1)/2).*1.1;
                xlim([0 limiteG]);
                ylim([0 limiteG]);

                ang=0:0.01:2*pi; 
                xp=(limiteG/2.2)*cos(ang);
                yp=(limiteG/2.2)*sin(ang);
                xpt=(limiteG/2.1)*cos(ang);
                ypt=(limiteG/2.1)*sin(ang);

                plot((limiteG/2)+xp,(limiteG/2)+yp,'LineWidth',1);

                AnoRa = round(size(xp,2)/size(matrixR,1));
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
                stIndx = 1;

                 for hg = 1:size(nEle,2)
                     str_nEle(hg) = stIndx;
                     stIndx = stIndx+nEle(hg);
                 end

                 cc=hsv(nOfNets);

                for i = 1:size(matrixR,1)

                    tester = i<=(str_nEle+nEle)-1;
                    networkNumb = find(tester==1,1);
                    colorV = cc(networkNumb,:);

                    if i<10
                        str = ['0' num2str(i)];
                    else
                        str = num2str(i);
                    end

                    text(xp2t(i),yp2t(i),str,'FontSize',9,'Color',colorV,'HorizontalAlignment','center','VerticalAlignment','middle');

                    vet = matrixR(:,i);
                    vet = find(vet);
                    for j = 1:size(vet)
                        line([xp2(i);xp2(vet(j))],[yp2(i);yp2(vet(j))],'Color',[0.2 0.5 1],'LineStyle','-','LineWidth',1);
                    end
                end
                set(gcf,'InvertHardCopy', 'off');
    %             saveas(circ,[Dir filesep 'FDRcorrected_Circ_Map'],'png')
                drawnow
                imgRR = getframe(gcf);
                imwrite(imgRR.cdata, [Dir filesep 'FDRcorrected_Circ_Map_Subj_' num2str(Sjo) '.png']);
                clear imgRR

                saveas(circ,[Dir filesep 'FDRcorrected_Circ_Map_Subj_' num2str(Sjo)],'fig')
                close(circ)

            end

            if get(handles.Graphos3D,'Value') %%%%% FDR corrected 3D Grapho
                GraphoFig3 = figure;
                set(GraphoFig3,'Name','Connectome3D',...
                    'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.6 ScreSize(1)*.4]),...
                                'Color',[1 .949 .867]);

                title('FDR Corrected 3D Connectome','FontSize',14);
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
                view([-1 -1 1])

                 stIndx = 1;
                 for hg = 1:size(nEle,2)
                     str_nEle(hg) = stIndx;
                     stIndx = stIndx+nEle(hg);
                 end
                FDRCorrLat_Right(Sjo) = 0;
                FDRCorrLat_Left(Sjo) = 0;
                
                 cc=hsv(nOfNets);
                 [rws,res,rdw] = sphere(15);
                 for huy = 1:size(CoordMS,1)
                    CoordQ = CoordMS(huy,:);
                    CoordList = CoordMS(huy+1:end,:);
                    CoordV4 = sum(FDRcor(:,huy)~=-1);
                    CoordV4 = 1+ (2*CoordV4/size(CoordMS,1));
                    for hgf = 1:size(CoordList,1)
                        isF = [CoordQ(1,1);CoordList(hgf,1)];
                        jsF = [CoordQ(1,2);CoordList(hgf,2)];
                        ksF = [CoordQ(1,3);CoordList(hgf,3)];
                        if FDRcor(hgf+huy,huy)>-1.*(str2num(get(handles.edit2,'String')))
                            [xre,yre,zre] = tubeplot([isF(1) isF(2);jsF(1) jsF(2);ksF(1) ksF(2)],0.4,10);
                            
                            if CoordQ(1)>=46
                                if CoordList(hgf,1)>=46
                                    FDRCorrLat_Right(Sjo) = FDRCorrLat_Right(Sjo) + 1;
                                end
                            end
                            if CoordQ(1)<=46
                                if CoordList(hgf,1)<=46
                                    FDRCorrLat_Left(Sjo) = FDRCorrLat_Left(Sjo) + 1;
                                end
                            end
                            
                            if (g1_mean(hgf+huy,huy)*g2_VAR(hgf+huy,huy,Sjo))>0;
                                if abs(g1_mean(hgf+huy,huy))>abs(g2_VAR(hgf+huy,huy,Sjo))
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

                teste = FDRcor;
                teste(teste==-1)=0;
                teste(teste~=0)=1;
                vetTTT = sum(teste,1);
                VetFF = vetTTT;
                VetFF(VetFF>1) = 1;
                binVet = (VetFF~=0);
                idxROIs = find(binVet);

                for i = 1:sum(VetFF~=0)

                    tester = idxROIs(i)<=(str_nEle+nEle)-1;
                    networkNumb = find(tester==1,1);
                    colorV = cc(networkNumb,:);
                    CoordQ = CoordMS(idxROIs(i),:);
                    CoordV4 = (vetTTT(idxROIs(i))/max(vetTTT))+1;
                    rNumb = find(VetFF~=0,i);
                    rNumb = rNumb(end);

                    if isequal(networkNumb,1)
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(3*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),['r' num2str(rNumb) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                    else
                        text(CoordQ(1)+(0.8*(3*CoordV4)),CoordQ(2)+(0.8*(3*CoordV4)),CoordQ(3)+(0.8*(3*CoordV4)),['r' num2str(rNumb-sum(nEle(1:networkNumb-1))) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',9);
                    end
                    hSurface = surf(CoordV4.*rws+CoordQ(1,1),CoordV4.*res+CoordQ(1,2),CoordV4.*rdw+CoordQ(1,3));
                    set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
                end


                 if exist('ex3_1') || exist('ex3_2')
                    if exist('ex3_1') && exist('ex3_2')
                        if exist('ex3_3')
                            legend([ySurface1,ySurface2,ySurface3],'Group 1 > Group 2','Group 2 > Group 1',['Group2 ',char(8596),' -(Group1)'])
                        else
                            legend([ySurface1,ySurface2],'Group 1 > Group 2','Group 2 > Group 1')
                        end
                    else if exist('ex3_1')
                            if exist('ex3_3')
                                legend([ySurface1,ySurface3],'Group 1 > Group 2',['Group2 ',char(8596),' -(Group1)'])
                            else
                                legend(ySurface1,'Group 1 > Group 2')
                            end
                        else
                             if exist('ex3_3')
                                 legend([ySurface2,ySurface3],'Group 2 > Group 1',['Group2 ',char(8596),' -(Group1)'])
                             else
                                legend(ySurface2,'Group 2 > Group 1')
                             end
                         end
                    end
                else if exist('ex3_3')
                         legend(ySurface3,['Group2 ',char(8596),' -(Group1)']) % Simbolos apenas copiando e colando...
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

                light('Style','infinite');%'Position',[1 0 1]);

                saveas(GraphoFig3,[Dir filesep 'FDR_Corrected_Graph_3D_Subj_' num2str(Sjo)],'fig')
                set(gcf, 'InvertHardCopy', 'off');

    %             saveas(GraphoFig3,[Dir filesep 'FDR_Corrected_Graph_3D'],'png')
                drawnow
                imgRR = getframe(gcf);
                imwrite(imgRR.cdata, [Dir filesep 'FDR_Corrected_Graph_3D_Subj_' num2str(Sjo) '.png']);
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
                    string = sprintf('Netowork %d color',uiC);
                    text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
                end
                set(gca,'xtick',[])
                set(gca,'xticklabel',[])
                set(gca,'ytick',[])
                set(gca,'yticklabel',[])


    %             saveas(CoLeg,[Dir filesep 'Color_Legend'],'png')
                drawnow
                imgRR = getframe(gcf);
                imwrite(imgRR.cdata, [Dir filesep 'Color_Legend' '.png']);
                clear imgRR

                saveas(CoLeg,[Dir filesep 'Color_Legend'],'fig')
                close(CoLeg)

            end

            if get(handles.Grapho2D,'Value')  %%%%% FDR corrected 2D Grapho
                CoordMS2 = 3.*CoordMS;
                GraphoFig4 = figure;
                set(GraphoFig4,'Name','Connectome2D',...
                    'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.45 ScreSize(1)*.45]),...
                                'Color',[1 .949 .867]);
                title('FDR Corrected 2D Connectome','FontSize',14);

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
                [rws,res,rdw] = sphere(15);

                for huy = 1:size(CoordMS2,1)
                    CoordQ2 = CoordMS2(huy,:);
                    CoordList = CoordMS2(huy+1:end,:);
                    CoordV4 = sum(FDRcor(:,huy)~=-1);
                    CoordV4 = 1+ (2*CoordV4/size(CoordMS2,1));
                    tester2 = huy<=(str_nEle+nEle)-1;
                    networkNumb2 = find(tester2==1,1);
                    colorV2 = cc(networkNumb2,:);
                    for hgf = 1:size(CoordList,1)
                        isF = [CoordQ2(1,1);CoordList(hgf,1)];
                        jsF = [CoordQ2(1,2);CoordList(hgf,2)];
                        ksF = [CoordQ2(1,3);CoordList(hgf,3)];
                        if FDRcor(hgf+huy,huy)>-1.*(str2num(get(handles.edit2,'String')))
                            if (g1_mean(hgf+huy,huy)*g2_VAR(hgf+huy,huy,Sjo))>0;
                                if abs(g1_mean(hgf+huy,huy))>abs(g2_VAR(hgf+huy,huy,Sjo))
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

                teste = FDRcor;
                teste(teste==-1)=0;
                teste(teste~=0)=1;
                vetTTT = sum(teste,1);
                VetFF = vetTTT;
                VetFF(VetFF>1) = 1;
                binVet = (VetFF~=0);
                idxROIs = find(binVet);

                for i = 1:sum(VetFF~=0)

                    tester = idxROIs(i)<=(str_nEle+nEle)-1;
                    networkNumb = find(tester==1,1);
                    colorV = cc(networkNumb,:);
                    CoordQ = CoordMS2(idxROIs(i),:);
                    CoordV4 = (vetTTT(idxROIs(i))/max(vetTTT))+1;
                    rNumb = find(VetFF~=0,i);
                    rNumb = rNumb(end);
                    if isequal(networkNumb,1)
                        text(CoordQ(1)+(0.9*(5*CoordV4)),CoordQ(2)+(0.9*(5*CoordV4)),CoordQ(3)+(0.9*(5*CoordV4)),['r' num2str(rNumb) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',11);
                    else
                        text(CoordQ(1)+(0.9*(5*CoordV4)),CoordQ(2)+(0.9*(5*CoordV4)),CoordQ(3)+(0.9*(5*CoordV4)),['r' num2str(rNumb-sum(nEle(1:networkNumb-1))) 'n' num2str(networkNumb)],'HorizontalAlignment','Right','FontSize',11);
                    end
                    hSurface = surf(CoordV4.*rws.*3+CoordQ(1,1),CoordV4.*3.*res+CoordQ(1,2),CoordV4.*3.*rdw+CoordQ(1,3));
                    set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if exist('ex4_1') || exist('ex4_2')
                    if exist('ex4_1') && exist('ex4_2')
                        if exist('ex4_3')
                            legend([pLine1,pLine2,pLine3],'Group 1 > Group 2','Group 2 > Group 1',['Group2 ',char(8596),' -(Group1)'])
                        else
                            legend([pLine1,pLine2],'Group 1 > Group 2','Group 2 > Group 1')
                        end
                    else if exist('ex4_1')
                            if exist('ex4_3')
                                legend([pLine1,pLine3],'Group 1 > Group 2',['Group2 ',char(8596),' -(Group1)'])
                            else
                                legend(pLine1,'Group 1 > Group 2')
                            end
                         else
                            if exist('ex4_3')
                                legend([pLine2,pLine3],'Group 1 > Group 2',['Group2 ',char(8596),' -(Group1)'])
                            else
                                legend(pLine2,'Group 1 > Group 2')
                            end
                        end
                     end
                else if exist('ex4_3')
                         legend(pLine3,['Group2 ',char(8596),' -(Group1)'])
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
                light('Style','infinite');

                set(gcf, 'InvertHardCopy', 'off');
                saveas(GraphoFig4,[Dir filesep 'FDR_Corrected_Graph_2D_Subj_' num2str(Sjo)],'fig')
    %             saveas(GraphoFig4,[Dir filesep 'FDR_Corrected_Graph_2D'],'png')
                drawnow
                imgRR = getframe(gcf);
                imwrite(imgRR.cdata, [Dir filesep 'FDR_Corrected_Graph_2D_Subj_' num2str(Sjo) '.png']);
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
                    string = sprintf('Netowork %d color',uiC);
                    text('Position',[0.1,uiC-0.5],'String',string,'FontSize',12)
                end
    %             saveas(CoLeg,[Dir filesep 'Color_Legend'],'png')
                set(gca,'xtick',[])
                set(gca,'xticklabel',[])
                set(gca,'ytick',[])
                set(gca,'yticklabel',[])

                drawnow
                imgRR = getframe(gcf);
                imwrite(imgRR.cdata, [Dir filesep 'Color_Legend' '.png']);
                clear imgRR

                saveas(CoLeg,[Dir filesep 'Color_Legend'],'fig')
                close(CoLeg)

            end
        else
            warndlg('No significant results with FDR correction','Ops!')
        end
    end
    clear Pmap FDRcor tPmap adj_p res unCorPmap matrixR
    clear('ex1_1','ex2_1','ex1_2','ex2_2','ex3_1','ex3_2','ex4_2','ex4_1')
end

NRight = sum(CoordMS(:,1)>=51);
NcombRight = size(combnk([1:NRight],2),1);

NLeft = sum(CoordMS(:,1)<=41);
NcombLeft = size(combnk([1:NLeft],2),1);

if get(handles.chec2,'Value') && get(handles.chec1,'Value')
    FDRCorrLat_Right  = (FDRCorrLat_Right'./NcombRight).*10;
    UnCorrLat_Right  = (UnCorrLat_Right'./NcombRight).*10;
    FDRCorrLat_Left2  = (FDRCorrLat_Left'./NcombLeft).*10;
    UnCorrLat_Left  = (UnCorrLat_Left'./NcombLeft).*10;
    save([Dir filesep 'FuncLateIdxs'],'FDRCorrLat_Right','FDRCorrLat_Left','UnCorrLat_Right','UnCorrLat_Left')
end

if get(handles.chec2,'Value') && get(handles.chec1,'Value')==0
    FDRCorrLat_Right  = (FDRCorrLat_Right'./NcombRight).*10;
    FDRCorrLat_Left  = (FDRCorrLat_Left'./NcombLeft).*10;
    save([Dir filesep 'FuncLateIdxs'],'FDRCorrLat_Right','FDRCorrLat_Left')
end

if get(handles.chec1,'Value') && get(handles.chec2,'Value')==0
    UnCorrLat_Right  = (UnCorrLat_Right'./NcombRight).*10;
    UnCorrLat_Left  = (UnCorrLat_Left'./NcombLeft).*10;
    save([Dir filesep 'FuncLateAlteValues'],'UnCorrLat_Right','UnCorrLat_Left')
end

set(handles.text13,'String','Done!')
disp('------------')
disp('You can open the *.fig files on the Matlab to interact (rotate, zoom, edit...) with the figures!')

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
