coordList = 0;
AddCovar = 1;

[fileG1,pathG1] = uigetfile({'*.mat','Matlab file'},...
    'Select the result file (All_Subjs-VAR.mat) for group 1','MultiSelect','off');

load([pathG1 fileG1])
g1_VAR = AllSubj3D;

[fileG2,pathG2] = uigetfile({'*.mat','Matlab file'},...
    'Select the result file (All_Subjs-VAR.mat) for group 2','MultiSelect','off');

% prompt={'Enter the uncorrected p-value (threshold)',...
%         'Enter the FDR-corrected p-value (threshold)'};
% name='Input thresholds';
% numlines=1;
% defaultanswer={'0.01','0.05'};
% answer=inputdlg(prompt,name,numlines,defaultanswer);

load([pathG2 fileG2])
g2_VAR = AllSubj3D;
Dir = uigetdir('','Define the output directory');
mkdir([Dir filesep 'UF2C_SecondLevel'])

if AddCovar
    options.Resize = 'on';
    options.WindowStyle = 'modal';
    options.Interpreter = 'none';
    numlines=1;
    defaultanswer = {''};
    COV = inputdlg('Ad a covariate','Covariate',[30 20],defaultanswer,options);
    COV = str2num(COV{1});
    if size(COV,1)<size(COV,2)
        COV = transpose(COV);
    end
end


if coordList
    [fileCor,pathCor] = uigetfile({'*.txt','Text file'},...
        'Select the cordinates list file','MultiSelect','off');
    
    seedListT = importdata([pathCor fileCor]);
    nOfNets = str2num(seedListT{end,1}(1:2));
    xs=1;
    nE = 1;
    while xs < size(seedListT,1)
        xStr = seedListT{xs,1}(1:3);
        net1 = strncmp(xStr,seedListT,3);
        nEle(nE) = sum(net1);
        xs = xs + sum(net1);
        nE = nE+1;
    end
    for i = 1:size(seedListT,1)
        seedList(i,:) = str2num(seedListT{i}(4:end));
    end
    for jh = 1:size(seedList,1)
        Vhdr = spm_vol(['sw',fileFunc]);
        TcM = Vhdr.mat;
        EPIcoord = [seedList(jh,1) seedList(jh,2) seedList(jh,3) 1]*(inv(TcM))';
        EPIcoord = round(EPIcoord);
        CoordMS(jh,:) = EPIcoord(1,:);
    end
else
    [fileROI,pathROI] = uigetfile({'*.nii','NIfTI files'},'Select 3D ROI masks','MultiSelect','on');
    if ~iscell(fileROI)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
        fileROI = {fileROI};
    end
    fileROI = sort(fileROI);
    xs=1;
    nE = 1;
    while xs < size(fileROI,2)
       xStr = fileROI{1,xs}(1:3);
       net1 = strncmp(xStr,fileROI,3);
       nEle(nE) = sum(net1);
       xs = xs + sum(net1);
       nE = nE+1;
    end
    nOfNets = size(nEle,2);
    loopSize = size(fileROI,2);
    
    for jh = 1:size(fileROI,2)
        struRoi = nifti([pathROI fileROI{jh}]);
        ROIMAt = struRoi.dat(:,:,:);
        ROIMAt(ROIMAt>0)=1;
        roi1 = ROIMAt;
        [iCo,jCo,kCo]= ind2sub(size(roi1), find(roi1>0));
        coords = [iCo,jCo,kCo];
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
    seedList = transpose(fileROI);
end

for i = 1:size(CoordMS,1)
    for j = 1:size(CoordMS,1)
        if i==j
            Pmap(i,j) = 1;
        else
        %[H,P,CI,tvalue] = ttest2(squeeze(g1_VAR(i,j,:)),squeeze(g2_VAR(i,j,:)));
            if AddCovar
                [T,P] = mancovan([squeeze(g1_VAR(i,j,:)); squeeze(g2_VAR(i,j,:))],...
                    [ones(size(g1_VAR,3),1); zeros(size(g2_VAR,3),1)],COV);
                Pmap(i,j) = P(1);
            else
                [T,P] = mancovan([squeeze(g1_VAR(i,j,:)); squeeze(g2_VAR(i,j,:))],...
                    [ones(size(g1_VAR,3),1); zeros(size(g2_VAR,3),1)],[],'group-group');
                Pmap(i,j) = P;
            end
        end
    end
end

% Uncorrected P Map
unCorPmap = Pmap;
unCorPmap(unCorPmap>str2num(answer{1})) = 1;
unCorPmap = -1.*unCorPmap;
fig1 = figure;
imagesc(unCorPmap)
title('Uncorrected P Map')
colorbar
caxis([-1.*(str2num(answer{1})) (max(max(unCorPmap))-(0.2*max(max(unCorPmap))))]) %reescala a volorbar para fugir dos extremos

saveas(fig1,[Dir filesep 'UF2C_SecondLevel' filesep 'Uncorrected_P_Map'],'png')
saveas(fig1,[Dir filesep 'UF2C_SecondLevel' filesep 'Uncorrected_P_Map'],'fig')
save([Dir filesep 'UF2C_SecondLevel' filesep 'unCorPmap-VAR'],'unCorPmap')
close(fig1)

GraphoFig = figure;
set(GraphoFig,'Name','Connecticome3D','Position', [130 130 1350 800],...
                'Color',[1 .949 .867]);
title('Connecticome3D','FontSize',14);
xlabel('X axis'  ,'FontSize',14);
ylabel('Y axis','FontSize',14);
zlabel('Z axis','FontSize',14);
scatter3(CoordMS(:,1),CoordMS(:,2),CoordMS(:,3),'Marker','square','MarkerFaceColor','red','MarkerEdgeColor','k');
set(gca, 'color', [0 0 0])
set(gca, 'gridcolor', [0.2 0.2 0.2])
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

 cc=hsv(12);
 [rws,res,rdw] = sphere(15);

for huy = 1:size(CoordMS,1)
    CoordQ = CoordMS(huy,:);
    CoordList = CoordMS(huy+1:end,:);
    CoordV4 = sum(unCorPmap(:,huy)~=-1);
    CoordV4 = 1+ (2*CoordV4/size(CoordMS,1));
    tester = huy<=(str_nEle+nEle)-1;
    networkNumb = length(find(tester==1));
    colorV = cc(networkNumb,:);
    for hgf = 1:size(CoordList,1)
        isF = [CoordQ(1,1);CoordList(hgf,1)];
        jsF = [CoordQ(1,2);CoordList(hgf,2)];
        ksF = [CoordQ(1,3);CoordList(hgf,3)];
        if unCorPmap(hgf+huy,huy)>-0.05
            hSurface = surf(CoordV4.*rws+CoordQ(1,1),CoordV4.*res+CoordQ(1,2),CoordV4.*rdw+CoordQ(1,3));
            set(hSurface, 'FaceColor',colorV, 'FaceAlpha',0.6, 'EdgeAlpha', 0,'LineSmoothing','on');
            [xre,yre,zre] = tubeplot([isF(1) isF(2);jsF(1) jsF(2);ksF(1) ksF(2)],0.4,10);
            ySurface = surf(xre,yre,zre);
            set(ySurface,'EdgeColor','none','LineStyle','none','FaceLighting',...
                'gouraud','FaceColor',[0.2 0.5 1],'FaceLighting','gouraud','LineSmoothing','on');
         else
             hSurface = surf(CoordV4.*rws+CoordQ(1,1),CoordV4.*res+CoordQ(1,2),CoordV4.*rdw+CoordQ(1,3));
             set(hSurface, 'FaceColor',colorV, 'FaceAlpha',0.6, 'EdgeAlpha', 0);
        end   
     end
end

hSurface = surf(CoordV4.*rws+CoordMS(end,1),CoordV4.*res+CoordMS(end,2),CoordV4.*rdw+CoordMS(end,3));
set(hSurface, 'FaceColor',colorV, 'FaceAlpha',0.7, 'EdgeAlpha', 0);

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
aphamap1(aphamap1>=300)=1;
set(ater(1),'alphadata',aphamap1,'facealpha','flat')

aphamap3 = aphamap3d(:,:,20);
aphamap3 = squeeze(aphamap3);
aphamap3(aphamap3<300)=0.0001;
aphamap3(aphamap3>=300)=1;
set(ater(2),'alphadata',aphamap3,'facealpha','flat')
light('Position',[1 3 2]);
set(ySurface,'LineSmoothing','off');
set(hSurface,'LineSmoothing','off');

saveas(GraphoFig,[Dir filesep 'UF2C_SecondLevel' filesep 'Uncorrected_Graph_3D_rot'],'fig')
tmpFig = myaa('publish');
saveas(tmpFig,[Dir filesep 'UF2C_SecondLevel' filesep 'Uncorrected_Graph_3D'],'tif')
saveas(tmpFig,[Dir filesep 'UF2C_SecondLevel' filesep 'Uncorrected_Graph_3D'],'fig')
close(tmpFig)

% saveas(GraphoFig,[Dir filesep 'UF2C_SecondLevel' filesep 'Uncorrected_Graph_3D'],'png')
% saveas(GraphoFig,[Dir filesep 'UF2C_SecondLevel' filesep 'Uncorrected_Graph_3D'],'fig')
% close(GraphoFig)

%%%%%%%%%% 2D Uncorrected
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

cc2=hsv(12);

for huy = 1:size(CoordMS,1)
    CoordQ2 = CoordMS(huy,:);
    CoordList = CoordMS(huy+1:end,:);
    CoordV4 = sum(unCorPmap(:,huy)~=-1);
    CoordV4 = 1+ (2*CoordV4/size(CoordMS,1));
    tester2 = huy<=(str_nEle+nEle)-1;
    networkNumb2 = length(find(tester2==1));
    colorV2 = cc2(networkNumb2,:);
    for hgf = 1:size(CoordList,1)
        isF = [CoordQ2(1,1);CoordList(hgf,1)];
        jsF = [CoordQ2(1,2);CoordList(hgf,2)];
        ksF = [CoordQ2(1,3);CoordList(hgf,3)];
        if unCorPmap(hgf+huy,huy)>-0.05
            hSurface = surf(CoordV4.*rws+CoordQ2(1,1),CoordV4.*res+CoordQ2(1,2),CoordV4.*rdw+CoordQ2(1,3));
            set(hSurface, 'FaceColor',colorV2, 'FaceAlpha',1, 'EdgeAlpha',0);
            yLine = line(isF,jsF,ksF,'Color',[0 0 0],'LineStyle','-',...
                'LineWidth',0.5,'LineSmoothing','on');
        else
            hSurface = surf(CoordV4.*rws+CoordQ2(1,1),CoordV4.*res+CoordQ2(1,2),CoordV4.*rdw+CoordQ2(1,3));
            set(hSurface, 'FaceColor',colorV2, 'FaceAlpha',1, 'EdgeAlpha', 0);
        end   
    end
end
hSurface = surf(CoordV4.*rws+CoordQ2(end,1),CoordV4.*res+CoordQ2(end,2),CoordV4.*rdw+CoordQ2(end,3));
set(hSurface, 'FaceColor',colorV2, 'FaceAlpha',1, 'EdgeAlpha', 0);

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
light('Position',[1 3 2]);

tmpFig = myaa('publish');
saveas(tmpFig,[Dir filesep 'UF2C_SecondLevel' filesep 'Uncorrected_Graph_2D'],'tif')
saveas(tmpFig,[Dir filesep 'UF2C_SecondLevel' filesep 'Uncorrected_Graph_2D'],'fig')
close(tmpFig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D FDR corrected

FDRcor = Pmap;
tPmap = tril(FDRcor,-1);
tPmap = tPmap(tPmap>0);

[h, crit_p, adj_p] = fdr_bh(tPmap,str2num(answer{2}),'pdep','yes');

res = zeros(size(CoordMS,1),size(CoordMS,1));
stt = 1;

for i = 1:size(CoordMS,1)
    res(1:i,i) = zeros(i,1);
    res(i+1:end,i) = adj_p(stt:stt+size(CoordMS,1)-i-1);
    stt = stt+(size(CoordMS,1)-i);
end

resF = res+transpose(res);
FDRcor = resF;

FDRcor(FDRcor>(str2num(answer{2}))) = 1;
FDRcor(FDRcor==0) = 1;
FDRcor = -1.*FDRcor;
fig2 = figure;
imagesc(FDRcor)
title('FDR corrected P Map')
colorbar
caxis([-1.*(str2num(answer{2})) (-1*(min(adj_p)-(0.2*min(adj_p))))])

saveas(fig2,[Dir filesep 'UF2C_SecondLevel' filesep 'FDR_corrected_P_Map'],'jpg')
saveas(fig2,[Dir filesep 'UF2C_SecondLevel' filesep 'FDR_corrected_P_Map'],'fig')
save([Dir filesep 'UF2C_SecondLevel' filesep 'unCorPmap-VAR'],'FDRcor')
close(fig2)

GraphoFig3 = figure;
set(GraphoFig3,'Name','Connecticome3D','Position', [130 130 1350 800],...
                'Color',[1 .949 .867]);
title('Connecticome3D','FontSize',14);
xlabel('X axis'  ,'FontSize',14);
ylabel('Y axis','FontSize',14);
zlabel('Z axis','FontSize',14);
scatter3(CoordMS(:,1),CoordMS(:,2),CoordMS(:,3),'Marker','square','MarkerFaceColor','red','MarkerEdgeColor','k');
set(gca, 'color', [0 0 0])
set(gca, 'gridcolor', [0.2 0.2 0.2])
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

 cc=hsv(12);
 [rws,res,rdw] = sphere(15);

 for huy = 1:size(CoordMS,1)
    CoordQ = CoordMS(huy,:);
    CoordList = CoordMS(huy+1:end,:);
    CoordV4 = sum(FDRcor(:,huy)~=-1);
    CoordV4 = 1+ (2*CoordV4/size(CoordMS,1));
    tester = huy<=(str_nEle+nEle)-1;
    networkNumb = length(find(tester==1));
    colorV = cc(networkNumb,:);
    for hgf = 1:size(CoordList,1)
        isF = [CoordQ(1,1);CoordList(hgf,1)];
        jsF = [CoordQ(1,2);CoordList(hgf,2)];
        ksF = [CoordQ(1,3);CoordList(hgf,3)];
        if FDRcor(hgf+huy,huy)>-0.05
            hSurface = surf(CoordV4.*rws+CoordQ(1,1),CoordV4.*res+CoordQ(1,2),CoordV4.*rdw+CoordQ(1,3));
            set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
            [xre,yre,zre] = tubeplot([isF(1) isF(2);jsF(1) jsF(2);ksF(1) ksF(2)],0.4,10);
            surf(xre,yre,zre,'EdgeColor','none','LineStyle','none','FaceLighting','gouraud','FaceColor',[0.2 0.5 1],'FaceLighting','gouraud')
        else
             hSurface = surf(CoordV4.*rws+CoordQ(1,1),CoordV4.*res+CoordQ(1,2),CoordV4.*rdw+CoordQ(1,3));
             set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);
        end   
    end
 end

hSurface = surf(CoordV4.*rws+CoordMS(end,1),CoordV4.*res+CoordMS(end,2),CoordV4.*rdw+CoordMS(end,3));
set(hSurface, 'FaceColor',colorV, 'FaceAlpha',1, 'EdgeAlpha', 0);

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
aphamap1(aphamap1>=300)=1;
set(ater(1),'alphadata',aphamap1,'facealpha','flat')

aphamap3 = aphamap3d(:,:,20);
aphamap3 = squeeze(aphamap3);
aphamap3(aphamap3<300)=0.0001;
aphamap3(aphamap3>=300)=1;
set(ater(2),'alphadata',aphamap3,'facealpha','flat')
light('Position',[1 3 2]);

saveas(GraphoFig3,[Dir filesep 'UF2C_SecondLevel' filesep 'FDR_Corrected_Graph_3D_rot'],'fig')
tmpFig = myaa('publish');
saveas(tmpFig,[Dir filesep 'UF2C_SecondLevel' filesep 'FDR_Corrected_Graph_3D'],'tiff')
saveas(tmpFig,[Dir filesep 'UF2C_SecondLevel' filesep 'FDR_Corrected_Graph_3D'],'fig')
close(tmpFig)

% saveas(GraphoFig3,[Dir filesep 'UF2C_SecondLevel' filesep 'FDR_Corrected_Graph_3D'],'png')
% saveas(GraphoFig3,[Dir filesep 'UF2C_SecondLevel' filesep 'FDR_Corrected_Graph_3D'],'fig')
% close(GraphoFig3)

%%%%%%%%%%
%%%%%%%%%%
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

cc2=hsv(12);

for huy = 1:size(CoordMS,1)
    CoordQ2 = CoordMS(huy,:);
    CoordList = CoordMS(huy+1:end,:);
    CoordV4 = sum(FDRcor(:,huy)~=-1);
    CoordV4 = 1+ (2*CoordV4/size(CoordMS,1));
    tester2 = huy<=(str_nEle+nEle)-1;
    networkNumb2 = length(find(tester2==1));
    colorV2 = cc2(networkNumb2,:);
    for hgf = 1:size(CoordList,1)
        isF = [CoordQ2(1,1);CoordList(hgf,1)];
        jsF = [CoordQ2(1,2);CoordList(hgf,2)];
        ksF = [CoordQ2(1,3);CoordList(hgf,3)];
        if FDRcor(hgf+huy,huy)>-0.05
            hSurface = surf(CoordV4.*rws+CoordQ2(1,1),CoordV4.*res+CoordQ2(1,2),CoordV4.*rdw+CoordQ2(1,3));
            set(hSurface, 'FaceColor',colorV2, 'FaceAlpha',1, 'EdgeAlpha',0);
            line(isF,jsF,ksF,'Color',[0 0 0],'LineStyle','-','LineWidth',0.5);
        else
            hSurface = surf(CoordV4.*rws+CoordQ2(1,1),CoordV4.*res+CoordQ2(1,2),CoordV4.*rdw+CoordQ2(1,3));
            set(hSurface, 'FaceColor',colorV2, 'FaceAlpha',1, 'EdgeAlpha', 0);
        end   
    end
end
hSurface = surf(CoordV4.*rws+CoordQ2(end,1),CoordV4.*res+CoordQ2(end,2),CoordV4.*rdw+CoordQ2(end,3));
set(hSurface, 'FaceColor',colorV2, 'FaceAlpha',1, 'EdgeAlpha', 0);

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
light('Position',[1 3 2]);

tmpFig = myaa('publish');
saveas(tmpFig,[Dir filesep 'UF2C_SecondLevel' filesep 'FDR_Corrected_Graph_2D'],'tif')
saveas(tmpFig,[Dir filesep 'UF2C_SecondLevel' filesep 'FDR_Corrected_Graph_2D'],'fig')
close(tmpFig)

% saveas(GraphoFig4,[Dir filesep 'UF2C_SecondLevel' filesep 'FDR_Corrected_Graph_2D'],'png')
% saveas(GraphoFig4,[Dir filesep 'UF2C_SecondLevel' filesep 'FDR_Corrected_Graph_2D'],'fig')
% close(GraphoFig4)

