function [OutPut4D] = uf2c_AnatDescrip(fileR,OutDir,OutName,type)
% inputs:
% fileR: Cell of strings with ALL images or masks to perform the anatomical
%       e.g.: {'Image1.nii';'Image2.nii';'Image3.nii'}
% OutDir: String with an output directoty
%       e.g.: 'c:\MyFoldrr\'
% OutName: String with and output filename
%       e.g.: 'AnateDescrip.txt'
% type: define if the image is:
% a mask/ROI with one connected component (one cluster): 'mask'
% a map/image with several voxel values and connected components
% (clusters): 'map'
 
if ~exist('OutDir','var')
    OutDir = uigetdir('','Define the output directory');
    if isequal(OutDir,0)
        OutDir = pwd;
    end
end
if ~exist('OutName','var')
    OutName = inputdlg('Define an output filename','Define an output filename',1,{'AnatResults.txt'});
    if isempty(OutName)
        OutName = 'AnatResults.txt';
    else
        OutName = OutName{1};
    end
end
if ~exist('fileR','var')
    [file,path] = uigetfile({'*.nii','NIfTI Files'},...
        'Add all images to analyze','MultiSelect','on');
    if isequal(file,0)
        warndlg('No Images added, process aborted')
        return
    else
        if ~iscell(file)   % CREATING A CELL VARIABLE FOR SINGULAR INPUTS
            file = {file};
        end
        for i = 1:numel(file)
            fileR{i,1} = fullfile(path,file{i});
        end
    end
end
if ~exist('type','var')
    type = questdlg('Define the input image(s) type ans mask (a mask/ROI with one connected component (one cluster)) or a map/image with several voxel values and connected components',...
        '','Mask','Map','Mask');
        if isempty(type)
            type = 'map';
        end
end

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

    fideGx = fopen([OutDir,filesep,OutName],'w+'); % CHARGE OUTPUT LOG FILE

    fprintf(fideGx,'UF²C: ROIs Anatomical Labelling and Quantifications\r\n\r\n');
    fprintf(fideGx,'These data are better visualized if imported to a MS Excel type software\r\n\r\n');
    fprintf(fideGx,'Credits:\r\n');
    fprintf(fideGx,'Maximum probability tissue labels derived from:\r\n');
    fprintf(fideGx,'''MICCAI 2012 Grand Challenge and Workshop on Multi-Atlas Labeling''\r\n');
    fprintf(fideGx,'These data were released under the Creative Commons Attribution-NonCommercial (CC BY-NC) with no end date.\r\n');
    fprintf(fideGx,'Labeled data provided by Neuromorphometrics, Inc. (http://Neuromorphometrics.com/) under academic subscription".\r\n\r\n\r\n');

    fprintf('\tPerforming Anatomical labeling\n');
    if ~exist('type','VAR')
        type = 'mask';
    end
    
    if isequal(type,'mask')
        for il = 1:numel(fileR)
            clear sFnl nOfvox Regions Posi roiFind centNs ncentr
            [~,bx,~] = fileparts(fileR{il});
            fprintf(fideGx,'Mask %d:\t%s\r\n', il,bx);
            anatDO = 1;
            ncentr = 1;
            roi11_Stru = nifti(fileR{il});

            roi11 = roi11_Stru.dat(:,:,:);
            
            roi11(roi11>0) = 1;
            roi11(roi11<0) = 0;
            roi11(isnan(roi11)) = 0;
            roi11 = double(roi11);
            OutPut4D(:,:,:,il) = roi11;

            nofVoxTOT = sum(sum(sum(roi11)));
            RoiVol = nofVoxTOT.*prod(roi11_Stru.hdr.pixdim(2:4));
            fprintf(fideGx,'\tNumber of Voxels:\t%d\r\n',nofVoxTOT);
            fprintf(fideGx,'\tTotal Volume:\t%d\tmm³\r\n',RoiVol);
            sFnl = zeros(size(roi11,3),3);
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
                end
            end
            sFnl(~any(sFnl,2), : ) = [];
            FinalCentrCoord(il,:) = round(median(sFnl,1));
            MultROILab = LabelImgMAT.*roi11;
            UniQ = unique(MultROILab);
            UniQ = UniQ(2:end);
            nofregions = numel(UniQ);
            if ~isequal(sum(FinalCentrCoord(il,:)),0)
                IDXcentroid = LabelImgMAT((size(roi11,1)-FinalCentrCoord(il,1)),FinalCentrCoord(il,2),FinalCentrCoord(il,3));
                if IDXcentroid == 0
                    CentroidRegions{il,1} = 'Extra Brain Region';
                else
                    Posi = find(NeuromorphometricsLabels.Index == IDXcentroid);
                    CentroidRegions{il,1} = NeuromorphometricsLabels.Label{Posi};
                end
                fprintf(fideGx,'\tCentroid Coordinate (voxel space):\t%dx%dx%d \r\n',FinalCentrCoord(il,1),FinalCentrCoord(il,2),FinalCentrCoord(il,3));
                fprintf(fideGx,'\tCentroid Exact Regions:\t%s \r\n\r\n',CentroidRegions{il,1});
                fprintf(fideGx,'\tRegions included on the ROI mask:\r\n');
                fprintf(fideGx,'\t In Brain Region Name\t Number of Voxels\t Percentage of the Total\r\n');

                if IDXcentroid
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
                end
            else
                fprintf(fideGx,'\tCentroid Coordinate (voxel space): No regions found');
            end

            fprintf(fideGx,'\r\n\r\n');
        end
    else
        for il = 1:numel(fileR)
            clear sFnl nOfvox Regions Posi roiFind centNs ncentr
            
            [~,bx,~] = fileparts(fileR{il});
            fprintf(fideGx,'Mask %d:\t%s\r\n', il,bx);
            anatDO = 1;
            ncentr = 1;
            roi11_Stru = nifti(fileR{il});

            roi11 = roi11_Stru.dat(:,:,:);
            roi11ori = roi11;
            
            roi11(roi11>0) = 1;
            roi11(roi11<0) = 0;
            roi11(isnan(roi11)) = 0;
            roi11 = double(roi11);
            OutPut4D(:,:,:,il) = roi11;
            
            if ~isequal(max(max(max(roi11ori))),0)
                [xGlo,yGlo,zGlo] = ind2sub(size(roi11ori),find(roi11ori==max(max(max(roi11ori)))));
                GlobMap = zeros(size(roi11ori));
                GlobMap(xGlo,yGlo,zGlo) = 1;

                GlobFind = GlobMap(:,:,zGlo);
                GlobFind = fliplr(GlobFind);
                GlobFind = permute(GlobFind,[2,1]);
                GlobFind = fliplr(GlobFind);
                GlobFind = flipud(GlobFind);

                GlobFind2 = logical(GlobFind);
                sg = regionprops(GlobFind2,'centroid');
                sFnlg(1,1:2) = sg.Centroid(1,:);
                sFnlg(1,3) = zGlo;  
                
                MultROILabG = LabelImgMAT.*GlobMap;
                UniQG = unique(MultROILabG);
                UniQG = UniQG(2:end);
                
                if isequal(UniQG,0) || isempty(UniQG)
                    CentroidRegionsG = 'Extra Brain Region';
                else
                    PosiG = find(NeuromorphometricsLabels.Index == UniQG);
                    CentroidRegionsG = NeuromorphometricsLabels.Label{PosiG};
                end

                nofVoxTOT = sum(sum(sum(roi11)));
                RoiVol = nofVoxTOT.*prod(roi11_Stru.hdr.pixdim(2:4));
                fprintf(fideGx,'\tNumber of Voxels:\t%d\r\n',nofVoxTOT);
                fprintf(fideGx,'\tTotal Volume:\t%d\tmm³\r\n',RoiVol);
                fprintf(fideGx,'\tGlobal Maxima Intensity:\t%.2f\r\n',max(max(max(roi11ori))));
                fprintf(fideGx,'\tGlobal Maxima Region:\t%s\r\n\r\n',CentroidRegionsG);

                L = bwlabeln(roi11);
                vet = unique(L);
                vet = vet(2:end);
                iidx = 0;

                for i = 1:size(vet,1)
                    clear CentroidRegions IDXcentroid UniQ nofregions FinalCentrCoord sFnl
                    TempRoiM = double(L==vet(i));
                    fprintf(fideGx,'\tSub Region %d:\t\r\n',i);
                    sFnl = zeros(size(roi11,3),3);
                    TempRoi = TempRoiM.*roi11ori;
                    TempRoi(isnan(TempRoi)) = 0;
                    SubRoiNvox = nnz(TempRoi);
                    MeanValue = sum(sum(sum(TempRoi)))/SubRoiNvox;
                    SubRoiVol = nnz(TempRoi).*prod(roi11_Stru.hdr.pixdim(2:4));


                    for yt = 1:size(TempRoiM,3)
                        roiFind = TempRoiM(:,:,yt);
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
                        end
                    end

                    sFnl(~any(sFnl,2), : ) = [];
                    FinalCentrCoord(il,:) = round(median(sFnl,1));
                    MultROILab = LabelImgMAT.*TempRoiM;

                    UniQ = unique(MultROILab);
                    UniQ = UniQ(2:end);
                    nofregions = numel(UniQ);
                    if ~isequal(sum(FinalCentrCoord(il,:)),0) && ~isequal(nofregions,0)
                        IDXcentroid = LabelImgMAT((size(TempRoiM,1)-FinalCentrCoord(il,1)),FinalCentrCoord(il,2),FinalCentrCoord(il,3));
                        if IDXcentroid == 0
                            CentroidRegions{il,1} = 'Extra Brain Region';
                        else
                            Posi = find(NeuromorphometricsLabels.Index == IDXcentroid);
                            CentroidRegions{il,1} = NeuromorphometricsLabels.Label{Posi};
                        end
                        fprintf(fideGx,'\t\tCentroid Coordinate (voxel space):\t%dx%dx%d \r\n',FinalCentrCoord(il,1),FinalCentrCoord(il,2),FinalCentrCoord(il,3));
                        fprintf(fideGx,'\t\tCentroid Exact Regions:\t%s \r\n',CentroidRegions{il,1});
                        fprintf(fideGx,'\t\tSub Region Mean Value:\t%.2f \r\n',MeanValue);
                        fprintf(fideGx,'\t\tSub Region Volume:\t%.2fmm³ \r\n\r\n',SubRoiVol);
                        fprintf(fideGx,'\t\tRegions included on the ROI mask:\r\n');
                        fprintf(fideGx,'\t\t In Brain Region Name\t Number of Voxels\t Percentage of the Total\r\n');

                        if IDXcentroid
                            for nr = 1:nofregions
                                Bint = double(MultROILab==UniQ(nr));
                                nOfvox(nr,1) = sum(sum(sum(Bint)));
                                Posi = find(NeuromorphometricsLabels.Index == UniQ(nr));
                                Regions{nr,1} = NeuromorphometricsLabels.Label{Posi};
                            end
                            [nOfvox,www] = sort(nOfvox,'descend');
                            Regions = Regions(www);

                            for nr = 1:nofregions
                                fprintf(fideGx,'\t\t%s\t%d\t%.2f\r\n',Regions{nr,1},nOfvox(nr,1),(100*(nOfvox(nr,1)/SubRoiNvox)));
                            end
                        end
                    else
                        fprintf(fideGx,'\t\tCentroid Coordinate (voxel space): No regions found');
                    end
                    fprintf(fideGx,'\r\n\r\n');
                end
            else
                fprintf(fideGx,'\t\t Empty Mask');
            end
        end
    end
    
    if anatDO
        fprintf('\tAnatomical labeling file saved\n');
    else
        fprintf('\tAttention! Your images resolutions do not enable anatomical descriptions.\r\n');
        fprintf(fideGx,'Attention! Your images resolutions do not enable anatomical descriptions.');
    end
    fclose(fideGx);
end
