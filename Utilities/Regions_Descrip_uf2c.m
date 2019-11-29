% Brunno Machado de Campos
% University of Campinas, 2014
%
% Copyright (c) 2014, Brunno Machado de Campos
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

[file,path] = uigetfile('*.nii','Add the subjects MAPS','MultiSelect','on');

[fileR,pathR] = uigetfile('*.nii','Add all masks','MultiSelect','on');

ThresholdMask = 0;

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
fideGx1 = fopen([path,'Results.txt'],'w+'); % CHARGE OUTPUT LOG FILE

fideGx = fopen([path,'Masks_Anat_Description.txt'],'w+'); % CHARGE OUTPUT LOG FILE

fprintf(fideGx,'UF²C: ROIs Anatomical Labelling and Quantifications\r\n\r\n');
fprintf(fideGx,'These data are better visualized if imported to a MS Excel type software\r\n\r\n');
fprintf(fideGx,'Credits:\r\n');
fprintf(fideGx,'Maximum probability tissue labels derived from:\r\n');
fprintf(fideGx,'''MICCAI 2012 Grand Challenge and Workshop on Multi-Atlas Labeling''\r\n');
fprintf(fideGx,'These data were released under the Creative Commons Attribution-NonCommercial (CC BY-NC) with no end date.\r\n');
fprintf(fideGx,'Labeled data provided by Neuromorphometrics, Inc. (http://Neuromorphometrics.com/) under academic subscription".\r\n\r\n\r\n');


for il = 1:numel(fileR)
    ncentr = 1;
    fprintf(fideGx,'Mask %d:\t%s\r\n', il,fileR{il});
    roi11_Stru = nifti([pathR,fileR{il}]);
    
    roi11 = roi11_Stru.dat(:,:,:);
    roi12(:,:,:,il) = roi11;
    
    roi11(roi11>0) = 1;
    roi11(roi11<0) = 0;
    roi11(isnan(roi11)) = 0;
    roi11 = double(roi11);
    
    nofVoxTOT = sum(sum(sum(roi11)));
    RoiVol = nofVoxTOT.*prod(roi11_Stru.hdr.pixdim(2:4));
    fprintf(fideGx,'\tNumber of Voxels:\t%d\r\n',nofVoxTOT);
    fprintf(fideGx,'\tTotal Volume:\t%d\tmm³\r\n',RoiVol);

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
    FinalCentrCoord(il,:) = round(median(sFnl,1));

    MultROILab = LabelImgMAT.*roi11;
    
    UniQ = unique(MultROILab);
    UniQ = UniQ(2:end);
    nofregions = numel(UniQ);

    IDXcentroid = LabelImgMAT((size(roi11,1)-FinalCentrCoord(il,1)),FinalCentrCoord(il,2),FinalCentrCoord(il,3));
    if IDXcentroid == 0
        CentroidRegions{il,1} = 'Extra Brain Region';
    else
        Posi = find(NeuromorphometricsLabels.Index == IDXcentroid);
        CentroidRegions{il,1} = NeuromorphometricsLabels.Label{Posi};
    end
    fprintf(fideGx,'\tCentroid Coordinate (voxel space):\t%dx%dx%d \r\n',FinalCentrCoord(il,1),FinalCentrCoord(il,2),FinalCentrCoord(il,3));
    fprintf(fideGx,'\tCentroid Exact Regions:\t%s \r\n\r\n',CentroidRegions{il,1});
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

    clear sFnl nOfvox Regions Posi roiFind centNs ncentr
end
fprintf('\tAnatomical labeling file saved\n');
fprintf('Done! ============== %s\r\n',spm('time'));
fclose(fideGx);

result = cell(size(file,2),size(fileR,2)+1);
result{1,1} = 'Name';
fprintf(fideGx1,'Name \t');
for i2 = 2:size(fileR,2)+1
    result{1,i2} = fileR{i2-1};
    fprintf(fideGx1,'%s \t ',fileR{i2-1});
end
fprintf(fideGx1,'\r\n',fileR{i2-1});

for i = 1:size(file,2)
    file{i} = [path file{i}];
end
for i = 1:size(fileR,2)
    fileR{i} = [pathR fileR{i}];
end

for i = 1:size(file,2)
   [a1,a2,a3] = fileparts(file{i});
   result{i+1,1} = a2;
   fprintf(fideGx1,'%s \t ',a2);
   
   stru = nifti(file{i});
   matI = double(stru.dat(:,:,:));
   matI(matI(:,:,:)<ThresholdMask) = 0;

    for j = 1:size(fileR,2)
        mat1 = roi12(:,:,:,j);
        try
            ResMat1 = matI.*mat1;
            catch
                warndlg('Maybe the  matrices sizes of the Maps and ROI do not match','Process aborted!')
                return
        end
        nume1 = sum(sum(sum(ResMat1)));
        nz1 = nnz(ResMat1);
        rela1 = nume1/nz1;
        result{i+1,j+1} = rela1;
        fprintf(fideGx1,'%.4f\t ',rela1);
    end
    fprintf(fideGx1,'\r\n',fileR{i2-1});
end
fclose(fideGx1);

