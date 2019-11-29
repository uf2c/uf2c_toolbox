fileFuncT1 = uipickfiles('Output','cell','Prompt','Add all images','REFilter','');

OnlyShape = 0; % use 1 if you want to consider only the shape (binary) or
% 0 if you want to accoun t to the voxels intensities also (normalized)
Theshold = 0.1; % Threshold on absolute values
OutLierMode = 'severe'; % 'minor', 'major' or 'severe'

for i = 1:size(fileFuncT1,2)
    stru = nifti(fileFuncT1{i});
    if numel(stru.dat.dim)>3
        mattmp = stru.dat(:,:,:,:);
        mattmp(mattmp<Theshold) = 0;
        matF(:,:,:,i) = median(mattmp,4);
    else
        mattmp = stru.dat(:,:,:);
        mattmp(mattmp<Theshold) = 0;
        matF(:,:,:,i) = mattmp;
    end
    [a,filename{i},c] = fileparts(fileFuncT1{i});
end

Dx = size(matF,1);
Dy = size(matF,2);
Dz = size(matF,3);
Dv = size(matF,4);

ReshMat = reshape(matF,prod([Dx,Dy,Dz]),Dv)';

if OnlyShape
    for u2 = 1:size(ReshMat,2)
        tmpVet = ReshMat(:,u2);
        tmpVet(tmpVet>0) = max(tmpVet);
        ReshMat(:,u2) = tmpVet;
    end
else
    for u = 1:size(ReshMat,1)
        ReshMat(u,:) = ReshMat(u,:)./max(ReshMat(u,:));
    end
end

finalOut = [];
MatOut = zeros(size(fileFuncT1,2),size(ReshMat,2));
for j = 1:size(ReshMat,2)
    if ~isequal(sum(ReshMat(:,j)),0)
        a = outlierdetec_uf2c(ReshMat(:,j),OutLierMode);
        MatOut(1:numel(a),j) = a;
        clear a
    end
end

volsoccur = unique(MatOut);
volsoccur = volsoccur(2:end);

for k = 1:size(volsoccur,1)
    s(k) = sum(sum(volsoccur(k)==MatOut));
end

[s2,inds] = sort(s,'descend');
filename2 = filename(inds);

OutImages = outlierdetec_uf2c(s2,OutLierMode,'upperside');

OutImages = inds(OutImages);

if ~isempty(OutImages)
    disp('Outliers:')
    for i3 = 1:numel(OutImages)
        fprintf('%s\n',fileFuncT1{OutImages(i3)})
    end
    fprintf('\n')
    
    switch OutLierMode
        case 'severe'
            whisker = 1;
        case 'minor'
            whisker = 1.5;
        case 'major'
            whisker = 3;
    end
else
    disp('No outliers found')
end

if ~isempty(OutImages)
    disp('Creating:')
    for l = 1:numel(OutImages)
        bin = MatOut==OutImages(l);
        bin = sum(bin,1);
        BinMat = reshape(bin',[Dx,Dy,Dz]);
        strutmp = stru;
        strutmp.dat.fname = [fileFuncT1{OutImages(l)}(1:end-4) '__OutilierVoxels.nii'];
        fprintf('%s\n',[fileFuncT1{OutImages(l)}(1:end-4) '__OutilierVoxels.nii'])
        strutmp.dat.dim = [Dx,Dy,Dz];
        strutmp.dat.dtype = 'FLOAT32-LE';
        strutmp.dat(:,:,:) = BinMat;
        create(strutmp)
    end
    
    figure;
    boxplot(s2,'Whisker',whisker);
    ylabel('HIGH <--- QUALITY ---> POOR','FontSize',14);
    mTextBox = uicontrol('style','text','HorizontalAlignment','left','position',round([1 1 450 30]));
    set(mTextBox,'String',{'Click on the outliers to check the filenames. Right Click removes the label.';'Press ''Enter'' to close this tool'});
    gname(filename2)
end
    
    
    
    
    
