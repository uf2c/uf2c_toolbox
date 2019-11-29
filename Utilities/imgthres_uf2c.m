function [finalEPI,histCut] = imgthres_uf2c(finalEPI,Dx,Dy,Dz,ScreSize,pathFunc,dirname)

    UDV = uf2c_defaults('imgthres_uf2c');
    
    fprintf('\r\n')
    fprintf('UF²C =============== %s\n',spm('time'));
    fprintf('Reshapeing and Thresholding\n');

    THREimg = mean(finalEPI,4);
    THREimgresh = reshape(THREimg,prod([Dx,Dy]),Dz)';
    
    THREimgresh(THREimgresh==0) = NaN;

    try %condition to old Matlab versions
        minVV = min(min(THREimgresh));
        maxVV = max(max(THREimgresh));
        nbins = round((maxVV-minVV)/5);
        [Nxs,edges,binxs] = histcounts(THREimgresh,nbins);
        [Mxs,Ixs] = max(Nxs);
        HigherFreq1 = edges(Ixs);
        HigherFreq2 = edges(Ixs+1);
        HigherFreq = (HigherFreq2+HigherFreq1)/2;
        CutoffValue = round(UDV.ITP*HigherFreq); %% Threshold: The cutoff intensity is 33% of the most frequent intensity
    catch
        disp('--- Old Matlab Version: Alternative code\r\n')
        THREimgresh = reshape(THREimgresh,prod([Dx,Dy,Dz]),1)';
        minVV = min(THREimgresh);
        maxVV = max(THREimgresh);
        nbins = ceil((maxVV-minVV)/5);
        [Nxs,edges] = hist(THREimgresh,nbins);
        [Mxs,Ixs] = max(Nxs);
        HigherFreq1 = edges(Ixs);
        HigherFreq2 = edges(Ixs+1);
        HigherFreq = (HigherFreq2+HigherFreq1)/2;
        CutoffValue = round(HigherFreq*UDV.ITP);
    end

    histCut = figure('Visible','off');
    set(histCut,'Name','Average Image Histogram',...
        'Position', round([ScreSize(1)*.1 ScreSize(2)*.1 ScreSize(1)*.4 ScreSize(1)*.3]),...
            'Color',[1 0.94 0.86]);

    try %condition to old Matlab versions
        histogram(THREimgresh,nbins);
    catch
        bar((minVV:5:maxVV),Nxs,'hist');
    end

    xlabel('Intensity'  ,'FontSize',14);
    ylabel('Frequency','FontSize',14);
    hold on
    ylimI = get(gca,'YLim');
    line([CutoffValue,CutoffValue],ylimI,'Color',[0 0 0],'LineStyle','-','LineWidth',0.5)
    str1 = ('\leftarrow');
    str2 = sprintf('Cutoff Intensity:\n %d',CutoffValue);
    strFi = [str1 str2];
    text(CutoffValue,ylimI(2)/2,strFi)
    title('Average Image Histogram','FontSize', 14);
    drawnow

    clear minVV maxVV str1 str2 strFi histCutoff THREimgresh Nxs edges THREimg

    finalEPI(finalEPI<CutoffValue) = 0;
    fprintf('Done! ============== %s\r\n',spm('time'));    
end