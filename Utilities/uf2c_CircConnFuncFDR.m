function uf2c_CircConnFuncFDR(matrixR,limiteG,ScreSize,cc,...
                Pva,nEle,LoopAte,fact,...
                str_nEle,tarNetLog,TarNetC,...
                DirF,FigNTpx,NetLabels,TranspMin)

            
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

mTextBox21 = uicontrol('style','text',...
 'position',round([ScreSize(1)*.001 ScreSize(1)*.4375 ScreSize(1)*.3  ScreSize(1)*.0125]),...
 'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
 'ForegroundColor',[1 .949 .867],'FontWeight','bold');
set(mTextBox21,'String','Circular Connectome - DRAWING........');

mTextBox22 = uicontrol('style','text',...
 'position',round([ScreSize(1)*.001 ScreSize(1)*.422 ScreSize(1)*.2  ScreSize(1)*.0125]),...
 'FontSize',16,'BackgroundColor',[0 0 0],'HorizontalAlignment','Left',...
 'ForegroundColor',[1 .949 .867],'FontWeight','bold');
set(mTextBox22,'String',['FDR-Corrected (alpha ' Pva ')']);
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

            if tarNetLog
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

set(mTextBox21,'String','Circular Connectome');
drawnow

set(gcf,'InvertHardCopy', 'off');
drawnow
imgRR = getframe(circ);
imwrite(imgRR.cdata, [DirF filesep FigNTpx '_FDR_corrected_Circ_Map' '.png']);
clear imgRR
%             saveas(circ,[DirF filesep FigNTpx '_FDR_corrected_Circ_Map'],'fig')
close(circ)
end