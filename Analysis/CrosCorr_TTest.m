
for i = 1:86
    for j = 1:86
        [H,P,CI,tvalue] = ttest2(AllSubj3D_Ctr(i,j,:),AllSubj3D_HS(i,j,:));
        Pmap(i,j) = P;
        Tmap(i,j) = tvalue.tstat;
        SDmap(i,j) = tvalue.sd;
    end
end
