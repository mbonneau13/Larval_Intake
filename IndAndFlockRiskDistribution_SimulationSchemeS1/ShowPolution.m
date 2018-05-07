function ShowPolution(Data,H,Idx,Title)

[A,R] = geotiffread('CarteBase.tif');
U = unique(Data.Idx);
U = U(Idx);


J = U(H>0)';
c = H(H>0);
x = [Data.Cx1(J) ; Data.Cx2(J) ; Data.Cx3(J) ; Data.Cx4(J)];
y = [Data.Cy1(J) ; Data.Cy2(J) ; Data.Cy3(J) ; Data.Cy4(J)];
h=figure('Name',Title);
hold on
% axis image
mapshow(A,R)
% patch(x,y,c,'EdgeColor','w')
scatterbar3(mean(x,1),mean(y,1),c,1)
C = colormap(jet(50));
colormap(C)
caxis([0 max(max(H))])
colorbar
% axis image
% truesize(h)