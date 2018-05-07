function ShowGrazingTime(Data,Time,Idx,Title)

[A,R] = geotiffread('CarteBase.tif');
U = unique(Data.Idx);
U = U(Idx);

for i = 1:size(Time,1)
    
    H = Time(i,:);
    
    J = U(H>0)';
    c = H(H>0);
    x = [Data.Cx1(J) ; Data.Cx2(J) ; Data.Cx3(J) ; Data.Cx4(J)];
    y = [Data.Cy1(J) ; Data.Cy2(J) ; Data.Cy3(J) ; Data.Cy4(J)];
    h= figure('Name',Title);
    hold on
    axis image
    mapshow(A,R)
    patch(x,y,c,'EdgeColor','w')
    C = rot90(rot90(colormap(gray(50))));
    colormap(C)
    caxis([0 max(max(Time))])
    colorbar
    axis image
    truesize(h)
end