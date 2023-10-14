function begin=Initialize(point,beimage,afimage,hght,wdth)

%point(1)---y coordinate;
%point(2)---x coordinate;
sub=beimage(point(1)+hght:-1:point(1)-hght,point(2)+wdth:-1:point(2)-wdth);
cf=conv2(afimage,sub,'same');

psub=ones(size(sub));
biaf=afimage.^2;
nm=conv2(biaf,psub,'same');
nm=((sum(sum(sub.^2)))^0.5)*sqrt(nm);


cf=cf./nm;
      
[y,x]=pm(cf);
begin=zeros(1,2);
begin(1)=y-point(1);
begin(2)=x-point(2);

function [y,x]=pm(A)
[u,y]=max(A);
[u,x]=max(u);
y=y(x);