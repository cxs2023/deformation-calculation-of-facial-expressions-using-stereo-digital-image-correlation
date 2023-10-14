function strain=strain_obtain(PtOI, stereotrans)
global wsp;
global hsp;
Num = size( stereotrans,1);    
for i=1:Num    
    clear coef_abc
    clear rp
    clear coef_XY
    clear stereotrans_uxy
    clear PtOI_xyz
    clear PtOI_s
    clear stereotrans_u
    clear stereotrans_s
    Point_near=NaN(Num,200);
    mn=0;
    for k=1:Num
        distance(i,k)=sqrt((PtOI(i,1)-PtOI(k,1))^2+(PtOI(i,2)-PtOI(k,2))^2+(PtOI(i,3)-PtOI(k,3))^2);
        if distance(i,k)<1.5*(wsp+hsp)/2
            mn=mn+1;
            coef_XY(mn,:)=PtOI(k,:);
            stereotrans_uxy(mn,:)=stereotrans(k,:);
            Point_near(i,mn)=k;
        end
    end
    %%%regress的核心也就是最小二乘了
    [coef_abc,~,rp] = regress(ones(size(coef_XY,1),1),coef_XY);
 %%   [seque,numb]=sort(distance(i,:));
 %%   if (numb(1)-numb(2)==1||-1)&&(numb(1)-numb(3)==1||-1)&&(numb(1)-numb(2)==yq||-yq)&&(numb(1)-numb(3)==yq||-yq)
 %%       coef_XY=[PtOI(numb(1),:);PtOI(numb(2),:);PtOI(numb(4),:)];
 %%       stereotrans_uxy=[stereotrans(numb(1),:);stereotrans(numb(2),:);stereotrans(numb(4),:)];
 %%   else
 %%       coef_XY=[PtOI(numb(1),:);PtOI(numb(2),:);PtOI(numb(3),:)];
 %%       stereotrans_uxy=[stereotrans(numb(1),:);stereotrans(numb(2),:);stereotrans(numb(3),:)];
 %%   end
 %%   coef_B=[1;1;1];
    %%%如果要控制精度，可以使用vpa函数;
    %%x0=[0.5;0.5;0.5];
    %%coef_abc=gauseidel(coef_XY,coef_B,x0,1e-4,200);
%%    coef_abc=linsolve(coef_XY,coef_B);%%%有问题
%%    coef_shishi(i,:)=coef_abc';
    
    coorn(i,:)=sqrt(1/((coef_abc(1))^2+(coef_abc(2))^2+(coef_abc(3))^2))*[coef_abc(1),coef_abc(2),coef_abc(3)];
    coort1(i,:)=sqrt(1/((coorn(i,1))^2+(coorn(i,3))^2))*[coorn(i,3),0,-coorn(i,1)];
    coort_2(i,:)=[coorn(i,2)*coort1(i,3)-coorn(i,3)*coort1(i,2),coorn(i,3)*coort1(i,1)-coorn(i,1)*coort1(i,3),coorn(i,1)*coort1(i,2)-coorn(i,2)*coort1(i,1)];
    coort2(i,:)=sqrt(1/((coort_2(i,1))^2+(coort_2(i,2))^2+(coort_2(i,3))^2))*[coort_2(i,1),coort_2(i,2),coort_2(i,3)];
    coortn=[coort1(i,:);coort2(i,:);coorn(i,:)];
    %%%坐标变换(只是转动，未平移)
    PtOI_xyz=coortn*(coef_XY');
    %%PtOI_s=PtOI_xyz';
    %%PtOI_s(:,3)=1;
    PtOI_s=PtOI_xyz';
    PtOI_s(:,3)=1;
    %%%位移坐标变换
    
    
    stereotrans_u=coortn*(stereotrans_uxy');
    stereotrans_s=stereotrans_u';
    stereotrans_s(:,3)=[];
    stereotrans_jubu(i,:)=coortn*(stereotrans(i,:)');
    
    %%%求系数
    %%%coef_s=linsolve(PtOI_s,stereotrans_s);%%%有问题,应该寻找更为适合的超静定方程组求解方法
    %%%用最小二乘法
    [coef_sa,~,rp] = regress(stereotrans_s(:,1),PtOI_s);
    [coef_sb,~,rp] = regress(stereotrans_s(:,2),PtOI_s);
    coef_s(:,1)=coef_sa;
    coef_s(:,2)=coef_sb;
    %%%求应变(小变形)
    strain_x(i)=coef_s(1,1);
    strain_y(i)=coef_s(2,2);
    strain_xy(i)=coef_s(1,2)+coef_s(2,1);
% % %     %%%求应变(大变形)
% % %     strain_x(i)=coef_s(1,1)+1/2*((coef_s(1,1))^2+(coef_s(2,1))^2);
% % %     strain_y(i)=coef_s(2,2)+1/2*((coef_s(1,2))^2+(coef_s(2,2))^2);
% % %     strain_xy(i)=coef_s(1,2)+coef_s(2,1)+coef_s(1,1)*coef_s(2,2)+coef_s(1,2)*coef_s(2,1);   
    
    
    %计算主应变
    strain_1(i)=(strain_x(i)+strain_y(i))/2+sqrt((strain_x(i)-strain_y(i)).^2+strain_xy(i).^2)/2;
    strain_2(i)=(strain_x(i)+strain_y(i))/2-sqrt((strain_x(i)-strain_y(i)).^2+strain_xy(i).^2)/2;
    if strain_x(i)<strain_y(i)
       arfa(i)=atand(strain_xy(i)/(strain_y(i)-strain_x(i)))/2;
    else
       arfa(i)=atand(strain_xy(i)/(strain_y(i)-strain_x(i)))/2+90; 
    end
    if arfa(i)<0
        arfa(i)=arfa(i)+180;
    end
end
strain=[strain_x, strain_y, strain_xy, strain_1, strain_2, arfa];
end