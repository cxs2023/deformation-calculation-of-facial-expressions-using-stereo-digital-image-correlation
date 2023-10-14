function [P,C]=intera(point,P0,deltax,deltay,hght,wdth,beimage,afimage)
    
    f=beimage(point(1)-hght:point(1)+hght,point(2)-wdth:point(2)+wdth);
    fm=mean(f(:));
    f_m=f-fm;
    funi=(f_m).^2;
    funi=sum(funi(:));
    
    c=0;

    for n=1:20

        u=P0(1);
        ux=P0(2);
        uy=P0(3);
        v=P0(4);
        vx=P0(5);
        vy=P0(6);

        %     位置
        posix=point(2)+deltax;
        posiy=point(1)+deltay;
        xafposi=posix+u+ux*deltax+uy*deltay;
        yafposi=posiy+v+vx*deltax+vy*deltay;

        %插值边界
        sx=floor(min(xafposi(:)))-5;
        ex=ceil(max(xafposi(:)))+5;
        sy=floor(min(yafposi(:)))-5;
        ey=ceil(max(yafposi(:)))+5;

        %插值
        interg=afimage(sy:ey,sx:ex);
        [A0,A1,A2,A3] = lrinter(interg);
        [A_0,A_1,A_2,A_3] = udinter(interg);
        g=zeros(41);
        for i=1:2*wdth+1
            for j=1:2*hght+1
                x=xafposi(j,i)-sx+1;
                y=yafposi(j,i)-sy+1;
                [a0,a1,a2,a3] = edinter(x,y,A0,A1,A2,A3);
                g(j,i)=a0+a1*(y-floor(y))+a2*(y-floor(y))^2+a3*(y-floor(y))^3;     %灰度插值
                dg_y(j,i)=a1+2*a2*(y-floor(y))+3*a3*(y-floor(y))^2;                      %y-方向灰度梯度插值
                [a_0,a_1,a_2,a_3] = ed_inter(x,y,A_0,A_1,A_2,A_3);
                dg_x(j,i)=a_1+2*a_2*(x-floor(x))+3*a_3*(x-floor(x))^2;                  %x-方向灰度梯度插值
            end
        end

        gm=mean(g(:));
        g_m=g-gm;
        guni=(g_m).^2;
        guni=sum(guni(:));

        C=((f_m)./sqrt(funi)).*((g_m)./sqrt(guni));
        C=-sum(C(:));
        
%         if abs(c-C)<0.0001
%             break
%         end

        DC=zeros(6,1);
        DDC_up=zeros(6,6);

        DC1=g_m.*dg_x;
        DC1=(sum(DC1(:)))/guni;
        DC2=f_m.*dg_x;
        DC2=(sum(DC2(:)))/sqrt(guni.*funi);
        DC_u=-C*DC1-DC2;
        DC(1)=DC_u;

        DC1=g_m.*dg_x.*deltax;
        DC1=(sum(DC1(:)))/guni;
        DC2=f_m.*dg_x.*deltax;
        DC2=(sum(DC2(:)))/sqrt(guni*funi);
        DC_ux=-C*DC1-DC2;
        DC(2)=DC_ux;

        DC1=g_m.*dg_x.*deltay;
        DC1=(sum(DC1(:)))/guni;
        DC2=f_m.*dg_x.*deltay;
        DC2=(sum(DC2(:)))/sqrt(guni*funi);
        DC_uy=-C*DC1-DC2;
        DC(3)=DC_uy;

        DC1=g_m.*dg_y;
        DC1=(sum(DC1(:)))/guni;
        DC2=f_m.*dg_y;
        DC2=(sum(DC2(:)))/sqrt(guni*funi);
        DC_v=-C*DC1-DC2;
        DC(4)=DC_v;

        DC1=g_m.*dg_y.*deltax;
        DC1=(sum(DC1(:)))/guni;
        DC2=f_m.*dg_y.*deltax;
        DC2=(sum(DC2(:)))/sqrt(guni*funi);
        DC_vx=-C*DC1-DC2;
        DC(5)=DC_vx;

        DC1=g_m.*dg_y.*deltay;
        DC1=(sum(DC1(:)))/guni;
        DC2=f_m.*dg_y.*deltay;
        DC2=(sum(DC2(:)))/sqrt(guni*funi);
        DC_vy=-C*DC1-DC2;
        DC(6)=DC_vy;

        DDC_u_u=dg_x.*dg_x;
        DDC_u_u=(sum(DDC_u_u(:)))/guni;
        DDC_u_u=C*DDC_u_u;
        DDC_up(1,1)=DDC_u_u;

        DDC_u_ux=dg_x.*dg_x.*deltax;
        DDC_u_ux=(sum(DDC_u_ux(:)))/guni;
        DDC_u_ux=C*DDC_u_ux;
        DDC_up(1,2)=DDC_u_ux;

        DDC_u_uy=dg_x.*dg_x.*deltay;
        DDC_u_uy=(sum(DDC_u_uy(:)))/guni;
        DDC_u_uy=C*DDC_u_uy;
        DDC_up(1,3)=DDC_u_uy;

        DDC_u_v=dg_x.*dg_y;
        DDC_u_v=(sum(DDC_u_v(:)))/guni;
        DDC_u_v=C*DDC_u_v;
        DDC_up(1,4)=DDC_u_v;

        DDC_u_vx=dg_x.*dg_y.*deltax;
        DDC_u_vx=(sum(DDC_u_vx(:)))/guni;
        DDC_u_vx=C*DDC_u_vx;
        DDC_up(1,5)=DDC_u_vx;

        DDC_u_vy=dg_x.*dg_y.*deltay;
        DDC_u_vy=(sum(DDC_u_vy(:)))/guni;
        DDC_u_vy=C*DDC_u_vy;
        DDC_up(1,6)=DDC_u_vy;

        DDC_ux_ux=dg_x.*deltax.*dg_x.*deltax;
        DDC_ux_ux=(sum(DDC_ux_ux(:)))/guni;
        DDC_ux_ux=C*DDC_ux_ux;
        DDC_up(2,2)=DDC_ux_ux;

        DDC_ux_uy=dg_x.*deltax.*dg_x.*deltay;
        DDC_ux_uy=(sum(DDC_ux_uy(:)))/guni;
        DDC_ux_uy=C*DDC_ux_uy;
        DDC_up(2,3)=DDC_ux_uy;

        DDC_ux_v=dg_x.*deltax.*dg_y;
        DDC_ux_v=(sum(DDC_ux_v(:)))/guni;
        DDC_ux_v=C*DDC_ux_v;
        DDC_up(2,4)=DDC_ux_v;

        DDC_ux_vx=dg_x.*deltax.*dg_y.*deltax;
        DDC_ux_vx=(sum(DDC_ux_vx(:)))/guni;
        DDC_ux_vx=C*DDC_ux_vx;
        DDC_up(2,5)=DDC_ux_vx;

        DDC_ux_vy=dg_x.*deltax.*dg_y.*deltay;
        DDC_ux_vy=(sum(DDC_ux_vy(:)))/guni;
        DDC_ux_vy=C*DDC_ux_vy;
        DDC_up(2,6)=DDC_ux_vy;

        DDC_uy_uy=dg_x.*deltay.*dg_x.*deltay;
        DDC_uy_uy=(sum(DDC_uy_uy(:)))/guni;
        DDC_uy_uy=C*DDC_uy_uy;
        DDC_up(3,3)=DDC_uy_uy;

        DDC_uy_v=dg_x.*deltay.*dg_y;
        DDC_uy_v=(sum(DDC_uy_v(:)))/guni;
        DDC_uy_v=C*DDC_uy_v;
        DDC_up(3,4)=DDC_uy_v;

        DDC_uy_vx=dg_x.*deltay.*dg_y.*deltax;
        DDC_uy_vx=(sum(DDC_uy_vx(:)))/guni;
        DDC_uy_vx=C*DDC_uy_vx;
        DDC_up(3,5)=DDC_uy_vx;

        DDC_uy_vy=dg_x.*deltay.*dg_y.*deltay;
        DDC_uy_vy=(sum(DDC_uy_vy(:)))/guni;
        DDC_uy_vy=C*DDC_uy_vy;
        DDC_up(3,6)=DDC_uy_vy;

        DDC_v_v=dg_y.*dg_y;
        DDC_v_v=(sum(DDC_v_v(:)))/guni;
        DDC_v_v=C*DDC_v_v;
        DDC_up(4,4)=DDC_v_v;

        DDC_v_vx=dg_y.*dg_y.*deltax;
        DDC_v_vx=(sum(DDC_v_vx(:)))/guni;
        DDC_v_vx=C*DDC_v_vx;
        DDC_up(4,5)=DDC_v_vx;

        DDC_v_vy=dg_y.*dg_y.*deltay;
        DDC_v_vy=(sum(DDC_v_vy(:)))/guni;
        DDC_v_vy=C*DDC_v_vy;
        DDC_up(4,6)=DDC_v_vy;

        DDC_vx_vx=dg_y.*deltax.*dg_y.*deltax;
        DDC_vx_vx=(sum(DDC_vx_vx(:)))/guni;
        DDC_vx_vx=C*DDC_vx_vx;
        DDC_up(5,5)=DDC_vx_vx;

        DDC_vx_vy=dg_y.*deltax.*dg_y.*deltay;
        DDC_vx_vy=(sum(DDC_vx_vy(:)))/guni;
        DDC_vx_vy=C*DDC_vx_vy;
        DDC_up(5,6)=DDC_vx_vy;

        DDC_vy_vy=dg_y.*deltay.*dg_y.*deltay;
        DDC_vy_vy=(sum(DDC_vy_vy(:)))/guni;
        DDC_vy_vy=C*DDC_vy_vy;
        DDC_up(6,6)=DDC_vy_vy;

        DDC=DDC_up+DDC_up';
        for i=1:6
            DDC(i,i)=DDC(i,i)/2;
        end
        DDC=-DDC;

        P=P0-(inv(DDC)*DC);

        if abs(P(1)-P0(1))<0.001&&abs(P(4)-P0(4))<0.001
            break
        end

        P0=P;
        c=C;
        n;
        
    end
            
end