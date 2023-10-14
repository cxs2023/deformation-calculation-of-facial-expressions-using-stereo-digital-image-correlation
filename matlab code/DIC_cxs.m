function [acu,acv,coef,pointx,pointy]=DIC_cxs(p,beimage,afimage)

dbstop if error
    
    N=size(p,1);
    startpt=p(1,:);

    hght=30;
    wdth=30;
    
    deltax=zeros(2*hght+1,2*wdth+1);
    deltay=zeros(2*hght+1,2*wdth+1);
    
    for i=1:2*wdth+1
        deltax(:,i)=i-(wdth+1);
    end
    for i=1:2*hght+1
        deltay(i,:)=i-(hght+1);
    end
    
    vu=Initialize(startpt,beimage,afimage,hght,wdth);
    v=vu(1);
    u=vu(2);
    ux=0;
    uy=0;
    vx=0;
    vy=0;
    
    P0=[u ux uy v vx vy]';

    acu=zeros(1,N);
    acv=zeros(1,N);
    coef=zeros(1,N);
    
    p_=p;
    p_(1,:)=[];
    p_re=zeros(N,2);
    p_re(1,:)=p(1,:);
    
    C=-1;
    D=0;
    point=startpt;
    
    for num=1:N 
               
            dis=[];
            if num>1
                for i=1:N-num+1
                    dis(i)=norm([point-p_(i,:)],2);
                end
                [D,d]=min(dis);
                point=p_(d,:);
                p_re(num,:)=p_(d,:);
                p_(d,:)=[];
            end
            
            if C>-0.9||D>500
                    VU=Initialize(point,beimage,afimage,hght,wdth);
                P0(1)=VU(2);
                P0(4)=VU(1);
            end
            
            P0(2)=0;
            P0(3)=0;
            P0(5)=0;
            P0(6)=0;
    
            try
               [P,C]=intera(point,P0,deltax,deltay,hght,wdth,beimage,afimage);
            catch
               C=0;
            end
            
            acu(num)=P(1);
            acv(num)=P(4);
            coef(num)=C;
            acux(num)=P(2);
            acuy(num)=P(3);
            acvx(num)=P(5);
            acvy(num)=P(6);
            
            num;
            
            P0=P;
            
    end
   
    acu=acu';
    acv=acv';
    
    pointx=p_re(:,1);
    pointy=p_re(:,2);
    

end

