function [a0,a1,a2,a3] = edinter(x,y,A0,A1,A2,A3)

    k=floor(x);
    j=floor(y);
   
    n=5;
    m=2*n+1;
    for i=1:n+1
        a=6-i;
        y1(i)=A0(j-a,k)+A1(j-a,k)*(x-k)+A2(j-a,k)*(x-k)^2+A3(j-a,k)*(x-k)^3;
    end
    for i=n+2:m
        a=i-6;
        y1(i)=A0(j+a,k)+A1(j+a,k)*(x-k)+A2(j+a,k)*(x-k)^2+A3(j+a,k)*(x-k)^3;
    end
    
    B(1,1)=2; B(1,2)=1; B(m,m-1)=2; B(m,m)=1;
    for i=2:m-1
        B(i,i-1)=1;
        B(i,i)=4;
        B(i,i+1)=1;
    end
    
    dy1(1)=y1(2)-y1(1);
    for i=2:m-1
        dy1(i)=y1(i+1)-y1(i-1);
    end
    dy1(m)=y1(m)-y1(m-1);
    dy1=3.*dy1;
    
    [L1,U1]=de(B);
    
   z1(1)=dy1(1);
   for i=2:m
       z1(i)=dy1(i)-L1(i,i-1)*z1(i-1);
   end

   diffy(m)=z1(m)/U1(m,m);
   for i=1:m-1
       diffy(m-i)=(z1(m-i)-U1(m-i,m-i+1)*diffy(m-i+1))/U1(m-i,m-i);
   end

   k=n+1;
   a0=y1(k); 
   a1=diffy(k);
   a2=3*(y1(k+1)-y1(k))-2*diffy(k)-diffy(k+1);
   a3=diffy(k)+diffy(k+1)-2*(y1(k+1)-y1(k));
   
end