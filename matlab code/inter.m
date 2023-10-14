function [a0,a1,a2,a3] = inter(y,L,U)
    n=length(y);
    dy(1)=y(2)-y(1);
    for i=2:n-1
        dy(i)=y(i+1)-y(i-1);
    end
    dy(n)=y(n)-y(n-1);
    dy=3.*dy;
    z(1)=dy(1);
    for i=2:n
        z(i)=dy(i)-L(i,i-1)*z(i-1);
    end
    diffy(n)=z(n)/U(n,n);
    for i=1:n-1
        diffy(n-i)=(z(n-i)-U(n-i,n-i+1)*diffy(n-i+1))/U(n-i,n-i);
    end
    for k=1:n-1
        a0(k)=y(k); a1(k)=diffy(k);
        a2(k)=3*(y(k+1)-y(k))-2*diffy(k)-diffy(k+1);
        a3(k)=diffy(k)+diffy(k+1)-2*(y(k+1)-y(k));
    end    
end