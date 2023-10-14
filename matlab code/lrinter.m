function [A0,A1,A2,A3] = lrinter(gradg)
    n=size(gradg,1);
    m=size(gradg,2);
    A(1,1)=2; A(1,2)=1; A(m,m-1)=2; A(m,m)=1;
    for i=2:m-1
        A(i,i-1)=1;
        A(i,i)=4;
        A(i,i+1)=1;
    end
    [L,U]=de(A);
    for j=1:n
        y=gradg(j,:); 
        [a0,a1,a2,a3]=inter(y,L,U);
        A0(j,:)=a0;
        A1(j,:)=a1;
        A2(j,:)=a2;
        A3(j,:)=a3;
    end  
end
