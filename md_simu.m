function x1=md_simu(x0,h,Beta,Delta,rows,cols)
M=8; Delta=reshape(Delta,[rows cols]); xb=0*x0;
NBR=[-1,  0,  1, -1, 1, -1, 0, 1;      % neighbor i-coordinate
     -1, -1, -1,  0, 0,  1, 1, 1];     % neighbor j-coordinate
for j=1:cols
    for i=1:rows
        k=(j-1)*rows+i; xn=zeros(M,1);
        nbi=NBR(1,:)+i; nbj=NBR(2,:)+j;
        for m=1:M
            ii=nbi(m); jj=nbj(m);
            if ii>=1 && ii<=rows && jj>=1 && jj<=cols
                xn(m)=x0(ii,jj);
            end
        end
        xb(i,j)=Beta(k,:)*xn;
    end
end
x1=x0+h*((1-x0).*xb-Delta.*x0);
end