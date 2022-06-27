function A=md_linear(Beta,Delta,rows,cols,h)
M=8;
NBR=[-1,  0,  1, -1, 1, -1, 0, 1;      % neighbor i-coordinate
     -1, -1, -1,  0, 0,  1, 1, 1];     % neighbor j-coordinate
n=rows*cols; beta=zeros(n,n);
for j=1:cols
    for i=1:rows
        k=(j-1)*rows+i; nbi=NBR(1,:)+i; nbj=NBR(2,:)+j;
        for m=1:M
            if Beta(k,m) > 0
                ii=nbi(m); jj=nbj(m); l=(jj-1)*rows+ii;
                beta(l,k)=Beta(k,m);
            end
        end
    end
end
ind=find(Delta>0); dd=Delta(ind); bb=beta(ind,ind);
A=diag(1-h*dd)+h*bb;
end