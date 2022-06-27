function Rm=rs_map(Veg,Beta)
M=8; [rows,cols]=size(Veg);
NBR  = [-1,  0,  1, -1, 1, -1, 0, 1;      % neighbor i-coordinate
        -1, -1, -1,  0, 0,  1, 1, 1];     % neighbor j-coordinate
for j=1:cols
    for i=1:rows
        k=(j-1)*rows+i; nbi=NBR(1,:)+i; nbj=NBR(2,:)+j;
        for m=1:M 
            ii=nbi(m); jj=nbj(m);
            if ii>=1 && ii<=rows && jj>=1 && jj<=cols && Veg(ii,jj) == 4
                Beta(k,m)=0;
            end
        end
    end
end
Rm=(Beta > 0);
end