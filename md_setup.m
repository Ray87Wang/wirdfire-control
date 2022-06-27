function [rows,cols,Veg,Beta,Delta] = md_setup()

% parameter setup, see above refs
c1=0.045; c2=0.131; sw=2*(sqrt(2)-1); sa=0.078; sb=sqrt(2); l=10; 
V=4;            % wind speed (m/s)
theta_w=-3*pi/4;% wind direction (northeast)
delta=0.5;      % baseline recovery rate
beta=0.5;       % baseline spreading rate

% veg: 1(forest), 2(grassland), 3(desert), 4(city), 5(water)
veg_table=[1.4, 1.0, 0.1, 0.5, 0.0];

% 8 nodes
M=8;
NBR  = [-1,  0,  1, -1, 1, -1, 0, 1;      % neighbor i-coordinate
        -1, -1, -1,  0, 0,  1, 1, 1];     % neighbor j-coordinate
Theta= [-1,  0,  1, -2, 2, -3, 4, 3]*pi/4; %fire propagation angle

Veg=csvread('Vegetation.csv');
E=csvread('Elevation.csv');
[rows,cols]=size(Veg); N=rows*cols;
Beta=zeros(N,M);    

for j=1:cols
    for i=1:rows
        k=(j-1)*rows+i; nbi=NBR(1,:)+i; nbj=NBR(2,:)+j;
        % spreading rate based on vegetation, wind and slope
        for m=1:M 
            pveg=0; pw=0; ps=0;
            ii=nbi(m); jj=nbj(m);
            if ii>=1 && ii<=rows && jj>=1 && jj<=cols && ...
                    Veg(ii,jj) < 5 && Veg(i,j) < 5
                pveg=veg_table(Veg(ii,jj));
                pw=exp(V*(c1+c2*(cos(Theta(m)-theta_w)-1)));
                slope=(E(i,j)-E(ii,jj))/l; 
                if ii~=i && jj~=j % diagonal
                    pw=pw*sw; slope=slope/sb;
                end
                ps=exp(sa*atan(slope));
            end
            Beta(k,m)=min(max(beta*pveg*pw*ps,0),1); % proj to [0,1]
        end
    end
end

Delta=delta*(sum(Beta,2) > 0);
end