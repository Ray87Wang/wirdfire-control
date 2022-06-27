clear; close all;

%% forward simulation
[rows,cols,Veg,Beta,Delta]=md_setup();
x0=zeros(rows,cols); x0(8:10,38:40)=1;
h=0.04; K=4; Gamma=10;
Rm=rs_map(Veg,Beta); Pm=(Veg<5);
A=md_linear(Beta,Delta,rows,cols,h);
rA=max(abs(eig(A)));
alpha=1/(0.05+rA);
C=0.001*ones(rows,cols); % cost weight map
% find city nodes and set higher cost
[r,c]=find(Veg == 4);    
for k=1:length(r)
    C(r(k),c(k))=1;
end

tic
[p,r,Jb,ts]=rs_alloc(x0,Beta,Delta,Pm,Rm,K,C,h,alpha,Gamma);
toc
