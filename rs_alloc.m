function [p,rs,Jb,ts]=rs_alloc(x,Beta,Delta,Pm,Rm,K,C,h,alpha,Gamma)
tic
beta_lb=1E-8; [rows,cols]=size(x); N=rows*cols; M=8; C=C(:);
ind=find(Delta>0); L=length(ind);
Br=Beta.*Rm; [o,u,s]=find(Br); S=length(s); rub=log(s/beta_lb);
sdpvar J
y=sdpvar(L,K); r=sdpvar(S,K); 
cons=[r>=0, sum(r,1)<=Gamma, sum(r,2)<=rub]; % resource bounds
for k=1:K
    kn=min(k+1,K);
    sr=sum(r(:,1:k),2); 
    R=full(sparse(o,u,sr,N,M)); 
    Y=full(sparse(ind,ones(L,1),y(:,k),N,1));
    Yn=full(sparse(ind,ones(L,1),y(:,kn),N,1));
    for l=1:L
        n=ind(l); nbr=neighbor(n,rows); m=find(Beta(n,:));
        q0=[log(C(n))-Y(n); Yn(n)-Y(n)+log(alpha*(1-Delta(n)))];
        q1=Yn(nbr(m))-Y(n)+log(alpha*h*Beta(n,m)')-R(n,m)';
        cons=[cons, logsumexp([q0; q1]) <= 0];
    end
end
xh=x(:); xh=xh(ind)+1e-6;
cons = [cons, logsumexp(log(xh)+y(:,1)-J) <= 0];
ts.ptime=toc;
res=optimize(cons,J,sdpsettings('debug',1,'convertconvexquad',0,'savesolveroutput',1));
ts.ytime=res.yalmiptime; %yalmip time
ts.stime=res.solvertime; %solver time
ts.mtime=res.solveroutput.res.info.MSK_DINF_OPTIMIZER_TIME; %mosek time
Jb=exp(double(J)); yv=double(y); rv=double(r); 
p=zeros(rows,cols,K); rs=zeros(N,M,K);
for k=1:K
    Yv=full(sparse(ind,ones(L,1),yv(:,k),N,1)); pv=exp(Yv);
    p(:,:,k)=reshape(pv,rows,cols).*Pm;
    rs(:,:,k)=full(sparse(o,u,rv(:,k),N,M));
end
end

function nbr=neighbor(n,rows)
NBR=[-1,  0,  1, -1, 1, -1, 0, 1;  % neighbor i-coordinate
     -1, -1, -1,  0, 0,  1, 1, 1]; % neighbor j-coordinate
j=floor(n/rows)+1; i=n-(j-1)*rows;
ii=NBR(1,:)+i; jj=NBR(2,:)+j;
nbr=(jj-1)*rows+ii;
end