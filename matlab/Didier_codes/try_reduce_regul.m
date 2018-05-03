% try qr decomposition on stacked derivatives

% 2 dimensions
ds1=get_l(60,2);
dt1=get_l(14,1);
ds2=kron(ds1,speye(14));
spy(ds2)
dt2=kron(speye(60),dt1);
spy(dt2)
dst2=[ds2;dt2];
spy(dst2)
spy(qr(dst2))
nnz(dst2),nnz(qr(dst2))
p=colamd(dst2);
nnz(dst2),nnz(qr(dst2)),nnz(qr(dst2(:,p)))
spy(qr(dst2(:,p)))
subplot(1,3,1),spy(dst2),title('original')
subplot(1,3,2),spy(qr(dst2)),title('qr')
subplot(1,3,3),spy(qr(dst2(:,p))),title('qr colamd')

% 3 dimensions: 
%Nphi=12,Ntheta=6,Nt=7
clear
Nphi=120,Ntheta=60,Nt=30
dtheta1=get_derivatives(Ntheta,2);
dphi1=get_derivatives(Nphi,2,1);%periodic in phi
dtheta2=kron(speye(Nphi),dtheta1);
spy(dtheta2)
dphi2=kron(dphi1,speye(Ntheta));
spy(dphi2)
d2=[dtheta2;dphi2];
spy(d2)

p=colamd(d2);
nnz(d2),nnz(qr(d2)),nnz(qr(d2(:,p)))
subplot(1,3,1),spy(d2),title('original')
subplot(1,3,2),spy(qr(d2)),title('qr')
subplot(1,3,3),spy(qr(d2(:,p))),title('qr colamd')

d2r = qr(d2);
%d2r = d2r(1:(Ntheta*Nphi-4),:); %remove half - 4 (nullspace is dim4)
d2r = d2r(1:(Ntheta*Nphi-2),:); %remove half - 2 (nullspace is dim2)
y=rand(Ntheta*Nphi,1);
x=d2'\y;
[c,R]=qr(d2',y);
x=R\c;

%remove p0=[1,2,Nphi+1,Nphi+2]
p=[3:Nphi,(Nphi+3):Nphi*Ntheta];
%p=[3:Nphi,(Nphi+3):Nphi*Ntheta,1,2,Nphi+1,Nphi+2];

d2p=d2(:,p);
y=rand(size(d2p,1),1);
x=d2p\y; % Lsolve (Lp y)

x=rand(size(d2p,2),1);
y=x'/d2p;y=y';% Ltsolve (Lp' x)


dt1=get_l(Nt,1);
dt3=kron(dt1,speye(Nphi*Ntheta));
d3=kron(speye(Nt),d2);
%d3=kron(speye(Nt),d2r); not interesting
d3=[d3;dt3];
clear dphi1 dphi2 dtheta1 dtheta2 dt1 d2 dt3
permamd=colamd(d3);
d3r=qr(d3(:,permamd));
d3r=d3r(1:(Ntheta*Nphi*Nt-2),:); %remove null rows
%remove 2 linearly dependant columns
%p=[3:Nphi,(Nphi+3):Nphi*Ntheta*Nt];
p=[3:Nphi*Ntheta*Nt];

d3p=d3(:,p);
y=rand(size(d3p,1),1);
x=d3p\y; % Lsolve (Lp y)

y=rand(size(d3r,1),1);
x=d3r\y; % Lsolve (Lp y)

x=rand(size(p),1);
y=x'/d3p;y=y';% Ltsolve (Lp' x)


