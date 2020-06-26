%Same code as in /Pulled/Parameter1/spec_wgt.m but with different values of alpha and gamma


clear all
close all
ind = 270;
sk = 2;

SOL = cell(ind,1);

for number=1:ind
    file_name = sprintf('AUTO_Existence/parallelsol_%d.dat',number);
    SOL{number} = load(file_name);
end

BIF = load('AUTO_Existence/parallelbif.dat');
C = BIF(:,1);
W = BIF(:,end);
CDAT = C(1:sk:ind*sk);
WDAT = W(sk:sk:ind*sk);

%%Concatenate two sides of solutions
SOLC = cell(ind,1);

for jj = 1:ind
    S = SOL{jj};
    M = length(S(:,1));
    SC = zeros(2*M,4);
    SC(:,1) = [S(:,1) - 1; S(:,1)];
    for ii = 2:4
        SC(:,ii) = [S(:,ii);S(:,ii+3)];
    end
    SOLC{jj} = SC;
end


%%%%%%%%Study spectra in the Beck et. al. coordinates
%%Read off length in c.het.1 = L = 500
%%Implies the domain is now [-500, 500]

%%%Overall process, interpolate solutions (z,r) onto even grid, convert into Beck variables,
%%%then calculate spectrum of discretized linearization with Neumann Boundary Conditions


%%Code to find spectrum at specific data point in CDAT
a = -0.1;
g = -0.9; %%%read from auto code,
LL = 500;  %interpolation width,   note, if not the same as AUTO domain size, may lose the zero eigenvalue from gauge symmetry
n = 20000;   %number of grid points
x = linspace(-LL,LL,n);
dx = x(2) - x(1);
mu = 2*(heaviside(-x)-0.5)'; %%inhomogeneity
count = 200;  %%how many eigenvalues to look for
eta = 0.6;  %exponential weight for spectral problem

clin = 2*sqrt(1+a^2);

Lamda = zeros(count,ind);
LambdaF = zeros(2*n,ind);
Eigfunc = cell(ind,1);

%%Form derivative matrices
E = ones(n,1);
%D = spdiags([-E/2 0*E E/2],[-1:1],n,n);D(1,2) = 0; D(n,n-1) = 0;
%D(1,2) = 0; D(n,n-1) = 0;
%%D(1,1) = -1/2; D(n,n) = 1/2;

D = spdiags([(1/12)*E (-2/3)*E 0*E (2/3)*E (-1/12)*E],-2:2,n,n);
D(1,2) = 0; D(1,3) = 0; D(2,4) = 0;
D(n,n-1) = 0;D(n,n-2) = 0;D(n-1,n-3) = 0;


D = D/dx;

%D2 = spdiags([E -2*E E],[-1:1],n,n);D2(1,2) = 2; D2(n,n-1) = 2;
%D2(1,1) = -1; D2(n,n) = -1;


D2 = spdiags([(-1/12)*E (4/3)*E (-5/2)*E (4/3)*E (-1/12)*E],-2:2,n,n);
D2(1,2) = 8/3;D2(1,3) = -2/12;D2(2,4) = -2/12;
D2(n,n-1) = 8/3; D2(n,n-2)=-2/12; D2(n-1,n-3) = -2/12;

%D2(1,1) =-dx; D2(1,2) = dx;
%D2(n,n-1) = -dx; D2(n,n) = dx;
D2 = D2/(dx^2);



for jj =1:ind%1:ind

IND = jj;
S = SOLC{IND};
c = CDAT(IND)
w = WDAT(IND)


%  XI = L*SOLC{IND}(:,1);
%  KAP = SOLC{IND}(:,2); %%real part of z
%  Q = SOLC{IND}(:,3);  %%imaginary part of z
%  R = SOLC{IND}(:,4); %%amplitude r
m1 = length(S(:,1));

SOLS = zeros(n,3);
%Interpolate the data onto a uniform grid on [-LL,LL]
for ii = 2:4
    [XX, index] = unique(LL*S(:,1));
    YY = interp1(XX, S(index,ii), x);

 %   F = interp1(L*S(:,1),S(:,ii),splines);
    SOLS(:,ii-1) = YY;%F(x);
end
kap = SOLS(:,1);
q = SOLS(:,2);
r = real(sqrt(SOLS(:,3)));

zz = kap + 1i*q;
zz_h = zz+c/2/(1+1i*a);

figure(3)
plot(x,r)


%%Calculate needed quatinities
%dq = D*q;
dq=-2*kap.*q + (w+a*c*kap-c*q+a*mu+(g-a)*r.^2)/(1+a^2);
%dkap = D*kap;
dkap=q.^2-kap.^2 +(a*w-a*c*q-c*kap-mu+(1+a*g)*r.^2)/(1+a^2);




figure(22)
plot(x,[dq./kap,dkap./kap])

%Calculate matrix multipliers
A2 = kron([1, -a; a, 1],D2);

A1 = kron(c*speye(2,2),D) - 2*kron([a, 1; -1, a],spdiags(q,0,n,n)*D) - 2*eta*kron([1, -a; a, 1],D);


A011 = spdiags([-q.^2-a*dq+mu-3*r.^2],0,n,n);
A012 = spdiags([a*(dkap+kap.^2)+2*q.*kap],0,n,n);
A021 = spdiags([-w-a*q.^2+dq-3*g*r.^2+c*q],0,n,n);
A022 = spdiags([-(dkap+kap.^2)+2*a*q.*kap-c*kap],0,n,n);
A0 = [A011, A012; A021, A022];

A0 = A0 + eta^2*kron([1, -a; a, 1],speye(n));
A0 = A0 - eta*(kron(c*speye(2,2),speye(n)) - 2*kron([a, 1; -1, a],spdiags(q,0,n,n)));


%L = A2 + A1 + A0;

L = A2 + A1 +A0;

condest(L)

%%Find evals close to zero
tic
[V2,LAM2] = eigs(L,count,0.1);%,'Tolerance',1e-17);
%[V2,LAM2] = eigs(L,count,'largestreal');
cpu_time_it = toc
DLAM2 = diag(LAM2);


[DLAM2,IS] = sort(DLAM2,'ComparisonMethod','real');
VS = V2(:,IS);
DLAM2 = DLAM2(end:-1:1);
VS = VS(:,end:-1:1);

Lambda(:,jj) = DLAM2;
Eigfunc{jj} = VS(:,1:25);


t = [0:0.01:10];
LABS = 1 - (t + CDAT(jj)^2/4)/(1+1i*a) - 1i*WDAT(jj);
figure(6)
plot(real(Lambda(:,jj)),imag(Lambda(:,jj)),'.')
hold on
plot(real(LABS),imag(LABS),'Color','Blue')
plot(real(LABS),-imag(LABS),'Color','Blue')
hold off
drawnow





figure(5)
subplot(4,1,1)
plot(x,real(VS(1:end/2,1))',x,real(VS(end/2+1:end,1))')
subplot(4,1,2)
plot(x,real(VS(1:end/2,2))',x,real(VS(end/2+1:end,2))')
subplot(4,1,3)
plot(x,real(VS(1:end/2,3))',x,real(VS(end/2+1:end,3))')
subplot(4,1,4)
plot(x,real(VS(1:end/2,4))',x,real(VS(end/2+1:end,4))')

figure(55)
subplot(4,1,1)
plot(x,imag(VS(1:end/2,1))',x,imag(VS(end/2+1:end,1))')
subplot(4,1,2)
plot(x,imag(VS(1:end/2,2))',x,imag(VS(end/2+1:end,2))')
subplot(4,1,3)
plot(x,imag(VS(1:end/2,3))',x,imag(VS(end/2+1:end,3))')
subplot(4,1,4)
plot(x,imag(VS(1:end/2,4))',x,imag(VS(end/2+1:end,4))')

% tic  %%%if all spectrum is desired
%gL = gpuArray(full(L));
% [DLAMF] = eig(gL);
% cpu_time = toc
% %DLAMF = diag(DLAMF);
% [DLAMF,IS] = sort(DLAMF);
% LambdaF(:,jj) = get(DLAMF);

figure(4)
plot(real(DLAM2),imag(DLAM2),'x')
%hold on
%plot(real(DLAMF),imag(DLAMF),'o')
%    hold off

end

%Plot real part of largest 6 eigenvalues

figure(2)
plot(CDAT,[real(Lambda(1:10,:))],'.-')
xlabel('c')
ylabel('Re \lambda')
drawnow

figure(3)
plot(CDAT,[imag(Lambda(1:10,:))],'.-')
xlabel('c')
ylabel('Im \lambda')
drawnow


% figure(4)
% plot(CDAT,[real(LambdaF(1:10,:))],'.-')
% xlabel('c')
% ylabel('Im \lambda')
% drawnow

%%If one wants to see spec for a specific speed
% PI = 11;
% figure(4)
% plot(real(Lambda(:,PI)),imag(Lambda(:,PI)),'x')

save('spec_data_wgtfull_newparam.mat', 'Lambda','CDAT','WDAT','a','g','x','n','LL')
