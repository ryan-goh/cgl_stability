%pushed_spectrum

%%Find spectrum of PDE traveling wave in a pushed triggered front


clear all
close all
gpuon = false; %use gpu to calculate eigenvalues,
%%% note that an NVIDIA GPU is required with both the NVIDIA drivers and the CUDA toolkit installed.




ind = 300;%850%565-17; %%index in data to stop at


%%%System parameters
a = 0.3;
g = -0.2;
b =  0.4;
dsh = 17;

%Other examples, requires data of run in AUTO_Existence with these parameters

%%%/pushedrho params
%a = 0.3;
%g = -0.2;
%b = 0.2;
%dsh = 15;
%%%/pushedS parameters
%a = 0.2;
%g = -0.3;
%g = -0.3;
%b = 0.4;
%dsh = 15;



eta = 0.2; %exponential weight value
mu = 4; %cubic nonlinearity value
LL = 100;  %width of spatial grid
n = 20000; %number of grid points
x = linspace(-LL,LL,n);
dx = x(2) - x(1);




%%Loads profile data from AUTO_Existence continuation
BIF = load('AUTO_Existence/pusheddata.dat');
sk = 2;
%sk = 5; %interval of profile saves in AUTO data

C = BIF(:,6); %Speed Continuation data variables
W = BIF(:,5); %Frequency Continuation data
CDAT = C(sk:sk:end); %C values for data points
WDAT = W(sk:sk:end); %W values for data points
RHO = BIF(:,1); %Trigger interface values
RHODAT = RHO(sk:sk:end); %Trigger interface values for data points
K = BIF(:,7); %Selected wavenumber
KDAT = K(sk:sk:end); %Selected wavenumber for data points



count = 200;  %%how many eigenvalues to look for in "eigs" calculations

Lamda = zeros(count,ind-dsh);
LambdaF = zeros(2*n,ind-dsh);
Eigfunc = cell(ind,1);

%%Form derivative matrices
E = ones(n,1);
%D = spdiags([-E/2 0*E E/2],[-1:1],n,n);
%D(1,2) = 0; D(n,n-1) = 0;
%%D(1,1) = -1/2; D(n,n) = 1/2;

D = spdiags([(1/12)*E (-2/3)*E 0*E (2/3)*E (-1/12)*E], -2:2,n,n); %Fourth-order with Neumann BC's
D(1,2:3) = 0; D(2,4) = 0;
D(n,n-2:n-1)=0; D(n-1,n-3) = 0;
D = D/dx;

%D2 = spdiags([E -2*E E],[-1:1],n,n);
%%D2(1,1) = -1; D2(n,n) = -1;
%D2(1,2) = 2; D2(n,n-1) = 2;
%%D2(1,1) =-dx; D(1,2) = dx;
%%D2(n,n-1) = -dx; D2(n,n) = dx;

D2 = spdiags([(-1/12)*E (4/3)*E (-5/2)*E (4/3)*E (-1/12)*E],-2:2,n,n); %Fourth-order with Neumann BC's
D2(1,2) = 8/3; D2(1,3) = -2/12; D2(2,4) = -2/12;
D2(n,n-1) = 8/3; D2(n,n-2)=-2/12; D2(n-1,n-3) = -2/12;
D2 = D2/(dx^2);


CHI =   1-5.5*(tanh(10000*(x-0.05))+1); %Set up inhomogeneity (same as in AUTO continnuation runs)

%%Main for loop, formulates linearization as a matrix and finds eigenvalues for each data point along the bifurcation curve.

for jj=dsh:ind%:100;%ind
	number = jj;
	IND = jj-dsh+1;
%file_name = sprintf('PushedS/psol_%d.dat',number);
	%file_name = sprintf('Pushed_rho/psol_%d.dat',number);
	file_name = sprintf('AUTO_Existence/psol_%d.dat',number);
    	S = load(file_name);;
	%S = SOL{IND};
	c = CDAT(IND)
	w = WDAT(IND)
	k = KDAT(IND)
	m1 = length(S(:,1));
	XI = S(:,5);
	SOLS = zeros(n,4);
	%Interpolate solutions onto a uniform grid
	for ii = 2:5
	    %[XX, index] = unique(LAT*S(:,1));
	    YY = interp1(XI, S(:,ii), x);
		%YY = YY1(end/2-n/2-1:end/2+n/2-2);
		SOLS(:,ii-1) = YY';
	 %   F = interp1(L*S(:,1),S(:,ii),splines);
	  %  SOLS(:,ii-1) = YY;%F(x);
	end
	%(z,R) values for this solution
	kap = SOLS(:,1);
	q = SOLS(:,2);
	r = real(sqrt(SOLS(:,3)));



	%%%
	%%Calculate needed quatinities

	dkap = q.^2-kap.^2 +(a*w-a*c*q-c*kap-CHI'-(r.^2).*(mu+a*g-r.^2*(1+a*b)))./(1+a^2);
%	dkap = D*kap;


	dq = -2*kap.*q + (w+a*c.*kap-c*q+a.*CHI'-(r.^2).*(g-a*mu-(r.^2)*(b-a)))./(1+a^2);
	%dq = D*q;

	%Calculate matrix multipliers and formulate the discretized linearized operator, obtain in a similar way as in the pulled case

	A2 = kron([1, -a; a, 1],D2);

	A1 = kron(c*speye(2,2),D) - 2*kron([a, 1; -1, a],spdiags(q,0,n,n)*D) - 2*eta*kron([1,-a;a,1],D);


	A011 = spdiags([-q.^2-a*dq+CHI'+3*mu*r.^2-5*r.^4],0,n,n);
	A012 = spdiags([a*(dkap+kap.^2)+2*q.*kap],0,n,n);
	A021 = spdiags([-w-a*q.^2+dq+3*g*r.^2-5*b*r.^4+c*q],0,n,n);
	A022 = spdiags([-(dkap+kap.^2)+2*a*q.*kap-c*kap],0,n,n);

	A0 = [A011, A012; A021, A022];
	%A0 = A0 + eta^2*kron([1,-a;a,1],speye(n)) - eta*(kron(c*speye(2,2),speye(n))) - 2*kron([a,1;-1,a],spdiags(q,0,n,n));
	A0 = A0 + eta^2*kron([1, -a; a, 1],speye(n));
	A0 = A0 - eta*(kron(c*speye(2,2),speye(n)) - 2*kron([a, 1; -1, a],spdiags(q,0,n,n)));

	JJ = 30;
	D0m = [A011(JJ,JJ), A012(JJ,JJ); A021(JJ,JJ), A022(JJ,JJ)]
	D0p = [A011(n-30,n-30), A012(n-30,n-30); A021(n-30,n-30), A022(n-30,n-30)]

	D1m = c*speye(2) - 2*q(JJ)*[a, 1; -1, a];
	D1p = c*speye(2) - 2*q(n-30)*[a, 1; -1, a];
	%l = [-10:0.1:10];
	l = -0.0;
	Lm = -(l).^2*[1, -a; a, 1] + D1m*(1i*l) + D0m;
	Lp = -(l).^2*[1, -a; a, 1] + D1p*(1i*l) + D0p;
	EGm = eigs(Lm)
	EGp = eigs(Lp)


	L = A2 + A1 + A0; %linear operator
	condnum = condest(L)  %print the condition number for a quick check

	if gpuon  %%formulates L, a sparse matrix, as a full matrix and uses "eig" to solve
	    fL = full(L);
	    gL = gpuArray(fL);
	    tic
	    [gV2,gLAM2] = eig(gL);
	    gpu_time_it = toc
	    LAM2 = gather(gLAM2);
	    V2 = gather(gV2);

	else
	    tic
	    [V2,LAM2] = eigs(L,count,0.4,'Tolerance',1e-17);
	    %[V2,LAM2] = eigs(L,count,'largestreal');
	    cpu_time_it = toc



	end

      DLAM2 = diag(LAM2);
    [DLAM2,IS] = sort(DLAM2,'ComparisonMethod','real'); %sort eigenvalues by real part
    DLAM2 = DLAM2(end:-1:1);
    VS = V2(:,IS(end:-1:1));

    Lambda(:,IND) = DLAM2;
    Eigfunc{IND} = VS(:,1:10);


	fig4=figure(4)  %plot spectrum
	plot(real(DLAM2),imag(DLAM2),'.')
	fig4.Position = [500 500 1000 1000];
	title(['c = ', num2str(c),' w = ',num2str(w),' rho = ',num2str(RHODAT(IND))])

	figure(5) %plot a few eigenfunctions (in the exponentially weighted space)
	subplot(2,1,1)
	plot(x,real(VS(1:end/2,1)),x,real(VS(end/2+1:end,1)))
	subplot(2,1,2)
	plot(x,real(VS(1:end/2,2)),x,real(VS(end/2+1:end,2)))

	figure(6) %plot speed-wavenumber curve, with dot labeling the current speed and wavenumber values
	plot(C,K)
	hold on
	plot(c,k,'.','MarkerSize',10)
	hold off
	xlabel('c')
	ylabel('\omega')




end

SIZ = size(Lambda);


%plot real parts of eigenvalues
figure(2)
plot(CDAT(1:SIZ(2)),[real(Lambda(1:10,1:end))],'.-')
xlabel('c')
ylabel('Re \lambda')
drawnow

%%plot imaginary part of eigennvalues
%figure(3)
%plot(CDAT(1:ind-dsh+1),[imag(Lambda(1:10,1:end))],'.-')
%xlabel('c')
%label('Im \lambda')
%drawnow


figure(3)
plot(CDAT(1:SIZ(2)),[imag(Lambda(1:10,1:end))],'.-')
xlabel('c')
ylabel('Im \lambda')
drawnow

%plot trigger front location agains the real part of eigenvalues
figure(44)
plot(abs(RHODAT(1:SIZ(2))),[real(Lambda(1:10,1:end))],'.-')
xlabel('\xi_{tf}')
ylabel('Re \lambda')
drawnow

% figure(4)
% plot(CDAT,[real(LambdaF(1:10,:))],'.-')
% xlabel('c')
% ylabel('Im \lambda')
% drawnow

%%If one wants to see spec for a specific speed
PI = 130;
figure(4)
plot(real(Lambda(:,PI)),imag(Lambda(:,PI)),'x')
%%Save spectral data for use in plot_spec.m
save('spec_dat4.mat', 'Lambda','CDAT','WDAT','RHODAT','a','g','b','x','n','LL')
