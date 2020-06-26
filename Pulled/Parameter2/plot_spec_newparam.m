ind = 280;
makespecvid = 0;
jf = 1;
a = -0.1;
g = -0.9; %%%read from auto code, 
LL = 500;  %interpolation width,   note, if not the same as AUTO domain size, may lose the zero eigenvalue from gauge symmetry
n = 20000;   %number of grid points
x = linspace(-LL,LL,n);
dx = x(2) - x(1);
clin = 2*sqrt(1+a^2);

SOL = cell(ind,1);
load('c-Lam.mat')
for number=1:ind
    file_name = sprintf('AUTO/parallelsol_%d.dat',number);
    SOL{number} = load(file_name);
end

BIF = load('AUTO/parallelbif.dat');
C = BIF(:,1);
W = BIF(:,end);
CDAT = C(1:5:ind*5);
WDAT = W(5:5:ind*5);
KDAT = (-CDAT+sqrt(CDAT.^2+4*(g+WDAT).*(g-a)))./(2*(g-a));

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




% %%%Plots of (c,k) curve and it's specific location
% figure(21)
% plot(clin-CDAT,abs(KDAT),'LineWidth',3)
% %     plot(clin-c,abs(k),'.','MarkerSize',20)
% %     hold off
%     xlim([0,clin-CDAT(90)])
%     xlabel('$\Delta c = c_{lin} - c$','Interpreter','latex')
%     ylabel('$|k_{tf}|$','Interpreter','latex')
%     f4 = gca;
%     f4.FontSize = 16;
%     
% %%%Plot trigger location versus speed
% load('trigger_loc.mat','ATRIG','TRIG','CH')
% 
% ATRIG = abs(TRIG);
% 
% % figure(21)
% % plot(CDAT,TRIG)
% 
% CH = clin-CDAT;
% figure(22)
% subplot(2,1,1)
% plot(CH,ATRIG,'LineWidth',3)
% xlabel('$\Delta c = c_{lin} - c$','Interpreter','latex')
% ylabel('$|\xi_{tf}|$','Interpreter','latex')
%  xlim([0,clin])
% ff = gca;
% ff.FontSize = 18
% subplot(2,1,2)
%     plot(clin-CDAT,abs(KDAT),'LineWidth',3)
%     hold on
% %     plot(clin-c,abs(k),'.','MarkerSize',20)
% %     hold off
%     xlim([0,clin])
%     xlabel('$\Delta c = c_{lin} - c$','Interpreter','latex')
%     ylabel('$|k_{tf}|$','Interpreter','latex')
%     f4 = gca;
%     f4.FontSize = 16;
%     drawnow
% 
% 
% 
% CS = CH.^(-1/2);
% figure(23)
% plot(CS(40:end), ATRIG(40:end),'LineWidth',3)
% xlabel('$\Delta c^{-1/2} = (c_{lin} - c)^{-1/2}$','Interpreter','latex')
% ylabel('$|\xi_{tf}|$','Interpreter','latex')
% xlim([min(CS), max(CS)])
% ff = gca;
% ff.FontSize = 18
%     




%%%Plot solution and spectrum from two specific speeds

%%%Plot solution and spectrum from one specific speed
IND = 200;
S = SOLC{IND};
c = CDAT(IND)
w = WDAT(IND)
k = KDAT(IND)


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

%%plots of solution components (z,R)
figure(200)
f = gcf
f.Position = [1309 765 560 420/2]
%ff.Position = [100 100 800 300]
idel1 = find(r<0.05);
idel = idel1(1);
z1 = (-c+sqrt(c^2-4*(1-1i*w)*(1+1i*a)))/(2+2i*a);
z2 = (-c-sqrt(c^2-4*(1-1i*w)*(1+1i*a)))/(2+2i*a);
plot(kap,q,'LineWidth',3)
hold on
plot(kap(idel),q(idel),'xk','MarkerSize',10,'LineWidth',2)
text(kap(idel)+0.05,q(idel)+0.003,'$ z(\xi_\delta)$','Interpreter','latex','FontSize',16)
plot(kap(1),q(1),'xb','MarkerSize',10,'LineWidth',2)
text(kap(1)-0.05,q(1)-0.005,'$-i k_\mathrm{tf}$','Interpreter','latex','FontSize',16)
plot(kap(end),q(end),'xr','MarkerSize',10,'LineWidth',2)
text(kap(end)+0.05,q(end),'$ z_+$','Interpreter','latex','FontSize',16)
plot(real(z1),imag(z1),'.k','MarkerSize',20)
text(real(z1)+0.03,imag(z1),'$ z_{0,+}$','Interpreter','latex','FontSize',16,'LineWidth',2)
plot(real(z2),imag(z2),'.k','MarkerSize',20)
text(real(z2)+0.03,imag(z2),'$ z_{0,-}$','Interpreter','latex','FontSize',16,'LineWidth',2)
%annotation('arrow', [kap(end/2),kap(end/2)+1],[q(end/2),q(end/2+1)]);
hold off
xlabel(['$\mathrm{Re} \,z_\mathrm{tf}$'],'Interpreter','latex')
ylabel(['$\mathrm{Im} \,z_\mathrm{tf}$'],'Interpreter','latex')
%title(['$z_{tf}\in \mathbb{C}$'])
f2 = gca
f2.FontSize = 16;
xlim([min(kap) - 0.3, max(kap+0.2)])
ylim([min(q) - 0.01,max(q) + 0.03])
title([' \,\, $\Delta c = c_{lin} - c = $', num2str(clin - c)],'Interpreter','latex')


figure(201)
f = gcf
f.Position = [1309 765 560 420/2]
plot(x,r.^2,'LineWidth',3)
ylim([-0.05,1])
ylabel(['$R_\mathrm{tf}$'],'Interpreter','latex')
xlabel(['$\xi$'],'Interpreter','latex')
%title(['$r_{tf}$'])
xlim([-200 20])
fx = gca
fx.FontSize = 16;
fx.YTick = [0 r(idel) 0.5 r(1)];
fx.YTickLabel = {'0','$\delta$','1/2','$ 1-k_\mathrm{tf}^2$'};
fx.XTick = [x(idel) 0 ];
fx.XTickLabel = {'$\xi_\mathrm{tf}$','0'};
fx.TickLabelInterpreter = 'latex';


%%%Plot of spectrum at this point
%%%Plot of spectrum at this point
 jp = IND-0;
        t = [0:0.01:10];
    LABS = 1 - (t + CDAT(jp)^2/4)/(1+1i*a) - 1i*WDAT(jp);
  figure(24)
  f4 = gcf;
   % f4.Position = [300 300 600 600];
    plot(real(Lambda(3:end,jp)),imag(Lambda(3:end,jp)),'.k','MarkerSize',14)
    hold on
    plot(real(Lambda(1,jp)),imag(Lambda(1,jp)),'.','MarkerSize',20,'Color',[0 0.4470 0.7410])
     plot(real(Lambda(2,jp)),imag(Lambda(2,jp)),'.','MarkerSize',20,'Color',[0.8500 0.3250 0.0980])
     plot(real(LABS),[imag(LABS);-imag(LABS)],'Color','Blue','LineWidth',3)
     hold off
     title(['c = ',num2str(CDAT(jp)), '  $ \Delta c = c_{lin} - c = $ ', num2str(clin-CDAT(jp))],'Interpreter','Latex')
     xlim([-0.15 real(LABS(1))+0.01])
     ylim([-0.02 0.02])
     xlabel(['Re \lambda'])
    ylabel(['Im \lambda'])
    f2 = gca;
    f2.FontSize = 16;
    axes('Position',[.64 .18 .22 .22])
    box on
    plot(real(Lambda(3:end,jp)),imag(Lambda(3:end,jp)),'.k','MarkerSize',14)
    hold on
    plot(real(Lambda(2,jp)),imag(Lambda(2,jp)),'x','MarkerSize',20,'Color',[0.8500 0.3250 0.0980],'LineWidth',3)
    plot(real(Lambda(1,jp)),imag(Lambda(1,jp)),'.','MarkerSize',20,'Color',[0 0.4470 0.7410])
     plot(real(LABS),[imag(LABS);-imag(LABS)],'Color','Blue','LineWidth',3)
     plot(real(LABS(1)),[imag(LABS(1));-imag(LABS(1))],'x','Color','Blue','MarkerSize',20,'LineWidth',3)
     hold off
     xlim([-0.0001 real(LABS(1))+0.0001])
     ylim([-0.0002 0.0002])
    drawnow

%%%%%

IND = 200;
S = SOLC{IND};
c = CDAT(IND)
w = WDAT(IND)
k = KDAT(IND)


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

%%plots of solution components (z,R)
figure(30)
subplot(2,1,1)
plot(kap,q,'LineWidth',3)
xlabel(['$\mathrm{Re} \,z$'],'Interpreter','latex')
ylabel(['$\mathrm{Im} \,z$'],'Interpreter','latex')
%title(['$z_{tf}\in \mathbb{C}$'])
f2 = gca
f2.FontSize = 16;
title(['$c =$ ', num2str(c), ', \,\, $\Delta c = c_{lin} - c = $', num2str(clin - c)],'Interpreter','latex')
subplot(2,1,2)
plot(x,r.^2,'LineWidth',3)
ylim([-0.05,1])
ylabel(['$R_{tf}$'],'Interpreter','latex')
xlabel(['$\xi$'],'Interpreter','latex')
%title(['$r_{tf}$'])
f2 = gca
f2.FontSize = 16;

%%Plot solution profile $r(\xi)exp(i \int_\xi0^xi q ds)
%%define phase by integrating q
phi = cumtrapz(x,q);
ATF = r.*exp(1i*phi);
figure(301)
plot(x,real(ATF),'LineWidth',2)
hold on
plot(x,imag(ATF),'LineWidth',2,'Color',[0.6350 0.0780 0.1840])
plot(x',(heaviside(-x)-0.5)*2,'-.','LineWidth',2)
hold off
xlabel('$\xi$','Interpreter','latex');
%ylabel('Re $A(\xi)$', 'Interpreter','latex');
legend(['Re $A(\xi)$'],['Im $A(\xi)$'],'Interpreter','latex','Location','best');
 f2 = gca;
    f2.FontSize = 16;
xticks([-120 0])
xticklabels({'\xi_{tf}','0'})
text(-400,-1.1,'diffusively stable','FontSize',14,'Interpreter','latex')
text(-120,-0.2,'abs. unstable','FontSize',14,'Interpreter','latex')
text(40,-0.2,'stable','FontSize',14,'Interpreter','latex')
text(20,1.0,'$-sign(x)$','FontSize',14,'Interpreter','latex')
%text
ylim([-1.2,1.2])
xlim([-500,140])
ax = gca;
ax.FontSize = 16;
% %subplot(2,1,2)
% %plot(x,imag(ATF),'LineWidth',3)
% %hold on
% %plot(x',(heaviside(-x)-0.5)*2,'.-','LineWidth',3)
% %hold off
% xlabel('$\xi$','Interpreter','latex');
% ylabel('Im $A(\xi)$', 'Interpreter','latex');
% % f2 = gca;
% %    f2.FontSize = 16;
% xticks([-120 0])
% xticklabels({'\xi_{tf}','0'})
% %text(-400,-1.1,'diffusively stable','FontSize',14,'Interpreter','latex')
% %text(-120,-0.2,'abs. unstable','FontSize',14,'Interpreter','latex')
% %text(40,-0.1,'stable','FontSize',14,'Interpreter','latex')
% %text(20,1.0,'$-sign(x)$','FontSize',14,'Interpreter','latex')
% %text
% ylim([-1.2,1.2])
% xlim([-500,100])
% ax = gca;
% ax.FontSize = 16;





%%%Plot of spectrum at this point
 jp = IND; 
        t = [0:0.01:10];
    LABS = 1 - (t + CDAT(jp)^2/4)/(1+1i*a) - 1i*WDAT(jp);
  f4 = figure(34)
    %f4.Position = [300 300 600 600]
    plot(real(Lambda(3:end,jp)),imag(Lambda(3:end,jp)),'.k','MarkerSize',14)
    hold on
    plot(real(Lambda(1,jp)),imag(Lambda(1,jp)),'x','MarkerSize',20,'linewidth',3)
     plot(real(Lambda(2,jp)),imag(Lambda(2,jp)),'.','MarkerSize',20)
     plot(real(LABS),[imag(LABS);-imag(LABS)],'Color','Blue','LineWidth',3)
     hold off
     title(['$c = $',num2str(CDAT(jp)), '  $ \Delta c = c_{lin} - c =  $', num2str(clin-CDAT(jp))],'Interpreter','latex')
     xlim([-0.15 real(LABS(1))+0.01])
     ylim([-0.02 0.02])
     xlabel(['$\mathrm{Re}\, \lambda$'],'Interpreter','latex')
    ylabel(['$\mathrm{Im}\, \lambda$'],'Interpreter','latex')
    f44 = gca;
    f44.FontSize = 16;
    
  
    
    
 %%%%%%%%%





% 
% figure(2)
% plot(CDAT,[real(Lambda(1:10,:))],'.-','LineWidth',3)
% xlabel('c')
% ylabel('Re \lambda')
% drawnow
% xlabel(['Re \lambda'])
% ylabel(['Im \lambda'])
% xlim([-1 0.1])
% f2 = gca
% f2.FontSize = 16;
% 
% figure(3)
% plot(clin - CDAT,[real(Lambda(1:10,:))],'.-','LineWidth',2)
% xlabel('\epsilon')
% ylabel('Re \lambda')
% f1 = gca
% f1.FontSize = 16;
% xlim([0,0.3])
% ylim([-1 0.1])

figure(50)
plot(CDAT,[real(Lambda(1:40,:))],'.-','LineWidth',2)
xlabel('c')
ylabel('Re \lambda')
f1 = gca
f1.FontSize = 16;
xlim([1.5 clin])
ylim([-0.2 0.01])

f51 = figure(51)
f51.Position = [300 300 1200 500]
subplot(2,4,1:3)
plot(CDAT,[real(Lambda(1:40,:))],'.-','LineWidth',2)
xlabel('$c$','Interpreter','Latex')
ylabel('$\mathrm{Re}\, \lambda$','Interpreter','Latex')
f1 = gca
f1.FontSize = 16;
xlim([0 clin])
ylim([-0.2 0.01])
subplot(2,4,4)
plot(CDAT,[real(Lambda(1:80,:))],'.-','LineWidth',3)
xlabel('$c$','Interpreter','Latex')
%ylabel('$\mathrm{Re}\, \lambda$','Interpreter','Latex')
xlim([clin-0.02 clin])
ylim([-0.1 0.01])
f1 = gca
f1.FontSize = 16

subplot(2,4,5)
I1 = 40;
 LABS = 1 - (t + CDAT(I1)^2/4)/(1+1i*a) - 1i*WDAT(I1);
plot(real(Lambda(:,I1)),imag(Lambda(:,I1)),'.k','MarkerSize',14)
hold on
plot(real(LABS),[imag(LABS);-imag(LABS)],'Color','Blue','LineWidth',2)
hold off
title(['$c = $', num2str(CDAT(I1))],'Interpreter','Latex')
ylabel('$\mathrm{Re}\, \lambda$','Interpreter','Latex')
xlabel('$\mathrm{Im}\, \lambda$','Interpreter','Latex')
f1 = gca;
f1.FontSize = 16;
xlim([-0.5 0.01])
ylim([-0.4 0.4])
subplot(2,4,6)
I2 = 90;
LABS = 1 - (t + CDAT(I2)^2/4)/(1+1i*a) - 1i*WDAT(I2);
plot(real(Lambda(:,I2)),imag(Lambda(:,I2)),'.k','MarkerSize',14)
hold on
plot(real(LABS),[imag(LABS);-imag(LABS)],'Color','Blue','LineWidth',2)
hold off
title(['$c = $', num2str(CDAT(I2))],'Interpreter','Latex')
f1 = gca;
f1.FontSize = 16;
xlim([-0.5 0.01])
ylim([-0.4 0.4])
subplot(2,4,7)
I3 = 130;
LABS = 1 - (t + CDAT(I3)^2/4)/(1+1i*a) - 1i*WDAT(I3);
plot(real(Lambda(:,I3)),imag(Lambda(:,I3)),'.k','MarkerSize',14)
hold on
plot(real(LABS),[imag(LABS);-imag(LABS)],'Color','Blue','LineWidth',2)
hold off
title(['$c = $', num2str(CDAT(I3))],'Interpreter','Latex')
f1 = gca;
f1.FontSize = 16;
xlim([-0.5 0.01])
ylim([-0.4 0.4])
subplot(2,4,8)
I4 = 160;
LABS = 1 - (t + CDAT(I4)^2/4)/(1+1i*a) - 1i*WDAT(I4);
plot(real(Lambda(:,I4)),imag(Lambda(:,I4)),'.k','MarkerSize',14)
hold on
plot(real(LABS),[imag(LABS);-imag(LABS)],'Color','Blue','LineWidth',2)
hold off
title(['$c = $', num2str(CDAT(I4))],'Interpreter','Latex')
f1 = gca;
f1.FontSize = 16;
xlim([-0.5 0.01])
ylim([-0.4 0.4])



figure(55)
plot(CDAT,[real(Lambda(1:80,:))],'.-','LineWidth',3)
xlabel('c')
ylabel('$\mathrm{Re}\, \lambda$','Interpreter','Latex')
f1 = gca
f1.FontSize = 16;
xlim([clin-0.1 clin])
ylim([-0.1 0.01])
axes('Position',[.4 .4 .35 .35])
box on
plot(CDAT,[real(Lambda(1:80,:))],'.-','LineWidth',3)
xlim([clin-0.01 clin])
ylim([-0.1 0.01])
drawnow


figure(555)
plot(CDAT,[real(Lambda(1:80,:))],'.-','LineWidth',3)
xlabel('c')
ylabel('$\mathrm{Re}\, \lambda$','Interpreter','Latex')
f1 = gca
f1.FontSize = 16;
xlim([clin-0.01 clin])
ylim([-0.1 0.01])

figure(56)
plot(CDAT,[real(Lambda(3:80,:))./(clin-CDAT')],'.-','LineWidth',3)
xlabel('c')
ylabel('$Re \lambda / (c_{lin}-c)$','Interpreter','Latex')
f1 = gca
f1.FontSize = 16;
ylim([-500 0])
xlim([clin-0.003,clin-0.0002])

figure(561)
plot(clin-CDAT,[real(Lambda(3:20,:))./(clin-CDAT')],'.-','LineWidth',3)
xlabel('$c_{lin}-c$','Interpreter','Latex')
ylabel('$Re \lambda / (c_{lin}-c)$','Interpreter','Latex')
f1 = gca
f1.FontSize = 16;
%ylim([-500 0])
xlim([0.0002,0.01])

figure(57)
plot((clin-CDAT).^(3/2),real(Lambda(2,:)),'.-','LineWidth',3,'Color',[0.8500, 0.3250, 0.0980])
xlabel('$(c_{lin} - c)^{3/2}$','Interpreter','latex')
ylabel('$\mathrm{Re}\, \lambda$','Interpreter','Latex')
title('Second Eigenvalue in $D_1$','Interpreter','latex')
f1 = gca
f1.FontSize = 16;
xlim([0 0.2])
%ylim([-0.115 0.01])



figure(60) %plot essential spectrum
tt = [-10:0.01:10];
DM1 = 1 - a^2*tt.^4 + 2*a^2*tt.^2 + 4i*a^2*(tt.^3 + tt);
DM2 = -tt.^2+2i*tt -1;
LP = (DM2 + sqrt(DM1))/(1+a^2);
LM = (DM2 - sqrt(DM1))/(1+a^2);
plot(real(LP),imag(LP),'.-','LineWidth',3);%,'Color',[0.8500, 0.3250, 0.0980])
hold on
plot(real(LM),imag(LM),'.-','LineWidth',3);
plot(real(LABS),[imag(LABS);-imag(LABS)],'Color','Blue','LineWidth',3)
hold off
xlabel('Re $\lambda_\pm(l)$','Interpreter','latex')
ylabel('Im $\lambda_\pm(l)$','Interpreter','Latex')
title('Essential spectrum curves','Interpreter','latex')
f1 = gca
f1.FontSize = 16;
xlim([-10 0.5])
ylim([-10 10])










if makespecvid
    
        fname = sprintf('spec_mov_new');
    v = VideoWriter(fname);
    v.Quality = 100;
    open(v)
    
    for jj = 90:length(CDAT)
        
        t = [0:0.01:10];
    LABS = 1 - (t + CDAT(jj)^2/4)/(1+1i*a) - 1i*WDAT(jj);
    c = CDAT(jj);
    w = WDAT(jj);
    k = KDAT(jj);
  f4 = figure(4)
    f4.Position = [300 300 600 600]
    subplot(3,1,1:2)
    plot(real(Lambda(3:end,jj)),imag(Lambda(3:end,jj)),'.k','MarkerSize',14)
    hold on
    plot(real(Lambda(1,jj)),imag(Lambda(1,jj)),'.','MarkerSize',20)
     plot(real(Lambda(2,jj)),imag(Lambda(2,jj)),'.','MarkerSize',20)
     plot(real(LABS),[imag(LABS);-imag(LABS)],'Color','Blue','LineWidth',3)
     hold off
     title(['c = ',num2str(CDAT(jj)), '   \epsilon = ', num2str(clin-CDAT(jj))])
     xlim([-0.15 real(LABS(1))+0.01])
     ylim([-0.05 0.05])
     xlabel(['Re \lambda'])
    ylabel(['Im \lambda'])
    f2 = gca;
    f2.FontSize = 16;
    
    subplot(3,1,3)
    plot(clin-CDAT,abs(KDAT),'LineWidth',3)
    hold on
    plot(clin-c,abs(k),'.','MarkerSize',20)
    hold off
    xlim([0,clin-CDAT(90)])
    xlabel('\epsilon = c_{lin} - c')
    ylabel('|k|')
    f4 = gca;
    f4.FontSize = 16;
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
     
    end
close(v)
end

