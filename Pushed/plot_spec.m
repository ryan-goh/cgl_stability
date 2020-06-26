%plot_spec.m
%%loads data from spectral solver and the AUTO_Existence continuation to plot various figures

load('spec_dat_save.mat');%, 'Lambda','CDAT','WDAT','RHODAT','a','g','b','x','n','LL')


SIZ = size(Lambda);

%locate the first two bifurcation points
IND1 = find(Lambda(1,:)>0.01);
bif1 = IND1(1)-1;
IND2 = find(Lambda(2,:)< -0.009);
iii = find(IND2(2:end) - IND2(1:end-1) >1);
bif2 = IND2(iii+1);

%Existence data, plot c vrs. frequency curve, c vrs. trigger location curve
figure(1);
subplot(2,1,1)
plot(CDAT,WDAT,'LineWidth',3)
hold on
plot(CDAT(bif1),WDAT(bif1),'o','MarkerSize',10,'LineWidth',1);
plot(CDAT(bif2),WDAT(bif2),'o','MarkerSize',10,'LineWidth',1);
hold off
xlabel('$c$','Interpreter','Latex')
ylabel('$\omega$','Interpreter','Latex')
ax = gca;
ax.FontSize = 16;
axes('Position',[.4 .8 .2 .1])
box on
plot(CDAT,WDAT,'LineWidth',3)
xlim([2.16,2.19])
ylim([-2.5,-2.44])

subplot(2,1,2)
plot(CDAT,RHODAT,'LineWidth',3)
hold on
plot(CDAT(bif1),RHODAT(bif1),'o','MarkerSize',10,'LineWidth',1);
plot(CDAT(bif2),RHODAT(bif2),'o','MarkerSize',10,'LineWidth',1);
hold off
xlabel('$c$','Interpreter','Latex')
ylabel('$\xi_\mathrm{tf}$','Interpreter','Latex')
ax = gca;
ax.FontSize = 16;
axes('Position',[.4 .2 .2 .2])
box on
plot(CDAT,RHODAT,'LineWidth',3)
xlim([2.18,2.185])
ylim([-7 0])



%%Plot trigger location vrs. real part of most unstable eigenvalues
figure(44)
plot(RHODAT(1:SIZ(2)),[real(Lambda(1:10,1:end))],'.-','Color',[0 0.4470 0.7410],'LineWidth',1)
hold on
plot(RHODAT(bif1),real(Lambda(1,bif1)),'o','MarkerSize',10,'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
plot(RHODAT(bif2),real(Lambda(1,bif2)),'o','MarkerSize',10,'LineWidth',1,'Color',[0.9290 0.6940 0.1250]);
hold off
xlabel('$\xi_{tf}$','Interpreter','Latex')
ylabel('$\mathrm{Re }\lambda$','Interpreter','Latex')
ax = gca;
ax.FontSize = 16;
ylim([-1,0.4])
xlim([-7,0])
axes('Position',[.2 .8 .4 .1])
box on
plot(RHODAT(1:SIZ(2)),[real(Lambda(1:10,1:end))],'-','Color',[0 0.4470 0.7410],'LineWidth',1)
xlim([-7,0])
ylim([-0.01,0.01])
drawnow

%Plot of c vrs. real part of most unstable eigenvalues
figure(2)
plot(CDAT(1:SIZ(2)),[real(Lambda(1:10,1:end))],'.-','Color',[0 0.4470 0.7410],'LineWidth',3)
hold on
plot(CDAT(bif1),real(Lambda(1,bif1)),'o','MarkerSize',10,'LineWidth',1);
plot(CDAT(bif2),real(Lambda(1,bif2)),'o','MarkerSize',10,'LineWidth',1);
hold off
xlabel('c','Interpreter','Latex')
ylabel('$\mathrn{Re} \lambda$','Interpreter','Latex')
ax = gca;
ax.FontSize = 16;
drawnow

%figure(3)
%plot(CDAT(1:ind-dsh+1),[imag(Lambda(1:10,1:end))],'.-')
%xlabel('c')
%label('Im \lambda')
%drawnow

%Plot of c vrs. real part of most unstable eigenvalues
figure(3)
plot(CDAT(1:SIZ(2)),[imag(Lambda(1:10,1:end))],'.-')
xlabel('c')
ylabel('$\mathrm{Im }\lambda$','Interpreter','Latex')
ax = gca;
ax.FontSize = 16;
drawnow

%%Plot of trigger location vrs. real part of most unstable eigenvalues
figure(44)
plot(abs(RHODAT(1:SIZ(2))),[real(Lambda(1:10,1:end))],'.-','Color',[0 0.4470 0.7410])
xlabel('$\xi_{tf}$','Interpreter','Latex')
ylabel('$\mathrm{Re }\lambda$','Interpreter','Latex')
ax = gca;
ax.FontSize = 16;
ylim([]
drawnow

%Pick a data point somewhere along the bifurcation curve and plot the associated spectrum.
PI = 100;
figure(4)
plot(real(Lambda(:,PI)),imag(Lambda(:,PI)),'x')
