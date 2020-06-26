%plot_spec.m


load('spec_dat4.mat');%, 'Lambda','CDAT','WDAT','RHODAT','a','g','b','x','n','LL')


SIZ = size(Lambda);

figure(1);
subplot(2,1,1)
plot(CDAT,WDAT)
xlabel('c')
ylabel('\omega')
subplot(2,1,2)
plot(CDAT,RHODAT)
xlabel('c')
ylabel('$\xi_\mathrm{tf}$','Interpreter','Latex')
ax = gca;
ax.FontSize = 16;


figure(2)
plot(CDAT(1:SIZ(2)),[real(Lambda(1:10,1:end))],'.-','Color',[0 0.4470 0.7410])
xlabel('c')
ylabel('Re \lambda')
ax = gca;
ax.FontSize = 16;
drawnow

%figure(3)
%plot(CDAT(1:ind-dsh+1),[imag(Lambda(1:10,1:end))],'.-')
%xlabel('c')
%label('Im \lambda')
%drawnow


figure(3)
plot(CDAT(1:SIZ(2)),[imag(Lambda(1:10,1:end))],'.-')
xlabel('c')
ylabel('Im \lambda')
ax = gca;
ax.FontSize = 16;
drawnow


figure(44)
plot(abs(RHODAT(1:SIZ(2))),[real(Lambda(1:10,1:end))],'.-','Color',[0 0.4470 0.7410])
xlabel('\xi_{tf}')
ylabel('Re \lambda')
ax = gca;
ax.FontSize = 16;
drawnow

figure(45)
plot(abs(RHODAT(1:SIZ(2))),[real(Lambda(1:2,1:end))],'.-','Color',[0 0.4470 0.7410])
xlabel('\xi_{tf}')
ylabel('Re \lambda')
xlim([0,7])
ylim([-0.1,0.1])
ax = gca;
ax.FontSize = 16;
drawnow



PI = 130;
figure(4)
plot(real(Lambda(:,PI)),imag(Lambda(:,PI)),'x')
