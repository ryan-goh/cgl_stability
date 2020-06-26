%%%%%Plots solution profiles from AUTO

ind = 300;
a = 0.3;
g = -0.2;
b =  0.4;
dsh = 17;
sk = 2;

SOL = cell(ind,1);

for number=1:ind
    file_name = sprintf('AUTO_Existence/psol_%d.dat',number);
    SOL{number} = load(file_name);
end




BIF = load('AUTO_Existence/pusheddata.dat');
C = BIF(:,6);
W = BIF(:,5);
CDAT = C(1:2:end);
%CDAT = [CDAT; C(164:5:end)];
WDAT = W(1:2:end);
%WDAT = [WDAT; W(164:5:end)];
RHO = BIF(:,1)



I = 50; %select index in bifurcation data to plot various things
figure(4)
plot3(SOL{I}(:,1),SOL{I}(:,2),SOL{I}(:,3))



M = length(S(:,1));



I = [1,50,100,130];
figure(1)
%Plot R vrs \xi
subplot(4,1,1)
plot(SOL{I(1)}(:,1),SOL{I(1)}(:,3),SOL{I(2)}(:,1),SOL{I(2)}(:,3),SOL{I(3)}(:,1),SOL{I(3)}(:,3),'-.',SOL{I(4)}(:,1),SOL{I(4)}(:,3),'--','LineWidth',1.005)
xx = xlabel('\xi');
yy= ylabel('q');
ll=legend(['c = ', num2str(CDAT(I(1)))], [' c=', num2str(CDAT(I(2)))],[' c =',num2str(CDAT(I(3)))] ,[' c = ', num2str(CDAT(I(4)))])
subplot(4,1,2)
plot(SOL{I(1)}(:,1),SOL{I(1)}(:,2),SOL{I(2)}(:,1),SOL{I(2)}(:,2),SOL{I(3)}(:,1),SOL{I(3)}(:,2),'-.',SOL{I(4)}(:,1),SOL{I(4)}(:,2),'--','LineWidth',1.005)
xx = xlabel('\xi');
yy= ylabel('\kappa');
ll=legend(['c = ', num2str(CDAT(I(1)))], [' c=', num2str(CDAT(I(2)))],[' c =',num2str(CDAT(I(3)))] ,[' c = ', num2str(CDAT(I(4)))])
subplot(4,1,3)
plot(SOL{I(1)}(:,1),SOL{I(1)}(:,4),SOL{I(2)}(:,1),SOL{I(2)}(:,4),SOL{I(3)}(:,1),SOL{I(3)}(:,4),'-.',SOL{I(4)}(:,1),SOL{I(4)}(:,4),'--','LineWidth',1.005)
xx = xlabel('\xi');
yy= ylabel('R');
ll=legend(['c = ', num2str(CDAT(I(1)))], [' c=', num2str(CDAT(I(2)))],[' c =',num2str(CDAT(I(3)))] ,[' c = ', num2str(CDAT(I(4)))])




% Here gca means get current figure.
%Plot z_i vrs \xi
set([ xx,yy,ll], ...
    'FontName'  , 'Times'            , ...
    'FontSize'  , 16          );

subplot(4,1,4)
plot(SOL{I(1)}(:,2),SOL{I(1)}(:,3),SOL{I(2)}(:,2),SOL{I(2)}(:,3),SOL{I(3)}(:,2),SOL{I(3)}(:,3),'-.',SOL{I(4)}(:,2),SOL{I(4)}(:,3),'--','LineWidth',1.005)
xx = xlabel('Re\{z\}')
yy= ylabel('Im\{z\}')
ll=legend(['c = ', num2str(CDAT(I(1)))], [' c=', num2str(CDAT(I(2)))],[' c =',num2str(CDAT(I(3)))] ,[' c = ', num2str(CDAT(I(4)))])

figure(3)
plot3(SOL{I(4)}(:,2),SOL{I(4)}(:,3),SOL{I(4)}(:,4))

figure(21)
plot(SOL{14}(:,1),SOL{14}(:,2).*sqrt(SOL{14}(:,4)), SOL{14}(:,1),-(SOL{14}(:,2).^2+CDAT(14)*SOL{14}(:,2)/(1+a^2)+1/(1+a^2)) + (1+a*g)*SOL{14}(:,4)/(1+a^2) + SOL{14}(:,3).^2 + (-CDAT(14)*a*SOL{14}(:,3) + a*WDAT(14))/(1+a^2) - SOL{14}(:,2).*sqrt(SOL{14}(:,4)) )


figure(20)
%Plot R vrs \xi
subplot(4,1,1)
plot(SOL{I(4)}(:,1),SOL{I(4)}(:,3),'--','LineWidth',1.005)
xx = xlabel('\xi');
yy= ylabel('q');
ll=legend([' c = ', num2str(CDAT(I(4)))])
subplot(4,1,2)
plot(SOL{I(4)}(:,1),SOL{I(4)}(:,2),'--','LineWidth',1.005)
xx = xlabel('\xi');
yy= ylabel('\kappa');
%ll=legend(['c = ', num2str(CDAT(I(1)))], [' c=', num2str(CDAT(I(2)))],[' c =',num2str(CDAT(I(3)))] ,[' c = ', num2str(CDAT(I(4)))])
subplot(4,1,3)
plot(SOL{I(4)}(:,1),SOL{I(4)}(:,4),'--','LineWidth',1.005)
xx = xlabel('\xi');
yy= ylabel('R');
%ll=legend(['c = ', num2str(CDAT(I(1)))], [' c=', num2str(CDAT(I(2)))],[' c =',num2str(CDAT(I(3)))] ,[' c = ', num2str(CDAT(I(4)))])
subplot(4,1,4)
plot(SOL{I(4)}(:,2),SOL{I(4)}(:,3),'--','LineWidth',1.005)
xx = xlabel('Re\{z\}');
yy= ylabel('Im\{z\}');
%ll=legend(['c = ', num2str(CDAT(I(1)))], [' c=', num2str(CDAT(I(2)))],[' c =',num2str(CDAT(I(3)))] ,[' c = ', num2str(CDAT(I(4)))])




% Here gca means get current figure.
set([ xx,yy,ll], ...
    'FontName'  , 'Times'            , ...
    'FontSize'  , 16          );


for ii = 1:ind
figure(10)
plot3(SOL{ii}(:,2),SOL{ii}(:,3),SOL{ii}(:,4))
%ll=legend(['c = ', num2str(CDAT(I(1)))], [' c=', num2str(CDAT(I(2)))],[' c =',num2str(CDAT(I(3)))] ,[' c = ', num2str(CDAT(I(4)))])
hold on

end
xx = xlabel('q')
yy= ylabel('kap')
zz = zlabel('r')
hold off
drawnow

for ii = 1:ind
  figure(11)
plot(SOL{ii}(:,1),SOL{ii}(:,4))
hold on

end
xlabel('\xi')
 ylabel('R')
hold off
drawnow
