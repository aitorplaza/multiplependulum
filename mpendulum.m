clear all
clc

%% Define constants

nPend=10; % Number of pendulums
nOscillations = 40% Number of oscillations to synchronize again

lengths = ComputeLengths(nPend,nOscillations)

theta_0 = 0.3;
dtheta_0 = 0;

t_final=nOscillations+0.24; % time adjusted due to non-linearities
delta_t=0.04;

%% Integrate

data.l = lengths';
data.g = 9.8;

state_init = [theta_0*ones(1,nPend),dtheta_0*ones(1,nPend)]';
tspan = [0:delta_t:t_final];

options = odeset('RelTol',1e-7);
[tSeries, stateSeries] = ode45 (@deriv, tspan, state_init, options, data);

%plot(stateSeries(:,1:nPend))

%% Plot
figure(1)
plot(stateSeries(:,1:nPend));

xSeries =  sin(stateSeries(:,1:nPend))*diag(lengths);
ySeries = -cos(stateSeries(:,1:nPend))*diag(lengths);

h=figure(2)
hold on;axis equal;grid on; 
title('Pendulum Waves','interpreter','latex','Fontsize', 18) 
yMin=min(min(ySeries));
xMax=max(max(xSeries));
axis([-1.2*xMax 1.2*xMax yMin*1.2 0.02])
set(gca,'xtick',[])
set(gca,'ytick',[])
hAxes = gca;     %Axis handle
hAxes.XRuler.Axle.LineStyle = 'none';  
hAxes.YRuler.Axle.LineStyle = 'none';

for i= 1:1:nPend  
    bar(i) = plot([0,xSeries(1,i)],[0,ySeries(1,i)],'k');
    ball(i)= scatter(xSeries(1,i),ySeries(1,i),100,'filled' );
end
plot(0,0,'.k');
hold off

% animate plot:
disp('Press a key to continue!')  
pause;
filename = 'PendulumWaves.gif';

dt=2*delta_t;
for i=1:length(tSeries)
    for j= 1:1:nPend
        x = xSeries(i,j);
        y = ySeries(i,j);
        set(bar(j),'Xdata',[0, x], 'YData',[0,y])
        set(ball(j),'Xdata', x, 'YData',y)
    end
    tSeries(i)
    
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    if tSeries(i) == 0;
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',dt);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',dt);
    end 
    
end

