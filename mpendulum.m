clear all
clc

%% Define constants

N=9; 
% dl=0.01;
% l=0.2;
% lengths = l:dl:l+dl*(N-1)
lengths = [0.205, 0.218, 0.232, 0.248, 0.266, 0.285 ,0.306, 0.330, 0.357];
theta_0 = 0.3;
dtheta_0 = 0;

t_final=30.0;
delta_t=0.01;

%% Integrate

data.l = lengths';
data.g = 9.8;

state_init = [theta_0*ones(1,N),dtheta_0*ones(1,N)]';
tspan = [0:delta_t:t_final];

options = odeset('RelTol',1e-7);
[tSeries, stateSeries] = ode45 (@deriv, tspan, state_init, options, data);

%% Plot
clrs=['r','g','b','c','m','y','k'];
xSeries =  sin(stateSeries(:,1:N))*diag(lengths);
ySeries = -cos(stateSeries(:,1:N))*diag(lengths);


h=figure(1)
hold on;axis equal;grid on; 
title('Multiple Pendulum','interpreter','latex','Fontsize', 12) 
yMin=min(min(ySeries));
xMax=max(max(xSeries));
axis([-1.1*xMax 1.1*xMax yMin*1.1 0.1])
set(gca,'xtick',[])
set(gca,'ytick',[])

for i= N:-1:1
    bar(i) = plot([0,xSeries(1,i)],[0,ySeries(1,i)],'k');
    ball(i)= scatter(xSeries(1,i),ySeries(1,i),100,'filled' );
end
hold off

% animate plot:
disp('Press a key to continue!')  
pause;
filename = 'multiPendulum.gif';

dt=delta_t*0.5;
for i=1:length(tSeries)
    for j= N:-1:1
        x = xSeries(i,j);
        y = ySeries(i,j);
        set(bar(j),'Xdata',[0, x], 'YData',[0,y])
        set(ball(j),'Xdata', x, 'YData',y)
    end
    
    
    pause(dt)
    
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    if tSeries(i) == 0;
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',dt);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',dt);
    end 
    
end

