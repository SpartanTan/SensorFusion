%This file draws a map, and allows us to manually
%draw the trajectory of the vehicle.

figure(1)
clf
hold on
plot([1+1i 1+9*1i 5+9*1i])
plot([7+9*1i 11+9*1i 11+1i 7+1i]);plot([5+1i 1+1i])
plot([2+5.2*1i 2+8.3*1i 4+8.3*1i 4+5.2*1i 2+5.2*1i])%House 1
plot([2+3.7*1i 2+4.4*1i 4+4.4*1i 4+3.7*1i 2+3.7*1i])%House 2
plot([2+2*1i 2+3.2*1i 4+3.2*1i 4+2*1i 2+2*1i])%House 3
plot([5+1i 5+2.2*1i 7+2.2*1i 7+1i])%House 4
plot([5+2.8*1i 5+5.5*1i 7+5.5*1i 7+2.8*1i 5+2.8*1i])%House 5
plot([5+6.2*1i 5+9*1i]);plot([7+9*1i 7+6.2*1i 5+6.2*1i])%House 6
plot([8+4.6*1i 8+8.4*1i 10+8.4*1i 10+4.6*1i 8+4.6*1i])%House 7
plot([8+2.4*1i 8+4*1i 10+4*1i 10+2.4*1i 8+2.4*1i])%House 8
plot([8+1.7*1i 8+1.8*1i 10+1.8*1i 10+1.7*1i 8+1.7*1i])%House 9

axis([0.8 11.2 0.8 9.2])
title('Map','FontSize',20)

disp('Start clicking in the graph to create a trajectory!')
disp('Press "Return" to finish.')

key = '1';
Xk = [];
clear h
while true
    [X,Y,key]=ginput(1);
    if isempty(key) || key ~= 1
        break
    end
    Xk = [Xk, [X; Y]];
    if exist('h', 'var')
        h.XData = Xk(1,:);
        h.YData = Xk(2,:);
    else
        h = plot(Xk(1,:),Xk(2,:),'k-*');
    end
    drawnow;
end


save Xk Xk
