function  ys = sensor(xs)
%Input
%   xs  Object position [n*1]

    range = sqrt(xs(1,:).^2+xs(2,:).^2);    
    angle = atan2(xs(2,:),xs(1,:));
    ys = [range;angle];
end