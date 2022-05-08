function  ys = sensor(xs)
    range = sqrt(xs(1,:).^2+xs(2,:).^2);    
    angle = atan2(xs(2,:),xs(1,:));
    ys = [range;angle];
end