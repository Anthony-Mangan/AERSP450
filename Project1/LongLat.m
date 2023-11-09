function [long,lat] = LongLat(x,y,z)
    long = zeros(size(x));
    lat = long;
    for i = 1:length(x)
        long(i) = atan2(y(i),x(i))*180/pi;
        lat(i) = atan2(z(i),sqrt(x(i)^2 + y(i)^2))*180/pi;
    end
end
