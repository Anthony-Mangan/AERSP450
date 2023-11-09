function [E1] = NR(e,M1)
error = 1;
absTol = 1E-8;
Eold = M1;
    while error > absTol
        y = M1 + e*sin(Eold) - Eold;
        z = e*cos(Eold) - 1;
        Enew = Eold - y/z;
        error = abs(Enew-Eold);
        Eold = Enew;
    end
E1 = Enew;
end
