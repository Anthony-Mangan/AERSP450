function[output] = DCM(angle,axis)
    switch axis
        case 1
            output = [1,         0,             0;
                      0,    cos(angle),    sin(angle);
                      0,    -sin(angle),   cos(angle)];
        case 2
            output =  [cos(angle),    0,       -sin(angle);
                           0,         1,           0;
                       sin(angle),    0,        cos(angle)];
         case 3
                   
            output =  [cos(angle),     sin(angle),       0;
                       -sin(angle),    cos(angle),       0;
                            0,              0,           1];
    end
     
end