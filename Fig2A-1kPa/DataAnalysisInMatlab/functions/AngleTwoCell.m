function [angle1,angle2,angle3] = AngleTwoCell(thetas);
%angle in degrees

%THIS CODE ASSUMES CELLS ARE PLACED NEXT TO EACHOTHER
theta1=thetas(1);
theta2=thetas(2);

if (theta1<0 && theta2 < 0)
    angle1 = 180 + min(thetas);
    angle2 = -max(thetas);
   angle3 =-min(thetas)+max(thetas);
end
if ( (theta1<=0 && theta2>=0) | (theta2<=0 && theta1>=0))
    angle1=abs(theta1);
    angle2=abs(theta2);
    angle3 = 180-max(thetas)+min(thetas); 
end
if (theta1>0 && theta2 > 0)
   angle1 = 180 - max(thetas);
   angle2 = min(thetas);
   angle3 = -min(thetas)+max(thetas);
end

end
