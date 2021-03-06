function I=getIntensity_QWLSI(xs,ys,delta,f,lambda,incidentAngle)
PI=3.14159265;
offset=f*tan(incidentAngle);
I=(cos(2*PI/lambda*sqrt((xs-delta/2).^2+(ys-offset).^2+f.^2))+...
    cos(2*PI/lambda*sqrt((xs+delta/2).^2+(ys-offset).^2+f.^2))+...
    cos(2*PI/lambda*sqrt((xs).^2+(ys-offset-delta/2).^2+f.^2))+...
    cos(2*PI/lambda*sqrt((xs).^2+(ys-offset+delta/2).^2+f.^2))+...
    4*cos(2*PI*ys/lambda*sin(incidentAngle))).^2+...
    (-sin(2*PI/lambda*sqrt((xs-delta/2).^2+(ys-offset).^2+f.^2))+...
    -sin(2*PI/lambda*sqrt((xs+delta/2).^2+(ys-offset).^2+f.^2))+...
    -sin(2*PI/lambda*sqrt((xs).^2+(ys-offset-delta/2).^2+f.^2))+...
    -sin(2*PI/lambda*sqrt((xs).^2+(ys-offset+delta/2).^2+f.^2))+...
    4*sin(2*PI*ys/lambda*sin(incidentAngle))).^2;