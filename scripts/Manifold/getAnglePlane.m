function angle = getAnglePlane(Vp1,Vp2)
angle = acosd(dot(Vp1,Vp2));
angle(angle > 90) = 180- angle(angle > 90);
