for t=theta
    inarc(320+10*cosd(t),240+10*sind(t),320,240,9,11,360-45,90)
    %th = cart2pol(10*cosd(t),10*sind(t)) * 180 / pi;
    %th = mod(th,360);
    %disp(th);
end