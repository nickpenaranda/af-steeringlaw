function in = inarc(x,y,cx,cy,ri,ro,ts,angle)
% function in = inarc(x,y,cx,cy,ri,ro,ts,angle)
%
% A big mess of a function that returns true if point (x,y) lies within
% an arc-shaped region defined by center point (cx,cy), inner bound radius
% ri, outer bound radius ro, starting angle ts (DEGREES) and arc angle te
% (DEGREES)
    [lTheta lRadius] = cart2pol(x - cx, y - cy);
    lTheta = mod(lTheta * 180 / pi,360);
    % Radius check
    %lRadius = sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
    if(lRadius > ro || lRadius < ri)
        in = false;
        return;
    else
        ts = mod(ts,360);
        te = ts + angle;
        if(te > 360)
            te = mod(te,360);
            if(lTheta < ts && lTheta > te)
                in = false;
                return;
            end
        else
            if(lTheta < ts || lTheta > te)
                in = false;
                return;
            end
        end
    end
    in = true;
end