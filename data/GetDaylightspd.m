function daylightspd = GetDaylightspd(CCT,wls)
load DSPD

%linearly interpolate DSPD
DSPD = [wls',interp1(DSPD(:,1),DSPD(:,[2 3 4]),wls,'linear')];

daylightspd = zeros(length(wls),2);

%calculate x_d,y_d based on CCT color temperature
if CCT <= 7000
    xd = .244063 + .09911*(1e3/CCT) + 2.9678*(1e6/(CCT^2)) - 4.6070*(1e9/(CCT^3));
else 
    xd = .237040 + .24748*(1e3/CCT) + 1.9018*(1e6/CCT^2) - 2.0064*(1e9/CCT^3);
end

yd = -3.000*xd^2 + 2.870*xd - 0.275;

%calculate relatative SPD
M = 0.0241 + 0.2562*xd - 0.7341*yd;
M1 = (-1.3515 - 1.7703*xd + 5.9114*yd)/M;
M2 = (0.03000 - 31.4424*xd + 30.0717*yd)/M;

daylightspd(:,1) = wls;
daylightspd(:,2) = DSPD(:,2) + M1*DSPD(:,3) + M2*DSPD(:,4);
end