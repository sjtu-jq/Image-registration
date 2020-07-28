function tanvalue = cp_atan(dy,dx)
%CP_ATAN calculate the Inclination angle of one line
tanvalue = zeros([length(dy) 1]);
for i = 1:length(dy)
    if  dy(i)==0 && dx(i)==0
        tanvalue(i) = 90;
        continue;
    else
        tanvalue(i) = round( 180/pi * atan(dy(i)/dx(i)) );
    end
    
end

