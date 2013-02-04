function mask = createMask(side,ratio)
rad2 = (side/2*ratio)^2;
%side = rad*2+1;
mask = zeros(side);
for i = 1:side
    for j = 1:side
        if ( (i-side/2)^2 + (j-side/2)^2 <= rad2 )
            mask(i,j) = 1;
        end
    end
end
        