function [longitude,latitude] = myConverter(radius,centX,centY,locX,locY)
% MYCONVERTER convert the location representation of a location from
% row-column coordinates to longitude-latitude coordinate (in degree)
%longitude = asin((locX - centX)/radius)*180/pi;
latitude = asin(  (locY - centY) / radius );
tmp_r =  sqrt(radius^2 - (locY - centY)^2);
longitude = asin( (locX - centX) / tmp_r );

latitude = latitude * 180 /pi;

latitude = -latitude;

longitude = longitude*180/pi;
