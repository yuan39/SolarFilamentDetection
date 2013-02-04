function surface = mysfit(img,mask,fitstr)
% MYSFIT Polynomial surface fitting a given image
%
% surface = mysfit(img,mask,fitstr) returns a double matrix which is the fitted 
% surface to a given double matrix img based on values specified by mask. The 
% degree of the fitting is specified by fitstr
%
% Example:
%	img = magic(40) 
%   mask = ones(40,40);
%   fitstr='poly44'
%   surface = mysfit(img,mask,fitstr)
%
% Yuan Yuan
% Copyright 2010-2011
%
%
% function mysfit is to approximate a give image img using
% polynomial surface fitting. 
% The input parameters are:
%	1. img : a double matrix which is the target image to be fitted
%   2. mask : a binary (0/1) matrix, where '1's represents the location where the value will be used for fitting.
%   3. fitstr is a string like 'poly44', 'poly43', 'poly23' to denote the degree on x and y direction for fitting.
% 
% The output paramter is:
%	surface : a double matrix reprsenting the fitted result.

     indices = find(mask);
     [rows,cols]=size(img);
     [X,Y] = meshgrid(1:cols,1:rows);
     xdata = reshape(X,numel(X),1);
     ydata = reshape(Y,numel(Y),1);
     zdata = reshape(img,numel(img),1);
     
     zdata = zdata(indices);
     xdata = xdata(indices);
     ydata = ydata(indices);
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 ft = fittype( fitstr );
	 opts = fitoptions( ft );
	 opts.Weights = zeros(1,0);
	 opts.Normalize = 'on';
	 [fitresult, gof] = fit( [xdata, ydata], zdata, ft, opts );

	
     figure( 'Name', 'Surface Fitting of Solar Disk Background' );
     h = plot( fitresult);
    
     shading interp
     colormap gray
     xlabel( 'xdata' );
     ylabel( 'ydata' );
     zlabel( 'zdata' );
     grid on
     view( -82.5, 32 );
     surface = fitresult(X,Y);  