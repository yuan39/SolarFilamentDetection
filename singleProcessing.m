%function [localStd,boxThreshs] = singleProcessing(...
%    filename,outputfolder,ori_xcent,ori_ycent,ori_rad,forfit,wsize,std_ratio,small_ratio,dilate_size,station,ori_date,p_angle,clip,steps)
 %clip = 0.5
 %wsize=800
 close all
 clear
 clc
 std_ratio = 0.00001;
 forfit = 0;
 steps = 10;
 filename = 'kanz_halph_fc_20100801_081431.fts'
 p_angle = fitsheader(filename,'SOLAR_P0')
 ori_rad = fitsheader(filename,'SOLAR_R')
 ori_ycent = fitsheader(filename,'CENTER_Y')
 ori_xcent = fitsheader(filename,'CENTER_X')
 ori_date = fitsheader(filename,'DATE-OBS')
 station = 'KANZ'
 %outputfolder = 'E:\yuan_backup\Dissertation\SolarPhysics\program';
 outputfolder = 'd:\Dissertation\SolarPhysics\program';
 %wsize = 400;
%filename = 'bbso_halph_fr_20020209_181538.fts';
%ori_xcent = 1021;
%ori_ycent = 1005;
%ori_rad = 915;
%std_ratio = 0.25;
small_ratio = 0.1;
dilate_size = 40;
img = fitsread(filename);
img = flipud(img);
%station = 'BBSO';
%ori_date = '20020209-181538';
%ori_date = fitsheader(filename,'DATE_OBS');
%if( length(ori_date) < 1)
%    ori_date = filename(end-22:end-8);
%end

%if ( ~isempty(strfind(filename,'mod')) )
%    ori_ycent = fitsheader(filename,'CENY');
%    ori_xcent = fitsheader(filename,'CENX');
%    ori_rad = fitsheader(filename,'RAD');
%    p_angle = 0;
%    station = filename(end-36:end-33);
%elseif ( ~isempty(strfind(filename,'kanz')) )
%    p_angle = fitsheader(filename,'SOLAR_P0');
%    ori_rad = fitsheader(filename,'SOLAR_R');
%    ori_ycent = fitsheader(filename,'CENTER_Y');
%    ori_xcent = fitsheader(filename,'CENTER_X');
%    station = 'kanz';
%end
ori_ycent = round(ori_ycent);
ori_xcent = round(ori_xcent);
ori_rad = round(ori_rad);
%make ori_rad a singular number
 if( mod(ori_rad,2) == 0)
     ori_rad= ori_rad - 1;
 end
croped_tmp = img(ori_ycent-ori_rad:ori_ycent+ori_rad,ori_xcent-ori_rad:ori_xcent+ori_rad);
croped = imresize(croped_tmp,[800 800]);
croped = imrotate(croped,-p_angle,'crop');

centerMask = createMask(800,sin(5/180*pi));
croped_Center = croped(centerMask>0);
centerMedian = median(croped_Center);

%if ( ~isempty(strfind(filename,'ynao')) )
%    croped = mat2gray(croped, [centerMedian-1000, centerMedian+1000]);
%elseif( ~isempty(strfind(filename,'kanz')) )
%    croped = mat2gray(croped, [centerMedian-0.5, centerMedian+0.5]);
%else
    %croped = mat2gray(croped, [centerMedian-2000, centerMedian+2000]);
    %croped = mat2gray(croped, [centerMedian-clip, centerMedian+clip]);

    %end
    

figure('Name','Croped')
%imshow(croped,[],'Border','tight','Colormap',colormap(hot),'InitialMagnification',70);
imshow(croped,[],'Border','tight');
hold on
text(2,20,station,'Color','w','FontSize',14);
text(2,785,ori_date,'Color','w','FontSize',13)
hold off
tim = getframe(gca); 
imwrite(tim.cdata,strcat('original','.png'));

%reply = input('Do you want more? Y/N [Y]: ', 's');
%imwrite(croped,strcat(outputfolder,'original.full','.png'));


mask = createMask(800, sin(80/180*pi) );

if forfit > 0
    surface = mysfit(croped,mask,'poly33');
else
    surface = zeros(size(croped));
end

balanced = (croped - surface) .* mask;
balanced(mask==0)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = balanced;
laplacian_filter = [1,1,1;1,-8,1;1,1,1];
g = imfilter(b, laplacian_filter);

b_onmask = b(mask==1);
b_std = std(b_onmask);
disp(['The global std is : ', num2str(b_std)]);
%steps = 30;
medianV = median(b_onmask);
disp(['The global median is : ', num2str(medianV)]);
minV = medianV - 3*std(b_onmask);
disp(['The global minV is : ', num2str(minV)]);
maxV = medianV - 1*std(b_onmask);
%maxV = medianV;
disp(['The global maxV is : ', num2str(maxV)]);
step_size = (maxV - minV)/steps;
t = minV : step_size : maxV;
sumGrad = zeros(1,length(t)-1);
for cc = 1:length(t)-1
    regions = (b < t(cc+1)) - (b < t(cc));
    size_diff = sum(sum( b < t(cc+1)))-sum(sum( b < t(cc)));
    sumGrad(cc) = sum(sum( regions .* g ))/size_diff;
end
% best_thresh = maxHeight(sumGrad);

% if show == 1
%     figure
%     plot(t(1:end-1),sumGrad,'-.or');
%     hold on
%     %plot(t(best_thresh),sumGrad(best_thresh),'*b')
%     %saveas(gcf, strcat(opath,'sumGrad','.png'), 'png');
%     hold off
% end
%value =  t(best_thresh);
%threshed = (b < value) ;
[maxG,maxGind]=max(sumGrad);
thresh = t(maxGind);
disp(['The global thresh is : ', num2str(thresh)]);

[rows,cols]=size(b)
threshed = zeros(rows,cols);
wsize=800
rowsI = floor(rows/wsize);
colsI = floor(cols/wsize);
%localThreshs = zeros(rowsI,colsI);
boxThreshs = zeros(rowsI,colsI);
localStd = zeros(rowsI,colsI);
for i = 1:rowsI
    for j = 1:colsI
        box = b( (i-1)*wsize+1 : wsize*i, (j-1)*wsize+1 : wsize*j);
        gbox = g( (i-1)*wsize+1 : wsize*i, (j-1)*wsize+1 : wsize*j);
        mbox = mask( (i-1)*wsize+1 : wsize*i, (j-1)*wsize+1 : wsize*j);
        if(sum(sum(mbox)) < 0.5*wsize*wsize) % if the effective area is less that half of the block , use global thresh instead
            threshed( (i-1)*wsize+1 : wsize*i, (j-1)*wsize+1 : wsize*j) = box<thresh;
            %localThreshs(i,j) = thresh;
            boxThreshs(i,j) = thresh;
            localStd(i,j) = 0;
            continue;
        end
        box_onmask = box(mbox==1);
        %localStd(i,j) = std(box(:));
        localStd(i,j) = std(box_onmask);
        if( localStd(i,j) < std_ratio*b_std)%if the block is too smooth, probably there is no filament.
        %if(std(box(:)) < 0.015)
            threshed( (i-1)*wsize+1 : wsize*i, (j-1)*wsize+1 : wsize*j) = 0;
            %localThreshs(i,j) = 0;
            boxThreshs(i,j) = 0;
            continue;
        end
        %box_onmask = box(mbox==1);
        %box_std = std(box_onmask);
        box_medianV = median(box_onmask);
        %box_maxV = box_medianV;
        box_maxV = box_medianV - 1*localStd(i,j);
        box_minV = box_medianV - 3*localStd(i,j);
        
        box_step_size = (box_maxV - box_minV)/steps;
        
        box_t = box_minV : box_step_size : box_maxV;
        box_sumGrad = zeros(1,length(box_t)-1);
        for ccc = 1:length(box_t)-1
            box_regions = (box < box_t(ccc+1)) - (box < box_t(ccc));
            box_size_diff = sum(sum( box < box_t(ccc+1)))-sum(sum( box < box_t(ccc)));
            box_sumGrad(ccc) = sum(sum( box_regions .* gbox ))/box_size_diff;
        end
        % best_thresh = maxHeight(sumGrad);
        
        %if show == 1
            %figure
            %plot(box_t(1:end-1),box_sumGrad,'-.or');
            %hold on
            %plot(t(best_thresh),sumGrad(best_thresh),'*b')
            %saveas(gcf, strcat(opath,'sumGrad','.png'), 'png');
            %hold off
        %end
        %value =  t(best_thresh);
        %threshed = (b < value) ;
        [box_maxG,box_maxGind]=max(box_sumGrad);
        box_thresh = box_t(box_maxGind);
        boxThreshs(i,j) = box_thresh;
        %threshed( (i-1)*wsize+1 : wsize*i, (j-1)*wsize+1 : wsize*j) = box< box_minV;%box_thresh;
        %if(box_thresh > medianV-1*b_std || box_thresh < medianV - 2*b_std)
        %    box_thresh = medianV-1.8*b_std;
        %end;
        %box_thresh
        %if(box_thresh<0.5)
            threshed( (i-1)*wsize+1 : wsize*i, (j-1)*wsize+1 : wsize*j) = box< box_thresh;
            %localThreshs(i,j) = box_thresh;
        %else
        %    threshed( (i-1)*wsize+1 : wsize*i, (j-1)*wsize+1 : wsize*j) = 0;
        %end
   
    end
end



% clean
figure( 'Name', 'Before Clean' );
imshow((threshed +(1-mask))>0,[],'Border','tight','InitialMagnification',70);
hold on
text(2,20,station,'FontSize',14);
text(2,785,ori_date,'FontSize',14)
hold off
tim = getframe(gca); 
imwrite(tim.cdata,strcat('before_clean','.png'));

[L,num] = bwlabel(threshed.*mask, 8);
cleaned = zeros(size(threshed));
for i = 1:num
    ind = find(L == i);
    if (length(ind) > small_ratio*size(threshed,1)/2)
        cleaned(ind) = 1;
    end
end
figure( 'Name', 'After Clean' );
imshow((cleaned +(1-mask))>0,[],'Border','tight','InitialMagnification',70);
hold on
text(2,20,station,'FontSize',14);
text(2,785,ori_date,'FontSize',14)
hold off
tim = getframe(gca); 
imwrite(tim.cdata,strcat('after_clean','.png'));


% apply morphology fill to fill the holes
filled = imfill(cleaned>0,'holes');

figure( 'Name', 'After Filling' );
imshow((filled+(1-mask))>0,[],'Border','tight','InitialMagnification',70);
hold on
text(2,20,station,'FontSize',14);
text(2,785,ori_date,'FontSize',14)
hold off
tim = getframe(gca); 
imwrite(tim.cdata,strcat('after_filling','.png'));
%se = strel('disk',dilate_size);
%filled = imdilate(filled,se);

% figure( 'Name', 'Filled Filaments' );
% imshow((filled+(1-mask))>0,[],'Border','tight','InitialMagnification',70);
% hold on
% text(2,785,ori_date,'FontSize',14)
% hold off


se = strel('disk',dilate_size);
dilated = imclose(filled,se);
figure( 'Name', 'Dilated Final Filaments' );
imshow((dilated+(1-mask))>0,[],'Border','tight','InitialMagnification',70);
hold on
text(2,20,station,'FontSize',14);
text(2,785,ori_date,'FontSize',14)
hold off

tim = getframe(gca); 
imwrite(tim.cdata,strcat('dilated_filaments','.png'));

%imwrite(dilated+(1-mask),strcat(outputfolder,'filaments.full','.png'));
%imwrite(dilated>0,strcat(outputfolder,'filaments.full2','.png'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sunspot removal
removed = zeros(size(dilated));
[B,L] = bwboundaries(dilated>0,'noholes');

% Display the label matrix and draw each boundary
figure
imshow(label2rgb(L, @jet, [.5 .5 .5]),'Border','tight')
hold on
for k = 1:length(B)
  boundary = B{k};
  plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end



stats = regionprops(L,'Area','Centroid');

threshold = 0.7;

% loop over the boundaries
for k = 1:length(B)

  % obtain (X,Y) boundary coordinates corresponding to label 'k'
  boundary = B{k};

  % compute a simple estimate of the object's perimeter
  delta_sq = diff(boundary).^2;
  perimeter = sum(sqrt(sum(delta_sq,2)));

  % obtain the area calculation corresponding to label 'k'
  area = stats(k).Area;

  % compute the roundness metric
  metric = 4*pi*area/perimeter^2;
  if (metric < 0.7)
      removed = removed + (L==k);
  end
  % display the results
  metric_string = sprintf('%2.2f',metric);

  % mark objects above the threshold with a black circle
  if metric > threshold
    centroid = stats(k).Centroid;
    %plot(centroid(1),centroid(2),'ro');
  end

  text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','y',...
       'FontSize',14,'FontWeight','bold');

end

title(['Metrics closer to 1 indicate that ',...
       'the object is approximately round']);

tim = getframe(gca); 
imwrite(tim.cdata,strcat('filaments_boundary','.png'));

figure,imshow(removed,[],'Border','tight');

tim = getframe(gca); 
imwrite(tim.cdata,strcat('filaments_sunspot_removal','.png'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%% post processing
% %imshow(dilated,[]);
% %reply = input('Do you want more? Y/N [Y]: ', 's');
% % characterizing
% %dilated = imfill(dilated,'holes');
% 
 [L,num] = bwlabel(removed,8);
 props = regionprops(L,{'Area','Centroid','BoundingBox'});
 % the diameter of the Sun is about 1392000 kilometers
 unit =  1392.000 / length(filled);
 fila_props = zeros(num,4);
% 
 skeletons = zeros(size(removed));
 main_skeletons = zeros(size(removed));
 
 figure( 'Name', 'Characterized Filaments' );
% imshow(croped,[],'Border','tight','InitialMagnification',70);
  imshow(croped,[],'Border','tight');
  hold on 
  %text(2,785,ori_date,'FontSize',14,'BackgroundColor',[1,1,1])
% 
 for i = 1:num
%     %figure,imshow(L==i);
     area = props(i).Area * unit^2;%unit in square kilometers
     %radius = ori_rad/2; centX = ori_rad/2; centY = ori_rad/2;
     %[longitude,latitude] = myConverter(radius,centX,centY,...
     %    props(i).Centroid(1),props(i).Centroid(2));
     [longitude, latitude] = myConverter(length(croped)/2,length(croped)/2,length(croped)/2,...
         props(i).Centroid(1),props(i).Centroid(2));
%     figure;
%     imshow(L==i);
     rectangle('Position', props(i).BoundingBox, 'EdgeColor','y');
     text(props(i).BoundingBox(1) + props(i).BoundingBox(3), ...
          props(i).BoundingBox(2) + props(i).BoundingBox(4),...
          num2str(i),'color','green','FontSize',20);
     hold off
     chared = getframe(gcf);
     imwrite(chared.cdata,strcat('char.png'))
%     [longitude, latitude] = myConverter(length(croped)/2,length(croped)/2,length(croped)/2,...
%         props(i).Centroid(1),props(i).Centroid(2));
%     sn = shapeNumber(L==i);
%     
     % get the skeleton of the i-th component
      skel = bwmorph(L==i,'thin',Inf);
      
%      figure,imshow(skel,'Border','tight','InitialMagnification',70);
      skeletons = skeletons+skel;
% %     if(sum(sum(skel))<3)
% %         continue;
% %     end
%         
%     figure,imshow(skel,[]);
     [main_skel,main_path] = graphCon(skel);
     main_skeletons = main_skeletons + main_skel;
%     figure,imshow(main_skel,[]);
     fila_len = length(main_path)*unit;
%     fila_slope = slope(ori_rad-main_path(1,2),main_path(1,1),...
%         ori_rad-main_path(length(main_path),2),main_path(length(main_path),1));
% %    fila_props(i,:) = [area,longitude,latitude,fila_len,sn,fila_slope];
      fila_props(i,:) = [area,longitude,latitude,fila_len];
%     %1st column store area for sorting
%     %2nd column longitude
%     %3rd column latitude
%     %4th column length
%     %5th column compactness
%     %6th column slope
 end
tim = getframe(gca); 
imwrite(tim.cdata,strcat('filaments_characterized','.png')); 
 
figure,imshow(skeletons +(1-mask),[],'Border','tight');
figure,imshow(main_skeletons +(1-mask),[],'Border','tight');
imwrite( (main_skeletons +(1-mask)).*255,'skeleton.png');

se = strel('disk',1);
main_skeletons_thick = imdilate(main_skeletons,se);
figure,imshow(main_skeletons_thick +(1-mask),[],'Border','tight');
hold on
text(2,20,station,'FontSize',14);
text(2,785,ori_date,'FontSize',14)
hold off
imwrite( (main_skeletons_thick +(1-mask)).*255,'skeleton_thick.png');