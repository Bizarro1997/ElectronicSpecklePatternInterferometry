function deltaX = correlation2D(idx1,idx2)
%clear all;
close all;

%%  Load and Crop Speckle Patterns

%pixel location of the crop
x1 = 61;
x2 = 121;

idex={num2str(idx1);num2str(idx2)};
img1 = imread(['Image00',idex{1}, '.bmp']); 
img2 = imread(['Image00',idex{2}, '.bmp']); 

img1= double(img1);
img2= double(img2);

img1 = img1(x1:x2,:);
img2 = img2(x1:x2,:);

figure, imshow(norm__(img1));
figure, imshow(norm__(img2));

%uncomment for saving the crops
%I = norm__(img1);
%imwrite(I,'.\fig1.png'); 
%I = norm__(img2);
%imwrite(I,'.\fig2.png'); 


%% Cross correlation
res = xcorr2(img1,img2);
res = norm__(res);
figure, imshow(res);
title(['2D Correlation between picture ',idex{1}, ' and ',idex{2}]);


% x position for analysis
mid = 59;

figure,
plot(res(mid,:))    
title(['Correlation at the middle']);



% filtering the data to return the local maximum - window length = 2*d+1
newRes = res(mid,:);
d = 1; 
for j = d+1:1:length(res(1,:))-d
    newRes(j) = max(res(mid,j-d:1:j+d));
end

figure, plot(newRes);
title(['Correlation at the middle - max among each three']);

% normalization of the cross correlation
t1 = 1279:-1:0;
t2 = 1:1:1279;
t = [t1,t2];

% elimination of possible outliers
M=20;
newNewRes = 1280*newRes(M:end-M)./(abs(1280-t(M:end-M)));
newNewRes = newNewRes/max(newNewRes);
figure, plot(newNewRes)

title(['Correlation at the middle - max among each three - normalized']);
[~,deltaX] = max(newNewRes);
deltaX=deltaX+M-1;


end

function [M] = norm__(M)
%% returns the normalized matrix M (0<=m_ij<=1) 
min_ = min(min(M));
max_ = max(max(M));
M = (M-min_)/(max_-min_);

end