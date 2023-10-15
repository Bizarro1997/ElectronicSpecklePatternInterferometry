clear all;
close all;



%%  Load Speckle Patterns
idex={'124';'121'};
img1 = imread(['Image0',idex{1}, '.bmp']); img1= double(img1);
img2 = imread(['Image0',idex{2}, '.bmp']); img2= double(img2);

size_ = size(img1);

%% Average over the set of images to get the systematic background fluctuations
indexRange = 106:1:129;
meanBackground = zeros(size_(1),size_(2));
for m = indexRange
    eval(['img = imread("Image0', num2str(m),'.bmp");']);
    meanBackground = meanBackground + double(img);
end
meanBackground = meanBackground/(129-106+1);

%% Interferrence Pattern
res0 = abs(img1-img2);
res_norm0 = norm__(res0);

figure, imshow(norm__(img1));
title(['Speckle Pattern ',idex{1}]);
figure, imshow(norm__(img2));
title(['Speckle Pattern ',idex{2}]);
figure, imshow(res_norm0);
title('Interferrence Pattern');
I = res_norm0;
imwrite(I,['.\image2\ip_' ,idex{1},'-',idex{2},'.png'  ]); 

%% Interferrence Pattern - Main Systematic Error mitigated
res1 = res0./meanBackground;
res_norm1 = norm__(res1);

figure, imshow(norm__(meanBackground));
title('Mean Background');
figure, imshow(norm__(img1./meanBackground));
title(['Speckle Pattern - Systematic-Error-Corrected ',idex{1}]);
figure, imshow(norm__(img2./meanBackground));
title(['Speckle Pattern - Systematic-Error-Corrected ',idex{2}]);
figure, imshow(res_norm1);
title('Interferrence Pattern - Systematic-Error-Corrected');
I = res_norm1;
imwrite(I,['.\image2\ipSystematicErrorCorrected_' ,idex{1},'-',idex{2}, '.png' ]); 

%% Modelling of the 2pi*nv regions (value approx 0) 
% Play with the parameters valleyThres, d and t
valleyThres = 0.15;
res2_valley_dot = ceil(-res_norm1+valleyThres);
figure, imshow(res2_valley_dot);
title(['Interferrence Pattern - Points at intensity <= ', num2str(valleyThres)] );

res2_valley_stripe = res2_valley_dot;
d = 10;
t = 0.776;

fid = fopen( ['image2\parameters_' ,idex{1},'-',idex{2}, '.txt'], 'a' );
fprintf( fid, 'valleyThres: %f\n', valleyThres);
fprintf( fid, 'd: %f\n', d);
fprintf( fid, 't: %f\n', t);
fclose(fid);

for i=1:1:960
    for j=1:1:1280  
        if res2_valley_dot(i,j) >=0
            meanBackground=0;
            for m=1:1:2*d
                for n=1:1:2*d
                    if i-d+m > 0 && i-d+m < 961 && j-d+n > 0 && j-d+n < 1281 && (-d+m)^2 + (-d+n)^2 <= d^2
                        meanBackground = meanBackground +res2_valley_dot(i-d+m,j-d+n);
                    end
                end
            end
            if meanBackground/(pi*d^2) > t
                res2_valley_stripe(i,j) = 1;
            else
                res2_valley_stripe(i,j) = 0;
                
            end
        end
    end
end
figure, imshow(res2_valley_stripe);
title('Interferrence Pattern - Valley Regions' );


%% Eliminating sigled out regions
res2_valley_stripe_cleaner = res2_valley_stripe;
d = 5;
for i=1:1:960
    for j=1:1:1280    
        if res2_valley_stripe(i,j) > 0
            v = 0;
            for m=1:1:2*d
                for n=1:1:2*d
                    if i-d+m > 0 && i-d+m < 961 && j-d+n > 0 && j-d+n < 1281 && (-d+m)^2 + (-d+n)^2 <= d^2
                        nv_ = res2_valley_stripe_cleaner(i-d+m,j-d+n);
                        if nv_ == 0
                            v = v+1;
                        end
                    end
                end
            end
            if v > 0.5*d^2*pi  
                res2_valley_stripe_cleaner(i,j) = 0; 
            end
        end
        
    end
end
figure, imshow(res2_valley_stripe_cleaner);
title('Interferrence Pattern - Valley Regions Cleaner' );

%% Smoothning Contour
res2_valley_stripe_cleanerNsmooth = res2_valley_stripe_cleaner;
d = 5;
for i=1:1:960
    for j=1:1:1280    
        if res2_valley_stripe_cleaner(i,j) < 1
            v = 0;
            for m=1:1:2*d
                for n=1:1:2*d
                    if i-d+m > 0 && i-d+m < 961 && j-d+n > 0 && j-d+n < 1281 && (-d+m)^2 + (-d+n)^2 <= d^2
                        nv_ = res2_valley_stripe_cleanerNsmooth(i-d+m,j-d+n);
                        if nv_ == 1
                            v = v+1;
                        end
                    end
                end
            end
            if v > 0.65*d^2*pi  
                res2_valley_stripe_cleanerNsmooth(i,j) = 1; 
            end
        end
        
    end
end
figure, imshow(res2_valley_stripe_cleanerNsmooth);
title('Interferrence Pattern - Valley Regions Cleaner and Smoother' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%% From Here on With Color %%%%%%%%%%%%%%%%%%%%%%



%% PEAK Separation of the 2pi*nv regions as a function of nv
x0 = 816;y0=134;
resAux = 1-res2_valley_stripe_cleanerNsmooth;
res3 = 1-res2_valley_stripe_cleanerNsmooth;
lmt1 = [y0, 1, size_(1); y0, 1, size_(1); y0, -1, 1; y0, -1, 1];
lmt2 = [y0, 1, -1; y0, 1, -1; y0,-1, 1; y0, -1, 1];
lmt3 = [x0, 1, size_(2); x0, -1, 1; x0, -1, 1; x0, 1, size_(2)];
for p = 1:1:4
    for i=lmt1(p,1):lmt1(p,2):lmt1(p,3)
        nv=0;
        if resAux(y0,x0) == 1
            nv=1;
        end
        b=0;
        for k=lmt2(p,1):lmt2(p,2):i+lmt2(p,3)
            if resAux(k,x0) > 0
                nv = nv+b;
                b=0;
            end
    
            if resAux(k,x0) < 1
                b=1;
            end
        end
        for j=lmt3(p,1):lmt3(p,2):lmt3(p,3)
            if resAux(i,j) > 0
                nv = nv+b;
                b=0;
                res3(i,j)=nv;
            end
    
            if resAux(i,j) < 1 
                res3(i,j)=0;
                b=1;
            end
        end
    
    end
end



figure,imagesc(res3);
title('Interferrence Pattern - Sepparated Peak Regions' );




res4 = res3;
d = 30;
for i=1:1:960
    for j=1:1:1280  
        if res3(i,j) > 0
            v = zeros(1,100);
            for m=1:1:2*d
                for n=1:1:2*d
                    if i-d+m > 0 && i-d+m < 961 && j-d+n > 0 && j-d+n < 1281 && (-d+m)^2 + (-d+n)^2 <= d^2
                        nv_ = res4(i-d+m,j-d+n);
                        if nv_ > 0
                            v(nv_) = v(nv_)+1;

                        end
                    end
                end
            end
            [~,nv] = max(v);
            res4(i,j) = nv;  
        end
    end
end
figure, imagesc(res4);
title('Interferrence Pattern - Peak Regions' );
res_peak = res4;





%% VALLEY Separation of the 2pi*nv regions as a function of nv
x0 = 816;y0=134;
resAux = res2_valley_stripe_cleanerNsmooth;
res3 = res2_valley_stripe_cleanerNsmooth;
lmt1 = [y0, 1, size_(1); y0, 1, size_(1); y0, -1, 1; y0, -1, 1];
lmt2 = [y0, 1, -1; y0, 1, -1; y0,-1, 1; y0, -1, 1];
lmt3 = [x0, 1, size_(2); x0, -1, 1; x0, -1, 1; x0, 1, size_(2)];
for p = 1:1:4
    for i=lmt1(p,1):lmt1(p,2):lmt1(p,3)
        nv=0;
        if resAux(y0,x0) == 1
            nv=1;
        end
        b=0;
        for k=lmt2(p,1):lmt2(p,2):i+lmt2(p,3)
            if resAux(k,x0) > 0
                nv = nv+b;
                b=0;
            end
    
            if resAux(k,x0) < 1
                b=1;
            end
        end
        for j=lmt3(p,1):lmt3(p,2):lmt3(p,3)
            if resAux(i,j) > 0
                nv = nv+b;
                b=0;
                res3(i,j)=nv;
            end
    
            if resAux(i,j) < 1 
                res3(i,j)=0;
                b=1;
            end
        end
    
    end
end



figure,imagesc(res3);
title('Interferrence Pattern - Sepparated Valley Regions' );




res4 = res3;
d = 30;
for i=1:1:960
    for j=1:1:1280  
        if res3(i,j) > 0
            v = zeros(1,100);
            for m=1:1:2*d
                for n=1:1:2*d
                    if i-d+m > 0 && i-d+m < 961 && j-d+n > 0 && j-d+n < 1281 && (-d+m)^2 + (-d+n)^2 <= d^2
                        nv_ = res4(i-d+m,j-d+n);
                        if nv_ > 0
                            v(nv_) = v(nv_)+1;

                        end
                    end
                end
            end
            [~,nv] = max(v);
            res4(i,j) = nv;  
        end
    end
end
figure, imagesc(res4);
title('Interferrence Pattern - Valley Regions' );
res_valley = res4;


valleyThres = 0.11;
peakThres = 0.25;
valleyThres_multiplier = ceil(-res_norm1+valleyThres);
peakThres_multiplier = ceil(-res_norm1+peakThres);
res_Total = res_peak+res_valley;
res_Total2 = res_peak+res_valley;
for i=1:1:960
    for j=1:1:1280  
        if res_peak(i,j) ~= 0 && res_valley(i,j) == 0
            res_Total2(i,j) = res_peak(i,j)+0.5*resAux(y0,x0);
            res_Total(i,j) = res_Total2(i,j)*peakThres_multiplier(i,j);
        elseif res_peak(i,j) == 0 && res_valley(i,j) ~= 0
            res_Total2(i,j) = res_valley(i,j)+0.5*(1-resAux(y0,x0));
            res_Total(i,j) = res_Total2(i,j)*valleyThres_multiplier(i,j);
        elseif res_peak(i,j) ~= 0 && res_valley(i,j) ~= 0
            res_Total2(i,j) = res_valley(i,j)+0.5*(1-resAux(y0,x0));
            res_Total(i,j) = res_Total2(i,j)*valleyThres_multiplier(i,j);
        else
            res_Total2(i,j) = 0;
            res_Total(i,j) = 0;
        end
    end
end

res_normTotal = norm__(-res_Total2);
figure, imshow(res_normTotal);
title('Interferrence Pattern - Peaks and Valleys');
figure, imagesc(res_Total2);
title('Interferrence Pattern - Peaks and Valleys');
I = res_normTotal;
imwrite(I,['.\image2\final_' ,idex{1},'-',idex{2},'.png' ]);

min_ = 0;
max_ = max(max(res_Total2))
I = (11-res_Total2)/(12);%(3,4.5,6,8,11)
imwrite(I,['.\image2\newNormfinal_' ,idex{1},'-',idex{2},'.png' ]);


function [M] = norm__(M)
%% returns the normalized matrix M (0<=m_ij<=1) 
min_ = min(min(M));
max_ = max(max(M));
M = (M-min_)/(max_-min_);
end