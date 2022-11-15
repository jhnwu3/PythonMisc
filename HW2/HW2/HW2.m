
% 1. perform gaussian smoothing on the grayscale image affleck_gray.png (or
% affleck_gray_flip.png). Try it for 5 values, sigma = 20, sigma = 10,
% sigma = 5, sigma = 1, sigma = 0.5
sigmas = [20, 10, 5, 1, 0.5];
%sigma=1.0; % use different values
for sigma = sigmas
    G = fspecial('gaussian', 2*ceil(3*sigma)+1, sigma); % create gaussian filter of 3sigma
    faceIm = double(imread('affleck_gray.png')); % read in image with double precision
    gIm = imfilter(faceIm, G, 'replicate'); % now let's filter it
    imshow(gIm/255); % double images need range of 0-1
    imwrite(uint8(gIm), 'gIm_Sigma'+ string(sigma) + '.bmp'); % convert back to integer
    pause;
end

%2.) Write a function to compute and display the 2D Gaussian derivative masks Gx and Gy for a given sigma.
%  Sample the Gaussian derivative/gradient (2D) equation directly (see class notes) at the appropriate x,y (col, row!) locations. Note: each mask is a square 2D matrix. Please ensure that the positive derivative lobe is on the side of the increasing direction of each axis for each mask (see notes regarding the negative sign in the equations). Plot each mask (either in 2D or 3D). [3 pts]
% function definitions below!!

sigma = 5;
[Gx, Gy] = gaussDeriv2D(sigma);
figure
contourf(Gx); % 2D / 3D representations, note all the axis labels to make things clear.
title('Gx');
xlabel('x');
ylabel('y');
saveas(gcf,'gX.png');
pause;
surf(Gx);
saveas(gcf,'gX3D.png');
title('Gx')
xlabel('x')
ylabel('y')

pause;
contourf(Gy)
saveas(gcf,'gY.png');
title('Gy')
xlabel('x')
ylabel('y')

pause;
surf(Gy)
saveas(gcf,'gY3D.png');
title('Gy')
xlabel('x')
ylabel('y')
pause;


% 3. Compute and display the gradient magnitude of an image on the web
myIm = imread('PinkGuy.png');
myIm = double(rgb2gray(myIm)); % convert to gray scale
imagesc(myIm);
colormap("gray");
pause;
gxIm = imfilter(myIm, Gx, 'replicate');
gyIm = imfilter(myIm, Gy, 'replicate');
magIm = sqrt(gxIm.^2 + gyIm.^2); % magnitude
imagesc(gxIm);
pause;
imagesc(gyIm);
pause;
% display magnitude image and save image
imagesc(magIm);
colormap(gray)
saveas(gcf, 'magnitudeImagePinkGuy.png')
pause;

% 4.) Threshold and display the "gradient magnitude” image with different threshold T
% levels. [2 pts] Please note that the 
maxThresh = max(max(magIm));
% pick 4 maxThreshs
Thresholds = [maxThresh / 10, maxThresh / 5, 2*maxThresh/5, 3*maxThresh/5, 4*maxThresh/5];
for T = Thresholds
    tIm = magIm > T;
    imagesc(tIm);
    saveas(gcf, 'pGuyThreshGauss' + string(T) + '.png');
    pause;
end

% 5) Compare the above results with the Sobel masks. [2 pts]
Im = myIm;
for T = Thresholds
    Fx = -fspecial('sobel')';
    fxIm = imfilter(Im,Fx);
    Fy = -fspecial('sobel');
    fyIm = imfilter(Im,Fy);
    magIm = sqrt(fxIm.^2 + fyIm.^2);
    tIm = magIm > T;
    imagesc(tIm);
    saveas(gcf, 'pGuySobelThresh' + string(T) + '.png');
    pause;
end

% 6.) Run the MATLAB canny edge detector, edge(Im,'canny'), on your image and display the default results. How does it compare? (Python: you can use the scikit-image package) [2 pts]
cannyIm = edge(Im, 'canny');
imagesc(cannyIm);
saveas(gcf, 'pGuyCannyEdge.png');
%----------------- FUNCTION DEFINITIONS ----------------- %



%2. Write a function to compute and display the 2D Gaussian derivative masks Gx and Gy for a given sigma.
%  Sample the Gaussian derivative/gradient (2D) equation directly (see class notes) at the appropriate x,y (col, row!) locations. Note: each mask is a square 2D matrix. Please ensure that the positive derivative lobe is on the side of the increasing direction of each axis for each mask (see notes regarding the negative sign in the equations). Plot each mask (either in 2D or 3D). [3 pts]
function dGc = dColumnGauss(r, c, centerColumn, centerRow, sigma) 
    mainTerm = (c - centerColumn) / (2 * pi * sigma^4);
    exponentialTerm = exp(-(((c - centerColumn)^2 + (r - centerRow)^2))/(2*sigma^2));
    dGc = (mainTerm * exponentialTerm);
end

function dGr = dRowGauss(r, c, centerColumn, centerRow, sigma)
    mainTerm = (r - centerRow) / (2 * pi * sigma^4);
    exponentialTerm = exp(-(((c - centerColumn)^2 + (r - centerRow)^2))/(2*sigma^2));
    dGr = mainTerm * exponentialTerm;
end

function [Gx, Gy] = gaussDeriv2D(sigma)
    % use 3 sigma mask size equation
    maskDim = ceil(3*sigma);
    maskSize = 2*maskDim + 1;
    centerRow =  maskDim + 1;
    centerColumn = maskDim + 1;
    Gx = zeros(-maskDim, maskDim); % x == column-wise
    Gy = zeros(-maskDim, maskDim); % y == row-wise
    for r=1:maskSize
        for c=1:maskSize
            Gx(r,c) = dColumnGauss(r, c, centerColumn, centerRow, sigma);
            Gy(r,c) = dRowGauss(r, c, centerColumn, centerRow, sigma);
        end
    end
    Sx = sum(abs(Gx), 'all'); % Probably same value, because square matrix
    Sy = sum(abs(Gy), 'all');
    Gx = Gx / Sx; % make sure mask adds up to 1.
    Gy = Gy / Sy;
end

