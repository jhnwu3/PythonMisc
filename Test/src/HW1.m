% John Wu
% CSE 5524
% 8/26/2022
%Test the MATLAB image functions to read, display, and write images. Use
%buckeyes_gray.bmp and buckeyes_rgb.bmp from the class materials webpage. [2 pts
grayIm = imread('buckeyes_gray.bmp');
imagesc(grayIm);
axis('image');
colormap('gray');
imwrite(grayIm, 'buckeyes_gray.jpg');
pause;
rgbIm = imread('buckeyes_rgb.bmp');
imagesc(rgbIm);
axis('image');
imwrite(rgbIm, 'buckeyes_rgb.jpg');
pause;
% Read and convert the rgb image to grayscale using NTSC conversion formula
% via rgb2gray also display it (saved just in case users wanted to see.)

grayIm = rgb2gray(rgbIm);
imagesc(grayIm);
axis('image');
imwrite(grayIm, 'buckeyes_rgb_to_gray.png')
pause;

% Test more fully by creating, writing, and reading a checker-board image. [2 pts]
zBlock = zeros(10,10);
oBlock = ones(10,10)*255;
pattern = [zBlock oBlock; oBlock zBlock];
checkerIm = repmat(pattern, 5, 5);
imwrite(uint8(checkerIm), 'checkerIm.bmp');
Im = imread('checkerIm.bmp');
imagesc(Im);
colormap('gray');
axis('image');

