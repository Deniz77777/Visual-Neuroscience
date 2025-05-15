close all
clc
A=ETdata.pos(:,4:5);
A=round(A);
img = imread('vy1024x768.png');
imshow(img)
hold on
scatter(A(:,1),A(:,2),'r')

figure



