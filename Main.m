clc
clear all
close all
I    = imread('2.jpg');
SI   = 5;                  
SX   = 6;                  
r    = 1.5;                
sNcut = 0.21;              
sArea = 80;
[Inc, Nnc]   = Ncut(I,SI,SX,r,sNcut,sArea);   
figure()
subplot(231); imshow(I);    title('Original');
subplot(232); imshow(rgb2gray(I)); title('Grayscale');
subplot(233); imshow(Inc);  title(['NormalizedCut',' : ',num2str(Nnc)]); 