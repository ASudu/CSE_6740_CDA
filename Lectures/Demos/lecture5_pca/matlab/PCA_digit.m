% In this example we reduce the digit data(16X16) into 2 dimension 
% using PCA.
clear; 
close all; 
% Loading the upsc digit dataset
load usps_all; 

pixelno = size(data, 1); 
digitno = size(data, 2); 
classno = size(data, 3); 
% Displaying the digit 1(data(:,:,1)) and 0(data(:,:,10)) in the dataset
H = 16; 
W = 16; 
figure; 
show_image([data(:,:,1), data(:,:,10)]', H, W); 
title('digit 1 and 0'); 
% figure; 
% show_image(data(:,:,10)', H, W); 
% title('digit 2'); 

% Extracting the digits 1 and 0 and converting into double
x0 = reshape(data(:,:,[1,10]), [pixelno, digitno*2]); 
x = double(x0); 
y = [ones(1,digitno), 2*ones(1,digitno)]; 

% number of data points to work with; 
m = size(x, 2); 

% Normalize the data
Anew = x'; 
stdA = std(Anew, 1, 1); 
Anew = Anew * diag(1./stdA); 
Anew = Anew'; 

% PCA
% Subtracting the mean of the dataset
mu=sum(Anew,2)./m;
xc = bsxfun(@minus, Anew, mu); 

% Finding the covariance 
C = xc * xc' ./ m; 

% Finding top 2 pricipal component(eigen vector of the covariance)
k = 2; 
[W, S] = eigs(C, k); 
diag(S);
%single(W' * W)

% projecting the 16X16 data on the 2 eigenvectors
dim1 = W(:,1)' * xc ./ sqrt(S(1,1));
dim2 = W(:,2)' * xc ./ sqrt(S(2,2));

% Displaying the the data along the 2 PCs
figure; 
hold on; 
plot(dim1(y==1), dim2(y==1), 'r.'); 
plot(dim1(y==2), dim2(y==2), 'b.'); 
hold off; 
