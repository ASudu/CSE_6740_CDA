% In this example we first generate a 2 dimensional dataset which is
% distributed along an axis which is 45 degrees to the x-axis. We then
% apply PCA to show that the axis it generates is alinged along this axis.
clear; 
close all; 

G = diag([2, 0.5]);
R = [cos(pi/4), -sin(pi/4); sin(pi/4), cos(pi/4)]; 
m = 1000; 
%Multiplying with G results in scaling of x-axis by 2 and y-axis by 0.5
%Multiplying with R results in rotation of the data by pi/4 rads
x = R * G * randn(2, m) + repmat([0.5;-0.5], 1, m); 
x = [x, R * G * randn(2, m) + repmat([-0.5;0.5], 1, m)]; 
y = [ones(1,m), 2*ones(1,m)]; 

iscolor = 1; 

%Plotting the orginal random data generated
figure; 
hold on; 
if iscolor == 1
  plot(x(1, y==1), x(2, y==1), 'r.'); 
  plot(x(1, y==2), x(2, y==2), 'b.'); 
else
  plot(x(1,:), x(2,:), 'b.'); 
end
axis equal 
hold off; 

drawnow; 

figure; 
subplot(1, 2, 1); 
hist(x(1,:), 100); 
subplot(1, 2, 2); 
hist(x(2,:), 100); 

% number of data points to work with; 
m = size(x, 2); 

input('press key to continue ...\n'); 

% Normalize the the data
Anew = x'; 
stdA = std(Anew, 1, 1); 
Anew = Anew * diag(1./stdA); 
Anew = Anew'; 

% PCA
% Subtracting by the mean
mu=sum(Anew,2)./m;
xc = bsxfun(@minus, Anew, mu); 

% Finding Covariance matrix
C = xc * xc' ./ m; 

% Finding eigenvectors of the Covarance matrix
k = 2; 
[W, S] = eigs(C, k); 
diag(S);
%single(W' * W)

%Pojecting the data along the principal component
dim1 = W(:,1)' * xc ./ sqrt(S(1,1));
dim2 = W(:,2)' * xc ./ sqrt(S(2,2));

% Plotting the data projected along the principal component
figure; 
hold on; 
if iscolor == 1 
  plot(dim1(y==1), dim2(y==1), 'r.'); 
  plot(dim1(y==2), dim2(y==2), 'b.'); 
else
  plot(dim1, dim2, 'b.'); 
end
hold off; 

figure; 
subplot(1, 2, 1); 
hist(dim1, 100); 
subplot(1, 2, 2); 
hist(dim2, 100); 
