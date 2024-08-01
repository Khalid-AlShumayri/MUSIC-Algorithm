clear all 
close all
clc 

% Comparison between three grid-search DoA approaches
% 1- Conventional Beamformer
% 2- Capon's Beamformer
% 3- spectrum-MUSIC algorithm

AntennaData = load('test');
DataMatrix=AntennaData.x;
N=15;    %Number of array elements
K=200;   %Number of data snapshots
d=0.5;   %Distance between elements in wavelengths
R=DataMatrix*DataMatrix'/K;
[U,S,V]=svd(R);
[Q ,D]=eig(R);
[D,I]=sort(diag(D),1,'descend');   %Find r largest eigenvalues
Q=Q (:,I); %Sort the eigenvectors to put signal eigenvectors first
%check the matrix to see how many column with zeros(or (i don't remember exactly and i don't have access to the data) maybe column with small values does it have and those are noises and the rest are signals in this code i had 3 signals )  
Qs=Q (:,1:3); %Get the signal eigenvectors
Qn=Q(:,4:N); %Get the noise eigenvectors
angles=(-90:0.1:90);
%Compute steering vectors corresponding values in angles
a1=exp(-i*2*pi*d*(0:N-1)'*sin([angles(:).']*pi/180));

for k=1:length(angles)
%Compute MUSIC ?spectrum?
music_spectrum(k)= 1./(a1(:,k)'*Qn*Qn'*a1(:,k));     % Music 
srch_bf(k)= (a1(:,k)'*R*a1(:,k))./(a1(:,k)'*a1(:,k));% convential beamformer
srch_cap(k)= 1./(a1(:,k)'*inv(R)*a1(:,k));        % capon's beamformer

end

subplot(311)
plot(angles,abs(music_spectrum)); grid on; 
title('Spatial (MUSIC) Spectrum')
xlabel('Angle in degrees')

subplot(312)
plot(angles,abs(srch_bf));    grid on; 
title('Spatial (conventional beamformer) Spectrum')
xlabel('Angle in degrees')

subplot(313)
plot(angles,abs(srch_cap));    grid on;
title('Spatial (Capon''s beamformer) Spectrum')
xlabel('Angle in degrees')

%% Remark

% we can note from the plot that MUSIC algorithm achieve the best
% performance (detecting the DoAs) in the range from -90 to 90 degrees
% Capon's beamformer is the second best approach
% Conventional beamformer is has the worst performance among the tree 
% teachniques
