clear all; close all; clc

%% Part 1 Averaging the Spectrum and Finding Frequency Signature
load subdata.mat % Imports the data as the 262144x49 (space by time) matrix called subdata

L = 10; % spatial domain
n = 64; % Fourier modes

x2 = linspace(-L,L,n+1); x = x2(1:n); y =x; z = x;
k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; ks = fftshift(k); % Creates frequency domain scaled by 2*pi/2L

[X,Y,Z]=meshgrid(x,y,z); % Creates a grid of the spatial domain
[Kx,Ky,Kz]=meshgrid(ks,ks,ks); % Creates a grid of frequency domain

Utn_ave = zeros(n,n,n);

for j=1:49
    Un(:,:,:)=reshape(subdata(:,j),n,n,n); % Reshapes raw data into 64x64x64 matrix at realization j 
    Utn = fftn(Un);  % Computes 3-D Fast Fourier Transform on data. 
    Utn_ave = Utn_ave + Utn; % Sums transformed data
end

Utn_ave = abs(fftshift(Utn_ave))/49; % fftshift becausee it would be difficult to center gaussian 
                                       % filter, the domain needs to be arranged appropiately. 

[Max_ave, Max_index] = max(abs(Utn_ave(:))); % Max_Ave = freq magnitude, Max_index = index of Max_Ave
[X_index, Y_index, Z_index] = ind2sub([n,n,n],Max_index); % Finds the indexes of the max 
                                                                                                                                 
%% Part 2 Applying the Filter and Finding Position of Submarine
tau = 0.2; % Width of filter
center_Kx = Kx(X_index,Y_index,Z_index); % Freq Signature with respect to Kx direction
center_Ky = Ky(X_index,Y_index,Z_index); % Freq Signature with respect to Ky direction
center_Kz = Kz(X_index,Y_index,Z_index); % Freq Signature with respect to Kz direction

filter = fftshift(exp(-tau*((Kx-center_Kx).^2 +(Ky-center_Ky).^2 +(Kz-center_Kz).^2))); 

for j=1:49
    Un(:,:,:)=reshape(subdata(:,j),n,n,n); % Reshapes data into 64x64x64 matrix
    Utn = fftn(Un); % 3-D Fourier Transform on data
    Unft = filter.*Utn; % Apply filter to transformed data
    Unf = ifftn(Unft); % 3-D Inverse Fourier Transform to go back to time domain
    
    [Max_Unf, pos_index] = max(abs(Unf(:))); 
    [X_index, Y_index, Z_index] = ind2sub([n,n,n],pos_index);
    
    x_pos(j) = X(X_index, Y_index, Z_index); % X_coordinate of Submarine
    y_pos(j) = Y(X_index, Y_index, Z_index); % Y_coordinate of Submarine
    z_pos(j) = Z(X_index, Y_index, Z_index); % Z_coordinate of Submarine
    
    plot3(x_pos, y_pos, z_pos,'-o','Color','b','MarkerSize',8)
    axis([-10 10 -10 10 -10 10]), grid on,
    xlabel('x cordinate')
    ylabel('y coordinate')
    zlabel('Depth')
    title('Submarine Trajectory')
    drawnow 
    hold on
end


   
