
% function displayModel(folder)
%
% displayModel plots the acceleration recorded in the modelling trials
% stored in the given folder. Acceleration data are decoded and filtered
% with median filtering.
%
% Input:
%   folder --> name of the folder containing the dataset to be displayed
%
% Output:
%   ---
%
% Example:
%%folder = 'Climb_stairs/';
%%displayModel(folder);

% READ THE ACCELEROMETER DATA FILES
folder = 'des_stairs/';
files = dir([folder,'*.txt']);
numFiles = length(files);
dataFiles = zeros(1,numFiles);
noisy_x = [];
noisy_y = [];
noisy_z = [];
for i=1:1:numFiles
    dataFiles(i) = fopen([folder files(i).name],'r');
    data = fscanf(dataFiles(i),'%d\t%d\t%d\n',[3,inf]);

    % Fix the array sizes for data vectors of differing lengths
    noisy_x = padarray(noisy_x, [0,max(size(data,2)-size(noisy_x,2),0)],0,'post');
    noisy_y = padarray(noisy_y, [0,max(size(data,2)-size(noisy_y,2),0)],0,'post');
    noisy_z = padarray(noisy_z, [0,max(size(data,2)-size(noisy_z,2),0)],0,'post');

    % CONVERT THE ACCELEROMETER DATA INTO REAL ACCELERATION VALUES
    % mapping from [0..63] to [-14.709..+14.709]
    noisy_x(i,1:size(data,2)) = -14.709 + (data(1,:)/63)*(2*14.709);
    noisy_y(i,1:size(data,2)) = -14.709 + (data(2,:)/63)*(2*14.709);
    noisy_z(i,1:size(data,2)) = -14.709 + (data(3,:)/63)*(2*14.709);
end
noisy_x = transpose(noisy_x);
noisy_y = transpose(noisy_y);
noisy_z = transpose(noisy_z);

% REDUCE THE NOISE ON THE SIGNALS BY MEDIAN FILTERING

%% Highpass Filter
% All frequency values are in Hz.
Fs = 32;  % Sampling Frequency

Fstop = 0;               % Stopband Frequency
Fpass = 2;               % Passband Frequency
Dstop = 0.0001;          % Stopband Attenuation
Dpass = 0.057501127785;  % Passband Ripple
dens  = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstop, Fpass]/(Fs/2), [0 1], [Dstop, Dpass]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);

filter_datax = filter(Hd,noisy_x);
filter_datay = filter(Hd,noisy_y);
filter_dataz = filter(Hd,noisy_z);

n = 3;      % order of the median filter
x_set = medfilt1(filter_datax,n);
y_set = medfilt1(filter_datay,n);
z_set = medfilt1(filter_dataz,n);
numSamples = length(x_set(: , 1));

A = [x_set;y_set;z_set];
fileID = fopen('filter_data.txt','w');
fprintf(fileID,'%6s %12s %18s\r\n','x_set','y_set','z_set');
fprintf(fileID,'%6.3f %12.3f %18.3f\r\n',A);
fclose(fileID);

% DISPLAY THE RESULTS
time = 1:1:numSamples;

%% Bandpass Filter

% All frequency values are in Hz.
Fs = 32;  % Sampling Frequency
Fstopbp1 = 0;               % First Stopband Frequency
Fpassbp1 = 1;               % First Passband Frequency
Fpassbp2 = 6;               % Second Passband Frequency
Fstopbp2 = 8;               % Second Stopband Frequency
Dstopbp1 = 0.001;           % First Stopband Attenuation
Dpassbp  = 0.057501127785;  % Passband Ripple
Dstop2bp = 0.0001;          % Second Stopband Attenuation
densbp   = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstopbp1 Fpassbp1 Fpassbp2 Fstopbp2]/(Fs/2), [0 1 ...
                         0], [Dstopbp1 Dpassbp Dstop2bp]);

% Calculate the coefficients using the FIRPM function.
bbp  = firpm(N, Fo, Ao, W, {densbp});
Hdbp = dfilt.dffir(bbp);

filter_databp_x = filter(Hdbp,x_set);
filter_databp_y = filter(Hdbp,y_set);
filter_databp_z = filter(Hdbp,z_set);


%% Plot the graphs
% noisy signal
figure,
    subplot(3,1,1)
    plot(time,noisy_x,'-')
    axis([0 numSamples -14.709 +14.709])
    xlabel('time [samples]')
    ylabel('acceleration [m/s^2]')
    title('Noisy accelerations along the x axis')
    subplot(3,1,2)
    plot(time,noisy_y,'-')
    axis([0 numSamples -14.709 +14.709])
    xlabel('time [samples]')
    ylabel('acceleration [m/s^2]')
    title('Noisy accelerations along the y axis')
    subplot(3,1,3)
    plot(time,noisy_z,'-')
    axis([0 numSamples -14.709 +14.709])
    xlabel('time [samples]')
    ylabel('acceleration [m/s^2]')
    title('Noisy accelerations along the z axis')
clean signal
    figure(2),
    subplot(3,1,1)
    plot(time,filter_databp_x,'-')
    axis([0 numSamples -14.709 +14.709])
    xlabel('time [samples]')
    ylabel('acceleration [m/s^2]')
    title('Bandpass Filtered accelerations along the x axis')
    subplot(3,1,2)
    plot(time,filter_databp_y,'-')
    axis([0 numSamples -14.709 +14.709])
    xlabel('time [samples]')
    ylabel('acceleration [m/s^2]')
    title('Bandpass Filtered accelerations along the y axis');
    subplot(3,1,3)
    plot(time,filter_databp_z,'-')
    axis([0 numSamples -14.709 +14.709])
    xlabel('time [samples]')
    ylabel('acceleration [m/s^2]')
    title('Bandpass Filtered accelerations along the z axis')
   

    figure(1),
    subplot(3,1,1)
    plot(time,noisy_x,'-')
    axis([0 numSamples -14.709 +14.709])
    title('Noisy accelerations along the x axis')
    subplot(3,1,2)
    plot(time,filter_data,'-')
    title('High pass data along the x axis')
    spectrogram(filter_datax,kaiser(281,1),0,64,32,'yaxis')
    title('Highpass Filter data set with x dataset')
    figure(2),
    spectrogram(filter_datay,kaiser(281,1),0,64,32,'yaxis')
    title('Highpass Filter data set with y dataset')
    figure(3),
    spectrogram(filter_dataz,kaiser(281,1),0,64,32,'yaxis')
    title('Highpass Filter data set with z dataset')


