function newdisplayModel(folder)
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
folder = 'Walk/';
files = dir([folder,'Accelerometer-2012-06-11-11-38-12-walk-m1.txt']);
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
n = 3;      % order of the median filter
x_set = medfilt1(noisy_x,n);
y_set = medfilt1(noisy_y,n);
z_set = medfilt1(noisy_z,n);
numSamples = length(x_set(:,1));

%% Spectrogram with kaiser window
% Increasing Beta for kaiser window results in decreased frequency resolution
% Increasing frequency points per segment decreases spectral resolution
% Increasing frequency increases temporal resolution
figure(1),
spectrogram(x_set,kaiser(281,1),0,64,32,'yaxis')  %% spectrogram(data,window,overlap,samples/segments,sampling frequency,yaxis)
title('spectrogram for xset')
figure(2),
spectrogram(y_set,kaiser(114,1),0,64,32,'yaxis')  %% spectrogram(data,window,overlap,samples/segments,sampling frequency,yaxis)
title('spectrogram for yset')
figure(3)
spectrogram(z_set,kaiser(114,1),0,64,32,'yaxis')  %% spectrogram(data,window,overlap,samples/segments,sampling frequency,yaxis)
title('spectrogram for zset')

%% Spectrogram with hamming window
%figure,
%spectrogram(x_set,hamming(114),0,64,32,'yaxis')  %% spectrogram(data,window,overlap,samples/segments,sampling frequency,yaxis)
%title('spectrogram at 114 frequency points per segment')
%figure,
%spectrogram(x_set,hamming(57),0,64,32,'yaxis')  %% spectrogram(data,window,overlap,samples/segments,sampling frequency,yaxis)
%title('spectrogram at 57 freuquency points per segment')

%% Spectrogram with flattop window
%figure,
% Periodic has higher spectral resoultion
% More frequency points per segment changes spectral resolution
%spectrogram(x_set,flattopwin(114,'symmetric'),0,64,32,'yaxis')  %% spectrogram(data,window,overlap,samples/segments,sampling frequency,yaxis)
%title('spectrogram at 114')
%figure,
%spectrogram(x_set,flattopwin(57,'symmetric'),0,64,32,'yaxis')  %% spectrogram(data,window,overlap,samples/segments,sampling frequency,yaxis)
%title('spectrogram at 57 ')  
   
end


