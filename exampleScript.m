% Sample script for MATLAB workshop
% 27 January 2016
% William Armstrong

clear all % clear all data in the workspace
close all % close all figure windows

%% Data entry
% Manual data entry
x = [0:10:100]; % creates a vector from 0 to 100, with each datapoint 10 greater than the previous
y1 = [10 2 9 100 3 0 19 28 90 50]; % manually enter values, separate each new value with a space or comma (or semicolon)
y2 = x*2; % multiply the vector x by 2
y3 = x.^2; % raise each entry in x to the power 2
y4 = y3 - y2; % create a new vector by a mathematical operation involving two other vectors

% Reading in simple delimited text files
tempData = dlmread('gglaTempData.csv',',',2,1); % read delimited (dlm) text file named 'gglaTempData.csv' using commas as a delimiter and beginning in the 2nd row, 1st column (to skip the header)

% Reading in delimited text files with numeric and string data
fid = fopen('gglaDataAll.csv'); % open the file and create a file id number
moreTempData = textscan(fid,'%s %s %f %f %f %f %f %f %f %f','HeaderLines',6,'delimiter',','); % read in the data, the %f corresponds to float-type numeric data. %s corresponds to string (text data)
fclose(fid)

% Rename fields from 'moreTempData' for ease of use later
dateString = moreTempData{2}; % date of observation as a string
dateNumber = moreTempData{4}; % date of observation as a matlab datenumber
daysSince = moreTempData{5}; % days since start of the record
temp2 = moreTempData{6}; % air temperature
relh = moreTempData{7}; % relative humidity

% Read in data from .mat file
load tempAndSnow.mat % read data in from mat file named tempAndSnow.mat

%% Processing data

% Make a quick plot to get a feel for data quality
figure(1)
clf
plot(date,snow)
hold on
datetick()

% Remove values greater than a threshold
dataThreshold = 175; % establish threshold
snowClean = snow; % make a new copy of the  original variable to work on
snowClean(snow>dataThreshold) = NaN; % make all values of snow greater than the threshold = no data

plot(date,snowClean,'r') % plot again to see what your cleaned data looks like
%%
% smoothing data
snowSmooth = smooth(date,snowClean,100,'rloess'); % smooth cleaned snow depths using 100 points (50 points on either side) and robust locally weighted smoothing
figure(1)
plot(date,snowSmooth,'k','linewidth',2) % plot the line

%% Plotting

% Plot meteorological data
figure(2) % create a new figure
clf % clear figure
orient landscape % make the plot with a landscape aspect ratio (wider than tall). The opposite of this is 'portrait'
subplot(2,1,1) % create a subplot. There will be 2 rows of figures, 1 column, and we are highlighting the first plot
plot(daysSince,temp2,'linewidth',2,'color','r','linestyle','-','marker','o') % plot temperature as a function of days since the start of the record. Plot a red solid line of width 2, using circles as data markers
text(1,12,'Station = GGLA2','fontsize',16) % annotate the plot at (x,y) = (1,12) with text and set the fontsize

subplot(2,1,2) % create a subplot. There will be 2 rows of figures, 1 column, and we are highlighting the second plot
plot(daysSince,relh,'linewidth',2,'color','b','linestyle','--','marker','None') % plot relative humidity as a function of days since the start of the record. Plot a blue dashed line of width 2, using no data markers

subplot(2,1,1) % highlight the first subplot
ylabel('Air temperature [^{\circ}C]','fontsize',16) % label the y-axis and set the fontsize
xlabel('Days since start of record [d]','fontsize',16) % label the x-axis and set the fontsize
axis([0 60 -20 15]) % set the axis limits. The first two values are the min/max of the x-axis. The second two values are the min/max for the y-axis
title('Air temperature and relative humidity record, beginning April 01, 2014','fontsize',20) % add a title to the plot and set the fontsize
grid on % turn grid lines on

subplot(2,1,2) % highlight the second subplot
ylabel('Relative humidity [%]','fontsize',16) % label the y-axis and set the fontsize
xlabel('Days since start of record [d]','fontsize',16) % label the x-axis and set the fontsize
axis([0 60 0 105]) % set the axis limits. The first two values are the min/max of the x-axis. The second two values are the min/max for the y-axis
grid on % turn grid lines on


%% imagesc() for plotting 3d data
load glacierVelocityData.mat % load data

figure(3) % create new figure
clf % clear figure
orient landscape % make the plot with a landscape aspect ratio (wider than tall). The opposite of this is 'portrait'
h = imagesc(y,z,U);  % create a plot with the handle 'h' (this allows editing later), and plot U as a function of y and z
hold on
caxis([0,0.25]) % set the limits of the colorscale -- 0 at the bottom, 0.25 at the top
v = 0:0.04:0.30; % create vector marking the values at which you want contours drawn
[c i] = contour(y,z,U,v); % draw contours, c and i are handles for future editing
clabel(c,i,v,'LabelSpacing',100) % label the contour lines and space labels out (I don't actually know what the 100 does, I just iteratively modify until it looks good).
set(i,'color',[1 1 1],'linewidth',2) % make the contours white and draw the lines with width 2 ([1,1,1] is the rgb value for white -- like in a rainbow, white actually has all the colors in it).
set(gca,'YDir','normal','fontsize',14) % 'gca' stands for 'get current axis'. This says, "get the current axis and make the y-axis direction normal (increasing upwards), and make the tick labels fontsize 14
set(h,'AlphaData',~isnan(U)) % where there is no velocity data, make the plot transparent (it defaults to the color of the zero value otherwise).
f = colorbar(); % draw a colorbar to show what the colors correspond to
set(get(f,'YLabel'),'string','Velocity [m d^{-1}]','fontsize',16) % give the colorbar a label and set its fontsize
xlabel('Cross-glacier distance [m]','fontsize',18) % label the x-axis and set the fontsize
ylabel('Elevation [m]','fontsize',18) % label the y-axis and set the fontsize
title('Cross-glacier transect of glacier velocity','fontsize',20) % add a title to the plot and set the fontsize





