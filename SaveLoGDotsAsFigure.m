function SaveLoGDotsAsFigure

%Generate a still image/picture of what the LoG dots look like for a publication.
%Based on the parameters of the Eye of Origin experiment dots stimulus.
%This is written to generate a single image (basically a single frame), however, multiple frames are generated so it could
%be used in principle to make an animation and/or sum several frames to generate a 'long exposure' type picture of the motion streaks.
%The aperture sizes can be changed if necessary by altering inrad & outrad
%
%Modified from moving_dots_ryan_fullfield.m
%R Maloney, June 2015


%Define some of the display parameters:
PPD = 46.4; %there are 46.39 pixels per degree on the PROpixx at 57cm viewing distance (on average over H/W)
screenid = max(Screen('Screens')); %get the screen number/identifier
screenRect = Screen('Rect',screenid); %get the screen resolution.
RefreshRate = Screen('NominalFrameRate', screenid);

ScaleBy = 5 %the amount by which you want to rescale/increase the dots for publication

% Define the dot texture, a square-shaped sheet of dots.
%Make the texture the same size as the height of the screen
imsize = screenRect(4); %* ScaleBy; don't scale this if you just want a stylised, enlarged picture of the dots

%Define dot density in dots/(deg^2):
dot_dens_per_deg2 = 3  / ScaleBy; %-> also scale the dot density, because as the dots get bigger, the space between them reduces.

%compute number of dots in the total available area
num_dots = round(dot_dens_per_deg2 * (imsize/PPD)^2) ; %this many dots in the full dot field

%specify dot size:
dot_sigma_in_degrees = 0.05 * ScaleBy; %size of SD of the dot profile in degs/vis angle
dot_sigma = dot_sigma_in_degrees * PPD; %sigma in pixels
dotsize = round(dot_sigma * 8); %make the dots some multiple of sigma
%NOTE: dotsize is simply half the length of the sides of a square patch of pixels that the dot profile is placed within.
%It is dot_sigma that really determines the size of the dots. Obviously, dotsize must be larger than dot_sigma.

%Specify dot speeds:
num_speeds = 1;
dot_speed_lower = 0.2;
dot_speed_upper = 8;
dot_speeds_degs_per_sec = logspace(log10(dot_speed_lower),log10(dot_speed_upper),num_speeds); %degrees travelled in 1 sec
dot_speeds_pixels_per_sec = dot_speeds_degs_per_sec * PPD; %pixels traveled in 1 sec
%Compute dot speed in pixels per frame.
%Because in 3D it is presented L eye, R eye, L eye, R eye...does that mean a single eye sees nothing while the other eye channel is open?
%Does this mean we need to double the dot speeds in pixels/frame, so that the dots can 'catch up' because they do nothing every 2nd frame?
dot_speeds_pixels_per_frame = dot_speeds_pixels_per_sec / RefreshRate; %distance in pixels travelled by a dot in a single video frame.

%Specify peak dot contrast:
pdc = 1; % peak dot contrast (used in the temporal windowing, if there is one).

dirn0 = 0; %full range of motion directions

frames = 150; % motion duration in frames
cos_smooth = dotsize / ScaleBy; %make it one dot size wide

cos_smooth = cos_smooth * 4; %Just double the width of the smoothing: looks better

inrad = PPD * 0.25 * ScaleBy;% inner radius of annulus (in pixels), for fixation spot
outrad = PPD * 1.5 * ScaleBy; %outer radius of annulus (in pixels)

motion_type = 'T'; % needed for call to update_pos.m

%for doing illustrations (eg in a paper):
% inrad=0;
% num_dots=num_dots/2;
% smooth=cos_smooth*10; %(optional)

coherence = ones(length(dirn0),1); %for translational:

J = ones(imsize);

x0 = (imsize+1)/2;
y0 = (imsize+1)/2;

% define colourmap for movie
x = 0:(1/255):1;
cm = [x;x;x]';

% define raised cosine annular window
%Note: pretty sure this is the CORRECT way to make the cosine; not sure why some code differs elsewhere
for ii=1:imsize
    for jj=1:imsize
        r2 = (ii-x0)^2 + (jj-y0)^2;
        if (r2 > (outrad)^2) || r2 < inrad^2
            J(ii,jj) = 0;
        elseif (r2 > (outrad -cos_smooth)^2)
            J(ii,jj) = cos(pi.*(sqrt(r2)-outrad+cos_smooth)/(2*cos_smooth))^2;
        elseif (r2 < (cos_smooth+inrad)^2)
            J(ii,jj) = cos(pi.*(sqrt(r2)-inrad-cos_smooth)/(2*cos_smooth))^2;
        end
    end
end

% define raised cosine temporal window over first & last 250ms
cont = pdc.*ones(1,frames); % initialize to peak dot contrast (defined above)
win_length = frames/10;     % define window length
cont(1:win_length) = pdc.*0.5.*(1-cos(pi*(1:win_length)./win_length)); % ramp up
cont(end:-1:end-win_length+1) = cont(1:win_length); % ... and down

% define Laplacian dot profile ...
G2 = zeros(2*dotsize+1);
var = dot_sigma.^2;
for x = -dotsize:1:dotsize
    for y = -dotsize:1:dotsize
        r = sqrt(x^2 + y^2);
        G0 = 2.5*exp(-(r.^2)./var);
        G2(x+dotsize+1,y+dotsize+1) = (2/var).*((2/var).*(r.^2)-1).*G0;
    end
end
G2 = G2./max(max(abs(G2))); % normalize peak to +/- 1

% define dot polarities (alternating black & white)
dot_col = ones(1,num_dots);
dot_col(1:2:end) = -1;

%allocate memory:
im_matrix = ones(imsize,imsize,frames,'uint8')*128; %the blank space holding the stimuli

for f = 1:frames
    
    % generate dot positions ...
    if (f==1) && motion_type == 'T' %=> for translational motion
        
        dot_pos = imsize.*rand(num_dots,2);
        
    elseif (f==1) && motion_type == 'C' %=> for complex motion
        
    else dot_pos = update_pos(imsize,dot_pos,dot_speeds_pixels_per_frame,dirn0,coherence,motion_type);
    end
    
    % put dots into matrix ...
    I = zeros(imsize);
    for X = 1:num_dots
        I(ceil(dot_pos(X,1)),ceil(dot_pos(X,2))) = dot_col(X);
    end
    
    % convolution to implement Laplacian dot profile...
    I = conv2(I,G2,'same');
    
    % implement spatial windowing: get rid of this if you just want a patch of dots
    I = I.*J;
    
    % transform & implement temporal windowing
    I = round(I.*(127*cont(f))+128);
    
    % clip to range 1-255
    I = min(max(I,1),255);
    
    im_matrix(:,:,f) =uint8(I);
    
end

%-------------------------------------------
%if you want to add a few frames & look at them for illustration, use this code:
%-------------------------------------------
%remember to set inrad=0 & num_dots=num_dots/2; cos_smooth=cos_smooth*10 (optional)
MM=(double(squeeze(im_matrix(:,:,frames/2)))/255) ;% +(double(squeeze(im_matrix(:,:,40)))/255); %+(double(squeeze(im_matrix(:,:,50)))/255); %must convert the unit8 back to double to perform the addition
%MM=MM/2; %normalise the 3 frames (only normalize when summing frames)
figure
imagesc(MM);
axis square
colormap(cm)
set(gca,'XTick',[]) %get rid of the tick marks
set(gca,'YTick',[])


