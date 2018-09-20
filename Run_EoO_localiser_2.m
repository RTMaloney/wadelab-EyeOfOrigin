
function Run_EoO_localiser_2(SubjCode, runNumber)

%The localiser for the EoO stimulus: used to define the area containing the stimulus
%and the region outside it!
%Radial checkerboard pattern that expands logarithmically.

%R Maloney, June 2015

%first, set default input values if not entered:
if nargin < 2
    runNumber = 1;
end

if nargin < 1
    SubjCode = 'test';
end

%%%%-------------------------------%%%%
%           Define parameters:
%%%%-------------------------------%%%%

TestFRs = false; %to test whether the screen(s) are at 120 Hz

PPD = 46.4; %there are 46.39 pixels per degree on the PROpixx at 57cm viewing distance (on average over H/W)

%If you're using the PROpixx or Viewpixx
UsingVP = false;
useHardwareStereo = false;

% work out output filenames:
dateString = datestr(now);
dateString(dateString == ' ') =  '_';
dateString(dateString == '-') =  '_';
dateString(dateString == ':') =  '_';
%Set aside information
subjectData.experimentDescriptor = 'EoO_localiser';
subjectData.subjectCode = SubjCode;
subjectData.runNumber = runNumber;

%File name for the data to be saved as:
subjectData.filename = [ SubjCode, '_', ...  %maybe put the /data prefix in later.
    subjectData.experimentDescriptor, ...
    ['_' , dateString, '_'] ...
    num2str( runNumber ), ...
    '_data.mat'];

%Just abort if the file already exists:
if exist(subjectData.filename,'file')
    userResponse = input('WARNING: Run file exists. Overwrite? Enter y or n: ','s');
    if ~strcmp( userResponse, 'y' )
        subjectData = [];
        error('Aborting function! File exists!');
    end
end

jheapcl; %clear the java heap space.

%%%%-------------------------------%%%%
%           Test display/s
%%%%-------------------------------%%%%

%First up, test the frame rate of the display/s if it's the first run.
%It is crucial for 3D/dichoptic presentations to be displayed at 120Hz on the Vpixx device.

if TestFRs && runNumber == 1
    NumScreens = Screen('Screens');
    for ii=1:length(NumScreens)
        ScreenFRs(ii) = Screen('NominalFrameRate', NumScreens(ii)) %print result to command window
    end
    %If any of the displays are not 120 Hz, provide a warning and check them again using the more accurate
    %'FrameRate' utility. Then abort the program.
    if any(ScreenFRs<120)
        beep
        sprintf('Nominal frame rate is not 120 Hz! \n Obtaining accurate frame rates. Please wait...')
        WaitSecs(1)
        for ii=1:length(NumScreens)
            ScreenFRs(ii) = FrameRate (NumScreens(ii)) %print result to command window
        end
        error ('Check screen display settings!')
    end
end

%Choose the screen: it is usually the max screen no. available.
%Frustratingly, the Shuttle XPC (purchased June 2015) always seems to make the Vpixx display == 1. Not sure why, & can't seem to change it.
%So if we're on that machine, need to -1 from the screen number:
[~, CompName] = system('hostname'); %find out the computer name
if strncmpi(CompName, 'pspcawshuttle', length('pspcawshuttle')) ... %normal strcmp not working here, can't figure out why...
        && (length(Screen('Screens'))>1) %and there is more than 1 display connected...
    WhichScreen = max( Screen( 'Screens' ) )-1;
else
    WhichScreen = max( Screen( 'Screens' ) ); %should be right for any other machine!
end

screenRect = Screen('Rect',WhichScreen); %get the screen resolution.
centreX = screenRect(3)/2;
centreY = screenRect(4)/2;
RefreshRate = Screen('NominalFrameRate', WhichScreen);

%Save the information about the PC hardware:
subjectData.HostPC = CompName;

%Unify the keyboard names in case we run this on a mac:
KbName('UnifyKeyNames')

%%%%----------------------------------------%%%%
%           Define timing parameters:
%%%%----------------------------------------%%%%

%Set a bunch of variables important in determining the timing of stimuli:
DummyTRs = 4 + 1; % dummy volumes at the beginning; we have one spare TR with the current timing parameters so make it one extra dummy I suppose
BlockLengthTRs = 3; %number of TRs in an block
TR = 3; %length of volume (TR), in sec
BlockLengthSec = BlockLengthTRs * TR;

%Set up the block conditions:
%These provide the given condition for each event of the scan, and are 'read in' when displaying the stimuli
LeftEye_Inner = 1;
LeftEye_Outer = 2;
RightEye_Inner = 3;
RightEye_Outer = 4;
Binocular_Inner = 5;
Binocular_Outer = 6;
Fixation = 7; %NOTE: We run twice as many fixation events just to balance out the design:
Dummy = -1; %dummy scans will be same stimulus as fixation/Nulls

%Work out a single cycle:
BlockOrder = [Fixation; LeftEye_Inner; Fixation; LeftEye_Outer; Fixation; RightEye_Inner; Fixation; ...
    RightEye_Outer; Fixation; Binocular_Inner; Fixation; Binocular_Outer];

%Replicate for 3 cycles, and place one final Fixation block at the end:
BlockOrder = [BlockOrder; BlockOrder; BlockOrder; Fixation];

%Insert the durations for each Block:
BlockOrder(:,2) = ones(length(BlockOrder),1) * BlockLengthSec;

%Insert the dummy scans into the beginning of 'BlockOrder':
BlockOrder = [Dummy, DummyTRs*TR; BlockOrder];
numBlocks = length(BlockOrder); %the total number of Blocks, including dummies at the beginning

%set aside for saving:
subjectData.BlockOrderDurn = BlockOrder;

%%%%-------------------------------%%%%
%           Define stimuli:
%%%%-------------------------------%%%%

% Define the dot texture, a square-shaped sheet of dots.
%Make the texture the same size as the height of the screen
imsize = screenRect(4);

x = (1-imsize)/2:(imsize-1)/2;
[x,y] = meshgrid(x,x);
%y = -y; %this doesn't seem to change much, apart from the phase...
[a,r] = cart2pol(x,y);
sf = 3.5; %seems like a good spatial frequency

%radially-expanding sinusoid
mask = square(2*sf*pi*log(r)); %taking the log of r makes the rings expand out logarithmically with the radius

%Plot the Checkerboard:
K=8; %K is the number of full cycles in the radial grating
mask2 = square(K*a); %Make the radial grating
mask3 = (mask.*mask2); %convolve the radial grating with the logarithmic rings

% *** contrast ***
%Specify peak contrast:
Contrast = 0.5; % max checkerboard contrast (used in the temporal windowing).

% define raised cosine temporal window over first & last 300ms of each stimulus event
cont = Contrast.*ones(1,BlockLengthSec * RefreshRate); % initialize to peak dot contrast (defined above)
win_length = BlockLengthSec * RefreshRate/20;     % define window length (5% of total event duration): ~450ms @ 9 sec block
cont(1:win_length) = Contrast.*0.5.*(1-cos(pi*(1:win_length)./win_length)); % ramp up
cont(end:-1:end-win_length+1) = cont(1:win_length); % ... and down
BlockContrast = [cont 0 0 0 0]; %add a few extra zeros in there, just in case it displays the stimulus past the final designated frame,
%before being re-synched to the next block (will crash otherwise). Hopefully these 0s will never be needed!

%set up the raised cosine annular windows. We need these to be at least as big as the screen,
%so that the texture doesn't poke out from behind it

%Do the INNER annulus first (determines inner aperture):
inrad = PPD * 0.25;% inner radius of annulus (in pixels) (for fixation spot)
outrad = PPD * 1.5; %outer radius of annulus (in pixels)
% define extent of spatial raised cosine at edge of aperture (in pixels)
cos_smooth = PPD/4; %make it 1/4 degree
%This should plonk the window in the middle of the matrix, which is what we want
imsize2 = imsize*2; %double the texture size
x0 = (imsize2+1)/2;
y0 = (imsize2+1)/2;
J = ones(imsize2);
for (ii=1:imsize2)
    for (jj=1:imsize2)
        r2 = (ii-x0)^2 + (jj-y0)^2;
        if (r2 > outrad^2)
            J(ii,jj) = 0;
        elseif (r2 < inrad^2)
            J(ii,jj) = 0;
        elseif (r2 > (outrad - cos_smooth)^2)
            J(ii,jj) = cos(pi.*(sqrt(r2)-outrad+cos_smooth)/(2*cos_smooth))^2;
        elseif (r2 < (inrad + cos_smooth)^2)
            J(ii,jj) = cos(pi.*(sqrt(r2)-inrad-cos_smooth)/(2*cos_smooth))^2;
        end
    end
end

%Now make the OUTER annulus:
inrad = PPD * 2;% inner radius of annulus: .5 deg outside the outer edge of the other annulus
outrad = (screenRect(4)/2); %outer radius: apparently as far as it can go. Will be covered by the fixn ring
% define extent of spatial raised cosine at edge of aperture (in pixels)
cos_smooth = PPD/8; %make it 1/8 degree
%This should plonk the window in the middle of the matrix, which is what we want
imsize2 = imsize*2; %double the texture size
x0 = (imsize2+1)/2;
y0 = (imsize2+1)/2;
J2 = ones(imsize2);
for (ii=1:imsize2)
    for (jj=1:imsize2)
        r2 = (ii-x0)^2 + (jj-y0)^2;
        if (r2 > outrad^2)
            J2(ii,jj) = 0;
        elseif (r2 < inrad^2)
            J2(ii,jj) = 0;
        elseif (r2 > (outrad - cos_smooth)^2)
            J2(ii,jj) = cos(pi.*(sqrt(r2)-outrad+cos_smooth)/(2*cos_smooth))^2;
        elseif (r2 < (inrad + cos_smooth)^2)
            J2(ii,jj) = cos(pi.*(sqrt(r2)-inrad-cos_smooth)/(2*cos_smooth))^2;
        end
    end
end

%%%%-------------------------------%%%%
%       Set up the fixation task:
%%%%-------------------------------%%%%

%Set up a 'shortcut' list of values to be scanned by KbCheck.
%These are the only keyboard/Current Designs fORP responses we are interested in.
%This should save a little bit of time with each call to 'KbCheck', because it only needs to scan these values, and ignore others
scanListVals = [KbName('1!'), KbName('2@'), KbName('3#'), KbName('4$'), KbName('5%'), KbName('q'), KbName('5')]; %5 is the trigger, q to quit.
scanList = zeros(1,256); scanList(scanListVals) = 1;

%Set up the fixation cross or spot:
%This is drawn directly to the screen using Screen('FillRect') or ('DrawDots')
%Define size of spot in pixels:
fix_spot_border = 10; %for a spot (if you use it instead of a cross): it's black border
fix_spot_centre = 6; %for the inner portion, that changes in lum with the task

%Define size of cross in pixels:
% crossWidth = 2;
% crossHeight = 10;
% fixationCross = fixation_cross(crossWidth,crossHeight,centreX,centreY);

%Make the fixation lock ring:
ringRadius = screenRect(4)/2; %outer edge (radius) of the ring: the edge of the screen
ringWidth = (screenRect(4)/2) - PPD/3; %1/3 of a degree thick
%Make the ring. It's in a 2*2 checkerboard pattern:
fixationRing = double(checkerboard(screenRect(4)/2,1) > 0.5);
%Define the ring:
xx = (1-imsize)/2:(imsize-1)/2;
[xx,yy] = meshgrid(xx,xx);
[~,r] = cart2pol(xx,yy);
%make the alpha mask for the ring.
ring_alpha = ((r>ringWidth+1) & (r<ringRadius-1)); %Make the alpha mask a tiny bit thinner than the ring itself.

%The fixation cross will randomly alternate between two different shades of grey (above 50% grey)
FixnGreyShift = [0.65 0.65 0.65; 0.8 0.8 0.8];
FixnCol = FixnGreyShift(1,:); %Parsed when fixation cross is drawn: just assign this here now

% generate frames when fixation cross dims / brightens
frameCount = 1;
totalFrames = sum(BlockOrder(:,2)) * RefreshRate; %all the frames for the whole experiment

% in frames: define parameters of the fixation task
earliestTimeBetweenTasks = 1.5 * RefreshRate; %minimum time of change
jitterRange = 7.5 * RefreshRate; %jitter the time between min time & jitterRange+min time

taskFrames = [];
while frameCount < totalFrames
    offset = round((rand(1,1)*jitterRange)+earliestTimeBetweenTasks);
    taskFrames(size(taskFrames,1)+1,1) = frameCount + offset;
    frameCount = frameCount + offset;
end
%To plot a distribution of the fixation task change times:
%FixnChangeDist = diff(taskFrames)
%hist(FixnChangeDist*1000/RefreshRate) %in ms

%%%%-------------------------------%%%%
%           Set up the screen:
%%%%-------------------------------%%%%

try %Start a try/catch statement, in case something goes awry with the PTB functions
    
    % initialization of the display
    AssertOpenGL;
    % Open PTB onscreen window: We request a 32 bit per colour component
    % floating point framebuffer if it supports alpha-blending. Otherwise
    % the system shall fall back to a 16 bit per colour component framebuffer:
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    %Set the color range to be normalised between 0 and 1 (rather than 0-255):
    PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange', 1);
    
    %Initialise the Vpixx device:
    
    if UsingVP        % Enable DATAPixx blueline support, and VIEWPixx scanning backlight for optimal 3D
        
        PsychImaging('AddTask', 'General', 'UseDataPixx');
        Datapixx('Open');
        %The following commands are included in demos that apparently work for both the Viewpixx AND the PROpixx, though they seem specific to the Viewpixx...
        Datapixx('DisableVideoScanningBacklight');    % optionally, turn it off first, in case the refresh rate has changed since startup
        Datapixx('EnableVideoScanningBacklight');     % Only required if a VIEWPixx.
        Datapixx('EnableVideoStereoBlueline');
        Datapixx('SetVideoStereoVesaWaveform', 2);    % If driving NVIDIA glasses
        
        if Datapixx('IsViewpixx3D') %If it's the Viewpixx3D
            
            Datapixx('EnableVideoLcd3D60Hz');
            subjectData.DisplayType = 'Viewpixx3D'; %set aside the device type for reference
            Datapixx('RegWr');
            
        elseif Datapixx('IsPropixx') %if it's the Propixx DLP projector
            
            subjectData.DisplayType = 'PROpixx'; %set aside the device type for reference
            Datapixx('SetPropixxDlpSequenceProgram',0); %set to normal RGB video processing for driving the LEDs & DLP MMDs
            %Datapixx('RegWr');
            
            %Modify the per-eye crosstalk on the PROpixx.
            %Apparently this cross-talk correction only works when using RB3D video mode,
            %where the red/blue channels contain the left/right eye greyscale images (which we are not using).
            %Datapixx('SetPropixx3DCrosstalkLR', 1);
            %Datapixx('SetPropixx3DCrosstalkRL', 1);
            Datapixx('RegWrRd'); %seem to need to do this after setting 'SetPropixxDlpSequenceProgram' to 0
        end
    end
    %No clue what the RegWr and RegWrRd commands are all about, but they are in the demos, so I've included them.
    
    % Open an on screen (grey) window and configure the imaging pipeline
    %Info about the 'blueline' mechanism for synching to the 3D glasses:
    % There seems to be a blueline generation bug on some OpenGL systems.
    % SetStereoBlueLineSyncParameters(windowPtr, windowRect(4)) corrects the
    % bug on some systems, but breaks on other systems.
    % We'll just disable automatic blueline, and manually draw our own bluelines!
    if useHardwareStereo
        [win, windowRect] = PsychImaging('OpenWindow', WhichScreen, 0.5, [], [], [], 1); %flag of 1 engages stereomode
        SetStereoBlueLineSyncParameters(win, windowRect(4)+10);
    else
        [win, windowRect] = PsychImaging('OpenWindow', WhichScreen, 0.5);
    end
    
    %Define the 'blue line' parameters
    blueRectLeftOn   = [0,                 windowRect(4)-1, windowRect(3)/4,   windowRect(4)];
    blueRectLeftOff  = [windowRect(3)/4,   windowRect(4)-1, windowRect(3),     windowRect(4)];
    blueRectRightOn  = [0,                 windowRect(4)-1, windowRect(3)*3/4, windowRect(4)];
    blueRectRightOff = [windowRect(3)*3/4, windowRect(4)-1, windowRect(3),     windowRect(4)];
    
    HideCursor;
    %Do the gamma correction. When using the Viewpixx/PROpixx through the PsychImaging pipeline, we shouldn't use
    %Screen(‘LoadNormalizedGamma’) (see http://www.jennyreadresearch.com/research/lab-set-up/datapixx/)
    %The Vpixx devices should have linear lUTS built in, but we will add this here for completeness.
    %R_gamma =
    %G_gamma =
    %B_gamma =
    %PsychColorCorrection('SetEncodingGamma', win, [1/R_gamma, 1/G_gamma, 1/B_gamma]);
    %raise priority level
    priorityLevel=MaxPriority(win); Priority(priorityLevel);
    %Query the screen refresh rate:
    ifi = Screen('GetFlipInterval',win); %in sec
    
    %Set the alpha-blending: note these settings are the same as with the EoO dots stimulus,
    %but probably don't need to be: yet they work so I won't change them...
    %For linear superposition:
    %Just like the Gabors in GarboriumDemo.m (see there for further info).
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE);
    % We also want alpha-blending for smooth (anti-aliased) textures...
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    %Make the checkerboard texture, and get rectangle coordinates
    CheckTex = Screen('MakeTexture',win,mask3,1,[],2); %set 'optimize for drawing angle' to 1, for fast rotation of animated textures.
    
    %Define a 'destinationRect' for the texture. This needs to be the same size as the 'sourceRect'
    %and hence the same size as the original texture image, because the area was squared during the convolution of the LoG dots into the texture
    TexRect = CenterRectOnPoint([0 0 imsize imsize], centreX, centreY)'; %...don't think we need this...
    
    %Generate the INNER annulus texture:
    AnnImages = 0.5.*ones(imsize2,imsize2,2); %specify RGB matrix of annulus, grey
    AnnImages(:,:,2) = 1.*J; %specify Alpha channel of annulus
    annulusInner = Screen('MakeTexture',win,AnnImages,[],[],2);
    
    %Generate the OUTER annulus texture:
    AnnImages = 0.5.*ones(imsize2,imsize2,2); %specify RGB matrix of annulus, grey
    AnnImages(:,:,2) = 1.*J2; %specify Alpha channel of annulus
    annulusOuter = Screen('MakeTexture',win,AnnImages,[],[],2);
    
    %Generate the ring/fixation lock texture:
    ringMat(:,:,1) = fixationRing;
    ringMat(:,:,2) = ring_alpha;
    fixationRingTexture = Screen('MakeTexture',win,ringMat,[],[],2);
    
    %Here we determine what the random phase of the checkerboard pattern is.
    %Phase will change at a rate of 1 Hz.
    %We will control this using the 'rotationAngle' argument in Screen('DrawTextures').
    %Technically, all blocks are the same length, so it could be a matrix.
    %But for the sake of making it easier, in case different-length trials are inserted later,
    %And because of the 15 sec of dummies at the start, we'll make this a cell array to be indexed on each screen draw.
    %And also, the dummies & fixations are blank so don't need rotation angles, but since they are invisible
    %anyway there is probably no harm in rotating something invisible!
    for ix=1:numBlocks %For each block
        R=[];
        for iz = 1:BlockOrder(ix,2) %for each second of the current block
            %make a vector of a random angle, R, in the range [0, 360] for each frame in that second
            R = [R, ones(1, RefreshRate) * round(rand(1) * 360)];
        end
        RotateTextureBy{ix} = [R, R(end), R(end), R(end), R(end)]; %add a few extras in there, in case it runs over time (& hence won't have an index out of bounds error)
    end
    subjectData.RotateTextureBy = RotateTextureBy; %set aside, just in case.
    
    %%%%-------------------------------%%%%
    %         Begin the experiment!
    %%%%-------------------------------%%%%
    
    %set some variables that increment over frames:
    missedFrames = 0;
    runRecord = [];
    tic; %needed for the keyboard checking sub-routine (see below)
    frame = 1;
    overallFrame = 1; %keep track of the total number of displayed frames.
    Block = 1;
    taskCount = 1; %keep track of change in fixation cross task
    dimmed = false;
    BlockTimes = nan(numBlocks,3);
    pulseTimes=[]; %this is where we will store the time of each scanner pulse
    fixationResponses = zeros(totalFrames,2); %Store the button presses & colour status of the fixation cross
    AnnTexIndex = annulusInner; %just for the first block (ie dummy scans)
    
    %wait for the trigger, obtain a starting time when received
    jheapcl; %clear the java heap space.
    KbCheck(); %take a quick KbCheck to load it now & flush any stored events
    startTime = GetSecs(); %load GetSecs and assign startTime into memory: doing this takes some time on the first occasion
    
    %%%%-------------------------------%%%%
    %         Wait for the trigger
    %%%%-------------------------------%%%%
    
    %Wait for either a TTL pulse from the scanner,
    %or the user presses the '5' key to begin.
    Screen('TextFont',win, 'Arial');
    Screen('TextSize',win, 20);
    %Display welcome screen to both eyes:
    if useHardwareStereo
        Screen('SelectStereoDrawBuffer', win, 0); %flag of 0= left eye
        DrawFormattedText(win,['Welcome ' SubjCode '. \n \nWaiting for trigger... \n \nPress '' 5 '' on keyboard to override, \n \nor '' q '' to quit at any time.'], ...
            'center', 'center', 1);
        Screen('SelectStereoDrawBuffer', win, 1); %flag of 1= right eye
        DrawFormattedText(win,['Welcome ' SubjCode '. \n \nWaiting for trigger... \n \nPress '' 5 '' on keyboard to override, \n \nor '' q '' to quit at any time.'], ...
            'center', 'center', 1);
        Screen('Flip', win); % , [], [], [], 1);
    else
        DrawFormattedText(win,['Welcome ' SubjCode '. \n \nWaiting for trigger... \n \nPress '' 5 '' on keyboard to override, \n \nor '' q '' to quit at any time.'], ...
            'center', 'center', 1);
        Screen('Flip', win); %, [], [], [], 1);
    end
    
    % if in the scanner, need to wait for the trigger before beginning
    start = false;
    while ~start
        StartResp = CheckForResponses (scanList);
        if StartResp ~= 0
            % If it's a five, it means we've either received the scanner trigger,
            % or someone has pressed '5' on the keyboard, so begin
            if StartResp == 5
                start = true; %exit this while loop & return 'startTime': continue playback to begin the program!
                %startTime = GetSecs();
            elseif StartResp == 81 %they've pressed 'q' to quit
                ExitGracefully;
                error( 'Run_EEO_Localiser_2: You quit the program!' );
            end
        end
    end
    
    %sync vbl to startTime
    vbl = Screen('Flip',win); %first flip after pulse
    startTime = vbl;
    %runTimer = GetSecs() - startTime; %stimuli are timed off this start point (should be zero, or very close to it!)
    
    %%%%-------------------------------%%%%
    %         Run through the blocks
    %%%%-------------------------------%%%%
    
    %The direction of the dots is controlled using the 'rotation Angle' parameter in DrawTextures
    while Block <= numBlocks %
        
        %determine correct block end time (in sec)
        if frame == 1
            BlockStartVBL = vbl; %from the last vbl
            %determine the correct length of block:
            BlockEndVBL = BlockStartVBL + (BlockOrder(Block,2)); %  - ifi); %don't subtract one ifi: don't need to! (runs on time on new machines)
            BlockTimes(Block,1) = BlockStartVBL;
            BlockTimes(Block,2) = BlockEndVBL;
        end
        
        %Determine the current block:
        %These can all be controlled by just 2 parameters:
        %contrast and the annulus texture drawn.
        %If a stimulus is absent (either from one eye or both), we simply set contrast to 0 for each frame
        switch BlockOrder(Block,1)
            
            case {Fixation, Dummy}
                
                %Fixation and null ISIs have no stimuli (just fixation)
                BlockContrastL = zeros(RefreshRate * BlockOrder(Block,2) + 3, 1); %add 3 just so there is a little bit extra on the end,
                BlockContrastR = BlockContrastL;                              %& we won't get any 'index out of bounds' errors if there is a presentation deadline or 2 missed...
                
            case Binocular_Inner
                
                BlockContrastL = BlockContrast;
                BlockContrastR = BlockContrast;
                AnnTexIndex = annulusInner;
                
            case Binocular_Outer
                
                BlockContrastL = BlockContrast;
                BlockContrastR = BlockContrast;
                AnnTexIndex = annulusOuter;
                
            case LeftEye_Inner
                
                BlockContrastL = BlockContrast;
                BlockContrastR = BlockContrast .* 0;
                AnnTexIndex = annulusInner;
                
            case LeftEye_Outer
                
                BlockContrastL = BlockContrast;
                BlockContrastR = BlockContrast .* 0;
                AnnTexIndex = annulusOuter;
                
            case RightEye_Inner
                
                BlockContrastL = BlockContrast .* 0;
                BlockContrastR = BlockContrast;
                AnnTexIndex = annulusInner;
                
            case RightEye_Outer
                
                BlockContrastL = BlockContrast .* 0;
                BlockContrastR = BlockContrast;
                AnnTexIndex = annulusOuter;
        end
        
        %Draw each stimulus event
        BlockEnd = false;
        while ~BlockEnd %should keep iterating across stim frames until vbl >= BlockEndVBL
            
            % if this is a frame where the dimming task should change
            if overallFrame == taskFrames(taskCount,1)
                % swap between true / false ie 0=light grey, 1=dark grey
                dimmed = mod(dimmed+1,2);
                taskCount = taskCount + 1;
                FixnCol = FixnGreyShift(dimmed+1, :); %will be either dark or light grey
            end
            
            % Keep record of what that fixation task status was
            % make the dimmed output either 1 or 2 (rather than 0 or 1)
            fixationResponses(overallFrame,1) = dimmed+1; %1=dark, 2=light grey
            
            % *** Draw the checkerboard! ***
            
            % Select left-eye image buffer for drawing:
            if useHardwareStereo
                Screen('SelectStereoDrawBuffer', win, 0);
            end
            
            %%%%------------------------------------------------%%%%
            %               Draw left eye stimulus:
            %%%%------------------------------------------------%%%%
            Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); %turn on alpha blending for lin superposition: see GarboriumDemo.m
            Screen('DrawTextures',win, ...
                CheckTex, ...%determines the texture to be drawn
                [], ... %sourceRect: the displayed subpart of the texture
                [], ... %determines location
                RotateTextureBy{Block}(frame), ... %Determines random checkerboard phase(by rotating the texture)
                [], BlockContrastL(frame)); %Final argument is for contrast (modulates global alpha)
            
            %Superimpose the annulus:
            Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
            Screen('DrawTexture',win,AnnTexIndex);
            Screen('Blendfunction', win, GL_ONE, GL_ZERO);
            
            %Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
            Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
            Screen('DrawTexture',win, fixationRingTexture);
            %Draw the fixation cross:
            %Screen('FillRect',win,FixnCol,fixationCross);
            %Or, if you want to draw a fixation spot (with a border) instead:
            Screen('DrawDots', win, [centreX centreY], fix_spot_border, [0 0 0], [], 2); %flag of 2 makes them rounded: high-quality anti-aliasing.
            Screen('DrawDots', win, [centreX centreY], fix_spot_centre, FixnCol, [], 2);
            Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); %Done drawing so flip alpha back again
            Screen('BlendFunction', win, GL_ONE, GL_ZERO);
            %Draw blue lines:
            Screen('FillRect', win, [0, 0, 1], blueRectLeftOn);
            Screen('FillRect', win, [0, 0, 0], blueRectLeftOff);
            
            % Select right-eye image buffer for drawing:
            if useHardwareStereo
                Screen('SelectStereoDrawBuffer', win, 1);
            else %Not sure what would happen if this 'else' was actually true. I suspect something would go wrong, as stim would be presented twice.
                %But this is more or less how the Vpixx people do it in their 'DatapixxImagingStereoDemo'
                Screen('DrawingFinished', win);
                [vbl , ~ , ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); %update display on next refresh (& provide deadline)
                frame = frame + 1;
                overallFrame = overallFrame + 1;
                %runTimer = GetSecs() - startTime %print out the timer, synched to clock
                runTimer = vbl - startTime; %print out the timer, synched to vbl
                if missed > 0
                    missedFrames = missedFrames + 1;
                end
            end
            
            %%%%------------------------------------------------%%%%
            %         Check for button presses/ scanner pulses
            %%%%------------------------------------------------%%%%
            resp = CheckForResponses (scanList);
            
            % if resp > 0, something has happened
            if resp ~= 0
                % if it was the trigger: store the time
                if resp == 5
                    pulseTimes(size(pulseTimes,1)+1,1) = GetSecs() - startTime; %could use round(GetSecs())
                    % if the 'q' key is pressed, quit the program
                elseif resp == KbName( 'q' )
                    ExitGracefully (UsingVP)
                    error( 'Run_EEO_localiser_2: You quit the program!' );
                else %any other response, should be a button press (or an accident/error)
                    fixationResponses( overallFrame, 2 ) = resp; %store the subject's response
                end
            end
            
            % if this is a frame where the fixation dimming task should change
            if overallFrame == taskFrames(taskCount,1)
                % swap between true / false ie 0=light grey, 1=dark grey
                dimmed = mod(dimmed+1,2);
                taskCount = taskCount + 1;
                FixnCol = FixnGreyShift(dimmed+1, :); %will be either dark or light grey
            end
            
            % Keep record of what that fixation task status was
            % make the dimmed output either 1 or 2 (rather than 0 or 1)
            fixationResponses(overallFrame,1) = dimmed+1; %1=dark, 2=light grey
            
            %%%%------------------------------------------------%%%%
            %               Draw right eye stimulus:
            %%%%------------------------------------------------%%%%
            Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); %turn on alpha blending for lin superposition: see GarboriumDemo.m
            Screen('DrawTextures',win, ...
                CheckTex, ...%determines the texture to be drawn
                [], ... %sourceRect: the displayed subpart of the texture
                [], ... %determines location
                RotateTextureBy{Block}(frame), ... %Determines random checkerboard phase(by rotating the texture)
                [], BlockContrastR(frame)); %Final argument is for contrast (modulates global alpha)
            
            %Superimpose the annulus:
            Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
            Screen('DrawTexture',win,AnnTexIndex);
            Screen('Blendfunction', win, GL_ONE, GL_ZERO);
            
            %Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
            Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
            Screen('DrawTexture',win, fixationRingTexture);
            %Draw the fixation cross:
            %Screen('FillRect',win,FixnCol,fixationCross);
            %Or, if you want to draw a fixation spot (with a border) instead:
            Screen('DrawDots', win, [centreX centreY], fix_spot_border, [0 0 0], [], 2); %flag of 2 makes them rounded: high-quality anti-aliasing.
            Screen('DrawDots', win, [centreX centreY], fix_spot_centre, FixnCol, [], 2);
            Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); %Done drawing so flip alpha back again
            Screen('BlendFunction', win, GL_ONE, GL_ZERO);
            %Draw blue lines:
            Screen('FillRect', win, [0, 0, 1], blueRectRightOn);
            Screen('FillRect', win, [0, 0, 0], blueRectRightOff);
            
            Screen('DrawingFinished', win);
            
            [vbl , ~ , ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); %update display on next refresh (& provide deadline)
            
            %keep record of any missed frames:
            if missed > 0
                missedFrames = missedFrames + 1;
            end
            
            %%%%------------------------------------------------%%%%
            %         Check for button presses/ scanner pulses
            %%%%------------------------------------------------%%%%
            resp = CheckForResponses (scanList);
            
            % if resp > 0, something has happened
            if resp ~= 0
                % if it was the trigger: store the time
                if resp == 5
                    pulseTimes(size(pulseTimes,1)+1,1) = GetSecs() - startTime; %could use round(GetSecs())
                    % if the 'q' key is pressed, quit the program
                elseif resp == KbName( 'q' )
                    ExitGracefully (UsingVP)
                    error( 'Run_EEO_localiser_2: You quit the program!' );
                else %any other response, should be a button press (or an accident/error)
                    fixationResponses( overallFrame, 2 ) = resp; %store the subject's response
                end
            end
            
            %%%%------------------------------------------------%%%%
            %               Prepare for next frame:
            %%%%------------------------------------------------%%%%
            
            %Increment frame and overallFrame
            frame = frame + 1;
            overallFrame = overallFrame + 1;
            %runTimer = GetSecs() - startTime %print out the timer, synched to clock
            runTimer = vbl - startTime; %print out the timer, synched to vbl
            
            % if we've reached the end of the BLOCK
            % synch to time rather than frames to overcome any slight
            % deviation from 120 Hz refresh accumulating over the run
            if vbl >= BlockEndVBL
                frame %print out the frame number it ends on, should be RefreshRate * block durn in sec, +1
                frame = 1;
                BlockTimes(Block,3) = (vbl - startTime); %store the time the event ended (could round this value)
                BlockEnd = true; %this should terminate current block execution
            end
            
        end %end of iterations across frames
        
        Block = Block + 1; %increment the event.
        
    end %end of iterations across events
    
catch
    
    % We throw the error again so the user sees the error description.
    psychrethrow(psychlasterror);
    ExitGracefully(UsingVP)
    error('Error!')
    
end % End of try/catch statement


% Done. Last flip to take end timestamp and for stimulus offset:
%And print diagnostic info to screen:
vbl = Screen('Flip', win, [], [], [], 1);
LastFrame = overallFrame
runTimer = vbl - startTime %should be about 1 frame greater than the correct run time

%Store relevant results, & write to disc:
subjectData.BlockTimes = BlockTimes;
subjectData.missedFrames = missedFrames;
subjectData.runTime = runTimer;
subjectData.LastFrame = LastFrame;
subjectData.pulseTimes = pulseTimes;
subjectData.responses = fixationResponses;
save(subjectData.filename, 'subjectData');

AverageFlipsPerS = overallFrame / (vbl - startTime) %report the average number of flips per sec. Should be close to the refresh rate.
missedFrames
AccuEr = runTimer-sum(BlockOrder(:,2))  %Accumulated error: how much total time the run was off by
ifi
%AccuEr/missedFrames %should be close to 1 frame, in sec

ExitGracefully (UsingVP)

end %end of main function

%Now for the sub-functions:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resp = CheckForResponses(scanList)
% resp is either 1,2,3, or 4 if coming from the Current Designs fibre optic response pad (fORP-932),
% or 5 from the scanner, or the keycode from the keyboard (all as numbers (doubles), not chars)
% Note that all inputs from the fORP-932 are received in the same way as the '1!', '2@', '3#', '4$' & '5%' American qwerty keys
% (not the numeric keypad).
% Find out the codes of different keys using 'KbName' eg KbName(53) = '5%'
% On the response pad, buttons are assigned in a clockwise direction beginning at the right (3 o'clock)

resp = 0;
[ keyIsDown, ~, keyCode ] = KbCheck([],scanList);
%if toc > 0.2 %wait at least 200 ms to re-check for response
if keyIsDown
    
    %tic;
    
    % return a trigger input (coded as a 5) or a key press of either of the 2 '5' keys to begin the program
    if find(keyCode,1) == 53; %KbName('5%') %=53 that's the trigger, or someone's pressed the qwerty 5/% key
        resp = 5;
    elseif find(keyCode,1) == 101; %KbName('5') %=101 someone's pressed 5 on the numeric keypad, that's ok
        resp = 5;
    elseif find(keyCode,1) == 49; %KbName('1!') %=49 button 1 = blue/right button has been pressed
        resp = 1;
    elseif find(keyCode,1) == 50; %KbName('2@') %=50 button 2 = yellow/bottom button has been pressed
        resp = 2;
    elseif find(keyCode,1) == 51; %KbName('3#') %=51 button 3 = green/left button has been pressed
        resp = 3;
    elseif find(keyCode,1) == 52; %KbName('4$') %=52 button 4 = red/top button has been pressed
        resp = 4;
    elseif find(keyCode,1) == 81; %KbName('q') %=81 someone has pressed 'q' key to quit the program!
        resp = 81; %to abort/exit
    else
        resp = 999; %if it's anything else, there's been some error/key pressed accidentally
    end
end
%end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExitGracefully (UsingVP)
%...need to shut everything down here...

%turn off the prioritisation:
Priority( 0 ); %restore priority

if UsingVP        % close down the ViewPixx or ProPixx
    %Datapixx('DisableVideoScanningBacklight');
    if Datapixx('IsViewpixx3D')
        Datapixx('DisableVideoLcd3D60Hz');
    end
    Datapixx('RegWr');
    %Datapixx('Close'); %closing it here might cause it to crash?
end

%Close down the screen:
Screen('CloseAll')

Datapixx('Close'); %closing the Datapixx here (after closing the screen) might stop it from crashing

%Bring back the mouse cursor:
ShowCursor();

end





