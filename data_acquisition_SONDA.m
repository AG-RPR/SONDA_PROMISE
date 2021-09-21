clc; sca; close all; clear all; commandwindow;
%% default parameters
if isunix % avoid hardcoding paths so that this trick should be enough to pass from one OS to another
    myPath = [pwd '/'];
else
    myPath = [pwd '\'];
end
cd(myPath); % go to path location
dummymode = 0; % flag to switch on / off dummymode (trial without eye-tracker) ---> SET BEFORE STARTING
subj ='demo'; % string variable for subject name
pursuit_cond = 2; %0: smooth | 1: saccadic | 2: both]
svfd_number = 0; % number of gaze contingent simulated visual field defect conditions
background = 127; % THIS VALUE MUST BE ADAPTED FROM THE CALIBRATION OF THE MONITOR
contrast = [10 100]; % contrast levels, in percentage
contrast_levels = length(contrast); % number of contrast levels
monitWidth = 540; % physical width of the display monitor (mm) ---> MEASURE BEFORE STARTING
viewDist = 600; % viewing distance of the observer (mm) ---> MEASURE BEFORE STARTING
%% custom parameters
default_params = input('Use default paramters? [0: NO| 1: YES]:...');
if default_params == 0
    pursuit_cond =  input('Which pursuit type do you want to test? [0: smooth | 1: saccadic | 2: both]:...');
    background = input('Background level? [greyscale 0-255 | default = 127]:...');
    contrast_levels = input('How many contrast levels? [Integer number]:...');
    contrast = [];
    for c = 1:contrast_levels
        contrast(c) = input(['Contrast #' num2str(c) ' [percentage 0-100]:...']);
    end
    svfd_number = input('How many simulated visual field conditions? [0 = no sVFD | any integer number]:...');
end
%% load trajectories
load('smooth_1.mat');
x_smooth(:,1) = x(1:2400); % short, just for trying it out
y_smooth(:,1) = y(1:2400);

for i = 1:3
load(['saccadic_' num2str(i) '.mat']);
x_sacc(:,i) = x(1:length(x_smooth)); % REMOVE THIS INDEX ONCE THE TIME SERIES ARE READY
y_sacc(:,i) = y(1:length(y_smooth));
end
%% Initialze PTB (PsychToolBox) screen & miscellaneous graphic options
Screen('Preference', 'SkipSyncTests', 1); % skip PTB sync-test
Screen('Preference', 'SkipSyncTests', 2);
Screen('Preference', 'VisualDebuglevel', 3); % skip PTB warning
screenNumber = max(Screen('Screens'));  % get the pointer to the screen (select secondary screen in case of multiple displays)
white = WhiteIndex(screenNumber); % define white (usually 255)
grey=GrayIndex(screenNumber); % define gray (usually 127)
[w, windowRect] = Screen('OpenWindow', 1, background); % open the window object and get a pointer to access that window, initialize the window with a grey background
Screen('Flip', w); % initial flip to clean the screen
ifi = Screen('GetFlipInterval', w); % get interframe interval
hz = 1/ifi; % get refresh rate
MaxPriority(w); % set top priority
[xCenter, yCenter] = RectCenter(windowRect); % get central coordinates (resolution-dependent)
[xScreen, yScreen]= Screen('WindowSize', w); % get screen coordinates (resolution-dependent)
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % enable alpha blending to smooth the edge of the stimulus and the VFD conditions
Screen('TextSize', w, 30); % set default text size to 30
AssertOpenGL; % ensure that OpenGL is running properly before running the data acquisition
%% Stimulus initialization (Gaussian luminance blob)
size_blob_deg = 0.86; % diameter of the stimulus in degrees. The resulting size of the stimulus correspond to a Goldmann III stimulus = 0.43 deg radius
size_blob=round(deg2pix(size_blob_deg,xScreen,monitWidth, viewDist)); % size of the stimulus in pixels
[xg,yg]=meshgrid(-size_blob/2:size_blob/2-1,-size_blob/2:size_blob/2-1); % create the meshgrid of coordinates to be used to create the PTB texture of the stimulus
gauss_sigma = size(xg)/4; gauss_sigma=gauss_sigma(1); % standard deviation of the gaussian aperture to mask the stimulus and giving it smooth edges
gaussian = exp(-((xg/gauss_sigma).^2)-((yg/gauss_sigma).^2)); % gaussian aperture
for c = 1:contrast_levels % define the stimulus value for each contrast level
    stimulus_val(c) = background + background*contrast(c)/100;
    gauss_blob{c}(:, :, 2)=  gaussian * white; % alpha layer of the texture containing the gaussian aperture. 255 = transparent, 0 = opaque
    gauss_blob{c}(:, :, 1) = ones(length(gauss_blob{c}(:,1,2)),length(gauss_blob{c}(:,1,2)))* stimulus_val(c); % background layer of the texture. 255 = white, 0 = black (for colored stimuli, create 4 layers, one for alpha and 3 for RGB)
    tex_blob(c)=Screen('MakeTexture', w, gauss_blob{c}); % create a PTB texture for the stimulus
    % OLD DEFINITION, MIGHT BE WRONG, CHECK gauss_blob(:, :, 1) = ones(length(gauss_blob(:,1,2)),length(gauss_blob(:,1,2)))* round(white*contrast(lum+1)) + grey; % adjust stimulus contrast depending on the current condition. "contrast" is a 1x2 vector
end
rect_blob = [0 0 size_blob size_blob]; % scale the stimulus to the appropriate size in pixel
rectim = [0 0 75 75]; % used for a fake cursor in dummy mode
x_mouse = xCenter; % position mouse coordinates at the center of the screen (useful only in dummy mode)
y_mouse = yCenter;
%% Offset trajectories to make them start from the center
x_smooth = x_smooth + xCenter;
y_smooth = y_smooth + yCenter;
x_sacc = x_sacc + xCenter;
y_sacc = y_sacc + yCenter;
%% Simulating visual field defects
% Empty for now, it will read a matrix calculated with the new HFA-> sVFD
% conversion method and transform it into a PTB texture.
%% Initialize Eyelink
el=EyelinkInitDefaults(w); % create a structure with graphic environment and tracker infos
EyelinkInit(dummymode,1) % initialization, if no eye-tracker is connected, start in dummy mode
[~, vs]=Eyelink('GetTrackerVersion'); % get eye tracker version
switch dummymode
    case 0
        fprintf('Running the test on a ''%s'' tracker.\n', vs ); % if not in dummy mode, print the actual version of the eye-tracker (sometimes it is useful to debug issues)
    case 1
        disp('Running the test in dummy mode [no eyetracker connected]'); % if in dummy mode, warn the user
end
Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA'); % make sure that we get gaze data from the Eyelink (this is some kind of data formatting)
EyelinkDoTrackerSetup(el); % calibrate the eye-tracker (standard 9-point calibration)
EyelinkDoDriftCorrection(el); % do a final calibration check using drift correction
Screen('FillRect',w,background); % clean the screen and fill it with the specified background
Screen('Flip', w);
edfFile=[subj '.edf']; %open file with subject name to record data to
Eyelink('Openfile', edfFile);
eye_used = -1; % -1 because we don't know yet which eye is available. The EyeLink will change this automatically once it picks-up either eye

%% Display
switch pursuit_cond
    case 0
        pursuit_trial = 0;
    case 1
        pursuit_trial = 1;
    case 2
        pursuit_trial = [0 1];
end
Eyelink('Message', 'START_SESSION'); % mark in edf file the global start
for p = 1:length(pursuit_trial)
    if pursuit_trial(p) == 0 % smooth trajectory, single trial
        n_trials = 1;
        pat_string = 'sm';
        x_pos = x_smooth;
        y_pos = y_smooth;
    elseif pursuit_trial(p) == 1 % saccadic trajectory, 3 trials
        n_trials = 3;
        x_pos = x_sacc;
        y_pos = y_sacc;
        pat_string = 'sc';
    end
    for c = 1:contrast_levels
        for t = 1:n_trials
            DrawFormattedText(w, '[LEFT CLICK]: START TRIAL', 'center', yCenter-50, white);
            DrawFormattedText(w, '[RIGHT CLICK]: QUIT (anytime)', 'center', yCenter+50, white);
            Screen(w,'Flip'); % flip to show the commands
            [clicks,~,~,whichButton] = GetClicks(w);
            switch whichButton
                case 1
                case 2
                    sca; return;  
            end
            rect_blob = CenterRectOnPoint(rect_blob,xCenter,yCenter);
            Eyelink('StartRecording'); % start recording eye position
            WaitSecs(0.1); % record a few samples before we actually start displaying
            SetMouse(xCenter,yCenter,w);
            HideCursor;
            Eyelink('Message', ['START_TRIAL_' num2str(t) '_' pat_string '_contrast' num2str(contrast(c))]); % mark zero-plot time in edf fil for current trial
            for i = 1:length(x_pos(:,t))
                rect_blob=CenterRectOnPoint(rect_blob,x_pos(i,t),y_pos(i,t));
                Screen('DrawTexture', w, tex_blob(c),[],rect_blob); % print stimulus on screen
                if dummymode % if NO Eyeleink is available
                    rectim = CenterRectOnPoint(rectim, x_mouse, y_mouse);
                    Screen('FrameOval', w, [0 0 0], rectim, 2, 2);
                    [x_mouse, y_mouse]=GetMouse(w);
                    x_resp(i,t)=x_mouse;
                    y_resp(i,t)=y_mouse;
                else % if the EyeLink is connected
                    if Eyelink( 'NewFloatSampleAvailable') > 0  % check for presence of a new sample update
                        evt = Eyelink( 'NewestFloatSample'); % get the sample in the form of an event structure
                        if eye_used ~= -1 % do we know which eye to use yet?
                            x = evt.gx(eye_used+1); % if we do, get current gaze position from sample
                            y = evt.gy(eye_used+1); % +1 as we're accessing MATLAB array
                            if x~=el.MISSING_DATA && y~=el.MISSING_DATA && evt.pa(eye_used+1)>0  % do we have valid data and is the pupil visible? if yes, store them
                                x_resp(i,t)=x; % store gaze coordinates
                                y_resp(i,t)=y;
                            else % if data is invalid (e.g. during a blink) draw scotoma in the last valid position
                                if i>1
                                    x_resp(i,t) = x_resp(i-1,t);
                                    y_resp(i,t) = y_resp(i-1,t);
                                end
                            end
                        else % if we don't know what eye is being tracked, find out
                            eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
                            if eye_used == el.BINOCULAR % if both eyes are tracked
                                eye_used = el.LEFT_EYE; % use left eye
                            end
                        end
                    end
                end
                vbl = Screen('Flip', w); % flip the screen, updating the position of the stimulus
                [~,~,buttons]=GetMouse(w);
                if any(buttons) && buttons(1)==0
                    sca;
                    return;
                end
            end
            Eyelink('Message', ['END_TRIAL_' num2str(t) '_' pat_string '_c' num2str(contrast(c))]); % print a message in the .edf file to mark the end of the current trial
            % save .mat file with data specific to the current trial
            filename = [myPath subj '_' pat_string '_c' num2str(contrast(c)) '_' num2str(t)]; % concatenate all the condition strings into a single filename
            x_resp = x_resp(:,t); % store only the current trial
            y_resp = y_resp(:,t); 
            x_stim = x_pos(:,t);
            y_stim = y_pos(:,t);
            save(filename, 'ifi','hz', 'pat_string','contrast','monitWidth', 'viewDist', 'myPath', 'subj', 'x_stim', 'x_resp', 'y_stim', 'y_resp', 'xScreen', 'yScreen', 'xCenter', 'yCenter');
        end
    end
end
Eyelink('Message', 'END_SESSION'); % mark in edf file the global end
Eyelink('StopRecording');
if dummymode == 0
    try % try to fetch edf file
        fprintf('Receiving data file ''%s''\n', edfFile);
        status=Eyelink('ReceiveFile'); % try to retrieve edf file (not necessary for the analysis, but useful as a backup)
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2==exist(edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd); % if successful, indicate where the edf can be found
        end
    catch
        fprintf('Problem receiving data file ''%s''\n', edfFile ); % if fail, print error
    end
end
ShowCursor;
ListenChar;
cd(myPath);
sca;