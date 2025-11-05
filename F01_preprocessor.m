
% run file that makes behavioral table
%analyzeWHY
function F01_preprocessing(arrayIndex)


cd /rds/projects/f/froemerr-accb-lab % go to right folder (RDS drive)

RESOURCEPATH = [pwd(), '/Resources']; % generate Resource path because that's where EEGLAB lives

addpath([RESOURCEPATH, '/EEG/eeglab13_6_5b']) % add EEGLAB to path

cd Experiments/WHY/ % enter main folder


BASEPATH= pwd();% save out current folder
fprintf('%s\n', BASEPATH) % show it to me

% define other folders relative to this one
RAWPATH = strcat(BASEPATH, '/data/EEG_data/raw/');
CLEANPATH = strcat(BASEPATH, '/data/EEG_data/clean/');
EXPORTPATH = strcat(BASEPATH, '/data/export/');

ETPATH = strcat(BASEPATH, '/data/ET_Results/402/');


%% preparation
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; %start EEG lab


% check if number of ET files match number of EEG files else stop right here and tell
% person to fix!
SOURCEFILES = dir(strcat(RAWPATH, '402_1*.bdf')); %all brain vision header files in the folder
ETFILES = dir(strcat(ETPATH, 'WHY1*.asc')); %all brain vision header files in the folder

if numel(SOURCEFILES) == numel(ETFILES)

    fprintf('Sweet! Same number of ET and EEG Files found... Proceeding to Analysis!\n')
    readytorun = 1;

else

    fprintf('Number of ET and EEG Files do not match! Convert ET files first and check EEG files!\n')

    readytorun = 0;

end


SUBJECTS = 1:numel(SOURCEFILES);




%%

if readytorun

    s= arrayIndex;

 %   for s =1:numel(SOURCEFILES)%SUBJECTS
        %% load Data

        fn   = SOURCEFILES(s).name(1:end-4);
        % grab this subject's data from the data table and put into a temporary
        % variable
        fprintf('Processing Participant with ID %s\n', fn);
        %% load EEG data
        EEG = pop_biosig(sprintf('%s%s',RAWPATH, SOURCEFILES(s).name));

        %% fix funny triggers

        EEGMarkers = [EEG.event(:).type];
        
        EEGMarkersFixed = num2cell([EEGMarkers-EEGMarkers(1)]');
        [EEG.event.type] = EEGMarkersFixed{:};

        %% get channel coordinates and remember chanlocs
        EEG=pop_chanedit(EEG, 'lookup',sprintf('%s/EEG/eeglab13_6_5b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp', RESOURCEPATH));

        chanlocs = EEG.chanlocs;
        %% find bad channels

        [EEG, indelec, measure, com] = pop_rejchan( EEG,'threshold',5,'norm','on','measure','prob');
        
        %% interpolate missing channels
        EEG = eeg_interp(EEG, chanlocs);


        %% added to detrend each channel to see if that fixes the analysis/plotting issue
        EEG = pop_rmbase( EEG, [] );

        %% Filtering
        % low-pass filter all raw data once at 100 Hz
        EEG = pop_eegfiltnew(EEG,[],100);


        %% do filtering at 0.0159 Hz (10 sec time constant)    
        % Olaf is right that this is ridiculously slow! We don't seem to
        % have the ERPLAB butterworth filter, though, so we are doing this
        % anyway. BEAR can handle!
        LOWCUTOFF = 0.015915494;
        EEG = pop_eegfiltnew(EEG,LOWCUTOFF,[]); 

        %% Rereference

        %load chanloc information
        EEG=pop_chanedit(EEG, 'lookup',sprintf('%s/EEG/eeglab13_6_5b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp', RESOURCEPATH));
        %re-reference to average
        EEG = pop_reref( EEG, []);


        %% match with eye-tracking data

        %%
        % load ET info & save as .mat
        % parseeyelink(ogFileNamewithLocation, NewFileWithLocation, Keyword)

% ==========================================GEMINI PART=================================================================================

        correct_asc_filename = ETFILES(s).name; %fcked past this point
        full_asc_filepath = [ETPATH, correct_asc_filename];
        output_mat_filepath = [ETPATH, correct_asc_filename(1:end-4), '.mat'];
        ET = parseeyelink(full_asc_filepath, output_mat_filepath, 'Trigger');

% ===========================================GEMINI PART================================================================================

% ET = parseeyelink([ETPATH, 'BASA',fn, '.asc'], [ETPATH, fn, '.mat'],'Trigger'); old code


        % load mat ET info & synchronize with EEG
        %pop_importeyetracker(EEG, matFileToLoad, startEndEvent, ...
        %   importColumns, newLabels, importEyeEvents, doRegression, ...
        %   filterEyetrack, plotFig)
        EEG = pop_importeyetracker(EEG,[ETPATH, fn, '.mat'],[50 60] ,[2:3] ,{'L-GAZE-X' 'L-GAZE-Y'},0,1,0,0);


        % append other trial information
        % requires running behavioral script first.
        ET_CHAN_L = [65 66];

        % after adding all necessary trial info, we can now reject data based
        % on bad eye-tracking
        % EEG = pop_rej_eyecontin(EEG,chns,minvals,maxvals,windowsize)
        EEG = pop_rej_eyecontin(EEG,[ET_CHAN_L] ,[1 1] ,[1920 1080] ,51,2);

        %% detect eye movments 
        % saccade detection parameters (Engbert & Mergenthaler)
        DEG_PER_PIXEL = 0.02;
        VELTHRESH   =   5; % velocity threshold (EXTRA SENSITIVE...)
        MINDUR_SMP  =   5; % minimum saccade duration (in smp)
        PLOTFIG_EM  =   1; % plot EM properties?     
        %EEG = pop_detecteyemovements(EEG,left_eye_xy,right_eye_xy,vfac,mindur,...
        %            degperpixel,smooth,globalthresh,clusterdist,clustermode,
        %            plotfig,writesac,writefix)
        EEG = pop_detecteyemovements(EEG,[ET_CHAN_L],[],VELTHRESH,MINDUR_SMP,DEG_PER_PIXEL,1,0,25,2,PLOTFIG_EM,1,1);
        EEG = eeg_checkset(EEG);


        %save set
        %pop_saveset(EEG, 'filename',sprintf('%s_synchFilt',CLEANPATH,fn));

        %% generate highly filtered copy of this for OC

        HIPASS           = 2;    % Filter's passband edge (in Hz)
        % Best results for scenes were obtained with values of 2 to 2.5 Hz
        % Possibly try even higher value for tasks like Reading
        OW_FACTOR        = 1;    % value for overweighting of SPs (1 = add spike potentials corresponding to 100% of original data length)
        REMOVE_EPOCHMEAN = true; % mean-center the appended peri-saccadic epochs? (strongly recommended)
        EEG_CHANNELS     = 1:64; % indices of all EEG channels (exclude any eye-tracking channels here)
        % I recommend to also include EOG channels (if they were recorded against the common reference)

        fprintf('\nCreating optimized ICA training data...')
        EEG_train = EEG;
        EEG_training = pop_eegfiltnew(EEG_train,HIPASS,[]);

        EEG_training = pop_epoch(EEG_training,{'1'},[-0.2 7]);


        %% Overweight spike potentials
        % Repeatedly append intervals around saccade onsets (-20 to +10 ms) to training data
        EEG_training = pop_overweightevents(EEG_training,'saccade',[-0.02 0.01],OW_FACTOR,REMOVE_EPOCHMEAN);

        %% Run ICA
        fprintf('\nRunning ICA on optimized training data...')
        EEG_training = pop_runica(EEG_training,'extended',1,'interupt','on','chanind',EEG_CHANNELS); % or use binary ICA for more speed


        %% Remember ICA weights & sphering matrix
        wts = EEG_training.icaweights;
        sph = EEG_training.icasphere;

        %% Remove any existing ICA solutions from your original dataset
        EEG.icaact      = [];
        EEG.icasphere   = [];
        EEG.icaweights  = [];
        EEG.icachansind = [];
        EEG.icawinv     = [];


        %% Transfer unmixing weights
        fprintf('\nTransfering ICA weights from training data to original data...')
        EEG.icasphere   = sph;
        EEG.icaweights  = wts;
        EEG.icachansind = EEG_CHANNELS;
        EEG = eeg_checkset(EEG); % let EEGLAB re-compute EEG.icaact & EEG.icawinv

        fprintf('\nIdentifying ocular ICs via saccade/fixation variance-ratio threshold...')

        %% Eye-tracker-guided selection of ICs
        IC_THRESHOLD     = 1.1;   % variance ratio threshold (determined as suitable in Dimigen, 2020)
        SACC_WINDOW      = [5 0]; % saccade window (in samples!) to compute variance ratios (see Dimigen, 2020)
        PLOTFIG          = true;  % plot a figure visualizing influence of threshold setting?
        ICPLOTMODE       = 2;     % plot component topographies (inverse weights)? (2 = only plot "bad" ocular ICs)
        FLAGMODE         = 3;     % overwrite existing rejection flags? (3 = yes)

        %% Automatically flag ocular ICs (Pl�chl et al., 2012)
        [EEG, varratiotable] = pop_eyetrackerica(EEG,'saccade','fixation',SACC_WINDOW,IC_THRESHOLD,FLAGMODE,PLOTFIG,ICPLOTMODE);
        close;
        
        % save out bad components (probably won't find blinks in our data)
        badcomps = EEG.reject.gcompreject;
        
        % run iclabel to identify those blinks
        EEG = pop_iclabel(EEG, 'default');
        % flag eye-artifacts
        EEG = pop_icflag(EEG, [NaN NaN;NaN NaN;0.9 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
        % combine the eye-tracker ICs with the blink ICs
        badcomps = badcomps' + EEG.reject.gcompreject;
        badcomps(badcomps >1) = 1; % in case both flagged the same components, set any 2s to 1
        %% Remove flagged ocular ICs
        EEG      = pop_subcomp(EEG,find(badcomps)); % remove them


        % low pass filter data at 40Hz
        EEG = applytochannels(EEG,[1:64],'pop_eegfiltnew(EEG,[],40);');

        % In interactive session: use this command to plot the cleaned data
        % and check what it looks like
        % pop_eegplot( EEG, 1, 1, 1);
        pop_saveset(EEG, 'filename',sprintf('%s%s_clean',CLEANPATH,fn));

   % end
    %%
    if exist(sprintf('%schanlocs.mat', EXPORTPATH), "file")==0
        chanlocs=EEG.chanlocs;
        save(sprintf('%schanlocs.mat', EXPORTPATH), 'chanlocs');
    end

end

end
