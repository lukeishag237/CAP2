%{
Author: NPelot
Created: June 2022
Revision by: EPena

One "dataset" here refers to one recruitment curve, so could have multiple
recruitment curves from a given experiment day

To-do
- Generalize number of CV ranges
- Automate loading of attenuation fits
- Fix attenuation fits to have zero intercept
- Generalize indexing into attenuation coefficients (rather than assuming
freq/10 is correct)
%}


%%
clear
close all
format compact

scrnsz          = get(0,'screensize');

% Add path to Jim Hokanson's MATLAB repo for opening LabChart files (ONLY
% WORKS IN WINDOWS)
addpath('external_dependencies\adinstruments_sdk_matlab-master');

%% Get raw_data_path from YAML file and add that to the 
local_config_file = 'local_config.yml';
fid = fopen(local_config_file);
line_i = fgetl(fid);
raw_data_path = ''; % initialize raw data path as empty
preprocessed_CAP_path = ''; % initialize preprocessed CAP path as empty
local_output_path = ''; % initialize local output path as empty
while isstr(line_i)
    target_starting_pattern1 = 'raw_data_path: ';
    target_pattern_len1 = length(target_starting_pattern1);
    if contains(line_i,target_starting_pattern1)
        assert(strcmp(raw_data_path,''),'%s defined more than once in %s',target_starting_pattern1,local_config_file);
        raw_data_path = line_i((target_pattern_len1+2):(end-1));
    end

    target_starting_pattern2 = 'preprocessed_CAP_path: ';
    target_pattern_len2 = length(target_starting_pattern2);
    if contains(line_i,target_starting_pattern2)
        assert(strcmp(preprocessed_CAP_path,''),'%s defined more than once in %s',target_starting_pattern2,local_config_file);
        preprocessed_CAP_path = line_i((target_pattern_len2+2):(end-1));
    end

    target_starting_pattern3 = 'local_output_path: ';
    target_pattern_len3 = length(target_starting_pattern3);
    if contains(line_i,target_starting_pattern3)
        assert(strcmp(local_output_path,''),'%s defined more than once in %s',target_starting_pattern3,local_config_file);
        local_output_path = line_i((target_pattern_len3+2):(end-1));
    end

    line_i = fgetl(fid);
end
fclose(fid);
assert(~strcmp(raw_data_path,''),'%s not found in %s',target_starting_pattern1,local_config_file);
assert(~strcmp(preprocessed_CAP_path,''),'%s not found in %s',target_starting_pattern2,local_config_file);
assert(~strcmp(local_output_path,''),'%s not found in %s',target_starting_pattern3,local_config_file);
clearvars -except local_config_file raw_data_path preprocessed_CAP_path local_output_path

%%
folder_output_figures  = [local_output_path,'/','OutputFigs'];
if ~exist(folder_output_figures,'dir')
    mkdir(folder_output_figures)
end

folder_output_metrics  = [local_output_path,'/','OutputMetrics'];
if ~exist(folder_output_metrics,'dir')
    mkdir(folder_output_metrics)
end

%% Define analysis parameters
CAP_avg_type            = 2;                    % 1 for mean; 2 for median
scaled_noise            = 1;
flag_attenuation        = 0;                    % Set to 0 to use nominal sinusoidal amplitudes; set to 1 to account for attenuation
subtractnoise           = 1;                    % Turn on/off subtraction of noise using AUC for t < 0
zerobaseline            = 1;                    % Turn on off zeroing of baseline by subtracting average "baseline" of stim-triggered responses from pre-KHF; why? AC-coupled force recordings without a HPF
digitalLPF              = 0;                    % Turn on/off digital LPF
digitalLPF_fc_bandpass  = 0;
digitaloffsetremoval    = 1;

% zerobaseline        = 1;                    

% CV_bounds


%% Define datasets for analysis
%{
1 dataset is 1 recruitment curve
requires ~25 s to analyze each dataset

dataset_subtag:             From each experiment day (expdate), give each recruitment curve a unique number
dataset_stim_polarity:      1 for orig polarity (stim cathode phase-first on rostral contact); 2 for flipped polarity
dataset_electrodetype:      Number of electrode contacts
dataset_monopolarrec:       1 for monopolar rec (each contact recorded on separate channel); 0 for bipolar/tripolar (only one channel with neural signal)
%}

%{
Datasets (all contain orig [1] & flip [2] polarity, unless otherwise
specified; PW means 'pulse width'; CD means 'conduction distance')
SN stim, SN rec
20220511, biphasic,   PW = 200 us
20220518, monophasic, PW = 50, 100, 150 us
20220525, monophasic, PW = 50, 100 us
20220531, monophasic, PW = 50, 100 us,							test-retest after removing/replacing rec cuff (retest only at 50 us; rec cuff same placement as distal SN stim + tibial rec trials)
20220720, monophasic, PW = 50 us,           CD = 13 mm
20221108, monophasic, PW = 50 us,           CD = 10.5 mm

SN stim, tibial rec 
20220525, monophasic, PW = 50 us,           CD = 10, 25 mm
20220531, monophasic, PW = 50 us,           CD = 10, 21.5 mm
20220720, monophasic, PW = 50 us,           CD = 11, 20 mm
20221108, monophasic, PW = 50 us,           CD = 13 mm

cVN stim, cVN rec (fast fibers)
20220616, biphasic,	  PW = 50, 100 us,      CD = 10, 12 mm,     test-retest at CD = 12 mm due to odd noise on the first attempt; CD = 10 mm only tested at 50 us
20220622, biphasic,   PW = 50 us,           CD = 11, 14 mm
20221004, biphasic,   PW = 50 us,           CD = 6, 11, 15 mm
20221103, biphasic,   PW = 50 us,           CD = 5.5, 9 mm

cVN stim, cVN rec (slow fibers)
20220622, biphasic,   PW = 200 us,          CD = 5, 11.5, 14 mm
20221004, biphasic,   PW = 200 us,          CD = 6, 11, 15 mm
20221103, biphasic,   PW = 200 us,          CD = 5.5, 9 mm

20220616
% proximal cVN stim, distal cVN rec
	a.	Trial 6-15, orig polarity, biphasic, 100 us
	22.	Flip it (recorded amps on Fluke)
	a.	Trial 16-25
	24.	Orig polarity, biphasic, 50 us
	a.	Trial 26-35
	25.	Flip it
	a.	Trial 36-45 (recorded some amps on Fluke)

	REPEAT (since we noticed weird offsets previously)
	27.	Repeat with original polarity (fluke not connected still))
	a.	Trial 57-67
	b.	High baseline noise (~15 kHz)
	c.	Offsets shifted again
	28.	Repeat with flip polarity (to run all trials in one go, so they have consistent offset & noise)
	a.	gallamine, flip stim connectiosn, adjust hot hands on belly; disconnected banans from bp probe to check impedance on blue impedance meter
	b.	Trial 68-78
	c.	Looks good! no offset, low noise
	29.	Try orig again, right away
	a.	Trial 79-89
	b.	Still looks good!

	39.	Orig polarity with shroter stim-rec distance; repeatfrom aabove, 50 us, biphasic, same list of amps

	distal cVN stim, distal cVN rec
	a.	Trial 90-100
	b.	Good nerve reseponse, still no overload, but PowerLab again is showing offset + more noise
	c.	Trial 101: repeat 0.05 mA; no change in offset or noise; did not touch ANYThING
	40.	Decide to leave issue still
	41.	Flip polarity without touching anything else
	a.	Trial 102-112
	b.	Same noise & offset as orig polarity above
	c.	seems that “touching�? something(?) makes the changes?

	Rat died before we could repeat above at 100 us; data quality TBD


20220622
% cVN stim, cVN rec, biphasic, PW = 50 us, orig & flip
	19.	Biphasic, 50 us/phase; still using G = 200
	20.	Trial 1-11
	a.	Looks good! Stopping at 0.2 mA for now. Going to flip polarity, and repeat list. Then we’ll change stim-rec distance. And then come back to higher amplitudes to get slower fibers.
	21.	Trial 12-22
	a.	Flip polarity; very careful not to touch anything
	b.	Add Fluke to monitor attenuation (put it on Edgar’s lap and gently pull out the Reseries cable)
	c.	Looks good
	d.	Larger stim artifact than in orig polarity?

% distal cVN stim, distal cVN rec, biphasic, PW = 50 us, orig & flip
	31.	Trial 23-33
	a.	Looks good! Although fastest signal is a bit in the artifact now
	32.	Trial 34-44
	a.	Flipped polarity

% (slow fibers) proximal cVN stim, distal cVN rec, biphasic, PW = 200 us
	46.	Trial 66-80
	a.	Trial 70: Nikki violently stood up and left – and caused ch1 to shift in offset and noise
	b.	Looks like there’s signal at ~1  m/s!!!
	c.	Theerfore, re-do 0.1:0.1:0.5, but also up to 1 mA
	d.	Looks like it saturated the slow signals!
	47.	Flip polarity, trial 81-90
	a.	Looks good!
	
% (slow fibers) distal cVN stim, distal cVN rec, biphasic, PW = 200 us
	48.	Slide stim cuff slightly caudal and repeat
	49.	Trial 91-	100
	a.	Looks good!
	50.	Flip polarity, Trial 101-110
	a.	Looks good!


% (slow fibers) *way* distal cVN stim, distal cVN rec, biphasic, PW = 200 us
	51.	Slide the stim cuff way down, just 2 mm from the rec cuff; also move the green needle to cage to between the cuffs again (see photo)
	52.	Trial 111-119
	a.	Only go to 0.7 mA since that seems to saturate the signal, to try to reduce the likelihood of overload
	b.	Looks great!
	c.	unit 1 had a couple flashes of overload for 0.9 mA
	53.	Drop gain to 100
	54.	Repeat, and then flip polarity
	55.	Trial 120-129
	a.	Looks good
	56.	Flip it, trial 130-139
	a.	Looks good! Artifact might be bleeding in? Maybe good to try blanking?


20220720
SN stim, SN rec
o	Monophasic, 50 us/phase, orig & flip polarity, G = 200, HPF = 100 Hz
	Stim proximal sciatic, rec distal sciatic
	Stim proximal sciatic, rec tibial
	Stim distal sciatic, rec tibial
63.	Start over the day’s plan, but with blanking.
b.	Unit 1 has more noise and big offset.
65.	Trial 77-91
a.	Flip it
b.	Added a couple amplitudes because thresh was higher
66.	Trial 92-106
a.	Repeat orig polarity with the amp list from trial 77-91
64.	Trial 66-76
a.	Looks good

68.	Remove sciatic tripolar recording cuff
69.	Measure stim-rec distance
81.	Trial 153-167
a.	Seems OK?
b.	Kinda boxy shape
82.	Trial 168-184
a.	Flipped
b.	Added some higher stim amps; higher thresholds; some anodic block?

83.	Slide stim cuff from proximal sciatic to distal sciatic; keep preamps on AC since they seem to need some time to settle
84.	Stim-rec distance
85.	Photo
86.	Add some saline, cover cavity
87.	Trial 185-204
a.	Orig polarity
b.	Looks good
88.	Trial 205-224
a.	Flip
b.	Looks good



20221004
cVN

20.	Biphasic, 50 us/phase; still using G = 200
21.	Trial 1-14
a.	Looks good! 
22.	Trial 15-28
a.	Repeat, with Fluke
23.	Trial 29-42
a.	Flip polarity
b.	Looks good!

slow fibers:
25.	Biphasic, 200 us/phase; still using G = 200
26.	Start with only 0.1 mA
28.	Trial 57-70
a.	Looks good!
29.	Trial 71-84
a.	Flip polarity
b.	Record on Fluke
c.	Looks good!

Sliding the stim cuff a bit toward the rec cuff. Re-measuring stim-to-recording distances
1.	Biphasic, 50 us/phase; still using G = 200
2.	Use same list of stim amps as used above for the fast fibers
3.	Trial 85-98
a.	Looks good!
4.	Trial 99-112
a.	Flipped polarity
b.	Looks good!

slow fibers:
5.	Biphasic, 200 us/phase; still using G = 200
6.	Had to declog trach; tried not to nudge cuffs, but in case a change is seen
7.	Use same list of stim amps as used above for the slow fibers
8.	Trial 113-126
a.	Orig polarity
b.	Looks good!
9.	Trial 127-140
a.	Flip polarity
b.	Looks good!



Sliding stim cuff down a bit more

1.	Biphasic, 50 us/phase; drop gain preemptively to G = 100
2.	Use same list of stim amps as used above for the fast fibers
3.	Trial 141-154
a.	Orig
b.	Looks good! (but signal in artifact, as expected)
4.	Trial 155-168
a.	Flip polarity
b.	Looks good! but different peaks for highest 4 amplitudes than those before (maybe artifact carryover)
5.	Biphasic, 200 us/phase; G = 100
6.	Use same list of stim amps as used above for the slow fibers
7.	Trial 169-182
a.	Looks good, from ~3 ms onward
8.	Trial 183-196
a.	Flip





20221103
cVN

20.	Trial 20-31
a.	50 us, biphasic
b.	stim-rec distance of 9 mm
c.	looks good for ch1, although ch2 and ch3 have a slowish signal; maybe the nerve is not touching those contacts very well?
21.	Trial 32-43
a.	flip polarity


slow fibers:
23.	Trial 56-65
a.	200 us, biphasic 
b.	flipped polarity
c.	recorded Fluke attenuation
d.	0.9 and 1 mA slow fiber CAPs seem to be sitting on an exponential artifact
24.	Trial 66-75
a.	orig polarity

fast fibers after sliding cuff caudally:
25.	We considered sliding stim cuff rostrally, and tried to but the stim to rec distance was not clearly larger than  original distance, so instead we just slid the stim cuff caudally only
26.	Short stim-rec distance
27.	Trial 76-87
a.	50 us, biphasic
b.	orig polarity
c.	same list of amps as above for 50 us stim
28.	Trial 88-99
a.	Flipped polarity

slow fibers:
30.	Drop gain to 100
31.	Trial 110-119
a.	Repeat; no overload
b.	200 us, orig poliarty
32.	Trial 120-129
a.	Flip polarity


20221108
SN / tibial

SN stim + SN rec
24.	Trial 34-53
a.	biphasic, 50 us/phase
b.	with Fluke
25.	Trial 54-73
a.	Flip
b.	Looks good!

SN stim + tibial rec
36.	Trial 74-93
a.	bipolar, 50 us, orig polarity
i.	Odd that stim artifact has positive phase first, but confirmed that black is connected to contact 2
37.	Trial 94-113
a.	flip polarity
b.	Notcied rat dead at end of trials



%}

%{
dataset_ind = 0; % Initialize

% ***** 20220511
% SN stim, SN rec, biphasic, PW = 200 us, orig polarity
dataset_ind = dataset_ind + 1; % 1
expdate                                     = '20220511';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 1;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 1; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

% ...                                   , flip polarity
dataset_ind = dataset_ind + 1; % 2
expdate                                     = '20220511';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 2;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 2;
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1;

% ***** 20220518
% SN stim, SN rec, monophasic, PW = 50 us, orig polarity
dataset_ind = dataset_ind + 1; % 3
expdate                                     = '20220518';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 1;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 1; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

% ...                                   , flip polarity
dataset_ind = dataset_ind + 1; % 4
expdate                                     = '20220518';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 2;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 2; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

% SN stim, SN rec, monophasic, PW = 100 us, orig polarity
dataset_ind = dataset_ind + 1; % 5
expdate                                     = '20220518';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 3;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 1; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

% ...                                   , flip polarity
dataset_ind = dataset_ind + 1; % 6
expdate                                     = '20220518';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 4;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 2; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

% SN stim, SN rec, monophasic, PW = 150 us, orig polarity
dataset_ind = dataset_ind + 1; % 7
expdate                                     = '20220518';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 5;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 1; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

% ...                                   , flip polarity
dataset_ind = dataset_ind + 1; % 8
expdate                                     = '20220518';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 6;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 2; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

% ***** 20220525
% SN stim, SN rec, monophasic, PW = 50 us, orig polarity
dataset_ind = dataset_ind + 1; % 9
expdate                                     = '20220525';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 1;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 1; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

% ...                                   , flip polarity
dataset_ind = dataset_ind + 1; % 10
expdate                                     = '20220525';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 2;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 2; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

% SN stim, SN rec, monophasic, PW = 100 us, orig polarity
dataset_ind = dataset_ind + 1; % 11
expdate                                     = '20220525';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 3;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 1; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

% ...                                   , flip polarity
dataset_ind = dataset_ind + 1; % 12
expdate                                     = '20220525';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 4;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 2; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

% Proximal SN stim, tibial rec, monophasic, PW = 50 us, orig polarity
dataset_ind = dataset_ind + 1; % 13
expdate                                     = '20220525';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 5;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 1; 
dataset_electrodetype(dataset_ind)          = 2;
dataset_monopolarrec(dataset_ind)           = 1; 

% ...                                   , flip polarity
dataset_ind = dataset_ind + 1; % 14
expdate                                     = '20220525';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 6;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 2; 
dataset_electrodetype(dataset_ind)          = 2;
dataset_monopolarrec(dataset_ind)           = 1; 

% Distal SN stim, tibial rec, monophasic, PW = 50 us, orig polarity
dataset_ind = dataset_ind + 1; % 15
expdate                                     = '20220525';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 7;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 1; 
dataset_electrodetype(dataset_ind)          = 2;
dataset_monopolarrec(dataset_ind)           = 1; 

% ...                                   , flip polarity
dataset_ind = dataset_ind + 1; % 16
expdate                                     = '20220525';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 8;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 2; 
dataset_electrodetype(dataset_ind)          = 2;
dataset_monopolarrec(dataset_ind)           = 1; 


% ***** 20220531
% SN stim, SN rec, monophasic, PW = 50 us, orig polarity
dataset_ind = dataset_ind + 1; % 17
expdate                                     = '20220531';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 1;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 1; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

% ...                                   , flip polarity
dataset_ind = dataset_ind + 1; % 18
expdate                                     = '20220531';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 2;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 2; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

% SN stim, SN rec, monophasic, PW = 100 us, orig polarity
dataset_ind = dataset_ind + 1; % 19
expdate                                     = '20220531';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 3;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 1; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

% ...                                   , flip polarity
dataset_ind = dataset_ind + 1; % 20
expdate                                     = '20220531';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 4;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 2; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

% Proximal SN stim, tibial rec, monophasic, PW = 50 us, orig polarity
dataset_ind = dataset_ind + 1; % 21
expdate                                     = '20220531';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 5;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 1; 
dataset_electrodetype(dataset_ind)          = 2;
dataset_monopolarrec(dataset_ind)           = 1; 

% ...                                   , flip polarity
dataset_ind = dataset_ind + 1; % 22
expdate                                     = '20220531';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 6;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 2; 
dataset_electrodetype(dataset_ind)          = 2;
dataset_monopolarrec(dataset_ind)           = 1; 

% Proximal SN stim, tibial rec, RETEST, monophasic, PW = 50 us, orig polarity
dataset_ind = dataset_ind + 1; % 23
expdate                                     = '20220531';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 7;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 1; 
dataset_electrodetype(dataset_ind)          = 2;
dataset_monopolarrec(dataset_ind)           = 1; 

% ...                                   , flip polarity
dataset_ind = dataset_ind + 1; % 24
expdate                                     = '20220531';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 8;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 2; 
dataset_electrodetype(dataset_ind)          = 2;
dataset_monopolarrec(dataset_ind)           = 1; 

% Distal SN stim, tibial rec, monophasic, PW = 50 us, orig polarity
dataset_ind = dataset_ind + 1; % 25
expdate                                     = '20220531';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 9;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 1; 
dataset_electrodetype(dataset_ind)          = 2;
dataset_monopolarrec(dataset_ind)           = 1; 

% ...                                   , flip polarity
dataset_ind = dataset_ind + 1; % 26
expdate                                     = '20220531';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 10;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 2; 
dataset_electrodetype(dataset_ind)          = 2;
dataset_monopolarrec(dataset_ind)           = 1; 

% % ***** 20220616
% % fast: 20220616, biphasic,	  PW = 50, 100 us,      CD = 10, 12 mm,     test-retest at CD = 12 mm due to odd noise on the first attempt; CD = 10 mm only tested at 50 us
% dataset_ind = dataset_ind + 1; % 27
% expdate                                     = '20220616';
% expdate_values{dataset_ind}                 = expdate;
% dataset_subtag                              = 1;
% file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
% dataset_stim_polarity(dataset_ind)          = 1; 
% dataset_electrodetype(dataset_ind)          = 3;
% dataset_monopolarrec(dataset_ind)           = 1; 
% 
% % ...                                   , flip polarity
% dataset_ind = dataset_ind + 1; % 28
% expdate                                     = '20220616';
% expdate_values{dataset_ind}                 = expdate;
% dataset_subtag                              = 2;
% file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
% dataset_stim_polarity(dataset_ind)          = 2; 
% dataset_electrodetype(dataset_ind)          = 3;
% dataset_monopolarrec(dataset_ind)           = 1; 


% ***** 20220622
% fast: 20220622, biphasic,   PW = 200 us,          CD = 5, 11.5, 14 mm
% slow: 20220622, biphasic,   PW = 50 us,           CD = 11, 14 mm
% (fast fibers) cVN stim, cVN rec, biphasic, PW = 50 us, CD = 14 mm, orig polarity
dataset_ind = dataset_ind + 1; % 27
expdate                                     = '20220622';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 1;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 1; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 

warning('skipping 20220622 dataset subtags 2-5...'); pause(5);
% ...                                   , flip polarity
dataset_ind = dataset_ind + 1; % 28
expdate                                     = '20220622';
expdate_values{dataset_ind}                 = expdate;
dataset_subtag                              = 6;
file_prefix_outputs_values{dataset_ind}     = [expdate '_' num2str(dataset_subtag)];
dataset_stim_polarity(dataset_ind)          = 2; 
dataset_electrodetype(dataset_ind)          = 3;
dataset_monopolarrec(dataset_ind)           = 1; 
%}

%{
% Make and save metadata table to summarize the datasets
experiment_metadata_table = table(...
    [1:dataset_ind]',expdate_values',file_prefix_outputs_values',dataset_stim_polarity',dataset_electrodetype',dataset_monopolarrec',...
    'VariableNames',...
    {'dataset_ind_all','expdate_values','file_prefix_outputs_values','dataset_stim_polarity','dataset_electrodetype','dataset_monopolarrec'});
save('experiment_metadata_table.mat','experiment_metadata_table');
%}

% Load the list of data to be excluded; these will not be preprocessed, and
% any existing preprocessed output from these datasets will be deleted 
exclusion_list = readtable('excluded_datasets.csv');

% Specify name of file containing the metadata for every valid
% experiemental recording
experiment_metadata_filename = [raw_data_path,'/AnalysisScript_CAPs_Batch_parameters.xlsx'];

% Load metadata
experiment_metadata_table = readtable(experiment_metadata_filename);

expdate_values = arrayfun(@(x) num2str(x), experiment_metadata_table.expdate_values, 'UniformOutput', false);
file_prefix_outputs_values = {};
for i = 1:length(experiment_metadata_table.file_prefix_outputs_values)
    file_prefix_outputs_values{i,1} = [preprocessed_CAP_path, '/', experiment_metadata_table.file_prefix_outputs_values{i}];
end
dataset_stim_polarity = experiment_metadata_table.dataset_stim_polarity;
dataset_electrodetype = experiment_metadata_table.dataset_electrodetype;
dataset_monopolarrec = experiment_metadata_table.dataset_monopolarrec;

dataset_ind = size(experiment_metadata_table,1);

disp(['You have ' num2str(dataset_ind) ' datasets!'])


%% Section 3: Print a table (to be pasted into a Table in the manuscript) summarizing the metadata on a per-nerve basis

%{

number of SN, cVN, tibial recordings?
with histology?
without?

We recorded from 


%}

% Specify the field name, the label for that field in the printed table,
% and how to print the field values assuming each nerve can have multiple
% unique and non-unique field values
%%% The entire summary table will summarize field values for each unique
%%% value of the specified unique_fields; if unique_fields has multiple
%%% elements, then table will be sorted by both (in that order) and
%%% summaries will be based on each unique *combination* of those two
%%% fields
unique_fields = {'dataset_stim_polarity','expdate_values','rec_nerve','target_fiber_types'}; 
unique_fields_print_info = {'%d','%d','%s','%s'};

%%% The unique fields will be printed with the labels specified here
unique_fields_labels = {'Polarity','Nerve ID','Nerve Target','Fiber Type Target'};

field_info = {
    'stim_rec_distance_center_to_center_mm','Conduction Distance (mm)', @(x) fprintf('%0.1f, ',x) % @(x) arrayfun(@(g,h) fprintf('%0.1f (%d)',g,length(x), unique(x), length()) %
    'dataset_electrodetype', 'Electrode Type', @(x) fprintf('%d ',unique(x))
    'flag_histology_available', 'Histology?', @(x) fprintf('%d ',unique(x))
    'gain_value', 'Amplifier Gain', @(x) fprintf('%d ',unique(x))
    'stim_PW', 'Pulse Width (us)', @(x) fprintf('%d ',unique(round(1e3*x)))
    'phase_type', 'Phases', @(x) cellfun(@(g) fprintf('%s, ',g), unique(x))
};

% For each unique combo, report summary of the fields in field_info
%%% Find all the rows with unique combos of unique_fields
[unique_field_combos,~,ic] = unique(experiment_metadata_table(:,...
    cellfun(@(x) find(strcmp(experiment_metadata_table.Properties.VariableNames,x)) ,unique_fields)),'rows');

% Print table headings
for j = 1:length(unique_fields)
    fprintf('%s\t',unique_fields_labels{j});
end

for j = 1:size(field_info,1)
    fprintf('%s\t',field_info{j,2});
end

% Print newline
fprintf('\n');

for i = 1:size(unique_field_combos,1)
    
    % Print the unique combo information in table format
    for j = 1:length(unique_fields)
        field_value = unique_field_combos(i,:).(unique_fields{j});
        if (iscell(field_value))
            assert(length(field_value)==1,'if the field value is a cell array type, it must have only one element');
            field_value = field_value{1};
        end
        fprintf([unique_fields_print_info{j},'\t'],field_value);
    end
    
    % Print the summary data from all rows with unique combo
    table_subset = experiment_metadata_table(ic==i,:);
    for j = 1:size(field_info,1)
        data_ij = table_subset.(field_info{j,1});
        field_info{j,3}(data_ij);
        fprintf('\t');
    end
    
    % Print newline
    fprintf('\n');
end

% Report the summary of # of unique animals cVN myel count (with histo),
% cVN unmyel count (with histo), # of unique animals SN (with histo), # of
% unique animals tibial (with histo)
% experiment_metadata_table.expdate_values
% unique(experiment_metadata_table(:,...
%     cellfun(@(x) find(strcmp(experiment_metadata_table.Properties.VariableNames,x)),{'rec_nerve','flag_histology_available','target_fiber_types','expdate_values'})),'rows');

% Report the total number of rows and columns
fprintf('\n%d rows (including header); %d columns (including unique_fields)\n',...
    size(unique_field_combos,1)+1, length(unique_fields)+size(field_info,1));

%% Analyze datasets
reanalyze_data = 1;

if reanalyze_data == 1
    disp('Consider if you need to go run empty_output_folders.m before re-analyzing data')
    tmp = input(['Are you sure you want to reanalyze ' num2str(dataset_ind) ' datasets? Enter y or n: ' ],'s');

    if ~strcmp(tmp,'y')
        error('I guess you want to set reanalyze_data = 0?')
    end

    warning('this version of code disables the CV bounds check that takes place in AnalysisScript_CAPs, since stim artifact is indeed sometimes overlapping with signal'); pause(2);
    % Reanalyze each dataset
    for dataset_ind = 1:numel(expdate_values)
        disp(['********** DATASET #' num2str(dataset_ind) ' **********'])
        expdate             = expdate_values{dataset_ind};
        file_prefix_outputs = file_prefix_outputs_values{dataset_ind};
        
        % If the dataset is listed for exclusion, delete its preprocessed
        % CAP data or skip preprocessing depending on whether the dataset
        % already exists
        trial_ID = split(file_prefix_outputs,'/');
        trial_ID = trial_ID{end};
        if (ismember(trial_ID,exclusion_list.trial_unique_ID))
            % If the preprocessed data exists for this excluded dataset,
            % delete it
            if exist(file_prefix_outputs,'dir')
                warning('\tDeleting existing preprocessed CAP for this dataset because it is listed for exclusion...')
                pause(5); % pause briefly to give an user brief opportunity to respond
                rmdir(file_prefix_outputs,'s');
            else
                fprintf('\tSkipping dataset because it is listed for exclusion\n\n');
            end
            continue;
        end

        % % If the dataset already exists, skip it
        % if exist(file_prefix_outputs,'dir')
        %     fprintf('\tSkipping dataset because its folder already exists\n\n');
        %     continue;
        % end

        % Preprocess CAPs
        mkdir(file_prefix_outputs)
        ch_num_Vinput           = 4;
        if (dataset_electrodetype(dataset_ind) == 1) || (dataset_monopolarrec(dataset_ind) == 0)
            ch_num_neural_values    = 1;
        elseif (dataset_electrodetype(dataset_ind) == 2) && (dataset_monopolarrec(dataset_ind) == 1)
            ch_num_neural_values    = [1 2];
        elseif (dataset_electrodetype(dataset_ind) == 3) && (dataset_monopolarrec(dataset_ind) == 1)
            ch_num_neural_values    = [1 2 3];
        else
            error('did not account for this possibility')
        end
        AnalysisScript_CAPs_Wrapper(expdate,file_prefix_outputs,CAP_avg_type,scaled_noise,flag_attenuation,...
            subtractnoise,zerobaseline,digitalLPF,digitalLPF_fc_bandpass,...
            digitaloffsetremoval,ch_num_Vinput,ch_num_neural_values,dataset_monopolarrec(dataset_ind),...
            experiment_metadata_table(dataset_ind,:),raw_data_path);
    end
end


%% Compile data: Calculate metrics
%{

Associated text (for reference):
we calculated the 
minimum
peak 
and time of minimum peak

that occurred after the stimulus pulses

for in vivo unmyelinated
we also
that were within the range of 0.4 to 2 m/s
and we verified visually that all the 
negative peaks were correctly within that range 



%}

% ***** For each signal, calculate metrics and store metrics as tables

% Define metrics
metric_functions = {
    @(t,y) find_CAP_amplitude_or_location(t,y,'peak_to_peak_value')
    @(t,y) find_CAP_amplitude_or_location(t,y,'negative_peak_time')
    };
metric_names = {
    'peak_to_peak_value_uV'
    'negative_peak_time_ms'
};

warning('adding title to the gut check figures... disable this eventually...'); pause(2);

% For each signal, calculate each metric
n_metrics = length(metric_functions);
output_metrics = nan(numel(expdate_values),n_metrics);
for dataset_ind = 1:numel(expdate_values)
    fprintf('processing %s\n',file_prefix_outputs_values{dataset_ind})
    
    % Load dataset
    expdate             = expdate_values{dataset_ind};
    file_prefix_outputs = file_prefix_outputs_values{dataset_ind};
    ch_num_neural = 1;

    % If the dataset is listed for exclusion, delete its preprocessed
    % CAP data or skip preprocessing depending on whether the dataset
    % already exists
    trial_ID = split(file_prefix_outputs,'/');
    trial_ID = trial_ID{end};
    if (ismember(trial_ID,exclusion_list.trial_unique_ID))
        fprintf('\tSkipping dataset because it is listed for exclusion\n\n');
        continue;
    end

    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg')
    
    % Calculate metrics and store in a table
    for metric_ind = 1:n_metrics
        % Only incldue valid time ranges (i.e., blank the stim pulse)
        t = t_window;
        y = ch_neural_avg(:,size(ch_neural_avg,2));
        if (strcmp(experiment_metadata_table.phase_type{dataset_ind},'monophasic'))
            stim_pulse_duration = experiment_metadata_table.stim_PW(dataset_ind);
        elseif (strcmp(experiment_metadata_table.phase_type{dataset_ind},'biphasic'))
            stim_pulse_duration = 2*experiment_metadata_table.stim_PW(dataset_ind);
        else
            error('Invalid value of phase type: %s',experiment_metadata_table.phase_type{dataset_ind});
        end
        y(t<stim_pulse_duration) = [];
        t(t<stim_pulse_duration) = [];
        
        % make both input vectors into colum vectors
        if (size(y,2)>1)
            y = y';
        end
        if (size(t,2)>1)
            t = t';
        end
        
        % Also, only include times that fall within the expected fiber
        % conduction velocity range per (Waxman 1978) "Physiology and
        % Pathobiology of Axons": 
        %   0.4-2 m/s (C fibers)
        %   3-120 m/s (A & B fibers) (but I reduced the upper bound to 100
        %   m/s since artifact was within the range for some of the data)
        if (strcmp(experiment_metadata_table.target_fiber_types{dataset_ind},'unmyel'))
            CV_bounds_m_per_s = [0.4 2]; % [m/s]
        elseif (strcmp(experiment_metadata_table.target_fiber_types{dataset_ind},'myel'))
            CV_bounds_m_per_s = [3 95]; % [m/s] (but I reduced the upper bound to 100  m/s since artifact was within the range for some of the data)
        else
            error('invalid value of experiment_metadata_table.target_fiber_types{dataset_ind}: %s',experiment_metadata_table.target_fiber_types{dataset_ind});
        end
        
        indices_within_bounds = find(t>experiment_metadata_table.stim_rec_distance_center_to_center_mm(dataset_ind)/CV_bounds_m_per_s(2) & ...
            t<experiment_metadata_table.stim_rec_distance_center_to_center_mm(dataset_ind)/CV_bounds_m_per_s(1));
        t = t(indices_within_bounds);
        y = y(indices_within_bounds);
        
        output_metrics(dataset_ind,metric_ind) = metric_functions{metric_ind}(t,y);
        
        if (metric_ind==1)
            title(file_prefix_outputs)
        end
    end
end

% Save metrics as a table
metrics_table = [file_prefix_outputs_values, array2table(output_metrics,'VariableNames',metric_names)];
metrics_table.Properties.VariableNames{1} = 'file_prefix_outputs_values';
save([folder_output_metrics,'/metrics_table.mat'],'metrics_table','metric_names');

%% Plot all the in vivo date vs. CV (color-coded per rat, separated by VN (myel), VN (unmyel), SN, and tibial) (all polarity 1)
% sort by polarity, then by fiber type, then by nerve, then by animal
% number
error('This line does not work currently. Need to add the .xlsx file to the data directory; the hardcoded path is obsolete')
plot_signal_comparison('C:\Grill_Lab_Files\Grill_Lab_Analysis_Code_and_Data\simulation_data_files\simulationA0103\look_at_data_variability.xlsx');

%% Compile data: Max CAP, SN stim & SN rec, monophasic, PW = 50 us

% ***** Plot CAPs
dataset_values  = [9 17 3];
ch_num_neural   = 1;

y_lim           = 4500;
lnwdth          = 3;

cmp_dataset      = parula;
Ncmp             = numel(dataset_values);
cmp_dataset      = cmp_dataset(1:floor(size(cmp_dataset,1)/Ncmp):end,:);

figure
hold on
for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
    
    plot(t_window,ch_neural_avg(:,end),'Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'DisplayName',[expdate ', ' num2str(stim_rec_dist) ' mm' ])
end

% Horizontal line at y = 0
lnwdth = 2;
x_tmp = t_window;
y_tmp = zeros(size(x_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Vertical line at x = PW
y_tmp = -y_lim:y_lim;
x_tmp = 0.05 .* ones(size(y_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Labels
xlabel('Time [ms]')
ylabel('Neural rec [\muV]')
axis square
axis([-0.1 1.5 -y_lim y_lim])
legend
title('Max CAP')
set(gca,'FontSize',14)

% Save figure
fname = [folder_output_figures '/MaxCAPs_SNrec_PW50us_OrigPolarity_ch1'];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);

% ***** Add flipped polarity
openfig(fname);

lnwdth = 2;

dataset_values = [9 17 3] + 1;

for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
    
    plot(t_window,ch_neural_avg(:,end),'--','Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'HandleVisibility','off')
end

title('Max CAP + flipped polarity')

% ***** Add other channels
openfig(fname);

dataset_values = [9 17 3];

ch_num_neural   = 2;
lnwdth          = 2;

for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
    
    plot(t_window,ch_neural_avg(:,end),'-','Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'HandleVisibility','off')
end

ch_num_neural   = 3;
lnwdth          = 1;

for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
    
    plot(t_window,ch_neural_avg(:,end),'-','Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'HandleVisibility','off')
end

title('Max CAP across 3 contacts of tripolar cuff')

% Save figure
fname = [folder_output_figures '/MaxCAPs_SNrec_PW50us_OrigPolarity_AllCh'];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);

% ***** Scale to same peak


%% Compile data: Max CAP vs. CV, SN stim & SN rec, monophasic, PW = 50 us
%{
% ***** Plot CAPs
dataset_values  = [9 17 3];
ch_num_neural   = 1;

y_lim           = 4500;
lnwdth          = 3;

cmp_dataset      = parula;
Ncmp             = numel(dataset_values);
cmp_dataset      = cmp_dataset(1:floor(size(cmp_dataset,1)/Ncmp):end,:);

figure
hold on
for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'CV_window','ch_neural_avg','stim_rec_dist')
    
    plot(CV_window,ch_neural_avg(:,end),'Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'DisplayName',[expdate ', ' num2str(stim_rec_dist) ' mm' ])
end

% Horizontal line at y = 0
lnwdth = 2;
x_tmp = CV_window;
y_tmp = zeros(size(x_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Vertical line at x = PW
y_tmp = -y_lim:y_lim;
x_tmp = 0.05 .* ones(size(y_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Labels
xlabel('CV [m/s]')
ylabel('Neural rec [\muV]')
axis square
axis([-0.1 150 -y_lim y_lim])
legend('Location','southeast')
title('Max CAP vs. conduction speed')
set(gca,'FontSize',14)

% Save figure
fname = [folder_outputs '/MaxCAPsVsCV_SNrec_PW50us_OrigPolarity_ch1'];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);
%}

%% Compile data: Max CAP, monophasic, PW = 50 us; compare across stim/rec locations

% ***** Plot CAPs, proximal SN stim, SN rec
dataset_values  = [9 17];
ch_num_neural   = 1;

y_lim           = 4500;
lnwdth          = 3;

cmp_dataset      = parula;
Ncmp             = 3; %numel(dataset_values); % 3 colors to match dataset colors with plots above
cmp_dataset      = cmp_dataset(1:floor(size(cmp_dataset,1)/Ncmp):end,:);

figure
hold on
for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
    
    plot(t_window,ch_neural_avg(:,end),'Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'DisplayName',[expdate ', ' num2str(stim_rec_dist) ' mm' ])
end

% Horizontal line at y = 0
lnwdth = 2;
x_tmp = t_window;
y_tmp = zeros(size(x_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Vertical line at x = PW
y_tmp = -y_lim:y_lim;
x_tmp = 0.05 .* ones(size(y_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Labels
xlabel('Time [ms]')
ylabel('Neural rec [\muV]')
axis square
axis([-0.1 1.5 -y_lim y_lim])
legend('Location','southeast')
title('Max CAP, proximal SN stim, distal SN rec')
set(gca,'FontSize',14)

% Save figure
fname = [folder_output_figures '/MaxCAPs_StimRecDists_Part1'];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);

% ***** Add tibial rec (still with proximal SN stim)
openfig(fname);

lnwdth = 2;

dataset_values = [13 21];

for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
    
    plot(t_window,ch_neural_avg(:,end),'-','Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'DisplayName',[expdate ', ' num2str(stim_rec_dist) ' mm' ])
end

title('Add tibial rec (still proximal SN stim)')

% Save figure
fname = [folder_output_figures '/MaxCAPs_StimRecDists_Part2'];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);

% ***** Add tibial rec with distal SN stim
openfig(fname);

lnwdth = 1;

dataset_values = [15 25];

for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
    
    plot(t_window,ch_neural_avg(:,end),'-','Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'DisplayName',[expdate ', ' num2str(stim_rec_dist) ' mm' ])
end

title('Keep tibial rec, but stim distal SN')

% Save figure
fname = [folder_output_figures '/MaxCAPs_StimRecDists_Part3'];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);


%% Compile data: Max CAP, test-retest
%{
% ***** Plot CAPs, proximal SN stim, SN rec
dataset_values  = [13 21];
ch_num_neural   = 1;

y_lim           = 2000;
lnwdth          = 3;

cmp_dataset      = parula;
Ncmp             = 3; %numel(dataset_values); % 3 colors to match dataset colors with plots above
cmp_dataset      = cmp_dataset(1:floor(size(cmp_dataset,1)/Ncmp):end,:);

figure
hold on
for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
    
    plot(t_window,ch_neural_avg(:,end),'Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'DisplayName',[expdate ', ' num2str(stim_rec_dist) ' mm' ])
end

% Horizontal line at y = 0
lnwdth = 2;
x_tmp = t_window;
y_tmp = zeros(size(x_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Vertical line at x = PW
y_tmp = -y_lim:y_lim;
x_tmp = 0.05 .* ones(size(y_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Labels
xlabel('Time [ms]')
ylabel('Neural rec [\muV]')
axis square
axis([-0.1 1.5 -y_lim y_lim])
legend('Location','southeast')
title('Max CAP, proximal SN stim, distal SN rec')
set(gca,'FontSize',14)

% ***** Add retest (note diff stim-rec distance, but more importantly, cuff
% was removed and replaced)
dataset_values = [23];
cmp_dataset_ind = 2;

ch_num_neural   = 1;
lnwdth          = 2;

for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
    
    plot(t_window,ch_neural_avg(:,end),'--','Color',cmp_dataset(cmp_dataset_ind,:),'LineWidth',lnwdth,'DisplayName',[expdate ', ' num2str(stim_rec_dist) ' mm' ])
end

title('Remove and replace rec cuff')
%}

%% Compile data: Max CAP, proximal SN stim & tibial rec, monophasic, PW = 50 us
%{
% ***** Plot CAPs (proximal stim)
dataset_values  = [13 21];
ch_num_neural   = 1;

y_lim           = 1200;
lnwdth          = 3;

cmp_dataset      = parula;
Ncmp             = 3; %numel(dataset_values); % 3 colors to match dataset colors with plots above
cmp_dataset      = cmp_dataset(1:floor(size(cmp_dataset,1)/Ncmp):end,:);

figure
hold on
for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
    
    plot(t_window,ch_neural_avg(:,end),'Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'DisplayName',[expdate ', ' num2str(stim_rec_dist) ' mm' ])
end

% Horizontal line at y = 0
lnwdth = 2;
x_tmp = t_window;
y_tmp = zeros(size(x_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Vertical line at x = PW
y_tmp = -y_lim:y_lim;
x_tmp = 0.05 .* ones(size(y_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Labels
xlabel('Time [ms]')
ylabel('Neural rec [\muV]')
axis square
axis([-0.1 1.5 -y_lim y_lim])
legend
title('Max CAP')
set(gca,'FontSize',14)

% Save figure
fname = [folder_outputs '/MaxCAPs_tibialrec_proximalstim'];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);

% ***** Add flipped polarity
openfig(fname);

lnwdth = 2;

dataset_values = [13 21] + 1;

for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
    
    plot(t_window,ch_neural_avg(:,end),'--','Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'HandleVisibility','off')
end

title('Max CAP + flipped polarity')

% ***** Add other channels
openfig(fname);

dataset_values = [13 21];

ch_num_neural   = 2;
lnwdth          = 1;

for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
    
    plot(t_window,ch_neural_avg(:,end),'-','Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'HandleVisibility','off')
end

% legend(gca,'off')

title('Max CAP across 2 contacts of bipolar cuff')
%}

%% Compile data: Max CAP, distal SN stim & tibial rec, monophasic, PW = 50 us
%{
% ***** Plot CAPs (proximal stim)
dataset_values  = [15 25];
ch_num_neural   = 1;

y_lim           = 2500;
lnwdth          = 3;

cmp_dataset      = parula;
Ncmp             = 3; %numel(dataset_values); % 3 colors to match dataset colors with plots above
cmp_dataset      = cmp_dataset(1:floor(size(cmp_dataset,1)/Ncmp):end,:);

figure
hold on
for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
    
    plot(t_window,ch_neural_avg(:,end),'Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'DisplayName',[expdate ', ' num2str(stim_rec_dist) ' mm' ])
end

% Horizontal line at y = 0
lnwdth = 2;
x_tmp = t_window;
y_tmp = zeros(size(x_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Vertical line at x = PW
y_tmp = -y_lim:y_lim;
x_tmp = 0.05 .* ones(size(y_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Labels
xlabel('Time [ms]')
ylabel('Neural rec [\muV]')
axis square
axis([-0.1 1.5 -y_lim y_lim])
legend
title('Max CAP')
set(gca,'FontSize',14)

% Save figure
fname = [folder_outputs '/MaxCAPs_tibialrec_distalstim'];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);

% ***** Add flipped polarity
openfig(fname);

lnwdth = 2;

dataset_values = [15 25] + 1;

for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
    
    plot(t_window,ch_neural_avg(:,end),'--','Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'HandleVisibility','off')
end

title('Max CAP + flipped polarity')

% ***** Add other channels
openfig(fname);

dataset_values = [15 25];

ch_num_neural   = 2;
lnwdth          = 1;

for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
    
    plot(t_window,ch_neural_avg(:,end),'-','Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'HandleVisibility','off')
end

% legend(gca,'off')

title('Max CAP across 2 contacts of bipolar cuff')
%}

%% Compile data: Max CAP vs. CV, proximal SN stim & tibial rec, monophasic, PW = 50 us
%{
% ***** Plot CAPs (proximal stim)
dataset_values  = [13 21];
ch_num_neural   = 1;

y_lim           = 1200;
lnwdth          = 3;

cmp_dataset      = parula;
Ncmp             = 3; %numel(dataset_values); % 3 colors to match dataset colors with plots above
cmp_dataset      = cmp_dataset(1:floor(size(cmp_dataset,1)/Ncmp):end,:);

figure
hold on
for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'CV_window','ch_neural_avg','stim_rec_dist')
    
    plot(CV_window,ch_neural_avg(:,end),'Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'DisplayName',[expdate ', ' num2str(stim_rec_dist) ' mm' ])
end

% Horizontal line at y = 0
lnwdth = 2;
x_tmp = CV_window;
y_tmp = zeros(size(x_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Vertical line at x = PW
y_tmp = -y_lim:y_lim;
x_tmp = 0.05 .* ones(size(y_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Labels
xlabel('CV [m/s]')
ylabel('Neural rec [\muV]')
axis square
axis([-0.1 150 -y_lim y_lim])
% legend
title('Max CAP vs. conduction speed')
set(gca,'FontSize',14)

% Save figure
fname = [folder_outputs '/MaxCAPsVsCV_tibialrec_proximalstim'];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);
%}

%% Compile data: Max CAP vs. CV, distal SN stim & tibial rec, monophasic, PW = 50 us
%{
% ***** Plot CAPs (proximal stim)
dataset_values  = [15 25];
ch_num_neural   = 1;

y_lim           = 2500;
lnwdth          = 3;

cmp_dataset      = parula;
Ncmp             = 3; %numel(dataset_values); % 3 colors to match dataset colors with plots above
cmp_dataset      = cmp_dataset(1:floor(size(cmp_dataset,1)/Ncmp):end,:);

figure
hold on
for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'CV_window','ch_neural_avg','stim_rec_dist')
    
    plot(CV_window,ch_neural_avg(:,end),'Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'DisplayName',[expdate ', ' num2str(stim_rec_dist) ' mm' ])
end

% Horizontal line at y = 0
lnwdth = 2;
x_tmp = CV_window;
y_tmp = zeros(size(x_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Vertical line at x = PW
y_tmp = -y_lim:y_lim;
x_tmp = 0.05 .* ones(size(y_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Labels
xlabel('CV [m/s]')
ylabel('Neural rec [\muV]')
axis square
axis([-0.1 150 -y_lim y_lim])
% legend
title('Max CAP vs. conduction speed')
set(gca,'FontSize',14)

% Save figure
fname = [folder_outputs '/MaxCAPsVsCV_tibialrec_distalstim'];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);
%}

%% Compile data: Max CAP, combine monopolar recordings into bipolar/tripolar recordings
%{
% ***** Plot CAPs (tripoles; proximal SN stim, distal SN rec)
dataset_values          = [9 17 3];
ch_num_neural_values    = [1 2 3];

y_lim           = 4500;
lnwdth          = 3;

cmp_dataset      = parula;
Ncmp             = numel(dataset_values);
cmp_dataset      = cmp_dataset(1:floor(size(cmp_dataset,1)/Ncmp):end,:);

figure
hold on
for dataset_ind = 1:numel(dataset_values)
    clear ch_neural_avg_all
    
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    
    for ch_num_neural_ind = 1:numel(ch_num_neural_values)
        ch_num_neural = ch_num_neural_values(ch_num_neural_ind);
        load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
        ch_neural_avg_all(ch_num_neural_ind,:,:) = ch_neural_avg;
    end
    
    % Max CAP
    ch_neural_avg_all = ch_neural_avg_all(:,:,end); % 3 monopoles x time x trials sorted by stim amp
        
    % Calculate tripole
    ch_neural_avg_all = ch_neural_avg_all(2,:) - ((ch_neural_avg_all(1,:) + ch_neural_avg_all(3,:))/2);
    
    plot(t_window,ch_neural_avg_all,'Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'DisplayName',[expdate ', ' num2str(stim_rec_dist) ' mm' ])
end

% Horizontal line at y = 0
lnwdth = 2;
x_tmp = t_window;
y_tmp = zeros(size(x_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Vertical line at x = PW
y_tmp = -y_lim:y_lim;
x_tmp = 0.05 .* ones(size(y_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Labels
xlabel('Time [ms]')
ylabel('Neural rec [\muV]')
axis square
axis([-0.1 1.5 -y_lim y_lim])
legend
title('Max CAP - Tripole')
set(gca,'FontSize',14)

% Save figure
fname = [folder_outputs '/MaxCAPs_SNrec_PW50us_OrigPolarity_tripole'];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);

% ***** Plot CAPs (bipoles; proximal SN stim, tibial rec)
dataset_values          = [13 21];
ch_num_neural_values    = [1 2];

y_lim           = 1200;
lnwdth          = 3;

cmp_dataset      = parula;
Ncmp             = 3; %numel(dataset_values); % 3 colors to match dataset colors with plots above
cmp_dataset      = cmp_dataset(1:floor(size(cmp_dataset,1)/Ncmp):end,:);

figure
hold on
for dataset_ind = 1:numel(dataset_values)
    clear ch_neural_avg_all
    
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    
    for ch_num_neural_ind = 1:numel(ch_num_neural_values)
        ch_num_neural = ch_num_neural_values(ch_num_neural_ind);
        load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
        ch_neural_avg_all(ch_num_neural_ind,:,:) = ch_neural_avg;
    end
    
    % Max CAP
    ch_neural_avg_all = ch_neural_avg_all(:,:,end); % 2 monopoles x time x trials sorted by stim amp
        
    % Calculate bipole
    ch_neural_avg_all = ch_neural_avg_all(1,:) - ch_neural_avg_all(2,:);
    
    plot(t_window,ch_neural_avg_all,'Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'DisplayName',[expdate ', ' num2str(stim_rec_dist) ' mm' ])
end

% Horizontal line at y = 0
lnwdth = 2;
x_tmp = t_window;
y_tmp = zeros(size(x_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Vertical line at x = PW
y_tmp = -y_lim:y_lim;
x_tmp = 0.05 .* ones(size(y_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Labels
xlabel('Time [ms]')
ylabel('Neural rec [\muV]')
axis square
axis([-0.1 1.5 -y_lim y_lim])
legend
title('Max CAP - Bipole')
set(gca,'FontSize',14)

% Save figure
fname = [folder_outputs '/MaxCAPs_tibialrec_PW50us_OrigPolarity_bipole_part1'];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);

% ***** Plot CAPs (bipoles; distal SN stim, tibial rec)
dataset_values          = [15 25];
ch_num_neural_values    = [1 2];

y_lim           = 2500;
lnwdth          = 3;

cmp_dataset      = parula;
Ncmp             = 3; %numel(dataset_values); % 3 colors to match dataset colors with plots above
cmp_dataset      = cmp_dataset(1:floor(size(cmp_dataset,1)/Ncmp):end,:);

figure
hold on
for dataset_ind = 1:numel(dataset_values)
    clear ch_neural_avg_all
    
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    
    for ch_num_neural_ind = 1:numel(ch_num_neural_values)
        ch_num_neural = ch_num_neural_values(ch_num_neural_ind);
        load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','stim_rec_dist')
        ch_neural_avg_all(ch_num_neural_ind,:,:) = ch_neural_avg;
    end
    
    % Max CAP
    ch_neural_avg_all = ch_neural_avg_all(:,:,end); % 2 monopoles x time x trials sorted by stim amp
        
    % Calculate bipole
    ch_neural_avg_all = ch_neural_avg_all(1,:) - ch_neural_avg_all(2,:);
    
    plot(t_window,ch_neural_avg_all,'Color',cmp_dataset(dataset_ind,:),'LineWidth',lnwdth,'DisplayName',[expdate ', ' num2str(stim_rec_dist) ' mm' ])
end

% Horizontal line at y = 0
lnwdth = 2;
x_tmp = t_window;
y_tmp = zeros(size(x_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Vertical line at x = PW
y_tmp = -y_lim:y_lim;
x_tmp = 0.05 .* ones(size(y_tmp));
plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')

% Labels
xlabel('Time [ms]')
ylabel('Neural rec [\muV]')
axis square
axis([-0.1 1.5 -y_lim y_lim])
legend
title('Max CAP - Bipole')
set(gca,'FontSize',14)

% Save figure
fname = [folder_outputs '/MaxCAPs_tibialrec_PW50us_OrigPolarity_bipole_part2'];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);
%}

%% Compile data: Graded CAPs, monophasic, PW = 50 us, orig polarity
dataset_values  = [9 17 3 13 21 15 25];
ch_num_neural   = 1;

clear dataset_legend
dataset_legend{1} = 'prox SN stim, SN rec';
dataset_legend{2} = 'prox SN stim, SN rec';
dataset_legend{3} = 'prox SN stim, SN rec';
dataset_legend{4} = 'prox SN stim, tibial rec';
dataset_legend{5} = 'prox SN stim, tibial rec';
dataset_legend{6} = 'dist SN stim, tibial rec';
dataset_legend{7} = 'dist SN stim, tibial rec';

clear y_lim_values
y_lim_values{1}   = 4500;
y_lim_values{2}   = 1500;
y_lim_values{3}   = 1300;
y_lim_values{4}   = 1200;
y_lim_values{5}   = 800;
y_lim_values{6}   = 2500;
y_lim_values{7}   = 2500;

% CV_bounds = [10 28 50 100];

% CV_bounds = [51.9	100	26.4	8;
%             59.5	100	27.8	8;
%             39.4	100	18.3	8;
%             58.1	100	31.6	8;
%             57.1	100	29.3	8;
%             38.5	100	20	    8;
%             43.5	100	20.4	8];

for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};
    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'t_window','ch_neural_avg','trial_values','amp_values_actual','stim_rec_dist','CV_bounds')
    
    figure
    hold on
    lnwdth = 1;
    y_lim = y_lim_values{dataset_ind};
    
    for trial_ind = 1:numel(trial_values)
        plot(t_window,ch_neural_avg(:,trial_ind),'LineWidth',lnwdth,'DisplayName',[num2str(amp_values_actual(trial_ind)) ' mA' ])
    end

    % Horizontal line at y = 0
    lnwdth = 2;
    x_tmp = t_window;
    y_tmp = zeros(size(x_tmp));
    plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')
    
    % Vertical line at x = PW
    lnwdth = 2;
    y_tmp = -y_lim:y_lim;
    x_tmp = 0.05 .* ones(size(y_tmp));
    plot(x_tmp,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')
    
    % Vertical line at each CV bound
    lnwdth = 1;
    y_tmp = -y_lim:y_lim;
    CV_bounds = unique(CV_bounds(:));
    
    for CV_ind = 1:numel(CV_bounds)
        x_tmp = stim_rec_dist / CV_bounds(CV_ind);
        x_tmp_values = x_tmp .* ones(size(y_tmp));
        plot(x_tmp_values,y_tmp,'Color',[0.5 0.5 0.5],'LineWidth',lnwdth,'HandleVisibility','off')
        text(x_tmp,-0.9*y_lim,num2str(CV_bounds(CV_ind)),'HorizontalAlignment','left','Rotation',90)
    end
        
    % Labels
    xlabel('Time [ms]')
    ylabel('Neural rec [\muV]')
    axis square
    axis([-0.1 2 -y_lim y_lim])
    legend('Location','eastoutside')
    title(['Graded CAP, ' expdate ', ' dataset_legend{dataset_ind}])
    set(gca,'FontSize',14)
    
    % Save figure
    trial_ID = split(file_prefix_outputs,'/');
    trial_ID = trial_ID{end};
    fname = [folder_output_figures '/GradedCAPs_PW50us_OrigPolarity_ch1_' trial_ID];
    saveas(gcf,[fname '.jpg']);
    saveas(gcf,[fname '.fig']);

end

%% Compile data: Recruitment curves, monophasic, PW = 50 us, orig polarity
dataset_values  = [9 17 3 13 21 15 25];
ch_num_neural   = 1;

clear dataset_legend
dataset_legend{1} = 'prox SN stim, SN rec';
dataset_legend{2} = 'prox SN stim, SN rec';
dataset_legend{3} = 'prox SN stim, SN rec';
dataset_legend{4} = 'prox SN stim, tibial rec';
dataset_legend{5} = 'prox SN stim, tibial rec';
dataset_legend{6} = 'dist SN stim, tibial rec';
dataset_legend{7} = 'dist SN stim, tibial rec';

title_strings{1}  = 'All CVs';
title_strings{2}  = 'First lobe';
title_strings{3}  = 'Second lobe';
title_strings{4}  = 'Third lobe';

% y_lim           = 2500;
lnwdth          = 3;

cmp_dataset      = parula;
Ncmp             = numel(dataset_values);
cmp_dataset      = cmp_dataset(1:floor(size(cmp_dataset,1)/Ncmp):end,:);

x_tmp           = 0:0.01:1; % [mA]

FIG_NUM_START = 100;

for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};

    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'amp_values_actual','AUC','ft_AUC','CV_bounds')
    for CVrange_ind = 1:size(CV_bounds,1)
        figure(FIG_NUM_START + CVrange_ind - 1)
%         set(gcf,'position',[scrnsz(3)*0.25 scrnsz(4)*0.2 scrnsz(3)*0.5 scrnsz(4)*0.6])
        hold on
        
        % Original data points
        y_tmp = AUC(CVrange_ind,:);
        plot(amp_values_actual,y_tmp./max(y_tmp),'.','MarkerSize',10,...
            'Color',cmp_dataset(dataset_ind,:),...
            'LineWidth',lnwdth,'HandleVisibility','off')
        
        % Sigmoid fit
        y_tmp = fitfunc_recruitmentcurve(x_tmp,ft_AUC{CVrange_ind}.a,...
                ft_AUC{CVrange_ind}.b,ft_AUC{CVrange_ind}.c);
        plot(x_tmp,y_tmp./max(y_tmp),'Color',cmp_dataset(dataset_ind,:),...
            'LineWidth',lnwdth,'DisplayName',[expdate ', ' dataset_legend{dataset_ind}])
    end
end

for CVrange_ind = 1:size(CV_bounds,1)
    figure(FIG_NUM_START + CVrange_ind - 1)
    
    % Labels
    xlabel('Stim amp [mA]')
    ylabel('Normalized AUC')
    axis square
    % axis([-0.1 1.5 -y_lim y_lim])
    ylim([-0.1 1.1])
%     legend('Location','eastoutside')
    title(['Dose-response curves - ' title_strings{CVrange_ind}])
    set(gca,'FontSize',14)

    % Save figure
    fname = [folder_output_figures '/DRC_CVrange' num2str(CVrange_ind)];
    saveas(gcf,[fname '.jpg']);
    saveas(gcf,[fname '.fig']);
end

%% Compile data: 50% thresholds vs. CV (ish...), monophasic, PW = 50 us, orig polarity
dataset_values  = [9 17 3 13 21 15 25];
ch_num_neural   = 1;

clear dataset_legend
dataset_legend{1} = 'prox SN stim, SN rec';
dataset_legend{2} = 'prox SN stim, SN rec';
dataset_legend{3} = 'prox SN stim, SN rec';
dataset_legend{4} = 'prox SN stim, tibial rec';
dataset_legend{5} = 'prox SN stim, tibial rec';
dataset_legend{6} = 'dist SN stim, tibial rec';
dataset_legend{7} = 'dist SN stim, tibial rec';

lnwdth          = 3;

cmp_dataset      = parula;
Ncmp             = numel(dataset_values);
cmp_dataset      = cmp_dataset(1:floor(size(cmp_dataset,1)/Ncmp):end,:);

thresh_target = 0.5;

clear thresh_values

figure
hold on

for dataset_ind = 1:numel(dataset_values)
    dataset_val         = dataset_values(dataset_ind);
    expdate             = expdate_values{dataset_val};
    file_prefix_outputs = file_prefix_outputs_values{dataset_val};

    load([file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'],'ft_AUC','CV_bounds')
    for CVrange_ind = 1:size(CV_bounds,1)
        y_tmp = ft_AUC{CVrange_ind}.a * thresh_target;
        thresh_values(dataset_ind,CVrange_ind) = fitfunc_recruitmentcurve_inverse(y_tmp,ft_AUC{CVrange_ind}.a,...
                                                                    ft_AUC{CVrange_ind}.b,ft_AUC{CVrange_ind}.c);
    end
    x_tmp = 1:size(CV_bounds,1);
    plot(x_tmp,thresh_values(dataset_ind,:),'.-','Color',cmp_dataset(dataset_ind,:),...
         'MarkerSize',20,'DisplayName',[expdate ', ' dataset_legend{dataset_ind}])
    
end

xticks(x_tmp)
xticklabels({'Full CV range','First lobe','Second lobe','Third lobe'})

% Labels
ylabel('50% thresh (mA)')
title('Threshold for different CV ranges')
axis square
% legend('Location','eastoutside')
set(gca,'FontSize',14)

%% Presumably obsolete, commented out code related to KHF; to be deleted eventually
%{
%% Compile data: KHF, normalized denoised AUC, sigmoid fit, thresh
% Exclude sigmoids as indicated by include_sigmoids = 0

for dataset_ind = 1:numel(expdate_values)
    expdate             = expdate_values{dataset_ind};
    file_prefix_outputs = file_prefix_outputs_values{dataset_ind};
    load([file_prefix_outputs '/' expdate '.mat'])
    
    Nfreq               = numel(freq_values);
    incl_idx            = include_sigmoids{dataset_ind};
    
    if (Nfreq ~= size(incl_idx,2))
        error('Inclusion criteria does not match number of sigmoids')
    end
   
    fibertype_values = dataset_fibertypes{dataset_ind};
    for fibertype_ind = 1:numel(fibertype_values)
        fibertype   = fibertype_values(fibertype_ind);
        thresh_tmp  = [];

        alldata_freqs{dataset_ind}{fibertype}      = freq_values(incl_idx(fibertype_ind,:) == 1);
        
        freq_count = 0;
        for freq_ind = find(incl_idx(fibertype_ind,:) == 1)
            freq_count = freq_count + 1;
            
            alldata_AUC{dataset_ind}{fibertype}{freq_count}       = AUC_all{freq_ind}{fibertype};
            alldata_sigmoidft{dataset_ind}{fibertype}{freq_count} = ft_all{freq_ind}{fibertype};
            thresh_tmp(end+1)                                     = thresh_all{freq_ind}{fibertype};
        end
        alldata_thresh{dataset_ind}{fibertype}     = thresh_tmp;
    end
end

%% Plot sigmoids
Nfreq_colors        = 80;
cmap                = parula(Nfreq_colors);

x_values_stimamp    = 0:0.01:20;    % [mA]
lnstyles{1}         = '-.';         % Monopolar (dataset_electrodetype == 1)
lnstyles{2}         = '-';          % Bipolar (dataset_electrodetype == 2)
lnstyles{3}         = '--';         % Tripolar (dataset_electrodetype == 3)
lnwdth              = 2;

y_thresh            = 0.5;

figure
set(gcf,'Position',[scrnsz(1) scrnsz(4)*0.3 scrnsz(3) scrnsz(4)*0.5])
subplot_rows = 1;
subplot_cols = 3;

for dataset_ind = 1:numel(expdate_values)

    fibertype_values = dataset_fibertypes{dataset_ind};
    for fibertype_ind = 1:numel(fibertype_values)
        fibertype       = fibertype_values(fibertype_ind);
        
        freqs_tmp       = alldata_freqs{dataset_ind}{fibertype};
        sigmoidft_tmp   = alldata_sigmoidft{dataset_ind}{fibertype};   
        thresh_tmp      = alldata_thresh{dataset_ind}{fibertype};

        lnstyle_tmp     = lnstyles{dataset_electrodetype(dataset_ind)};

        if max(freqs_tmp) > Nfreq_colors
            error('Need to expand freq cmap')
        end
        
        % VNS
        if dataset_type(dataset_ind) == 1
            if fibertype == 1         % B-fibers
                subplot_idx = 2;
            elseif fibertype == 2     % C-fibers
                subplot_idx = 3;
            end

        % SNS with force rec
        elseif dataset_type(dataset_ind) == 2
            subplot_idx = 1;
        end

        subplot(subplot_rows,subplot_cols,subplot_idx); hold on

        for freq_ind = 1:numel(freqs_tmp)
            freq_val = freqs_tmp(freq_ind);

            ft_tmp = sigmoidft_tmp{freq_ind};
            plot(x_values_stimamp,KHFAUC_fitfunc(x_values_stimamp,ft_tmp.b,ft_tmp.c),...
                'LineStyle',lnstyle_tmp,'LineWidth',lnwdth,'Color',cmap(round(freq_val),:),'HandleVisibility','off')

            th_tmp = thresh_tmp(freq_ind);
            plot(th_tmp,y_thresh,'.','MarkerSize',20,'Color',cmap(round(freq_val),:))
        end
    end
end

for subplot_idx = 1:3
    subplot(subplot_rows,subplot_cols,subplot_idx)
    set(gca,'FontSize',14)
    xlabel('KHF amp (mA)')
    ylabel('Normalized denoised AUC')
    axis square
    
    c = colorbar;
    c_min = 1;
    c_max = Nfreq_colors;
    tick_labels = 10:10:Nfreq_colors;
    c.Ticks = (tick_labels - c_min) / (c_max - c_min);%linspace(0,1,Nfreq_colors/10);
    c.TickLabels = tick_labels;
    set(get(c,'title'),'string','kHz');
end

subplot_idx = 1;
subplot(subplot_rows,subplot_cols,subplot_idx)
title('SNS, Force rec')
xmax = 15;
xlim([0 xmax])

% Manual legend for linestyles corresponding to # contacts of KHF electrode
xtext_factor = 0.6;
ytext = 0.95; ytext_delta = 0.05;
text(xmax*xtext_factor,ytext,'-.-  Monopolar','FontSize',14)
ytext = ytext - ytext_delta;
text(xmax*xtext_factor,ytext,'—  Bipolar','FontSize',14)
ytext = ytext - ytext_delta;
text(xmax*xtext_factor,ytext,'- -  Tripolar','FontSize',14)

subplot_idx = 2;
subplot(subplot_rows,subplot_cols,subplot_idx)
title('VNS, B-fibers')
xmax = 15;
xlim([0 xmax])

subplot_idx = 3;
subplot(subplot_rows,subplot_cols,subplot_idx)
title('VNS, C-fibers')
xmax = 15;
xlim([0 xmax])


%% Plot block thresh vs. freq
electrodetype_markers       = '.^s';        % Monopolar, bipolar, tripolar
electrodetype_marker_sizes  = [10 5 8];
fibertype_colors            = 'kbr';        % VNS, B-fibers; VNS, C-fibers; SNS, force rec

figure; hold on
set(gcf,'Position',[scrnsz(3)*0.2 scrnsz(4)*0.3 scrnsz(3)*0.5 scrnsz(4)*0.5])

for dataset_ind = 1:numel(expdate_values)
    
    fibertype_values = dataset_fibertypes{dataset_ind};
    for fibertype_ind = 1:numel(fibertype_values)
        fibertype       = fibertype_values(fibertype_ind);
        
        freqs_tmp       = alldata_freqs{dataset_ind}{fibertype};
        thresh_tmp      = alldata_thresh{dataset_ind}{fibertype};
        
        % VNS, B-fibers
        if (dataset_type(dataset_ind) == 1) && (fibertype == 1)
            color_idx = 1;
            
            % VNS, C-fibers
        elseif (dataset_type(dataset_ind) == 1) && (fibertype == 2)
            color_idx = 2;
            
            % SNS, force rec
        elseif (dataset_type(dataset_ind) == 2) && (fibertype == 1)
            color_idx = 3;
            
        else
            error('No color to account for this dataset type and fiber type in block vs. freq plot')
        end
        color_tmp = fibertype_colors(color_idx);
        
        plot(freqs_tmp,thresh_tmp,'LineStyle','none',...
            'Marker',electrodetype_markers(dataset_electrodetype(dataset_ind)),...
            'MarkerSize',electrodetype_marker_sizes(dataset_electrodetype(dataset_ind)),...
            'Color',color_tmp)
    end
end

axis square
set(gca,'FontSize',20)
xlabel('Freq (kHz)')
ylabel('50% block thresh (mA)')
xlim([0 90])

% Manual legend
xtext = 100;
ytext = 10; ytext_delta = 0.5;

color_idx = 3;
text(xtext,ytext,'SNS, force rec','FontSize',14,'Color',fibertype_colors(color_idx))
ytext = ytext - ytext_delta;

color_idx = 1;
text(xtext,ytext,'VNS, B-fibers','FontSize',14,'Color',fibertype_colors(color_idx))
ytext = ytext - ytext_delta;

color_idx = 2;
text(xtext,ytext,'VNS, C-fibers','FontSize',14,'Color',fibertype_colors(color_idx))
ytext = ytext - ytext_delta*2;

text(xtext,ytext,'\bullet   Monopolar','FontSize',14)
ytext = ytext - ytext_delta;

text_tmp = texlabel('Delta');
text(xtext,ytext,[text_tmp '  Bipolar'],'FontSize',14)
ytext = ytext - ytext_delta;


%}