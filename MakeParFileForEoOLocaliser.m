
function MakeParFileForEoOLocaliser

%Generate a .par file describing the timing in the Localiser for the Eye of Origin experiment.
%This is in tab-delimited text format, as expected by mrVista when setting up GLM parameters.
%This includes the labelling of all fixation/baseline blocks as 'Fix'
%Note that this file does not include the dummy TRs at the beginning: they must be removed manually in the mrInit process.
%Only 1 file should be needed, since the localiser (should) be the same for all subjects.
%R Maloney 22 June 2015

%The file name to be saved:
parFileName = 'EoO_localiser_design_vista.par';

%Most of below comes directly from Run_EoO_localiser where the timing is determined
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
Fix = 0; %NOTE: In MrVista conventions, 'fixation' periods are always labelled as 'Fix' and given code 0

%Work out a single cycle:
BlockOrder = [Fix; LeftEye_Inner; Fix; LeftEye_Outer; Fix; RightEye_Inner; Fix; ...
    RightEye_Outer; Fix; Binocular_Inner; Fix; Binocular_Outer];

%Replicate for 3 cycles, and place one final Fixation block at the end:
BlockOrder = [BlockOrder; BlockOrder; BlockOrder; Fix];

%Determine the onset times for each BLOCK. These onset times ignore the 5 dummy scans at the beginning
%Ie the first is set as time 0: this is the same as with the event-related scans as determined using Optseq, where the dummy scans are also inserted manually (ie not by Optseq)
%The dummy scans must be manually removed in mrVista at the start of the analysis.
%We want the onset times for each block only, so subtract the length of the dummy scans, and one block, so the final value is the ONSET time of the FINAL block
BlockOnsets = 0:BlockLengthSec:348-DummyTRs*TR-BlockLengthSec; %(348 is total length, in sec of scan, including the dummies).

%Now make a cell array containing the condition names as strings.
%This seems optional, but we will include it because the mrVista files seem to have it:
for ii = 1: length(BlockOrder)
    
    switch BlockOrder(ii)
        case Fix
            BlockName{ii} = 'Fix'
        case LeftEye_Inner
            BlockName{ii} = 'LeftEye_Inner'
        case LeftEye_Outer
            BlockName{ii} = 'LeftEye_Outer'
        case RightEye_Inner
            BlockName{ii} = 'RightEye_Inner'
        case RightEye_Outer
            BlockName{ii} = 'RightEye_Outer'
        case Binocular_Inner
            BlockName{ii} = 'Binocular_Inner'
        case Binocular_Outer
            BlockName{ii} = 'Binocular_Outer'
    end
end


%Now open and save the appropriate text file:

fileID = fopen(parFileName,'wt'); %open the file to be saved, 'w' for writing, 't' for text

%Set up the formatting for the file, this is the order of each item in each row of the file.
% 3.2f means fixed point notation, so up to 3 integer parts and 2 decimal places (for onset time)
% d means a signed integer  (for block code)
% s means a string          (for block/condition name)
% \t means separate by a horizontal tab 
% \n means go to a new line.
formatSpec = '%3.2f\t %d\t %s\t\n'; 

%loop across each line and add it to the file:
for ii = 1:length(blockOrder)
    
    fprintf(fileID, formatSpec, BlockOnsets(ii), blockOrder(ii), BlockName{ii});
    
end

fclose(fileID); %close the file when done. It should be saved in the pwd


    






