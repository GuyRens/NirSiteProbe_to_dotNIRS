%% Step 1: loading in the data
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IMPORTANT COMMENTS

%the wavelengths are hardcoded so they need to be changed for our specific
%ones. I don't think this matter too much for sensitivity modelling as it
%already contains many assumptions such as the propogation vector (?) =>
%how easy the wavelength goes to certain tissue.

% but to be correct, we might wanna change them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load in probe. 
folder_probe = uigetdir;
reg_probe = true; %this registers the probe to the Colin Mesh. %Not sure if this is necessary but I leave it in for now
probe = nirs.io.loadNIRxProbe(folder_probe,reg_probe);


%% Step 2: Ensuring that the link data (channels) and types (wavelengths) are correct
%When loading in a probe using loadNIRxProbe, it does not contain
%wavelengths. To run the simulate date, these wavelengths have to be
%included. Wavelengths are stored in probe.types. However probe.types is
%dependent on probe.link which contains the source detector pairs and their
% wavelengths. Accordingly, the probe.link variable has to adjusted.

%to see how a probe.link variable looks like run:
%exampleprobe = nirs.testing.simData;

%CHECK IF NIRSTORM CAN RUN MORE THAN TWO WAVELENGTHS AND OTHERS THAN THESE
%ONES
wls = [690;830]; %changes to the correct wavelengths for our hardware
wls_amount = length(wls);

%the steps I use to make the probe.link variable are kind of tedious due to the
%reshape function. There might be a better way in doing so. But with these
%steps I am sure that the generated probe.link variable is in the exact
%same format as nirstoolbox would show them...
%If there is a better way then I call this coding style a guyism
temp_link = table2array(probe.link); % arrays > tables
temp_link = repmat(temp_link,[1,wls_amount]); %repeat the matrix for amount of wavelengths
temp_link = temp_link'; %flip before reshaping
temp_link = reshape(temp_link,[2,numel(temp_link)/2]); %reshape
temp_link = temp_link'; % flip again to have col_one sources and col_two dets

%add the wavelengths to the array
temp_link(:,3) = repmat(wls,[length(temp_link)/wls_amount,1]); 

%back to table and providing expected header names
temp_link = array2table(temp_link); 
temp_link.Properties.VariableNames = {'source','detector','type'}; %this naming has to be exact

probe.link = temp_link; % attach the newly made link variable to probe again.



%% Step3: Generate random data for the probe (montage)
% nirstorm does not allow to import montages without data. So we create
% fake/random data for our montage

%simulate noise for the random data
noise = nirs.testing.simARNoise(probe);

% hard coded values are arbitrary. They should not matter for the forward
% model. Hard coded values and the steps below are from the the nirstoolbox
% examples within the simData function.
stim  = nirs.testing.randStimDesign(noise.time, 2,7,3);
beta  = [3 2 1]';

sd = unique([noise.probe.link.source noise.probe.link.detector], 'rows');
channels = sd(1:round(end/2),:);

[data, truth] = nirs.testing.simData(noise, stim, beta, channels );

%This is redundant but has been included to generate non-empty auxillary data
%when writing as .snirf 
data.auxillary.keys = {'ones'};
data.auxillary.values = {ones(length(data.time),1)};

%% Step4: Final check whether the montage is correct

data.probe = data.probe.SetFiducialsVisibility(false); %hide fiducials (on by default - gets in the way)
data.probe.defaultdrawfcn='3D Mesh';
data.probe.draw %finally draw to check whether things are ok

%% Step5: Save the data.
%Nirstorm can read .snirf files, however it seems there are some
%differences between the .snirf format of nirstoolbox and nirstorm. It
%seems that saving to .nirs (Homer format) works perfectly fine. So we safe
%as that.

dir_out = [folder_probe filesep 'format_dotnirs' filesep];

if ~exist(dir_out,'dir')
    mkdir(dir_out);
end

%we save back in the original folder but create a subfolder
nirs.io.saveDotNirs(data,dir_out)

%this would be to save as snirf but there seem to be some formatting issues
%between nirstorm and analyzIR
%nirs.io.saveSNIRF(data,'C:\Users\guy_r\Documents\nirstorm_montage_tester\test_registered')

%% step6: Add additional files
% nirstorm requires a .nirs file but also an optodes.txt and fiducials.txt
% file. These need to be included in the folder to succesfully import data
% into nirstorm and have the probe registered to the scalp.

%--------------------------------------------
%generating the fiducials.txt file.
%The fiducials are 'hardcoded'. These values are taken from nirstorm and
%represent the fiducials from the Colin template that is used for forward
%modeling.
%To make a long story short: NIRx uses ICBM and nirstorm uses Colin. The
%fiducial coordinates are a bit different and the montage starts floating
%in nirstorm when using the ICBM fiducials.

%some of the info in the fiducials file might be redundant, but it works
%fine this way
fids_headers = {'#' 'Sample Name'	'Session Name'	'Index'	'Loc. X'	'Loc. Y'	'Loc. Z'	'Offset'};

fids_nasion = {'Nasion'	'Session' '1'	'2'	'0.75'	'81.22'	'-42.64' 	'0.000'};
fids_LE = {'LeftEar'	'Session' '1'	'3'	'-79.37'	'-27.52'	'-48.07'	'0.000'};
fids_RE = {'RightEar'	'Session' '1'	'4'	'80.64'	'-25.90'	'-46.60'	'0.000'};

% in case when the ICBM template might be needed, I expect you can use the
% coordinates below (these are the original ones for ICBM in nirsite)
% I think atlas viewer uses these ones
% Nz 0.400 85.900 -47.600 
% RPA 83.900 -16.600 -56.700 
% LPA -83.800 -18.600 -57.200 

fiducials = [fids_headers;fids_nasion;fids_LE;fids_RE];

%because the fiducials.txt uses both white space and tabs as delimiters, we
%need to write the output in a speficic way
fids_headers_delims = '%s %s\t%s\t%s\t%s\t%s\t%s\t%s\n'; %only for the headers
fids_delims = '%s\t%s %s\t%s\t%s\t%s\t%s\t%s\n'; %for the landmarks

%guyism for transparancy and because it works well this way
fileID = fopen([dir_out 'fiducials.txt'],'w');
fprintf(fileID,fids_headers_delims,string(fids_headers));
fprintf(fileID,fids_delims,string(fids_nasion));
fprintf(fileID,fids_delims,string(fids_LE));
fprintf(fileID,fids_delims,string(fids_RE));
fclose(fileID);

%----------------------------------------

%--------------------------------------------
%generating the optodes.txt file.

if ~exist([folder_probe filesep 'Standard_Optodes.txt'],'file')
%if standard optodes do no exist extract optode locs from
%standard_probeinfo file

    NS_probe = load([folder_probe filesep 'Standard_probeInfo.mat']); % NS = nirsiteprobe
    NS_probe = NS_probe.probeInfo.probes;

    NS_optode_locs = [NS_probe.coords_s3;NS_probe.coords_d3];
    NS_optode_locs = NS_optode_locs*10; %for some reason it is in cm and has to to be changed to mm

    num_s = (1:nirsite_probe.nSource0)';
    num_d = (1:nirsite_probe.nDetector0)';

    sources = [repmat('s',[length(num_s),1]),char(string(num_s))];
    dets = [repmat('d',[length(num_d),1]),char(string(num_d))];

    optodes = [sources;dets];
    optodes = [cellstr(optodes),num2cell(NS_optode_locs)];
    
else %if standard file optode file exist extract from the text file
    
    fid = fopen([folder_probe filesep 'Standard_Optodes.txt']);
    optodes = textscan(fid,'%s','delimiter',',');
    fclose(fid);
    %yes these are the famous guyisms of coding again
    optodes = optodes{1,1};
    optodes = reshape(optodes,[4,length(optodes)/4]);
    optodes = optodes';

    optodes = optodes(contains(optodes(:,1),{'S', 'D'}),:); 
    %some anatomical points are included which start with letter O, we yeet
    %them 
    
end


%getting the right format
optodes = [optodes(:,1),cellstr(repmat('placeholder',[size(optodes,1),1])),num2cell(zeros([size(optodes,1),1])),optodes(:,2:end),num2cell(zeros([size(optodes,1),1]))];

optode_headers = {'#' 'Sample Name'	'Session Name'	'Index'	'Loc. X'	'Loc. Y'	'Loc. Z'	'Offset'};
optode_header_delims = '%s %s\t%s\t%s\t%s\t%s\t%s\t%s\n'; %only for the headers

optode_delims = '%s\t%s\t%s\t%s\t%s\t%s\t%s\n'; %only for the optodes


fileID = fopen([dir_out 'optodes.txt'],'w');
printtxt = fprintf(fileID,optode_header_delims,string(optode_headers));


for L_opt = 1:size(optodes,1) %saving the optodes to the textfile
    %needs to be looped because this amount can vary
    
    printtxt = fprintf(fileID,optode_delims,string(optodes(L_opt,:)));
    
end
fclose(fileID);
%--------------------------------------------


%% step7: generate the file that can be used for atlasviewer
% Atlasviewer requires a digpts.txt file but and the 1.nirs file
%. These need to be included in the folder to succesfully import data
% into atlasviewer and have the probe registered to the scalp.

%--------------------------------------------
%generating the digpts.txt file
%The fiducials are 'hardcoded'. These values are taken from nirstorm and
%represent the fiducials from the Colin template that is used for forward
%modeling.

% in case when the ICBM template might be needed, I expect you can use the
% coordinates below (these are the original ones for ICBM in nirsite)
% I think atlas viewer uses these ones
% RPA: 83.900 -16.600 -56.700 
% Nz: 0.400 85.900 -47.600 
% Cz: -0.461 -8.416 101.365 
% LPA: -83.800 -18.600 -57.200 
% Iz: 0.200 -120.500 -25.800 

%formatting the data in a correct manner
AV_NZ = {'NZ:' '0.400' '85.900' '-47.600'};
AV_LP = {'LPA:' '-83.800' '-18.600' '-57.200'};
AV_RP = {'RPA:' '83.900' '-16.600' '-56.700'};
AP_CZ = {'CZ:' '-0.461' '-8.416' '101.365'};
AP_IZ = {'IZ:' '0.200' '-120.500' '-25.800'} ;

AV_fiducials = [AV_NZ;AV_LP;AV_RP;AP_CZ;AP_IZ];
AV_optodes = optodes(:,[1,4,5,6]);

AV_optodes_name = char(AV_optodes(:,1));
AV_optodes_name = cellstr([AV_optodes_name,repmat(':',size(AV_optodes_name,1),1)]);
AV_optodes(:,1) = AV_optodes_name;

AV_all = [AV_fiducials;AV_optodes];

%prep for saving
optode_delims = '%s\t%s\t%s\t%s\n'; %only for the optodes

fileID = fopen([dir_out 'digpts.txt'],'w');

for L_opt = 1:size(AV_all,1) %saving the optodes to the textfile
    %needs to be looped because this amount can vary
    
    printtxt = fprintf(fileID,optode_delims,string(AV_all(L_opt,:)));
    
end
fclose(fileID);


