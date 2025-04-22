function column_function_axon_diameter_dog_2024(idx_x,Delta,delta)
addpath(genpath('/home/gustavo/Downloads/misst_emmanuelle'))
addpath(genpath('/home/gustavo/Downloads/misst_emmanuelle/matlab_scripts'))
mat_files_path = '/home/gustavo/dogSpinalCord_125um_matfiles';
% 
% addpath(genpath('/home/users/gchau/misst_emmanuelle'))
% addpath(genpath('/home/users/gchau/misst_emmanuelle/matlab_scripts'))
% mat_files_path = '/home/users/gchau/dogSpinalCord_125um_matfiles';

HEIGHT = 41;

%%  -----------  protocol ------------------------------------------------
% waveform ---------------------------------------------
param.delta = Delta; %0.015; % duration in s
param.b = [3000,6000] ;%100:100:8000; % s/mm2
param.smalldel = delta; %0.01;  % duration of the waveform in s (of the gradient)
param.tau = 1e-4 ; % sampling interval
param.Npoints = param.smalldel / param.tau ; % number of points
param.p180 = param.delta-param.smalldel; % time interval between the first and second waveforms (duration of rf pulse)
param = gradient_amplitude(param); % maximum gradient strength in T/m depending on the b value
param.sequence = 'PGSE_no_ramping';

% ---- output file

out = '../result_sim_diffusion/';
mkdir(out);

sim_number = 1;
model_location = ['cylinder_prueba/sim_' num2str(sim_number) '/'];

% output folder
out_folder = [out model_location];
test_directory(out_folder);

% initialization
sim = 1;

disp('generate tensor orientations');
disp(['sim ' num2str(sim) ' starts'])
%number_directions = 40;
%dir = uniform_sphere(number_directions,200);
%dir = [0 0 0 ; dir];%load('tensor_40.mat');
load('/home/users/gchau/axon_code/dir_bruker.mat','dir')
% load('/home/gustavo/Gdrive/Stanford/Lab/ODF_prediction/sherlock_bash/dir_bruker.mat','dir')

param.dir = dir;
% ------ design of the waveform
disp('create waveform') 
[protocol_init,protocolGEN, M] = protocol_design(param); % M is how many diffusion volumes you will obtain
save([out_folder 'protocol_init.mat' ],'protocol_init'); % for corresponding b values, directions 


%% tissue model
disp('generate tissue model and simulate the signal')
model.name = 'Cylinder';
% signal_stack = zeros(M,numel(tissue_model.theta),numel(tissue_model.phi));

%di = 1.7E-9; % intrinsic diffusivity % TODO: check the boundaries of diff for water in brain in white matter
di = 0.35E-9; % parallel diffusivity
theta = 0; %15*pi/180;
phi = 0; %15*pi/180;

all_signal = zeros(HEIGHT,M);
for idx_y=1:HEIGHT
    
    disp(['voxel ' num2str(idx_y) ' of ' num2str(HEIGHT)])
    
    stats_file = [mat_files_path '/tile_' num2str(idx_x) '_' num2str(idx_y) '_statistics.mat'];
    if ~isfile(stats_file)
        all_signal(idx_y,:) = nan;
        disp('tile not found')
        continue
    end
    
    tic
    % get list of axon diameters
    load(stats_file,'valid_statistics','avf','ecf')
    axon_diameters_um = valid_statistics(:,2);
    axon_areas = pi*(axon_diameters_um.^2);
    total_area = sum(axon_areas);
    num_axons = length(axon_diameters_um);
    disp('simulating')
    % simulate signal from cylinders
    signal_cylinders = 0;
    for idx_ax=1:num_axons
        rad = axon_diameters_um(idx_ax)*1E-9/2;
        model.params = [di rad phi theta];
        protocolGEN = MMConstants(model,protocolGEN);
        signal_cylinders = signal_cylinders + (axon_areas(idx_ax)/total_area)*SynthMeas(model,protocolGEN);
    end
    
    model_ec.name = 'Zeppelin';
    d_ort = di*(1-avf);
    model_ec.params = [di d_ort 0 0];
    iso_signal = SynthMeas(model_ec,protocolGEN);
    
    
    signal = avf*signal_cylinders + ecf*iso_signal;
    all_signal(idx_y,:) = signal;
    
    
    
    toc
     
end

save(['diffusion_column' num2str(idx_x) '_Delta_' num2str(param.delta) '_delta_' num2str(param.smalldel) '.mat'],'all_signal','model'); % parameters
        
end



