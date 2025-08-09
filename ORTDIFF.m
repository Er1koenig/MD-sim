clear; close all hidden;

%% Initialization

tic
MAP = load('C:\Users\Tao Yuhao\Documents\MATLAB\MD SIM\AnalysisORT(2022.10.28)\GRAD\output_mod.dat'); % MAP(ID, grain size, grain ID)
MAP = sortrows(MAP,1);
toc

%% Compute the orientation of each grain at 0% tensile strain

tic
% load files
DUMP_0 = load('C:\Users\Tao Yuhao\Documents\MATLAB\MD SIM\AnalysisORT(2022.10.28)\GRAD\Dump\dump50000.dat');
ID_0 = DUMP_0(:,1);
STR_0 = DUMP_0(:,6);

QORT_0 = load('C:\Users\Tao Yuhao\Documents\MATLAB\MD SIM\AnalysisORT(2022.10.28)\GRAD\ORT\ort000.dat'); % orientation vector in quaternion form
if(length(DUMP_0) ~= length(QORT_0))
    errordlg('Length difference between DUMP_0 and QORT_0!'); % check length
end

% convert quaternion to angle
ORT_0 = NaN(length(QORT_0), 1); % orientation for z-axis in degree
alpha = NaN;
beta = NaN;
gamma = NaN;
counter_weird = zeros(max(MAP(:,3)), 1);
for i = 1: length(QORT_0)
    alpha = atan2(2 * (QORT_0(i, 4) * QORT_0(i, 1) + QORT_0(i, 2) * QORT_0(i, 3)), 1 - 2 * (QORT_0(i, 1) ^ 2 + QORT_0(i, 2) ^ 2)) * 180 / pi;
    beta = asin(2 * (QORT_0(i, 4) * QORT_0(i, 2) - QORT_0(i, 3) * QORT_0(i, 1))) * 180 / pi;
    gamma = atan2(2 * (QORT_0(i, 4) * QORT_0(i, 3) + QORT_0(i, 2) * QORT_0(i, 1)), 1 - 2 * (QORT_0(i, 2) ^ 2 + QORT_0(i, 3) ^ 2)) * 180 / pi;
    if (abs(beta) > 35) && (abs(beta) < 55) && (abs(alpha) > -10) && (abs(alpha) < 10) && (gamma > 40) && (gamma < 50)
        ORT_0(i) = gamma - 90;
    else
        if (abs(alpha) > 35) && (abs(alpha) < 55) && (abs(beta) < 10)
            ORT_0(i) = gamma;
        elseif (abs(beta) > 35) && (abs(beta) < 55) && (abs(alpha) < 10)
            if (gamma + 90 > -45) && (gamma + 90 < 135)
                ORT_0(i) = gamma + 90;
            elseif (gamma - 90 > -45) && (gamma - 90 < 135)
                ORT_0(i) = gamma - 90;
            end
        else
            ORT_0(i) = -100;
        end
    end
end

if(length(DUMP_0) ~= length(QORT_0))
    errordlg('Length difference between QORT_0 and ORT_0!'); % check length
end
for i = 1: length(QORT_0)
    if isnan(ORT_0(i))
        errordlg('NaN found in ORT_0!'); % check NaN
    end
end

SORTER_0 = [ID_0, ORT_0, STR_0];
SORTER_0 = sortrows(SORTER_0, 1);

ORTG_0 = [MAP(:,3), SORTER_0(:,2), SORTER_0(:,3)]; % ORTG_0(grain ID, lattice orientation, structure type)

disp('Rotation angle are processed for every atom.')
toc
disp(' ')

tic
ort_init = zeros(max(MAP(:,3)), 1);
counter = zeros(max(MAP(:,3)), 1);
for i = 1:length(ORTG_0(:,1))
    if ORTG_0(i,2) == -100
        counter_weird(ORTG_0(i,1)) = counter_weird(ORTG_0(i,1)) + 1;
    else
        if ORTG_0(i,3) == 1
            counter(ORTG_0(i,1)) = counter(ORTG_0(i,1)) + 1;
            ort_init(ORTG_0(i,1)) = ort_init(ORTG_0(i,1)) + ORTG_0(i,2);
        end
    end
end
ort_init = ort_init ./ counter;

disp('The initial frame is ready.')
toc
disp(' ')