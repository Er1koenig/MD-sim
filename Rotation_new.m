clear; close all hidden;

%% Initialization

tic
MAP = load('C:\Users\Tao Yuhao\Documents\MATLAB\MD SIM\AnalysisORT(2022.10.28)\GRAD\output_mod.dat');
MAP = sortrows(MAP,1);
meas = 10; % to modify

counter = zeros(meas,7);
ROT_REL = zeros(meas,7);
ROT = zeros(meas, 7);
INFO = cell(meas, 1); % data output

noneFCC = load('C:\Users\Tao Yuhao\Documents\MATLAB\MD SIM\AnalysisORT(2022.10.28)\GRAD\noneFCC.dat');
toc

%% Check Out the Frames during Tensile Deformation

for frame = 1:meas
    % load files
    tic
    DUMP_X = load(['C:\Users\Tao Yuhao\Documents\MATLAB\MD SIM\AnalysisORT(2022.10.28)\GRAD\Dump\dump' num2str(frame + 5) '0000.dat']);
    ID_X = DUMP_X(:,1);
    STR_X = DUMP_X(:,6);
    if frame < 10
        QROT_X = load(['C:\Users\Tao Yuhao\Documents\MATLAB\MD SIM\Analysis(2022.05.29)\Rotation\ROT\r0' num2str(frame) '0.dat']);
    else
        QROT_X = load(['C:\Users\Tao Yuhao\Documents\MATLAB\MD SIM\Analysis(2022.05.29)\Rotation\ROT\r' num2str(frame) '0.dat']);
    end

    if(length(DUMP_X) ~= length(QROT_X))
        errordlg('Length difference between DUMP and QROT!'); % check length
    end

    disp(['The frame with applied tensile strain of ' num2str(frame) '% is loaded.'])
    toc
    disp(' ')

    % convert quaternion to angle
    tic
    ROT_X = NaN(length(QROT_X), 1); % orientation for z-axis in degree
    for i = 1: length(QROT_X)
        ROT_X(i) = atan2(2 * (QROT_X(i,4) * QROT_X(i,3) + QROT_X(i,2) * QROT_X(i,1)), 1 - 2 * (QROT_X(i,2) ^ 2 + QROT_X(i,3) ^ 2)) * 180 / pi;
    end
    if(length(DUMP_X) ~= length(QROT_X))
        errordlg('Length difference between QROT_X and ROT_X!'); % check length
    end
    for i = 1: length(QROT_X)
        if isnan(ROT_X(i))
            errordlg('NaN found in ROT_X!'); % check NaN
        end
    end

    SORTER_X = [ID_X, ROT_X, STR_X];
    SORTER_X = sortrows(SORTER_X, 1);
    OUT_X = [MAP(:,2), SORTER_X(:,2), SORTER_X(:,3), MAP(:,3)]; % OUT(grain size, rotation, structure type, grain ID)

    disp('Quaternions are converted to angles.')
    toc
    disp(' ')

    % verify the rotation angles of every atom: only consider fcc atoms, process the rotation angles till rot ∈ [-45, 45]
    tic
    VER = [OUT_X, (1:length(OUT_X))']; % VER(grain size, rotation, structure type, grain ID, atom ID)
    if length(noneFCC) ~= length(VER)
        disp('Lengths of VER and noneFCC are not the same.')
    end
    for i = 1:length(VER)
        if (noneFCC(i) == 0) % only consider fcc atoms
            if (VER(i,2) >= -45) && (VER(i,2) <= 45)
                continue
            elseif (VER(i,2) >= -90) && (VER(i,2) < -45)
                VER(i,2) = VER(i,2) + 90;
            elseif (VER(i,2) > 45) && (VER(i,2) <= 90)
                VER(i,2) = VER(i,2) - 90;
            else
                VER(i,2) = NaN;
            end
        else
            VER(i,2) = NaN;
        end
    end
    VER(isnan(VER(:,2)),:) = [];

    disp('Rotation angle are processed for every atom.')
    toc
    disp(' ')

    % detect twinning areas and compute absolute rotation for every grain
    tic
    FCC = VER;
    FCC(:,3) = []; % FCC(grain size, rotation, grain ID, atom ID)

    Grain = cell(max(FCC(:,3)), 1); % cell for features of every grain, including atoms, rotation angle distribution and valid atoms
    counter_grain = zeros(max(FCC(:,3)), 1);

    for i = 1:length(FCC) % initialize atom list
        counter_grain(FCC(i,3)) = counter_grain(FCC(i,3)) + 1;
    end
    for i = 1:length(Grain)
        Grain{i}.atoms = NaN(counter_grain(i), length(FCC(1,:))); % atoms(grain size, rotation, grain ID, atom ID)
    end

    for i = 1:length(Grain) % initialize atom distribution for every grain
        Grain{i}.DTB = zeros(18,4); % DTB(rot_min, rot_max, atom quantity, sum of rotation)
        for j = 1:9
            Grain{i}.DTB(2 * j - 1, 1) = 5 * (j - 1);
            Grain{i}.DTB(2 * j - 1, 2) = 5 * j;

            Grain{i}.DTB(2 * j, 1) = -5 * j;
            Grain{i}.DTB(2 * j, 2) = -5 * (j - 1);
        end
    end

    for i = 1:length(Grain) % initialize valid atom list
        Grain{i}.vatoms = NaN(counter_grain(i), length(FCC(1,:))); % vatoms(grain size, rotation, grain ID, atom ID)
    end

    rot_grain = NaN(length(Grain), 4); % initialize final rotation without twinning areas: rot_grain(valid atom percentage, valid rotation angle, rot_min, rot_max)

    counter_grain = zeros(max(FCC(:,3)), 1);
    for i = 1:length(FCC) % build atom list for every grain
        counter_grain(FCC(i,3)) = counter_grain(FCC(i,3)) + 1;
        Grain{FCC(i,3)}.atoms(counter_grain(FCC(i,3)), :) = FCC(i,:);
    end

    for i = 1:length(Grain)
        if ~isempty(Grain{i}.atoms)
            for j = 1:counter_grain(i) % build list of rotation angle distribution for every grain
                for k = 1:length(Grain{i}.DTB(:,1))
                    if (Grain{i}.atoms(j,2) >= Grain{i}.DTB(k,1)) && (Grain{i}.atoms(j,2) < Grain{i}.DTB(k,2))
                        Grain{i}.DTB(k,3) = Grain{i}.DTB(k,3) + 1;
                        Grain{i}.DTB(k,4) = Grain{i}.DTB(k,4) + Grain{i}.atoms(j,2);
                        break
                    end
                end
            end
            Grain{i}.DTB = sortrows(Grain{i}.DTB, 3, 'descend');
            if (max([Grain{i}.DTB(1,1), Grain{i}.DTB(2,1), Grain{i}.DTB(3,1)]) - min([Grain{i}.DTB(1,1), Grain{i}.DTB(2,1), Grain{i}.DTB(3,1)]) > 9.99)... % top 3 distribution are neighbors
                    && (max([Grain{i}.DTB(1,1), Grain{i}.DTB(2,1), Grain{i}.DTB(3,1)]) - min([Grain{i}.DTB(1,1), Grain{i}.DTB(2,1), Grain{i}.DTB(3,1)]) < 10.01)
                rot_grain(i,1) = (Grain{i}.DTB(1,3) + Grain{i}.DTB(2,3) + Grain{i}.DTB(3,3)) / counter_grain(i);
                rot_grain(i,2) = (Grain{i}.DTB(1,4) + Grain{i}.DTB(2,4) + Grain{i}.DTB(3,4)) / (Grain{i}.DTB(1,3) + Grain{i}.DTB(2,3) + Grain{i}.DTB(3,3));
                rot_grain(i,3) = min([Grain{i}.DTB(1,1), Grain{i}.DTB(2,1), Grain{i}.DTB(3,1)]);
                rot_grain(i,4) = max([Grain{i}.DTB(1,2), Grain{i}.DTB(2,2), Grain{i}.DTB(3,2)]);
            else
                quant = 0;
                left = NaN;
                right = NaN;
                quant = quant + Grain{i}.DTB(1,3);
                for j = 2:length(Grain{i}.DTB(:,1))
                    if (Grain{i}.DTB(1,1) - Grain{i}.DTB(j,1) < 5.01) && (Grain{i}.DTB(1,1) - Grain{i}.DTB(j,1) > 4.99)
                        left = j;
                        quant = quant + Grain{i}.DTB(j,3);
                    end
                    if (Grain{i}.DTB(j,1) - Grain{i}.DTB(1,1) < 5.01) && (Grain{i}.DTB(j,1) - Grain{i}.DTB(1,1) > 4.99)
                        right = j;
                        quant = quant + Grain{i}.DTB(j,3);
                    end
                end
                if isnan(left)||isnan(right)
                    errordlg('One of the neighbor is NaN.') % verify neighbors
                end
                if quant / counter_grain(i) > 0.6 % the top 1 distribution stands for the majority
                    rot_grain(i,1) = quant / counter_grain(i);
                    rot_grain(i,2) = (Grain{i}.DTB(1,4) + Grain{i}.DTB(left,4) + Grain{i}.DTB(right,4)) / (Grain{i}.DTB(1,3) + Grain{i}.DTB(left,3) + Grain{i}.DTB(right,3));
                    rot_grain(i,3) = min([Grain{i}.DTB(1,1), Grain{i}.DTB(left,1), Grain{i}.DTB(right,1)]);
                    rot_grain(i,4) = max([Grain{i}.DTB(1,2), Grain{i}.DTB(left,2), Grain{i}.DTB(right,2)]);
                else
                    if (abs(Grain{i}.DTB(2,1) - Grain{i}.DTB(3,1)) < 5.1) && ((Grain{i}.DTB(2,3) + Grain{i}.DTB(3,3)) / counter_grain(i) > 0.55)
                        rot_grain(i,1) = (Grain{i}.DTB(2,3) + Grain{i}.DTB(3,3)) / counter_grain(i);
                        rot_grain(i,2) = (Grain{i}.DTB(2,4) + Grain{i}.DTB(3,4)) / (Grain{i}.DTB(2,3) + Grain{i}.DTB(3,3));
                        rot_grain(i,3) = min([Grain{i}.DTB(2,1), Grain{i}.DTB(3,1)]);
                        rot_grain(i,4) = max([Grain{i}.DTB(2,2), Grain{i}.DTB(3,2)]);
                    else
                        if (abs(Grain{i}.DTB(1,1) - Grain{i}.DTB(2,1)) < 5.1)
                            rot_grain(i,1) = (Grain{i}.DTB(1,3) + Grain{i}.DTB(2,3)) / counter_grain(i);
                            rot_grain(i,2) = (Grain{i}.DTB(1,4) + Grain{i}.DTB(2,4)) / (Grain{i}.DTB(1,3) + Grain{i}.DTB(2,3));
                            rot_grain(i,3) = min([Grain{i}.DTB(1,1), Grain{i}.DTB(2,1)]);
                            rot_grain(i,4) = max([Grain{i}.DTB(1,2), Grain{i}.DTB(2,2)]);
                        elseif (abs(Grain{i}.DTB(2,1) - Grain{i}.DTB(3,1)) < 5.1)
                            rot_grain(i,1) = (Grain{i}.DTB(2,3) + Grain{i}.DTB(3,3)) / counter_grain(i);
                            rot_grain(i,2) = (Grain{i}.DTB(2,4) + Grain{i}.DTB(3,4)) / (Grain{i}.DTB(2,3) + Grain{i}.DTB(3,3));
                            rot_grain(i,3) = min([Grain{i}.DTB(2,1), Grain{i}.DTB(3,1)]);
                            rot_grain(i,4) = max([Grain{i}.DTB(2,2), Grain{i}.DTB(3,2)]);
                        elseif (abs(Grain{i}.DTB(1,1) - Grain{i}.DTB(3,1)) < 5.1)
                            rot_grain(i,1) = (Grain{i}.DTB(1,3) + Grain{i}.DTB(3,3)) / counter_grain(i);
                            rot_grain(i,2) = (Grain{i}.DTB(1,4) + Grain{i}.DTB(3,4)) / (Grain{i}.DTB(1,3) + Grain{i}.DTB(3,3));
                            rot_grain(i,3) = min([Grain{i}.DTB(1,1), Grain{i}.DTB(3,1)]);
                            rot_grain(i,4) = max([Grain{i}.DTB(1,2), Grain{i}.DTB(3,2)]);
                        else
                            rot_grain(i,1) = quant / counter_grain(i);
                            rot_grain(i,2) = (Grain{i}.DTB(1,4) + Grain{i}.DTB(left,4) + Grain{i}.DTB(right,4)) / (Grain{i}.DTB(1,3) + Grain{i}.DTB(left,3) + Grain{i}.DTB(right,3));
                            rot_grain(i,3) = min([Grain{i}.DTB(1,1), Grain{i}.DTB(left,1), Grain{i}.DTB(right,1)]);
                            rot_grain(i,4) = max([Grain{i}.DTB(1,2), Grain{i}.DTB(left,2), Grain{i}.DTB(right,2)]);
                            errordlg(['Probelm! The major distribution only account for ' num2str(100 * quant / counter_grain(i)) '% of the atoms in the grain ' num2str(i) '.'])
                        end
                    end
                end
            end
        end
    end

    disp('Twinning areas are trimmed, rotation angles are calculated.')
    toc
    disp(' ')

    % compute the relative rotation angle of every grain
    load('Neighborlist.mat');
    rel_rot_grain = NaN(length(rot_grain), 1); % initialize relative rotation vector: rel_rot_grain(grain ID)
    for i = 1:length(rot_grain)
        if ~isnan(rot_grain(i,2))
            adder = 0;
            for j = 1:length(Neighb{i}.grain)
                adder = adder + rot_grain(i,2) - rot_grain(Neighb{i}.grain(j), 2);
            end
            rel_rot_grain(i) = adder / length(Neighb{i}.grain);
        end
    end

    % compute average rotation intensity every grain size area (gradient structure)
    tic
    for i = 1:length(Grain)
        if ~isempty(Grain{i}.atoms) && ~isnan(rel_rot_grain(i))
            counter(frame, Grain{i}.atoms(1,1)) = counter(frame, Grain{i}.atoms(1,1)) + 1;
            ROT_REL(frame, Grain{i}.atoms(1,1)) = ROT_REL(frame, Grain{i}.atoms(1,1)) + abs(rel_rot_grain(i));
            ROT(frame, Grain{i}.atoms(1,1)) = ROT(frame, Grain{i}.atoms(1,1)) + abs(rot_grain(i,2));
        end
    end
    ROT_REL(frame,:) = ROT_REL(frame,:) ./ counter(frame,:);
    ROT(frame,:) = ROT(frame,:) ./ counter(frame,:);

    disp(['The rotation angle of tensile frame ' num2str(frame) '% is ready.'])
    disp(ROT_REL(frame, :))
    disp(ROT(frame, :))
    toc
    disp(' ')

    % save data into output file
    INFO{frame}.Grain = Grain;
    INFO{frame}.rot_grain = rot_grain;
end

save('INFO_new.mat', 'INFO')