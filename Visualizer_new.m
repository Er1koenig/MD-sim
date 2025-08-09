clear; close all hidden;

tic
load('INFO_new.mat');
disp('Data is loaded.')
toc

tic
strain = 10; % to modify
C = fopen('ROTVerifier_10%_new.dat','w'); % to modify
DUMP = load('C:\Users\Tao Yuhao\Documents\MATLAB\MD SIM\AnalysisORT(2022.10.28)\GRAD\Dump\dump150000.dat'); % to modify

atoms = 0; % quantity of atoms before cropping the twinning areas
for i = 1:length(INFO{strain}.Grain)
    if ~isempty(INFO{strain}.Grain{i}.atoms)
        atoms = atoms + length(INFO{strain}.Grain{i}.atoms);
    end
end

counter = 0;
Verifier = NaN(atoms, 2);
for i = 1:length(INFO{strain}.Grain)
    if ~isempty(INFO{strain}.Grain{i}.atoms)
        for j = 1:length(INFO{strain}.Grain{i}.atoms)
            if (INFO{strain}.Grain{i}.atoms(j,2) >= INFO{strain}.rot_grain(i,3)) && (INFO{strain}.Grain{i}.atoms(j,2) <= INFO{strain}.rot_grain(i,4))
                counter = counter + 1;
                Verifier(counter,1) = INFO{strain}.Grain{i}.atoms(j,4); % atom ID
                Verifier(counter,2) = INFO{strain}.Grain{i}.atoms(j,2); % atom rotation
            end
        end
    end
end
Verifier(counter + 1:end, :) = [];
Verifier = sortrows(Verifier,1);
toc

tic
DUMP(:, 6:end) = [];
DUMP(:, 2) = [];
DUMP = sortrows(DUMP, 1);

for i = 1:length(Verifier)
    Verifier(i,3) = DUMP(Verifier(i,1), 2);
    Verifier(i,4) = DUMP(Verifier(i,1), 3);
    Verifier(i,5) = DUMP(Verifier(i,1), 4);
end
toc

tic
for i = 1:length(Verifier)
     if i == 1
         fprintf(C,'%d %f %f %f %f',Verifier(i,1), Verifier(i,2), Verifier(i,3), Verifier(i,4), Verifier(i,5));
     else
         fprintf(C,'\n%d %f %f %f %f',Verifier(i,1), Verifier(i,2), Verifier(i,3), Verifier(i,4), Verifier(i,5));
     end
end
fclose(C);
toc

%{

ITEM: TIMESTEP
150000
ITEM: NUMBER OF ATOMS
16474723
ITEM: BOX BOUNDS pp pp pp
-1.4939350827549626e+02 3.0293935082752982e+03
2.4119276072746334e+02 5.9941902392723223e+03
5.4419394590024872e-02 1.4934830605408148e+01
ITEM: ATOMS id rot x y z

%}
