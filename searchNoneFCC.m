clear; close all hidden;

tic
strain = 10; % to modify
C = fopen('noneFCC.dat','w'); % to modify
DUMP_0 = load('C:\Users\Tao Yuhao\Documents\MATLAB\MD SIM\AnalysisORT(2022.10.28)\GRAD\Dump\dump0.dat');
l = length(DUMP_0);
toc
disp('Initialization completed.')
disp('')

noneFCC = NaN(l,1);
counter_defect = 0; % counter of defect atoms in the whole history
for frame = 1:strain
    tic
    DUMP = load(['C:\Users\Tao Yuhao\Documents\MATLAB\MD SIM\AnalysisORT(2022.10.28)\GRAD\Dump\dump' num2str(frame + 5) '0000.dat']);
    UZFL = [DUMP(:,1), DUMP(:,6)]; % UZFL(ID, structure type)
    UZFL = sortrows(UZFL, 1);

    if(length(UZFL) ~= l) 
        errordlg('Length deviation between DUMP files!'); % check length
    end
    for j = 1:l
        if (UZFL(j,2) ~= 1) && isnan(noneFCC(j))
            noneFCC(j) = 1;
            counter_defect = counter_defect + 1;
        end
    end
    toc
    disp(['Defects in the frame with ' num2str(frame) '% tensile strain are considered'])
    disp('')
end

counter_FCC = 0;
for i = 1:l
    if isnan(noneFCC(i))
        noneFCC(i) = 0;
        counter_FCC = counter_FCC + 1;
    end
end

if (counter_FCC + counter_defect ~= l)
    disp('Length not match between FCC and defects') % check defect quantity
end

tic
for i = 1:l
     if i == 1
         fprintf(C, '%d', noneFCC(i));
     else
         fprintf(C, '\n%d', noneFCC(i));
     end
end
fclose(C);
toc
disp('Printer is done.')
disp('')