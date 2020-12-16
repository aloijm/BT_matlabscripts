function [mparams, v1, v2, rpe1, rpe2] = noveltyRLFit()

% cd('C:\Bruno\Papers\In preparation\Blair novelty\NoveltyData_9.25');
currentFolder = pwd
cd(currentFolder);

% data = xlsread('subj_Novelty6.xlsm', 8);
data = xlsread('subj_Novelty6.xlsm', 8);


options = optimset('MaxFunEvals', 100000, 'TolFun', 0.001);

minll = Inf;

for ipi = 1 : 10
    
    %%% LR
    iparams(1) = rand(1,1);
    %%% Beta
    iparams(2) = 4 + sqrt(2)*rand(1,1);
    %%% Novelty bonuses for 2
    %iparams(3) = 1 + sqrt(0.5)*rand(1,1);
        
    [mparams, lla] = fminsearch(@(iparams) fitNoveltyRL(iparams, data), iparams, options);
    
    if lla < minll
        minll = lla;
        minParams = mparams;
    end
    
end

mparams = minParams;   

[ll, d1, d2, v1, v2, rpe1, rpe2] = f_fitNoveltyRL(mparams, data);

quit

% plot(v1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ll] = fitNoveltyRL(params, data)

[ll] = f_fitNoveltyRL(params, data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ll, d1, d2, v1, v2, rpe1, rpe2] = f_fitNoveltyRL(params, data)

lr    = params(1);
beta  = params(2);
%bonus1 = params(3);
%bonus2 = params(3);
%bonus3 = params(4);

ll = 0;
for runi = 1 : 2
    
    if runi == 1
        choiceData = data(:, 1:6);
    else
        choiceData = data(:, 7:12);
    end
    
    nTrials = size(choiceData, 1);
    
    v   = zeros(nTrials+1, 3);
    %v(1,1:3) = 0.15;
    d   = zeros(nTrials+1, 3);
    imageIndex = choiceData(1, 1:3);
    imageIndices = zeros(nTrials+1, 3);
    
    
    vchosen   = zeros(nTrials+1,1);
    rpechosen = zeros(nTrials+1, 1);
    dchosen   = zeros(nTrials+1, 1);
    NovelNumbers = zeros(nTrials+1,1);
    NovelChosen = zeros(nTrials+1,1);
    novelImageNumber = 999;
    
    %Create a string array for ssid (for exporting output later)
    ssid = string("subj");
    ssidtable = strings(nTrials+1,1);
    ssidtable(:,:) = ssid;
    ssidtable2 = [ssidtable;ssidtable];
    
    for triali = 1 : nTrials
        
        %% Check if cells are empty and set values to 0 and quit if they are empty
        if isnan(choiceData(triali, 1:4))
        	v1 = zeros(nTrials+1,1);
        	d1 = zeros(nTrials+1,1);
        	v2 = zeros(nTrials+1,1);
        	d2 = zeros(nTrials+1,1);
        	lr = 0;
        	beta = 0;
        	%bonus2 = 0;
        	%bonus3 = 0;
        	
        	%a = [lr beta bonus2];
        	a = [lr beta];
			T = table(a);
			writetable(T, 'results_temp.txt')
			
			run1and2 = [v1 d1; v2 d2];
			U = table(ssidtable2, run1and2);
			writetable(U, 'decisiondata_temp.txt')
        	quit
        end
                         
        %%% Need to assign new option when it is introduced and some value
        %%% for novel option.  Will make separate parameter for now.
        if choiceData(triali, 6) == 1
            
            novelImageIndex = ~ismember(imageIndex, choiceData(triali, 1:3));
            
            novelImageLoc = ~ismember(choiceData(triali, 1:3), imageIndex);
            
            novelImageNumber = choiceData(triali, novelImageLoc);
            
            NovelNumbers(triali, 1) = novelImageNumber;
            
            z = find(novelImageIndex == 1);
            
            imageIndex(z) = novelImageNumber;
            
            %v(triali, z) = bonus1;
            
        elseif choiceData(triali, 6) == 0
        
        	NovelNumbers(triali, 1) = novelImageNumber;
            
        end 
        
        %%% Bonus shows up on second trial after option introduced
       % if triali > 1 && choiceData(triali-1, 6) == 1

        %    v(triali, z) = bonus2;
            
        %end
        
        %%% Bonus shows up on third trial after option introduced
        %if triali > 2 && choiceData(triali-2, 6) == 1

         %   v(triali, z) = bonus3;
            
        %end
        
        d(triali, :) = exp(beta*v(triali, :))/sum(exp(beta*v(triali, :)));
        
        if choiceData(triali, 4) == 0            
            v(triali+1, :) = v(triali, :);
            rpechosen(triali, 1) = -Inf;
            continue;
        end
        
        
        choiceNumber = choiceData(triali, choiceData(triali, 4));
        
        choiceIndex = find(imageIndex == choiceNumber);
        
        vchosen(triali, 1)   = v(triali, choiceIndex);
        rpechosen(triali, 1) = choiceData(triali, 5) - v(triali, choiceIndex);
        dchosen(triali, 1) = d(triali, choiceIndex);
        imageIndices(triali, :) = imageIndex(1, :);
        
        if choiceNumber == novelImageNumber
            NovelChosen(triali) = 1;
        end
            
        v(triali+1, :) = v(triali, :);
        
        v(triali+1, choiceIndex) = v(triali, choiceIndex) + lr*(choiceData(triali, 5) - v(triali, choiceIndex));
        
                
                
        ll = ll -log(d(triali, choiceIndex));
        
%         [d(triali, :) v(triali, :)]
        
        fprintf('');
        
    end
    
    if runi == 1
        v1   = vchosen;
        vrun1 = v;
        rpe1 = rpechosen;
        d1chosen   = dchosen;
        d1 = d;
        imageIndices1 = imageIndices;
        choiceData1 = [choiceData; 0 0 0 0 0 0];
        NovelNumbers1 = NovelNumbers;
        NovelChosen1 = NovelChosen;
       % Run1 = ones(141,1);
    else
        v2   = vchosen;
        vrun2 = v;
        rpe2 = rpechosen;
        d2chosen   = dchosen;
        d2 = d;
        imageIndices2 = imageIndices;
        choiceData2 = [choiceData; 0 0 0 0 0 0];
        NovelNumbers2 = NovelNumbers;
        NovelChosen2 = NovelChosen;
       % Run2 = 2*ones(141,1);
    end
    
end

%prints out ssid learning rate, inverse temp, bonus 2 and bonus 3 to a text file in the file
%a = [lr beta bonus2];
a = [lr beta];
T = table(a);
writetable(T, 'results_temp.txt')

run1and2 = [v1 vrun1 d1 imageIndices1 d1chosen choiceData1 NovelNumbers1 NovelChosen1; v2 vrun2 d2 imageIndices2 d2chosen choiceData2 NovelNumbers2 NovelChosen1];
U = table(ssidtable2, run1and2);
writetable(U, 'decisiondata_temp.txt')
