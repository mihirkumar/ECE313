close all;
clear all;

%%Task 0
%Load in data sets and make it readable.
dat1_const = load('patient_data/1_a41178.mat');
dat2_const = load('patient_data/2_a42126.mat');
dat3_const = load('patient_data/3_a40076.mat');
dat4_const = load('patient_data/4_a40050.mat');
dat5_const = load('patient_data/5_a41287.mat');
dat6_const = load('patient_data/6_a41846.mat');
dat7_const = load('patient_data/7_a41846.mat');
dat8_const = load('patient_data/8_a42008.mat');
dat9_const = load('patient_data/9_a41846.mat');

%put in array for simplification of data management
dat_array = [dat1_const, dat2_const, dat3_const, dat4_const, dat5_const, dat6_const, dat7_const, dat8_const, dat9_const];

%X_length: stores number of data points per patient for data, training,
%test
all_length = zeros(1, 9);
training_length = zeros(1,9);
testing_length = zeros(1, 9);
for k = 1:9
    %store the number of elements to perform operations with
    all_length(k) = length(dat_array(k).all_data);
    training_length(k) = floor(all_length(k) * 2/3);
    testing_length(k) = all_length(k) - training_length(k);
    dat_array(k).all_data = floor(dat_array(k).all_data);
   
end

%initialize all data sets for permuatations
for k = 1:9
    train(k).all_data = dat_array(k).all_data(:,1:training_length(k));
    train(k).all_labels = dat_array(k).all_labels(:,1:training_length(k));
    test(k).all_data = dat_array(k).all_data(:, training_length(k)+1:all_length(k));
    test(k).all_labels = dat_array(k).all_labels(:, training_length(k)+1:all_length(k));
end
for k = 1:9
   for j = 1:7
      index1 = 1;
      index0 = 1;
      for p = 1: training_length(k)
          if(train(k).all_labels(1, p) == 1)
              train(k).goldens(j, index1) = train(k).all_data(j,p);
              index1 = index1 + 1;
          end
          if(train(k).all_labels(1, p) == 0)
              train(k).nongoldens(j,index0) = train(k).all_data(j,p);
              index0 = index0 + 1;
          end
      end
   end
end
%open a result file for printing things
fid = fopen('ECE313_Final_group5', 'w');
%completed: 5:55 AM
%update: 8:08 AM
%%Task 1
%a: Calculate prior probabilities of P(H1) an dP(H0)
for k = 1:9
    prior_H1(k) = sum(train(k).all_labels)/training_length(k);
    prior_H0(k) = 1 - prior_H1(k);
end
%b: construct likelihood matrices for each of the seven features

feature_labels = {'Mean Area under the Heart Beat','Mean R-to-R peak interval','Number of beats per minute (Heart Rate)','Peak to peak interval for Blood Pressure','Systolic Blood Pressure','Diastolic Blood Pressure','Pulse Pressure'};
feature_max_length = [17, 220, 220, 220, 115, 86, 74];

for k = 1:9
    figure;
    for j = 1:7
        train(k).H1{j,1} = tabulate(train(k).goldens(j,:))';
        train(k).H1{j,1}(3,:) = train(k).H1{j,1}(3,:) / 100;
        
        train(k).H0{j,1} = tabulate(train(k).nongoldens(j,:))';
        train(k).H0{j,1}(3,:) = train(k).H0{j,1}(3,:) / 100;        
        
        subplot(7, 1, j); 
        plot(train(k).H0{j,1}(3,:));
        
        % Add feature titles and set axis for each subplot
        title(feature_labels(j));
        axis([0 feature_max_length(j) 0 1]);
        
        hold on; 
        plot(train(k).H1{j,1}(3,:));
    end
    legend('H0 pmf', 'H1 pmf');
end

%c: show results by generating a seperate figure for each patient 
%c Executed in B
for k = 1:9
    for j = 1:7
        %make H1, H0 same length
        delta = length(train(k).H1{j,1}(3,:)) - length(train(k).H0{j,1}(3,:));
        if(delta < 0)
            old_start = length(train(k).H1{j,1}(3,:));
            last_val = train(k).H1{j,1}(1,old_start);
            first_val = train(k).H1{j,1}(1,1);
            if(train(k).H1{j,1}(1,1) < train(k).H0{j,1}(1,1))
                train(k).H1{j,1} = horzcat(train(k).H1{j,1}(:,:), zeros(3,abs(delta)));
                q = 0;
                for p = old_start:(old_start + abs(delta))
                    q = q+1;
                    train(k).H1{j,1}(1,p) = last_val + q;
                end
            else
                train(k).H1{j,1} = horzcat(zeros(3,abs(delta)), train(k).H1{j,1}(:,:));
                for p = 1:abs(delta)
                    train(k).H1{j,1}(1,p) = first_val - abs(delta) + p - 1;
                end
            end
        end
        if (delta > 0)
            if(train(k).H0{j,1}(1,1) < train(k).H1{j,1}(1,1))
                train(k).H0{j,1} = horzcat(train(k).H0{j,1}(:,:), zeros(3,abs(delta)));
                q = 0;
                for p = old_start:(old_start + abs(delta))
                    q = q+1;
                    train(k).H1{j,1}(1,p) = last_val + q;
                end
            else
                train(k).H0{j,1} = horzcat(zeros(3,abs(delta)),train(k).H0{j,1}(:,:));
                for p = 1:abs(delta)
                    train(k).H1{j,1}(1,p) = first_val - abs(delta) + p - 1;
                end
            end
        end
    end
end
%d: Calculate ML and MAP decision rule vectors
for k = 1:9
    for j = 1:7
        for p = 1 : length(train(k).H1{j,1}); %MK 3:34pm no need of the
       % loop
            for q = 1:7
                %p = length(train(k).H1{j,1}(1,:));
                if( train(k).H1{j,1}(3,p) >= train(k).H0{j,1}(3,p))
                    train(k).ML{j,1}(q,p) = 1;
                else 
                    train(k).ML{j,1}(q,p) = 0;
                end
                if((train(k).H1{j,1}(3,p) * prior_H1(j)) >= (train(k).H0{j,1}(3,p)*prior_H0(j)))
                    train(k).MAP{j,1}(q,p) = 1;
                else
                    train(k).MAP{j,1}(q,p) = 0;
                end
            end
        end
    end
end

for k = 1:9
    for j = 1:7
       
        %declare arrays for tabulation
        Xi = union(train(k).H1{j,1}(1,:), train(k).H0{j,1}(1,:)); %Xi the same size as max(H1,H0)
        H1 = train(k).H1{j,1}(3,:);
        H0 = train(k).H0{j,1}(3,:);
        ML = train(k).ML(j,:);
        MAP = train(k).MAP(j,:);
        
        HT_table_array{k, j} = table(Xi, H1, H0, ML, MAP);
          
    end
end

%e: save the results in a 9 by 7 cell array called HT_table_array
HT_table_array = cell(9,7);

%% Task 1.2
%Use HT_table_array to generate alarms based on ML/MAP
for k = 1:9
    for j = 1:7
        for p = 1:testing_length(k)
            test(k).ML(j,p) = HT_table_array{k,j}.ML(4);
            test(k).MAP(j,p) = HT_table_array{k,j}.MAP(4);
            for q = 1:length(HT_table_array{k,j}.Xi)
                value = floor(test(k).all_data(j,p)); Office hour code for building ML MAP arrays
                [~, idx] = find(HT_table_array{k,j}.Xi == value);
                if isnan(idx) == 0
                    h0_val = HT_table_array{k,j}.H0(idx);
                    h1_val = HT_table_array{k,j}.H1(idx);
                    h0_val_with_prior = h0_val*prior_H0(k);
                    h1_val_with_prior = h1_val*prior_H1(k);
                    if h0_val > h1_val
                        test(k).ML(j,p) = 0;
                    else
                        test(k).ML(j,p) = 1;
                    end
                
                    if h0_val_with_prior > h1_val_with_prior
                        test(k).MAP(j,p) = h0_val_with_prior;
                    else
                        test(k).MAP(j,p) = h1_val_with_prior;
                    end
                
                else
                    test(k).ML(j,p) = 1;
                    test(k).MAP(j,p) = 1;
                end
        end
%           for q = 1:4 MK 2:23 p      
%                 %if the data fall within the range given, assign the ML/MAP
%                 %value
%                 %if(test(k).all_data(j,p) <= HT_table_array{k,j}.Max_Value(5-q))
%                     if value == HT_table_array{k,j}.Xi(q);
%                         test(k).ML(j,p) = HT_table_array{k,j}.ML(j,k);
%                 end
%                 test(k).ML(j,p) = HT_table_array{k,j}.ML(5-q);
%                 test(k).MAP(j,p) = HT_table_array{k,j}.MAP(5-q);
%            end
        end %testing_length    
    end %7
end %9 

%b: evaluate ML/MAP rules
for k = 1:9 
    for j = 1:7
      missed_detect_MAP = 0;
      false_alarm_MAP = 0;
      error_MAP = 0;

      missed_detect_ML = 0;
      false_alarm_ML = 0;
      error_ML = 0;
      for i = 1:testing_length(k)
          %MAP
          if(test(k).MAP(j,i) == 0 & test(k).all_labels(i) == 1)
              missed_detect_MAP = missed_detect_MAP +1;
          end
          if (test(k).MAP(j,i) == 1 & test(k).all_labels == 0)
              false_alarm_MAP = false_alarm_MAP +1;
          end
          if((test(k).MAP(j,i) == 0 & test(k).all_labels == 1) | (test(k).MAP(j,i) == 1 & test(k).all_labels == 0))
              error_MAP = error_MAP +1;
          end
          %ML
          if(test(k).ML(j,i) == 0 & test(k).all_labels(i) == 1)
              missed_detect_ML = missed_detect_ML +1;
          end
          if (test(k).ML(j,i) == 1 & test(k).all_labels == 0)
              false_alarm_ML = false_alarm_ML +1;
          end
          if((test(k).ML(j,i) == 0 & test(k).all_labels == 1) | (test(k).ML(j,i) == 1 & test(k).all_labels == 0))
               error_ML = error_ML +1;
          end
      end
      num_nongoldens = testing_length(j) - sum(test(k).all_labels);
      num_goldens = sum(test(k).all_labels);
      num_totals = length(test(k).all_labels);
   
      test(k).MAP_MD(j) = missed_detect_MAP / num_goldens;
      test(k).ML_MD(j) = missed_detect_ML / num_goldens;
      test(k).MAP_FA(j) = false_alarm_MAP / num_nongoldens;
      test(k).ML_FA(j) = false_alarm_ML / num_nongoldens;
      test(k).MAP_E(j) = error_MAP / num_totals;
      test(k).ML_E(j) = error_ML / num_totals;
    end
end


Error_table_array = cell(9,7);

for k = 1:9
    for j = 1:7
           falseAlarm(:,1) = [test(k).ML_FA(j), test(k).MAP_FA(j)];
           missDetect(:,1) = [test(k).ML_MD(j), test(k).MAP_MD(j)];
           error(:,1) = [test(k).ML_E(j), test(k).MAP_E(j)];
           Error_table_array{k, j} = table(falseAlarm, missDetect, error);
          
    end
end
                
