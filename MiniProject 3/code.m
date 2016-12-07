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
% the above function creates non-golden data and golden data
% Task 1
% a: Calculate prior probabilities of P(H1) an dP(H0)
for k = 1:9
    prior_H1(k) = sum(train(k).all_labels)/training_length(k);
    prior_H0(k) = 1 - prior_H1(k);
end

%b: construct likelihood matrices for each of the seven features
feature_labels = {'Mean Area under the Heart Beat','Mean R-to-R peak interval','Number of beats per minute (Heart Rate)','Peak to peak interval for Blood Pressure','Systolic Blood Pressure','Diastolic Blood Pressure','Pulse Pressure'};
feature_max_length = [17, 220, 220, 220, 115, 86, 74];

likelihood_matrix_all_features = cell(9,7);
    
for k = 1:9
    for j = 1:7
        goldens_tabulated = tabulate(train(k).goldens(j,:))';
        nongoldens_tabulated = tabulate(train(k).nongoldens(j,:))';
        Xi = union(goldens_tabulated(1,:), nongoldens_tabulated(1,:));
        likelihood_matrix_all_features{k,j} = zeros(5,length(Xi));
        
        for idx = 1:length(Xi) %idx iterates over Xi
            likelihood_matrix_all_features{k,j}(1,idx) = Xi(idx);
            %populate h1 row in likelihood matrix of jth feature
            if ismember(Xi(idx),goldens_tabulated(1,:))
                golden_index_of_Xi_value = find(goldens_tabulated(1,:)==Xi(idx),1);
                likelihood_matrix_all_features{k,j}(2,idx) = goldens_tabulated(3, golden_index_of_Xi_value)/100;
            else
                likelihood_matrix_all_features{k,j}(2,idx) = 0;
            end
            
            %populate h0 row in likelihood matrix of jth feature
            if ismember(Xi(idx),nongoldens_tabulated(1,:))
                nongolden_index_of_Xi_value = find(nongoldens_tabulated(1,:)==Xi(idx),1);
                likelihood_matrix_all_features{k,j}(3,idx) = nongoldens_tabulated(3, nongolden_index_of_Xi_value)/100;
            else
                likelihood_matrix_all_features{k,j}(3,idx) = 0;
            end
            
        end
    end
end
        

for k = 1:9
    figure;
    for j = 1:7
       subplot(7, 1, j); 
        plot(likelihood_matrix_all_features{k,j}(3,:));
        
        % Add feature titles and set axis for each subplot
        title(feature_labels(j));
        axis([0 feature_max_length(j) 0 1]);
        
        hold on; 
        plot(likelihood_matrix_all_features{k,j}(2,:));
    end
    legend('H0 pmf', 'H1 pmf');
end

for k = 1:9
    for j = 1:7
        for i = 1:length(likelihood_matrix_all_features{k,j}(1,:))
                if(likelihood_matrix_all_features{k,j}(2,i) >= likelihood_matrix_all_features{k,j}(3,i))
                    likelihood_matrix_all_features{k,j}(4,i) = 1;
                else 
                    likelihood_matrix_all_features{k,j}(4,i) = 0;
                end
                if(likelihood_matrix_all_features{k,j}(2,i) * prior_H1(k) >= likelihood_matrix_all_features{k,j}(3,i)*prior_H0(k))
                    likelihood_matrix_all_features{k,j}(5,i) = 1;
                else
                    likelihood_matrix_all_features{k,j}(5,i) = 0;
                end
        end
    end
end

HT_table_array = likelihood_matrix_all_features;
for k = 1:9
    for j = 1:7
        HT_table_array{k,j} = HT_table_array{k,j}';
    end
end

Error_table_array = cell(9,7);

for k = 1:9
    prior_H1_test(k) = sum(test(k).all_labels);
    prior_H0_test(k) = testing_length(k) - prior_H1_test(k);
end

for k = 1:9
    for j = 1:7
        Error_table_array{k,j} = zeros(2,3);
        
        missed_detect_MAP = 0;
        false_alarm_MAP = 0;
        error_MAP = 0;
        
        missed_detect_ML = 0;
        false_alarm_ML = 0;
        error_ML = 0;
        
        for i = 1:testing_length(k)
            if ismember(test(k).all_data(j,i), HT_table_array{k,j}(:,1))
            %if that measurement value was in HT_table_array
                measurement_index = find(HT_table_array{k,j}(:,1)==test(k).all_data(j,i),1);
                ML = HT_table_array{k,j}(measurement_index,4); %ML value
                MAP = HT_table_array{k,j}(measurement_index,5); %MAP value
            else
            %if not then rule in favor of H1
                ML = 1;
                MAP = 1;
            end
            
            %physician alarm for this measurement
            golden_value = test(k).all_labels(i);
            
            %ML counter increments
            if (ML == 1 && golden_value == 0)
                false_alarm_ML = false_alarm_ML + 1;
            end
            
            if (ML == 0 && golden_value == 1)
                missed_detect_ML = missed_detect_ML + 1;
            end
            
            if ((ML == 1 && golden_value == 0) || (ML == 0 && golden_value == 1))
                error_ML = error_ML + 1;
            end
            
            %MAP counter increments
            if (MAP == 1 && golden_value == 0)
                false_alarm_MAP = false_alarm_MAP + 1;
            end
            
            if (MAP == 0 && golden_value == 1)
                missed_detect_MAP = missed_detect_MAP + 1;
            end
            
            if ((MAP == 1 && golden_value == 0) || (MAP == 0 && golden_value == 1))
                error_MAP = error_MAP + 1;
            end
        end
        Error_table_array{k,j}(1,1) = false_alarm_ML/prior_H0_test(k);
        Error_table_array{k,j}(1,2) = missed_detect_ML/prior_H1_test(k);
        Error_table_array{k,j}(1,3) = error_ML/testing_length(k);
        
        Error_table_array{k,j}(2,1) = false_alarm_MAP/prior_H0_test(k);
        Error_table_array{k,j}(2,2) = missed_detect_MAP/prior_H1_test(k);
        Error_table_array{k,j}(2,3) = error_MAP/testing_length(k);
    end
end

%% Task 2
%For every patient, find the lowest correlations
    %Make data sets that are all the same length (1412, shortest data 
    % length) for comparison
    feature_array = {'Mean Area under the heart beat','Mean R-to-R peak interval','Number of beats per minute (Heart Rate)','Peak to peak interval for Blood pressure','Systolic Blood Pressure','Diastolic Blood Pressure','Pulse Pressure'};
for k = 1:9
    for j = 1:7
        for h = 1:1412
        train(k).short(j, h) = train(k).all_data(j,h);
        end
    end
end

for k = 1:9
    for j = 1:7
        for h = 1:7
            correlation = corrcoef(train(k).short(j,:), train(k).short(h,:));
            train(k).corr(j,h) = correlation(2,1);
        end
    end
end

%Two features with lowest MAP errors (MAP error lower than ML error)
for k = 1:9
    for j = 1:7
        error_table(k,j) = Error_table_array{k,j}(2,3);
    end
end

%combine these features to find the best two features to use
for k = 1:9
    for j = 1:7
        for h = 1:7
            error_array(h,j) = abs(error_table(k,j)*error_table(h,j));
            patient_pref(k,j,h) = abs(error_table(k,j)*error_table(k,h)*train(k).corr(h,j));
        end
    end
end

for k = 1:9
    best_error(k) = min(min(patient_pref(k,:,:)));
    for j = 1:7
        for h = 1:7
            if patient_pref(k,j,h) == best_error(k)
                features(k,1) = h;
                features(k,2) = j;
            end
        end
    end
end

