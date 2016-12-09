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
% a: Calculate prior probabilities of P(H1) and P(H0)
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
for k = 1:9
    for j = 1:7
        for h = 1:7
            %get a correlation matrix
            correlation = corrcoef(train(k).all_data(j,:), train(k).all_data(h,:));
            train(k).corr(j,h) = correlation(2,1);
        end
    end
end
%Make data sets that are all the same length for comparison



%Minimize error
%error_table: 
%feature|MAP e|MAP md||MAP fa|ML e|ML md|ML fa|
testing_col = 2; %MAP error
for k = 1:9
    for j = 1:7
        train(k).error_table(j,1) = j;
        train(k).error_table(j,2) = Error_table_array{k,j}(2,3);
        train(k).error_table(j,3) = Error_table_array{k,j}(2,2);
        train(k).error_table(j,4) = Error_table_array{k,j}(2,1);
        train(k).error_table(j,5) = Error_table_array{k,j}(1,3);
        train(k).error_table(j,6) = Error_table_array{k,j}(1,2);
        train(k).error_table(j,7) = Error_table_array{k,j}(1,1);
    end
    train(k).error_table = sortrows(train(k).error_table,testing_col); 
end
%better (lower) errors will sort to the top. Use
% sortrows(train(k).error_table, COL_INDEX) to get best features for each
% decision facet.

%use error_table to get the lowest 3 errors, and pick the two with the
%lowest correlation
%FEATURES:
%|Feature 1|Feature2|Error of 1|Error of 2|Correlation of 1,2|
for k = 1:8 % Leave off 9 because 9 and 6 are identical
    corr12 = train(k).corr(train(k).error_table(1,1), train(k).error_table(2,1));
    corr13 = train(k).corr(train(k).error_table(1,1), train(k).error_table(3,1));
    corr23 = train(k).corr(train(k).error_table(1,1), train(k).error_table(3,1));
    lowest_corr = min([corr12, corr13, corr23]);
    if corr12 == lowest_corr
        features(k,1) = train(k).error_table(1,1);
        features(k,2) = train(k).error_table(2,1);
        features(k,3) = train(k).error_table(1,testing_col);
        features(k,4) = train(k).error_table(2,testing_col);
        features(k,5) = corr12;
    elseif corr13 == lowest_corr
        features(k,1) = train(k).error_table(1,1);
        features(k,2) = train(k).error_table(3,1);
        features(k,3) = train(k).error_table(1,testing_col);
        features(k,4) = train(k).error_table(3,testing_col);
        features(k,5) = corr13;
    else
        features(k,1) = train(k).error_table(2,1);
        features(k,2) = train(k).error_table(3,1);
        features(k,3) = train(k).error_table(2,testing_col);
        features(k,4) = train(k).error_table(3,testing_col);
        features(k,5) = corr23;
    end
    features(k,6) = features(k,3)*features(k,4)*features(k,5);
    features(k,7) = k;
end

for k = 1:9
    for j = 1:7
        for h = 1:7
            error_array(j,h) = Error_table_array{k,j}(2,3)*Error_table_array{k,h}(2,3);
        end
    end
    figure;
    surf(abs(train(k).corr));
    title(['Correlation, patient ', num2str(k)]);
    xlabel('Feature 1');
    ylabel('Feature 2');
    figure;
    surf(abs(error_array));
    title(['Error, patient ', num2str(k)]);
    xlabel('Feature 1');
    ylabel('Feature 2');
end
 features = sortrows(features,6);
%% Task 3.1
%3.1a - Generate likelihood matrices from feature pairs
%Assuming independent features, so P(X=k,Y=j) = P(X=k)P(Y=j)
patient_itt = [features(1,7), features(2,7), features(3,7)]; % use this to itterate over all 3 patients
feature_itt = [features(1,1) features(1,2); features(2,1) features(2,2); features(3,1) features(3,2);];
%patient_itt = [1, 3, 5];
%feature_itt = [2 3; 2 3; 2 3;]
for p = 1:3
    pos = 0;
    k_max = length(HT_table_array{patient_itt(p),feature_itt(p,1)}(:,1));
    j_max = length(HT_table_array{patient_itt(p),feature_itt(p,2)}(:,1));
    for k = 1:k_max; % Iterate over all of the first feature
        for j = 1:j_max % Iterate over all of the second feature
            pos = pos + 1;
            Joint_HT_table{p}(pos,1) = HT_table_array{patient_itt(p),feature_itt(p,1)}(k,1);
            Joint_HT_table{p}(pos,2) = HT_table_array{patient_itt(p),feature_itt(p,2)}(j,1);
            % Calculate P(X=k,Y=j | H1)
            Joint_HT_table{p}(pos,3) = HT_table_array{patient_itt(p),feature_itt(p,1)}(k,2) * HT_table_array{patient_itt(p),feature_itt(p,2)}(j,2);
            % Calculate P(X=k,Y=j | H0)
            Joint_HT_table{p}(pos,4) = HT_table_array{patient_itt(p),feature_itt(p,1)}(k,3) * HT_table_array{patient_itt(p),feature_itt(p,2)}(j,3);


            %3.1b - Generate MAP and ML decision rule vectors
            % ML first
            if(Joint_HT_table{p}(pos,3) >= Joint_HT_table{p}(pos,4))
                Joint_HT_table{p}(pos,5) = 1;
            else
                Joint_HT_table{p}(pos,5) = 0;
            end
            % MAP with H1 favorable for breaking ties
            if(Joint_HT_table{p}(pos,3) * prior_H1(patient_itt(p)) >= Joint_HT_table{p}(pos,4) * prior_H0(patient_itt(p)))
                Joint_HT_table{p}(pos,6) = 1;
            else
                Joint_HT_table{p}(pos,6) = 0;
            end

            %prep for 3.1d
            mesh_table{p,1}(k,j) = Joint_HT_table{p}(pos,3);
            mesh_table{p,2}(k,j) = Joint_HT_table{p}(pos,4);
        end
    end
end


%3.1d - Plot conditional joint PMFs
for p = 1:3 % TODO: Rename the below titles to be clearer
    title(['H0 Conditional Joint PMF for Patient ', num2str(patient_itt(p))]);
    figure;
    mesh(mesh_table{p,1}); % H0 hypothesis
    title(['H1 Conditional Joint PMF for Patient ', num2str(patient_itt(p))]);
    figure;
    mesh(mesh_table{p,2}); % H1 hypothesis
end

%% Task 3.2
%3.2a Calculate alarms for the testing data set based on
for p = 1:3 % All 3 patients
    for sample_itt = 1:testing_length(patient_itt(p));
        f1_itt = test(patient_itt(p)).all_data(feature_itt(p,1),sample_itt);
        f2_itt = test(patient_itt(p)).all_data(feature_itt(p,2),sample_itt);
        %sample_val is the index where we have the requested values for
        %feature 1 and feature 2
        sample_val = find((Joint_HT_table{p}(:,1) == f1_itt) & (Joint_HT_table{p}(:,2) == f2_itt));
        if(isempty(sample_val))
            ML_test_alarms{p}(sample_itt) = 1;
            MAP_test_alarms{p}(sample_itt) = 1; 
        else
            ML_test_alarms{p}(sample_itt) = Joint_HT_table{p}(sample_val,5);
            MAP_test_alarms{p}(sample_itt) = Joint_HT_table{p}(sample_val,6);
        end
        majority_test_alarms{p}(sample_itt) = ML_test_alarms{p}(sample_itt) && MAP_test_alarms{p}(sample_itt);
    end
end



%3.2b Calculate probabilities
for p = 1:3
        Test_data_error_table_array{p} = zeros(2,3);

        missed_detect_MAP = 0;
        false_alarm_MAP = 0;
        error_MAP = 0;

        missed_detect_ML = 0;
        false_alarm_ML = 0;
        error_ML = 0;

        for j = 1:testing_length(patient_itt(p))
            ML = ML_test_alarms{p}(j);
            MAP = MAP_test_alarms{p}(j);

            %physician alarm for this measurement
            golden_value = test(patient_itt(p)).all_labels(j);

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
        Test_data_error_table_array{p}(1,1) = false_alarm_ML/prior_H0_test(patient_itt(p));
        Test_data_error_table_array{p}(1,2) = missed_detect_ML/prior_H1_test(patient_itt(p));
        Test_data_error_table_array{p}(1,3) = error_ML/testing_length(patient_itt(p));

        Test_data_error_table_array{p}(2,1) = false_alarm_MAP/prior_H0_test(patient_itt(p));
        Test_data_error_table_array{p}(2,2) = missed_detect_MAP/prior_H1_test(patient_itt(p));
        Test_data_error_table_array{p}(2,3) = error_MAP/testing_length(patient_itt(p));
end

%3.2c - Plot generated alarms
for p = 1:3
    figure;

    % MAP
    subplot(4, 1, 1);
    bar(MAP_test_alarms{p});
    title(['MAP for Patient ', num2str(patient_itt(p))]);

    %ML
    subplot(4, 1, 2);
    bar(ML_test_alarms{p});
    title(['ML for Patient ', num2str(patient_itt(p))]);
    
    %Majority Voter
    subplot(4, 1, 3);
    bar(majority_test_alarms{p});
    title(['Majority Voter for Patient ', num2str(patient_itt(p))]);    
    

    % Golden
    subplot(4, 1, 4);
    bar(test(patient_itt(p)).all_labels);
    title(['Golden Alarms for Patient ', num2str(patient_itt(p))]);
end

%% Task 3.3

% Calculate average P(error) across all 3 patients for MAP and ML
%average_ML_error = (Test_data_error_table_array{1}(1,3) + Test_data_error_table_array{2}(1,3) + Test_data_error_table_array{3}(1,3))/3;
%average_MAP_error = (Test_data_error_table_array{1}(1,4) + Test_data_error_table_array{2}(1,4) + Test_data_error_table_array{3}(1,4))/3;
