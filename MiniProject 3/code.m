close all;
clear all;

%%Task 0
%Load in data sets and make it readable.
dat1_const = load('1_a41178.mat');
dat2_const  = load('2_a42126.mat');
dat3_const = load('3_a40076.mat');
dat4_const = load('4_a40050.mat');
dat5_const = load('5_a41287.mat');
dat6_const = load('6_a41846.mat');
dat7_const = load('7_a41846.mat');
dat8_const = load('8_a42008.mat');
dat9_const = load('9_a41846.mat');

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
    test(k).all_data = dat_array(k).all_data(:, training_length(k):all_length(k));
    test(k).all_labels = dat_array(k).all_labels(:, training_length(k):all_length(k));
end
for k = 1:9
   for j = 1:7
      index1 = 1;
      index0 = 1;
      for i = 1: training_length(k)
          if(train(k).all_labels(1, i) == 1)
              train(k).goldens(j, index1) = train(k).all_data(j,i);
              index1 = index1 + 1
          end
          if(train(k).all_labels(1, i) == 0)
              train(k).nongoldens(j,index0) = train(k).all_data(j,i);
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
    prior_H1(k) = sum(train(1).all_labels)/training_length(1);
    prior_H0(k) = 1 - prior_H1(k);
end
%b: construct likelihood matrices for each of the seven features
for k = 1:9
    for j = 1:7
        train(k).likelyH1= tabulate(train(k).goldens(j,:));
            train(k).likelyH1(:,3) = train(k).likelyH1(:,3) / 3;
        train(k).likelyH0 = tabulate(train(k).nongoldens(j,:));
            train(k).likelyH0(:,3) = train(k).likelyH0(:,3) / 3;
    end
end
labels = {'MEAN AREA UNDER HEART BEAT', 'MEAN R TO R PEAK INTERVAL', 'NUMBER OF BEATS PER MINUTE', 'PEAK TO PEAK INTERVAL FOR BP','SYSTOLIC BP', 'DIASTOLIC BP', 'PULSE PRESSURE'}; 
for k = 1:9
    for j = 1:7
        percent_1st = 0;
        percent_2nd = 0;
        percent_3rd = 0;
        percent_4th = 0;
        index = 0;
        max1 = max(train(k).likelyH1(:,1)); 
        max2 = max(train(k).likelyH0(:,1));
        max_val = max(max1,max2);
        min1 = min(train(k).likelyH1(:,1)); 
        min2 = min(train(k).likelyH0(:,1));
        min_val = min(min1,min2);
        middle = floor((max_val - min_val)/2);
        first_quarter = floor((middle - min_val)/2);
        third_quarter = floor((max_val - middle)/2);
        while(train(k).likelyH1(:,1) < first_quarter)
            percent_1st = train(k).likelyH1(:,1) + percent_first;
        end
    end
end
%c: show results by generating a seperate figure for each patient
for k = 1:9
 subplot(7, 1, k); 
 plot(H0_pmf(:)); 
 hold on; 
 plot(H1_pmf(:));
 %legend(’H0 pmf', ’H1 pmf’);
end
%%Task 2
