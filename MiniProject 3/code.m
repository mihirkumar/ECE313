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
    test(k).all_data = dat_array(k).all_data(:, training_length(k)+1:all_length(k));
    test(k).all_labels = dat_array(k).all_labels(:, training_length(k)+1:all_length(k));
end
for k = 1:9
   for j = 1:7
      index1 = 1;
      index0 = 1;
      for i = 1: training_length(k)
          if(train(k).all_labels(1, i) == 1)
              train(k).goldens(j, index1) = train(k).all_data(j,i);
              index1 = index1 + 1;
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
    figure;
    for j = 1:7
        likelyH1= tabulate(train(k).goldens(j,:));
            likelyH1(:,3) = likelyH1(:,3) / 100;
        likelyH0 = tabulate(train(k).nongoldens(j,:));
            likelyH0(:,3) = likelyH0(:,3) / 100;
        subplot(7, 1, j); 
        plot(likelyH0(:,3)); 
        hold on; 
        plot(likelyH1(:,3));
    %figure out quartiles and save them
        max1 = max(likelyH1(:,1)); 
        max2 = max(likelyH0(:,1));
        max_val = max(max1,max2);
        min1 = min(likelyH1(:,1)); 
        min2 = min(likelyH0(:,1));
        min_val = min(min1,min2);
        train(k).middle = floor((max_val - min_val)/2);
        train(k).first_quarter = floor((train(k).middle - min_val)/2);
        train(k).third_quarter = floor((max_val - train(k).middle)/2) + train(k).middle;
        train(k).max = max_val;
        train(k).H1(j,1:4) = zeros();
        train(k).H0(j,1:4) = zeros();
        for i = 1:size(likelyH1) 
           if (likelyH1(i,1) < train(k).first_quarter)
                percent = likelyH1(i,3);
                p = 1;
           elseif (likelyH1(i,1) < train(k).middle)
                percent = likelyH1(i,3);
                p = 2;
           elseif (likelyH1(i,1) < train(k).third_quarter)
                percent = likelyH1(i,3);
                p = 3;
           else
                percent = likelyH1(i,3);
                p = 4;
           end
           train(k).H1(j,p) = train(k).H1(j,p)+ percent;
        end
        
        for i = 1:size(likelyH0) 
           if (likelyH0(i,1) < train(k).first_quarter) 
                percent = likelyH0(i,3);
                p = 1;
           elseif (likelyH0(i,1) < train(k).middle)
                percent = likelyH0(i,3);
                p = 2;
           elseif (likelyH0(i,1) < train(k).third_quarter)
                percent = likelyH0(i,3);
                p = 3;
           else
                percent = likelyH0(i,3);
                p = 4;
           end
           train(k).H0(j,p) = train(k).H0(j,p)+ percent;
        end
    end
    
    legend('H0 pmf', 'H1 pmf');
    
end

%c: show results by generating a seperate figure for each patient 
%c Executed in B

%d: Calculate ML and MAP decision rule vectors
for k = 1:9
    for j = 1:7
        for p = 1:4
            if( train(k).H1(j,p) >= train(k).H0(j,p))
                train(k).ML(j,p) = 1;
            else 
                train(k).ML(j,p) = 0;
            end
            if(train(k).H1(j,p) * prior_H1 >= train(k).H0(j,p)*prior_H0)
                train(k).MAP(j,p) = 1;
            else
                train(k).MAP(j,p) = 0;
            end
        end
    end
end
%e: save the results in a 9 by 7 cell array called HT_table_array
%we'll figure this out later today before the report is due but it is
%cosmetic so I'm moving on
HT_table_array = cell(9,7);
for k = 1:9
    for j = 1:7
        HT_table_array(k, j) = mat2cell(HT_table(j));
           HT_table(j) = [train(k).first_quarter, train(k).H1(j,1), train(k).H0(j,1), train(k).ML(j,1), train(k).MAP(j, 1);
                          train(k).middle, train(k).H1(j,2), train(k).H0(j,2), train(k).ML(j,2), train(k).MAP(j, 2);
                          train(k).third_quarter, train(k).H1(j,3), train(k).H0(j,3), train(k).ML(j,3), train(k).MAP(j, 3);
                          train(k).max, train(k).H1(j,4), train(k).H0(j,4), train(k).ML(j,4), train(k).MAP(j,4)];
    end
end
%Task 1.2

                
