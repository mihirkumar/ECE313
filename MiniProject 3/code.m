close all;
clear all;

patient1 = load('1_a41178.mat');
patient2 = load('2_a42126.mat');
patient3 = load('3_a40076.mat');
patient4 = load('4_a40050.mat');
patient5 = load('5_a41287.mat');
patient6 = load('6_a41846.mat');
patient7 = load('7_a41846.mat');

patient_array = [patient1, patient2, patient3, patient4, patient5, patient6, patient7];

patient_signals = {
	'Mean Area under the heart beat',
	'Mean R-to-R peak interval',
	'Number of beats per minute (Heart Rate)',
	'Peak to peak interval for Blood pressure',
	'Systolic Blood Pressure',
	'Diastolic Blood Pressure',
	'Pulse Pressure'
};

for k = 1:9
	patient_array(k).all_data = floor(patient_array(k).all_data);
	splitpoint = floor(length(patient_array(k).all_data)*2/3);
	full_length = length(patient_array(k).all_data);
	patient_array(k).training = patient_array(k).all_data(:,1:splitpoint);
	patient_array(k).label_training = patient_array(k).all_labels(:,1:splitpoint);

	patient_array(k).testing = patient_array(k).all_data(:,splitpoint+1:full_length);
	patient_array(k).label_testing = patient_array(k).all_labels(:,splitpoint+1:full_length);
end

table = cell(9,7);

%old code underneath
all_data = floor(all_data);

training = all_data(:,1:2866)
label_training = all_labels(1:2866)

testing = all_data(:, 2867:4299)
label_testing = all_labels(2867:4299)

h0_occurrences = 0;
h1_occurrences = 0;

for k=1:2866
	if label_training(k) == 0
		h0_occurrences = h0_occurrences + 1;
	end

	if label_training(k) == 1
		h1_occurrences = h1_occurrences + 1;
	end
end

h0_probability = h0_occurrences/2866;
h1_probability = h1_occurrences/2866;

% Minimum value of each signal
var1_min = 0; % 1. Mean Area under the heart beat
var1_max = 0;

var2_min = 0; % 2. Mean R-to-R peak interval
var2_max = 0;

var3_min = 0; % 3. Number of beats per minute (Heart Rate)
var3_max = 0;

var4_min = 0; % 4. Peak to peak interval for Blood pressure
var4_max = 0;

var5_min = 0; % 5. Systolic Blood Pressure
var5_max = 0;

var6_min = 0; % 6. Diastolic Blood Pressure
var6_max = 0;

var7_min = 0; % 7. Pulse Pressure
var7_max = 0;

var1_min = min(training(1, :));
var2_min = min(training(2, :));
var3_min = min(training(3, :));
var4_min = min(training(4, :));
var5_min = min(training(5, :));
var6_min = min(training(6, :));
var7_min = min(training(7, :));

var1_max = max(training(1, :));
var2_max = max(training(2, :));
var3_max = max(training(3, :));
var4_max = max(training(4, :));
var5_max = max(training(5, :));
var6_max = max(training(6, :));
var7_max = max(training(7, :));
%for var1
h0_counter = 0;
h1_counter = 0;

var1_min_iterator = var1_min
var1_max_iterator = var1_max
index = 1;
for i=var1_min_iterator:var1_max_iterator
	for k=1:2866
		if training(1, k) == var1_min_iterator
			if label_training(1, k) == 0
				h0_counter = h0_counter+1;
			end

			if label_training(1, k) == 1
				h1_counter = h1_counter+1;
			end
		end
	end
	var1_h0_probability = h0_counter/h0_occurrences;
	var1_h1_probability = h1_counter/h1_occurrences;

	var1_likelihood_matrix(1,index) = var1_h1_probability;
	var1_likelihood_matrix(2,index) = var1_h0_probability;
	index = index+1;
end

%I'm repeating this code for each of the signal types, all the varX's
%for var2
h0_counter = 0;
h1_counter = 0;

var2_min_iterator = var2_min
var2_max_iterator = var2_max
index = 1;
for i=var2_min_iterator:var2_max_iterator
	for k=1:2866
		if training(1, k) == var2_min_iterator
			if label_training(1, k) == 0
				h0_counter = h0_counter+1;
			end

			if label_training(1, k) == 1
				h1_counter = h1_counter+1;
			end
		end
	end
	var2_h0_probability = h0_counter/h0_occurrences;
	var2_h1_probability = h1_counter/h1_occurrences;

	var2_likelihood_matrix(1,index) = var1_h1_probability;
	var2_likelihood_matrix(2,index) = var1_h0_probability;
	index = index+1;
end

%for var3
h0_counter = 0;
h1_counter = 0;

var3_min_iterator = var3_min
var3_max_iterator = var3_max
index = 1;
for i=var3_min_iterator:var3_max_iterator
	for k=1:2866
		if training(1, k) == var3_min_iterator
			if label_training(1, k) == 0
				h0_counter = h0_counter+1;
			end

			if label_training(1, k) == 1
				h1_counter = h1_counter+1;
			end
		end
	end
	var3_h0_probability = h0_counter/h0_occurrences;
	var3_h1_probability = h1_counter/h1_occurrences;

	var3_likelihood_matrix(1,index) = var3_h1_probability;
	var3_likelihood_matrix(2,index) = var3_h0_probability;
	index = index+1;
end

%for var4
h0_counter = 0;
h1_counter = 0;

var4_min_iterator = var4_min
var4_max_iterator = var4_max
index = 1;
for i=var4_min_iterator:var4_max_iterator
	for k=1:2866
		if training(1, k) == var4_min_iterator
			if label_training(1, k) == 0
				h0_counter = h0_counter+1;
			end

			if label_training(1, k) == 1
				h1_counter = h1_counter+1;
			end
		end
	end
	var4_h0_probability = h0_counter/h0_occurrences;
	var4_h1_probability = h1_counter/h1_occurrences;

	var4_likelihood_matrix(1,index) = var4_h1_probability;
	var4_likelihood_matrix(2,index) = var4_h0_probability;
	index = index+1;
end

%for var5
h0_counter = 0;
h1_counter = 0;

var5_min_iterator = var5_min
var5_max_iterator = var5_max
index = 1;
for i=var5_min_iterator:var5_max_iterator
	for k=1:2866
		if training(1, k) == var5_min_iterator
			if label_training(1, k) == 0
				h0_counter = h0_counter+1;
			end

			if label_training(1, k) == 1
				h1_counter = h1_counter+1;
			end
		end
	end
	var5_h0_probability = h0_counter/h0_occurrences;
	var5_h1_probability = h1_counter/h1_occurrences;

	var5_likelihood_matrix(1,index) = var5_h1_probability;
	var5_likelihood_matrix(2,index) = var5_h0_probability;
	index = index+1;
end


%for var6
h0_counter = 0;
h1_counter = 0;

var6_min_iterator = var6_min
var6_max_iterator = var6_max
index = 1;
for i=var6_min_iterator:var6_max_iterator
	for k=1:2866
		if training(1, k) == var1_min_iterator
			if label_training(1, k) == 0
				h0_counter = h0_counter+1;
			end

			if label_training(1, k) == 1
				h1_counter = h1_counter+1;
			end
		end
	end
	var6_h0_probability = h0_counter/h0_occurrences;
	var6_h1_probability = h1_counter/h1_occurrences;

	var6_likelihood_matrix(1,index) = var6_h1_probability;
	var6_likelihood_matrix(2,index) = var6_h0_probability;
	index = index+1;
end


%for var7
h0_counter = 0;
h1_counter = 0;

var7_min_iterator = var7_min
var7_max_iterator = var7_max
index = 1;
for i=var7_min_iterator:var7_max_iterator
	for k=1:2866
		if training(1, k) == var7_min_iterator
			if label_training(1, k) == 0
				h0_counter = h0_counter+1;
			end

			if label_training(1, k) == 1
				h1_counter = h1_counter+1;
			end
		end
	end
	var7_h0_probability = h0_counter/h0_occurrences;
	var7_h1_probability = h1_counter/h1_occurrences;

	var7_likelihood_matrix(1,index) = var7_h1_probability;
	var7_likelihood_matrix(2,index) = var7_h0_probability;
	index = index+1;
end

