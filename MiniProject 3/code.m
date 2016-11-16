close all;
clear all;

patient1 = load('1_a41178.mat');
patient2 = load('2_a42126.mat');
patient3 = load('3_a40076.mat');
patient4 = load('4_a40050.mat');
patient5 = load('5_a41287.mat');
patient6 = load('6_a41846.mat');
patient7 = load('7_a41846.mat');
patient8 = load('8_a42008.mat');
patient9 = load('9_a41846.mat');

patient_array = [patient1, patient2, patient3, patient4, patient5, patient6, patient7, patient8, patient9];

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


%old code underneath

%training = all_data(:,1:2866)
%label_training = all_labels(1:2866)

%testing = all_data(:, 2867:4299)
%label_testing = all_labels(2867:4299)

for k = 1:9
	temp = tabulate(patient_array(k).label_training);
	% 1,k is for H0
	h0h1_data(1,k) = temp(1,3)/100;
	
	% 2,k is for H1
	h0h1_data(2,k) = temp(2,3)/100;

	full_length = length(patient_array(k).training);

	for i = 1:7 %for each signal type
		H0_vector = [];
		H1_vector = [];

		for j = 1:full_length
			if patient_array(k).label_training(j) == 1
				H1_vector = [H1_vector patient_array(k).training(i, j)];
			end
			if patient_array(k).label_training(j) == 0
				H0_vector = [H0_vector patient_array(k).training(i, j)];
			end
		end

		H0_pmf = transpose(tabulate(H0_vector));
		H1_pmf = transpose(tabulate(H1_vector));

		likelihoodMatrixLength = max(length(H0_pmf),length(H1_pmf));
		likelihoodMatrix = zeros(3, likelihoodMatrixLength);

		for x = 1:likelihoodMatrixLength
			if H1_pmf
				likelihoodMatrix(1, x) = H1_pmf(1,x);
			else
				likelihoodMatrix(1,x) = H0_pmf(1,x);
			end
		end

		for x = 1:likelihoodMatrixLength
			if H1_pmf
				likelihoodMatrix(2,x) = H1_pmf(3,x)/100;
			end
			if H0_pmf
				likelihoodMatrix(3,x) = H0_pmf(3,x)/100;
			end
		end
	end
end
