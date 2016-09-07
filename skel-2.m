%   FILL IN / MODIFY THE CODE WITH "" or comments with !!

% import monitor_alarms.mat and put it in a 
load('\\ad.uillinois.edu\engr\Instructional\calbers2\documents\MATLAB\monitor_alarms.mat');
data = MonitorAlarms;
% not that the data is in a table format,
% try the following in the command window:
%   data
%   data(1)
%   data.StartTime
%   data.StartTime(1)
% now you have learned how to access data in the table

% remove single quotes from the time data
%data.StartTime = strrep(data.StartTime, 0, ); 
%data.StopTime = strrep(data.StopTime, 2, 4);

% open the result file
% !! replace # with your own groupID
fid = fopen('ECE313_Mini1_group5', 'w');

% T1.1
% !! subset your data for each alarm_type
data_SYSTEM = data(ismember(data.Alarm_Type, 'SYSTEM_ALARM'), :);
data_ADVISORY = data(ismember(data.Alarm_Type, 'ADVISORY_ALARM'), :);
data_WARNING = data(ismember(data.Alarm_Type, 'WARNING_ALARM'), :);
data_CRISIS = data(ismember(data.Alarm_Type, 'CRISIS_ALARM'), :);

% !! count the number of alarms for each alarm_type
numSYSTEM = height(data_SYSTEM);
numADVISORY = height(data_ADVISORY);
numWARNING = height(data_WARNING);
numCRISIS = height(data_CRISIS);

% !! calaculate the probability for each alarm_type

numTOTAL = numSYSTEM + numADVISORY + numWARNING + numCRISIS;
fprintf(fid, 'Task 1.1 1\n\n');
fprintf(fid, 'P(SYSTEM) = %f\n',  numSYSTEM/numTOTAL);
fprintf(fid, 'P(ADVISORY) = %f\n', numADVISORY/numTOTAL);
fprintf(fid, 'P(WARNING) = %f\n', numWARNING/numTOTAL);
fprintf(fid, 'P(CRISIS) = %f\n',  numCRISIS/numTOTAL);

% T1.2. 
causeSYSTEM = {'APP_ERR', 'SIG_ARTIFACT', 'LEADS_FAILURE', 'NW_ERR'};
causeWARNING = {'LOW_OXY_SAT', 'HA_BRADY', 'SLEEP_DISORDER', 'PAUSE'};

fprintf(fid, '\n\nTask 1.2\n\n');
fprintf(fid, 'Probability of causes for SYSTEM alarms\n');
for i=1:4,
    % !! subset your data
    %data_SYSTEM_curr = data_SYSTEM(ismember(data_SYSTEM.Cause, causeSYSTEM{i}), :);
    % !! do the counts to derive your answers
    data_SYSTEM_curr = height(data_SYSTEM(ismember(data_SYSTEM.Cause, causeSYSTEM{i}), :));
    prob_i = data_SYSTEM_curr/numSYSTEM;
    fprintf(fid, 'P(%s) = %f\n', cell2mat(causeSYSTEM(i)), prob_i);
end
fprintf(fid, 'Probability of causes for WARNING alarms\n');
for i=1:4,
    % !! subset your data
    % !! do the counts to derive your answers
    % fprintf(fid, 'P(%s) = %f\n', cell2mat(causeWARNING(i)), "PROBABILITY FOR CAUSE(i) and WARNING");
    data_WARNING_curr = height(data_WARNING(ismember(data_WARNING.Cause, causeWARNING{i}), :));
    prob_i = data_WARNING_curr/numWARNING;
    fprintf(fid, 'P(%s) = %f\n', cell2mat(causeWARNING(i)), prob_i);
end
    
% T1.3.

% !! Split the data in terms of hours of the start time
% Please note that the time format is 'HH:MM:SS.FFF'
fprintf(fid, '\n\nTask 1.3\n\n');
% Print the result
fprintf(fid, '\t\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t22\t23\n');
fprintf(fid, 'count\t');
hourSYSTEM = {1:24};
hourWARNING = {1:24};
hourADVISORY = {1:24};
hourCRISIS = {1:24};
formatIn = 'HH:MM:SS.FFF';
zero_str = '00:00:00.000';
time_str = '00:17:51.000';
one_constant = datenum(zero_str, formatIn) - datenum(zero_str, formatIn);
fprintf('%s ', (daten
for i=1:24,
    % !! Count the number of alarms for the given hour(i)
um(time_str, formatIn) - datenum(zero_str, formatIn)));
    %time_freq_SYSTEM = height(data_SYSTEM(ismember(data_SYSTEM.StartTime, date)
    %fprintf(fid, '%d\t', "COUNT_RESULT");
end
fprintf(fid, '\n');
figure;
% !! plot the histograms of alarms per hour

% Labeling the figure
title('Histogram of Alarms per hour');
hours = {'00h', '01h', '02h','03h','04h','05h','06h','07h','08h','09h','10h','11h','12h','13h','14h','15h','16h','17h','18h','19h','20h','21h','22h','23h'};
set(gca, 'XTick', 0:23);
set(gca,'XTickLabel',hours);
ylabel('number of alarms');

% T1.4.

% !! convert the time into units of minutes
% e.g. 01:35:01.432 is the 95th minute of a day

% !! count how many minutes contain at least one alarm

fprintf(fid, '\n\nTask 1.4\n\n');
%fprintf(fid, 'P(Alarm In a minute) = %f\n', "PROBABILITY OF AN ALARM IN A RANDOMLY CHOSEN MINUTE");

% T2.1 

% !! Using the results from Task 1, derive the probability

fprintf(fid, '\n\nTask 2.1\n\n');
%fprintf(fid, 'P(LOW_OXY_SAT and SYSTEM) = %f\n', "PROBABILITY OF LOW_OXY_SAT and SYSTEM");

% T2.2

% !! Calculate the duration for each alarms
% make sure that the durations are in seconds

% T2.2.a 

fprintf(fid, '\n\nTask 2.2.a\n\n');
%for i=1:"NUMBER OF ALARM TYPES",
    % !! Split the data interms of alarm type
    
    % !! Count the number of alarms for each alarm type
    
    % !! Calculate the average duration of each alarm type
    
    %fprintf(fid, 'Average duration for alarm type %s: %f\n', "ALARMTYPE", "AVERAGE_DURATION");
%end

% T2.2.b

fprintf(fid, '\n\nTask 2.2.b\n\n');
for i=1:24,
    % Please not that i loop from 1 to 24
    % The hours in the data are from 0 to 23

    % !! Split the data in terms of hours
    
    % !! Calculate the average duration for each hour(i)
    
    %fprintf(fid, 'Average duration for hh=%d = %f\n', "HOUR", "AVERAGE DURATION PER HOUR");
end
figure;
% !! Draw a bar chart to plot the average duration per hour

% label the plot
title('Average Duration for each hour of the day');
ylabel('avg duration');
hours = {'00h', '01h', '02h','03h','04h','05h','06h','07h','08h','09h','10h','11h','12h','13h','14h','15h','16h','17h','18h','19h','20h','21h','22h','23h'};
set(gca, 'XTick', 0:23);
set(gca,'XTickLabel',hours);

% BONUS Find the two patient beds with the largest number of alarms,
fprintf(fid, '\n\nBONUS\n\n');
% !! No codes provided (Bonus Question)


% T3.
fprintf(fid, '\n\nTask3\n\n');
% !! Extract (split) the data that applies to your group

% !! Write your own code for analysis (ref. codes for Task1 and Task2)

fclose(fid);