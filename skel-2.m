%   FILL IN / MODIFY THE CODE WITH "" or comments with !!

% import monitor_alarms.mat and put it in a 
load('\\ad.uillinois.edu\engr\Instructional\alitz2\documents\MATLAB\monitor_alarms.mat');
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

% !! calculate the probability for each alarm_type
%divide overall by hieght of each
%calculate overall
numTOTAL = numSYSTEM + numADVISORY + numWARNING + numCRISIS;
fprintf(fid, 'Task 1.1 1\n\n');
%print averages in out file
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
    % !! do the counts to derive your answers
    %find height of each cause
    data_SYSTEM_curr = height(data_SYSTEM(ismember(data_SYSTEM.Cause, causeSYSTEM{i}), :));
    %divide each cause by total number in system
    prob_i = data_SYSTEM_curr/numSYSTEM;
    fprintf(fid, 'P(%s) = %f\n', cell2mat(causeSYSTEM(i)), prob_i);
end
fprintf(fid, 'Probability of causes for WARNING alarms\n');
for i=1:4,
    % !! subset your data
    % !! do the counts to derive your answers
    %find height of each cause
    data_WARNING_curr = height(data_WARNING(ismember(data_WARNING.Cause, causeWARNING{i}), :));
    %divide each cause by total number in system
    prob_i = data_WARNING_curr/numWARNING;
    fprintf(fid, 'P(%s) = %f\n', cell2mat(causeWARNING(i)), prob_i);
end
    
% T1.3.

% !! Split the data in terms of hours of the start time
% Please note that the time format is 'HH:MM:SS.FFF'
t = cell2mat(data.StartTime);
t_hours = hour(t);
fprintf(fid, '\n\nTask 1.3\n\n');
% Print the result
fprintf(fid, '\t\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t22\t23\n');
fprintf(fid, 'count\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t', sum(t_hours(:) == 0),sum(t_hours(:) == 1),sum(t_hours(:) == 2),sum(t_hours(:) == 3),sum(t_hours(:) == 4),sum(t_hours(:) == 5),sum(t_hours(:) == 6),sum(t_hours(:) == 7),sum(t_hours(:) == 8),sum(t_hours(:) == 9),sum(t_hours(:) == 10),sum(t_hours(:) == 11),sum(t_hours(:) == 12),sum(t_hours(:) == 13),sum(t_hours(:) ==14),sum(t_hours(:) == 15),sum(t_hours(:) == 16),sum(t_hours(:) == 17),sum(t_hours(:) == 18),sum(t_hours(:) == 19),sum(t_hours(:) == 20),sum(t_hours(:) == 21),sum(t_hours(:) == 22),sum(t_hours(:) == 23));

fprintf(fid, '\n');
figure;
% !! plot the histograms of alarms per hour
histogram(t_hours);
% Labeling the figure
title('Histogram of Alarms per hour');
hours = {'00h', '01h', '02h','03h','04h','05h','06h','07h','08h','09h','10h','11h','12h','13h','14h','15h','16h','17h','18h','19h','20h','21h','22h','23h'};
set(gca, 'XTick', 0:23);
set(gca,'XTickLabel',hours);
ylabel('number of alarms');

% T1.4.
% !! convert the time into units of minutes
% e.g. 01:35:01.432 is the 95th minute of a day
t_minutes = minute(t);
t_overall_minute = t_hours*60 + t_minutes;
% !! count how many minutes contain at least one alarm
minutes_w_alarm = length(unique(t_overall_minute));
fprintf(fid, '\n\nTask 1.4\n\n');
%divide total minutes in which an alarms start by 1440 (minutes in a day)
fprintf(fid, 'P(Alarm In a minute) = %f\n', minutes_w_alarm/1440);

% T2.1 
% !! Using the results from Task 1, derive the probability
 data_SYSTEM_curr = height(data_SYSTEM(ismember(data_SYSTEM.Cause,'APP_ERR'), :));
 prob_app_given_sys = data_SYSTEM_curr/numSYSTEM;
 fprintf(fid, '\n\nTask 2.1\n\n');
 %Baye's Law: P(A|B) * P(B) = P(A)
 fprintf(fid, 'P(APP_ERR and SYSTEM) = %f\n', prob_app_given_sys* numSYSTEM/numTOTAL);

% T2.2
% !! Calculate the duration for each alarms
% make sure that the durations are in seconds
f = cell2mat(data.StopTime);
%translate start time and stop time to datenum
t = datenum(t);
f = datenum(f);
%use datenum to calculate length of alarms, then convert to seconds using
%*24*3600
alarm_duration = 24*3600 * (f - t);
%add alarm duration to data table
alarm_duration = array2table(alarm_duration);
data = [data alarm_duration];
% T2.2.a 
fprintf(fid, '\n\nTask 2.2.a\n\n');
ALARMTYPE = {'SYSTEM_ALARM', 'ADVISORY_ALARM', 'WARNING_ALARM', 'CRISIS_ALARM'};
for i=1:4,
    % !! Split the data in terms of alarm type
    data_i = data(ismember(data.Alarm_Type, ALARMTYPE{i}), :);
    % !! Calculate the average duration of each alarm type
    AVERAGE_DURATION = mean(data_i.alarm_duration);
    fprintf(fid, 'Average duration for alarm type %s: %f\n', ALARMTYPE{i}, AVERAGE_DURATION);
end

% T2.2.b

fprintf(fid, '\n\nTask 2.2.b\n\n');
  t_hours_tab = hour(t);
  t_hours_tab = array2table(t_hours_tab);
  data = [data t_hours_tab];
  avg_table = {1:24};
for i=1:24,
    % Please note that i loop from 1 to 24
    % The hours in the data are from 0 to 23
    % !! Split the data in terms of hours (t_hours(:) == 0)
    data_i = data(ismember(data.t_hours_tab, (i-1)), :);
    % !! Calculate the average duration for each hour(i)
    avg_i = mean(data_i.alarm_duration);
    fprintf(fid, 'Average duration for hh=%d = %f\n', i-1, avg_i);
    avg_table{i} = avg_i;
end
figure;
y = [avg_table{1:24}];
% !! Draw a bar chart to plot the average duration per hour
bar(y);
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

%   Our group was assigned hours 4 and 5.
%   We will be checking changes in:
%       - Number of Alarms
%       - Alarm Types
%       - Alarm Causes (for System Alarms)
%       - Alarms per Bed

%   Check for number of alarms
hour_4_sum = sum(t_hours(:) == 4);
hour_5_sum = sum(t_hours(:) == 5);
fprintf(fid, 'Number of Alerts in Hour 4:\t%d\n', hour_4_sum);
fprintf(fid, 'Number of Alerts in Hour 5:\t%d\n', hour_5_sum);
fprintf(fid, 'Percent Change from Hour 4 to Hour 5:\t%f%%\n', 100*(hour_5_sum-hour_4_sum)/hour_4_sum);
fprintf(fid, '\n');

%   Check alarms by type
hour_4_alarms = data(ismember(data.t_hours_tab, 4), :);
hour_5_alarms = data(ismember(data.t_hours_tab, 5), :);
fprintf(fid, 'In Hour 4, there were:\n');
fprintf(fid, '\t => %5d  System Alarms\n', height(hour_4_alarms(ismember(hour_4_alarms.Alarm_Type, 'SYSTEM_ALARM'), :)));
fprintf(fid, '\t\t\t -> %5d  APP_ERR Instances\n', height(hour_4_alarms(ismember(hour_4_alarms.Cause, 'APP_ERR'), :)));
fprintf(fid, '\t\t\t -> %5d  SIG_ARTIFACT Instances\n', height(hour_4_alarms(ismember(hour_4_alarms.Cause, 'SIG_ARTIFACT'), :)));
fprintf(fid, '\t\t\t -> %5d  LEADS_FAILURE Instances\n', height(hour_4_alarms(ismember(hour_4_alarms.Cause, 'LEADS_FAILURE'), :)));
fprintf(fid, '\t\t\t -> %5d  NW_ERR Instances\n', height(hour_4_alarms(ismember(hour_4_alarms.Cause, 'NW_ERR'), :)));
fprintf(fid, '\t => %5d  Advisory Alarms\n', height(hour_4_alarms(ismember(hour_4_alarms.Alarm_Type, 'ADVISORY_ALARM'), :)));
fprintf(fid, '\t => %5d  Warning Alarms\n', height(hour_4_alarms(ismember(hour_4_alarms.Alarm_Type, 'WARNING_ALARM'), :)));
fprintf(fid, '\t => %5d  Crisis Alarms\n', height(hour_4_alarms(ismember(hour_4_alarms.Alarm_Type, 'CRISIS_ALARM'), :)));
fprintf(fid, 'In Hour 5, there were:\n');
fprintf(fid, '\t => %5d  System Alarms\n', height(hour_5_alarms(ismember(hour_5_alarms.Alarm_Type, 'SYSTEM_ALARM'), :)));
fprintf(fid, '\t\t\t -> %5d  APP_ERR Instances\n', height(hour_5_alarms(ismember(hour_5_alarms.Cause, 'APP_ERR'), :)));
fprintf(fid, '\t\t\t -> %5d  SIG_ARTIFACT Instances\n', height(hour_5_alarms(ismember(hour_5_alarms.Cause, 'SIG_ARTIFACT'), :)));
fprintf(fid, '\t\t\t -> %5d  LEADS_FAILURE Instances\n', height(hour_5_alarms(ismember(hour_5_alarms.Cause, 'LEADS_FAILURE'), :)));
fprintf(fid, '\t\t\t -> %5d  NW_ERR Instances\n', height(hour_5_alarms(ismember(hour_5_alarms.Cause, 'NW_ERR'), :)));
fprintf(fid, '\t => %5d  Advisory Alarms\n', height(hour_5_alarms(ismember(hour_5_alarms.Alarm_Type, 'ADVISORY_ALARM'), :)));
fprintf(fid, '\t => %5d  Warning Alarms\n', height(hour_5_alarms(ismember(hour_5_alarms.Alarm_Type, 'WARNING_ALARM'), :)));
fprintf(fid, '\t => %5d  Crisis Alarms\n', height(hour_5_alarms(ismember(hour_5_alarms.Alarm_Type, 'CRISIS_ALARM'), :)));
fprintf(fid, '\n');

%   Most Problematic Bed
hour_4_problem_bed = mode(hour_4_alarms.Bed_No);
hour_4_problem_bed_subset = hour_4_alarms(ismember(hour_4_alarms.Bed_No, hour_4_problem_bed), :);
hour_5_problem_bed = mode(hour_5_alarms.Bed_No);
hour_5_problem_bed_subset = hour_5_alarms(ismember(hour_5_alarms.Bed_No, hour_5_problem_bed), :);
fprintf(fid, 'The Most Problematic Bed during Hour 4 was Bed #%d', hour_4_problem_bed);
fprintf(fid, ' with %d alarms (%.5f%% of Hour 4''s total)\n', height(hour_4_problem_bed_subset), 100*height(hour_4_problem_bed_subset)/height(hour_4_alarms)); 
fprintf(fid, '\t => %5d  System Alarms\n', height(hour_4_problem_bed_subset(ismember(hour_4_problem_bed_subset.Alarm_Type, 'SYSTEM_ALARM'), :)));
fprintf(fid, '\t => %5d  Advisory Alarms\n', height(hour_4_problem_bed_subset(ismember(hour_4_problem_bed_subset.Alarm_Type, 'ADVISORY_ALARM'), :)));
fprintf(fid, '\t => %5d  Warning Alarms\n', height(hour_4_problem_bed_subset(ismember(hour_4_problem_bed_subset.Alarm_Type, 'WARNING_ALARM'), :)));
fprintf(fid, '\t => %5d  Crisis Alarms\n', height(hour_4_problem_bed_subset(ismember(hour_4_problem_bed_subset.Alarm_Type, 'CRISIS_ALARM'), :)));
fprintf(fid, 'The Most Problematic Bed during Hour 5 was Bed #%d', mode(hour_5_alarms.Bed_No));
fprintf(fid, ' with %d alarms (%.5f%% of Hour 5''s total)\n', height(hour_5_problem_bed_subset), 100*height(hour_5_problem_bed_subset)/height(hour_5_alarms)); 
fprintf(fid, '\t => %5d  System Alarms\n', height(hour_5_problem_bed_subset(ismember(hour_5_problem_bed_subset.Alarm_Type, 'SYSTEM_ALARM'), :)));
fprintf(fid, '\t => %5d  Advisory Alarms\n', height(hour_5_problem_bed_subset(ismember(hour_5_problem_bed_subset.Alarm_Type, 'ADVISORY_ALARM'), :)));
fprintf(fid, '\t => %5d  Warning Alarms\n', height(hour_5_problem_bed_subset(ismember(hour_5_problem_bed_subset.Alarm_Type, 'WARNING_ALARM'), :)));
fprintf(fid, '\t => %5d  Crisis Alarms\n', height(hour_5_problem_bed_subset(ismember(hour_5_problem_bed_subset.Alarm_Type, 'CRISIS_ALARM'), :)));
fprintf(fid, '\n');

% !! Write your own code for analysis (ref. codes for Task1 and Task2)

fclose(fid);

