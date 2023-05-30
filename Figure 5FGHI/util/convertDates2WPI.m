% Weeks post-surgery
function weekArr = convertDates2WPI(dates,animalIDX)

surgeryDates = {'17-Sep-2021','13-Dec-2021','03-Dec-2021','17-Dec-2021','17-Feb-2022'};
surgeryDate = surgeryDates{animalIDX};
nDATES = numel(dates);
weekArr = zeros(nDATES,1);

for i = 1:nDATES
    date_i = dates{i};
    D = datetime(date_i(1:11),'InputFormat','dd-MMM-yyyy');
    % compute day with respect to surgery date
    week_i = between(surgeryDate, D,'weeks');
    week_i = split(week_i,'WEEKS');
    weekArr(i) = week_i;
end

