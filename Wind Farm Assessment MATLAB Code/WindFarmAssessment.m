% Design and Computing Summative - Orkland Islands Wind Farm Assessment

% This code generates a wind farm assessment for the Orkland Islands wind
% farm based on the user's input of the specific site and month the user
% would like to conduct the assessment.

% This code was written and functions normally in MATLAB r2020b

% The assessment is split into 4 parts: 
% Part 1: Generating a layout of the wind farm based on the desired site
% Part 2: Plotting the wind speed throughout the farm based on the desired
%         month
% Part 3: Plotting a graph of power output throughout the farm of the farm 
%         for the desired site from the average, upper and lower quartile 
%         wind speeds
% Part 4: Plotting a graph of varying power outputs for varying wind
%         directions across a specific site during a specified month

% Imran Rizki Putranto
% 16 April 2021
% Department of Mechanical Engineering
% University of Bristol

clc
clear all

Repeat = 1; 
while Repeat == 1; % Ask user if they want to repeat for another location after the script ends
    
    % Computing Task 1 - Layout of Wind Turbine
    disp('____________________________________________________________________________');
    disp('Computing Task 1 - Generating a grid layout for turbine based on location')
    disp('____________________________________________________________________________');

    fprintf('\nThis part of the assessment generates a layout for the user-specified location');
    fprintf('\nNOTE THAT THIS LAYOUT WILL BE USED THROUGHOUT THE ASSESSMENT.')
    
    % Ask user to specify location of site
    Proceed = 1;
    while Proceed  == 1 % Loop active when user chooses KerryGold
        GEOLOC = input('\n\n     Please specify which location you would like to conduct the assessment in: \n     (Deltling, Collinfirth, Nestling, KerryGold): ','s'); 
        GEOLOC = lower(GEOLOC); % Convert input to lowercase
        
        % Loop to prevent user from entering invalid inputs
        while GEOLOC ~= "deltling" & GEOLOC ~= "collinfirth" & GEOLOC ~= "kerrygold" & GEOLOC ~= "nestling" | isempty(GEOLOC);
            disp('     Invalid Input');
            GEOLOC = input('     Please specify one of the four given locations (Deltling, Collinfirth, Nestling, KerryGold): ','s');
            GEOLOC = lower(GEOLOC);
        end
        
        % Defining number of turbines for each location
        if GEOLOC == "deltling";
                c1 = 9; % c1 gives Number of turbines in y direction
                c2 = 0:550:(c1-1)*550; % Generating TURBINE_CENTRES variable for FLORIS in y direction
                d1 = 10; % d1 gives number of turbines in x direction
                d2 = 0:770:(d1-1)*770; % Generating TURBINE_CENTRES variable for FLORIS in x direction
                Proceed  = 0; % Switch off loop due to valid location
            elseif GEOLOC == "collinfirth"
                c1 = 6;
                c2 = 0:550:(c1-1)*550;
                d1 = 6;
                d2 = 0:770:(d1-1)*770;
                Proceed  = 0;
            elseif GEOLOC == "nestling"
                c1 = 4;
                c2 = 0:550:(c1-1)*550;
                d1 = 5;
                d2 = 0:770:(d1-1)*770;
                Proceed  = 0;
            elseif GEOLOC == "kerrygold"
                fprintf('\n     This site will not have any turbines as it is near a RSBP nature reserve and contains important flora and fauna');
                Proceed  = 1;
        end
    end

    [C2,D2] = meshgrid(c2,d2); 

    % Plotting scatter graph wind farm layout
    figure(1)
    c2 = C2(:);
    d2 = D2(:);
    scatter(c2,d2, "filled");
    xlabel('Longitudinal Position Along x Direction (m)');
    ylabel('Lateral Position Along y Direction (m)');
    title('Wind Farm Layout Visualisation');


    % Computing Task 2 - Wind Speeds throughout the farm for ONE oncoming speed
    disp('________________________________________________________________________________');
    disp('Computing Task 2 - Wind Speeds throughout farm for a given oncoming wind speed. '); 
    disp('________________________________________________________________________________');
    
    fprintf('\nThis part of the assessment will generate the wind speeds across the farm for the user-specified month');
    fprintf('\nUser will be asked to input the number of the month and the script will automatically generate the required plot');
    fprintf('\nNote that it may take a few minutes to generate the plot due to the number of turbines being modelled');

    % Entering data for different months from task3data excel file
    % All wind speeds are at 10m above sea level
    MonthData = readmatrix("task3data.xlsx");
    MONTHS = MonthData(:,1);
    WIND_SPEED10 = MonthData(:,2); % Average wind speed 
    WIND_SPEED10 = WIND_SPEED10.*(1000/3600);
    Density = MonthData(:,3);
    Wind_Direction = MonthData(:,4);
    WIND_SPEEDLQ = MonthData(:,5); % Lower quartile wind speed
    WIND_SPEEDLQ = WIND_SPEEDLQ.*(1000/3600); % Converting from kph to m/s
    WIND_SPEEDUQ = MonthData(:,6);% Upper quartile wind speed
    WIND_SPEEDUQ = WIND_SPEEDUQ.*(1000/3600); % Converting from kph to m/s
    

    % Hellman Exponent and Height above sea level of location
    if GEOLOC == "deltling";
        HEIGHT = 200;  % Includes 100m hub height of turbine
        ALPHA = 0.11; % Hellman exponent
    elseif GEOLOC == "collinfirth"
        HEIGHT = 100;
        ALPHA = 0.27;
    elseif GEOLOC == "nestling"
        HEIGHT = 450;
        ALPHA = 0.06;
    elseif GEOLOC == "kerrygold"
        HEIGHT = 100;
        ALPHA = 0.4;
    end

    NewMonth = 1;
    while NewMonth == 1; % Loop to ask user if they would like to repeat PART 2 for a different month

        WIND_SPEED = WIND_SPEED10.*((HEIGHT/10)^ALPHA); % Converting wind speed at 10m to hub height

        % User input for month to determine wind speed
        MONTH = input('\n\n     Please input the number of the month in which the assessment will be carried out(i.e, January = 1): ', 's');
        MONTH = str2double(MONTH);
        % Loop to prevent user from entering invalid inputs
        while isempty(MONTH)| MONTH < 1 | MONTH > 12 | MONTH ~= floor(MONTH);
            disp('     Invalid Input.');
            MONTH = input('     Please enter an integer between 1 and 12: ', 's');
            MONTH = str2double(MONTH);
        end
        
        % User input for wind direction
        WIND_DIRECTION = input('     Please enter the desired wind direction in degrees between 0 and 360(0 for North, 90 for East, etc): ','s');
        WIND_DIRECTION = str2double(WIND_DIRECTION);
        % Loop to prevent user from entering invalid inputs
        while isempty(WIND_DIRECTION)| WIND_DIRECTION < 0 | WIND_DIRECTION > 360 | WIND_DIRECTION ~= floor(WIND_DIRECTION);
            disp('     Invalid Input.');
            WIND_DIRECTION = input('     Please enter an integer between 0 and 360: ','s');
            WIND_DIRECTION = str2double(WIND_DIRECTION);
        end

        % Choosing air density and Wind speed for different months
        for i = MONTH
            DENSITY = Density(i);
            WIND_SPEED = WIND_SPEED(i);
        end

        % Defining Cut in and Cut off wind speeds
        CUT_IN = 4;
        CUT_OFF = 25;

        % Creating array for Turbine centres from location chosen in PART 1
        xCENTRE = c2;
        yCENTRE = d2;
        zCENTRE = HEIGHT*ones(length(xCENTRE),1); % Assume constant height for all turbines
        TURBINE_CENTRES= [xCENTRE, yCENTRE, zCENTRE]; % Turbine centres vartiable for FLORIS

        YAW_ANGLES = zeros(length(xCENTRE),1); % Assume yaw angles remain 0
        DIAMETERS = 110*ones(length(xCENTRE),1);
        POWER_CURVE = readmatrix("farmpowercurve.xlsx"); % Reads power curve data from excel file
        
        % Generating location variable for FLORIS 
        xLOCATION = [-550:50:(c1)*550];
        yLOCATION = [-770:50:(d1)*770];
        [XLOCATION, YLOCATION] = meshgrid(xLOCATION, yLOCATION);

        FINALSPEED = zeros(length(yLOCATION),length(xLOCATION)); % Creating empty array for speed across farm
        
        % If statement for when Cp will be beyond the tolerable limits
        if WIND_SPEED < CUT_IN | WIND_SPEED > CUT_OFF
            disp('     WARNING!');
            disp('     Wind Speed is either below cut in or above cut off speed. ');
        else
            for i = 1:length(xLOCATION)
                for j = 1:length(yLOCATION)
                    SPEED = 0;
                    LOCATION = [XLOCATION(j,i), YLOCATION(j,i), HEIGHT];
                    [POWER,SPEED]=floris(WIND_SPEED, DENSITY, WIND_DIRECTION, TURBINE_CENTRES, YAW_ANGLES, DIAMETERS, POWER_CURVE, LOCATION);
                    FINALSPEED(j,i) = SPEED;
                end
            end
            % Plotting wind speed across the farm
            figure(2)
            pcolor(XLOCATION,YLOCATION,FINALSPEED);
            xlabel('Horizontal position along length of farm (m)');
            ylabel('Lateral position along width of farm (m)');
            title('Wind Speed Across Farm for a Given Month');
            shading interp;
            colorbar
            ylabel(colorbar, 'Wind Speed (ms^{-1})', 'Fontsize',11);
        end
        
        % Asking user if they want to repeat for another month
        disp('     Would you like to perform a new assessment for a different month?');
        NewMonth = input('     Enter 1 for YES and 0 for NO: ','s');
        NewMonth = str2double(NewMonth);
        % While loop to prevent invalid inputs
        while isempty(NewMonth)| NewMonth > 1 | NewMonth < 0 | NewMonth ~= floor(NewMonth);
            disp('     Invalid input');
            NewMonth = input('     Please enter an 1 for YES and 0 for NO: ','s');
            NewMonth = str2double(NewMonth);
        end
    end


    % Computing Task 3 - Power Output through the Year and Expected Variations
    disp('__________________________________________________________________________');
    disp('Computing Task 3 - Power output through the year and expected variations. '); 
    disp('__________________________________________________________________________');
    
    fprintf('\nThis part of the assessment will generate the power output of the farm across the year');
    fprintf('\nVariations across the lower and upper quartile speeds will also be generated');
    fprintf('\nNote that this plot is automatically generated according to the specificed location (i.e. no user input needed)');
    
    % Resetting Hellman Exponent and Height above sea level of location
    if GEOLOC == "deltling";
        HEIGHT = 200;
        ALPHA = 0.11;
    elseif GEOLOC == "collinfirth"
        HEIGHT = 100;
        ALPHA = 0.27;
    elseif GEOLOC == "nestling"
        HEIGHT = 450;
        ALPHA = 0.06;
    elseif GEOLOC == "kerrygold"
        HEIGHT = 100;
        ALPHA = 0.4;
    end
  
    %MONTHSCHAR = ["JAN"; "FEB"; "MAR"; "APR"; "MAY"; "JUN"; "JUL"; "AUG"; "SEP"; "OCT"; "NOV"; "DEC"];
  
    % Entering Variables for FLORIS
    LOCATION  = [0 ,0 ,0]; % 0 since there is no need to determine speed
    FINALPOWER1 = zeros(12,1); % Empty array for average power throughtout year
    AVG_POWER1 = zeros(13,1); % Empty array for overall average power 
    
    % Looping FLORIS to get power for each month
    for n = 1:length(WIND_SPEED10) % Average wind speed to give average power output
        WIND_SPEED = WIND_SPEED10(n,1)*((HEIGHT)/10)^ALPHA; % Convert wind speed at 10m to hub height
        if WIND_SPEED < CUT_IN | WIND_SPEED > CUT_OFF % If statement to give warning if Cp limit is reached
            disp('     WARNING');   
            str1 = "     Average wind speed for ";
            str2 = " is below the cut in speed or above the cut off speed";
            str = str1 + MONTHSCHAR(n) + str2;
            disp(str);
            disp('     Based on the data for this month, this site (on average) will not produce any power');
        else % Looping FLORIS if no error
            DENSITY = Density(n,1);
            WIND_DIRECTION = Wind_Direction(n,1);
            [POWER,SPEED]=floris(WIND_SPEED, DENSITY, WIND_DIRECTION, TURBINE_CENTRES, YAW_ANGLES, DIAMETERS, POWER_CURVE, LOCATION);
            for i = 1:length(POWER)
                if POWER(i) > 4500000;
                    POWER(i) = 4500000;
                end
            end
            SUMPOWER = sum(POWER);
            FINALPOWER1(n,1) = SUMPOWER./1000000; % Finding total average power per month in MW
        end
    end

    AVG_POWER1(13,1) = mean(FINALPOWER1); % Calculating average power for the year

    % Plotting Power throughout year for lower quartile speeds
    FINALPOWER2 = zeros(12,1); % Empty array for upper quartile power throughtout year
    AVG_POWER2 = zeros(13,1); % Empty array for  overall upper quartile power 

    for n = 1:length(WIND_SPEEDLQ)
        WIND_SPEED = WIND_SPEEDLQ(n,1)*((HEIGHT)/10)^ALPHA;
        if WIND_SPEED < CUT_IN | WIND_SPEED > CUT_OFF
            disp('     WARNING');   
            str1 = "     Average wind speed for ";
            str2 = " is below the cut in speed or above the cut off speed";
            str = str1 + MONTHSCHAR(n) + str2;
            disp(str);
            disp('     Based on the data for this month, this site (on average) will not produce any power');
        else
            DENSITY = Density(n,1);
            WIND_DIRECTION = Wind_Direction(n,1);
            [POWER,SPEED]=floris(WIND_SPEED, DENSITY, WIND_DIRECTION, TURBINE_CENTRES, YAW_ANGLES, DIAMETERS, POWER_CURVE, LOCATION);
            for i = 1:length(POWER)
                if POWER(i) > 4500000;
                    POWER(i) = 4500000;
                end
            end
            SUMPOWER = sum(POWER);
            FINALPOWER2(n,1) = SUMPOWER./1000000;
        end
    end

    AVG_POWER2(13,1) = mean(FINALPOWER2);
    
    % Plotting Power Output for Upper Quartile Speeds
    FINALPOWER3 = zeros(12,1); % Empty array for lower quartile power throughtout year
    AVG_POWER3 = zeros(13,1); % Empty array for overall lower quartile 

    for n = 1:length(WIND_SPEEDUQ)
        WIND_SPEED = WIND_SPEEDUQ(n,1)*((HEIGHT)/10)^ALPHA;
        if WIND_SPEED < CUT_IN | WIND_SPEED > CUT_OFF
            disp('     WARNING');   
            str1 = "     Average wind speed for ";
            str2 = " is below the cut in speed or above the cut off speed";
            str = str1 + MONTHSCHAR(n) + str2;
            disp(str);
            disp('     Based on the data for this month, this site (on average) will not produce any power');
        else
            DENSITY = Density(n,1);
            WIND_DIRECTION = Wind_Direction(n,1);
            [POWER,SPEED]=floris(WIND_SPEED, DENSITY, WIND_DIRECTION, TURBINE_CENTRES, YAW_ANGLES, DIAMETERS, POWER_CURVE, LOCATION);
            for i = 1:length(POWER)
                if POWER(i) > 4500000;
                    POWER(i) = 4500000;
                end
            end
            SUMPOWER = sum(POWER);
            FINALPOWER3(n,1) = SUMPOWER./1000000;
        end
    end

    AVG_POWER3(13,1) = mean(FINALPOWER3);
    
    % Plotting Power output throughout year along with variations
    FINALPOWER_ALL = [FINALPOWER1, FINALPOWER2, FINALPOWER3]; % Concatenating all powers through year 
    AVGPOWER_ALL = [AVG_POWER1, AVG_POWER2, AVG_POWER3]; % Concatenating all overall powers
    
    % Plotting data on bar chart
    figure(3);
    bar(FINALPOWER_ALL);
    title('Power Throughout Year With Expected Variations');
    hold on
    bar(AVGPOWER_ALL);
    set(gca, 'xtick', 1:13);
    set(gca, 'xticklabel', ["JAN"; "FEB"; "MAR"; "APR"; "MAY"; "JUN"; "JUL"; "AUG"; "SEP"; "OCT"; "NOV"; "DEC"; "AVG"]);
    xlabel('Month');
    ylabel('Total Power Output in Site (MW)');
    lgd = legend('Average Power','Lower Quartile Power','Upper Quartile Power', 'Annual Mean', 'Annual Lower Quartile', 'Annual Upper Quartile');
    lgd.NumColumns = 3;
    grid on
    hold off

    % Computing Task 4 - Showing effect of changing oncoming flow direction to
    % the power output
    
    fprintf('\n');
    disp('_________________________________________________________________________________________');
    disp('Computing Task 4 - Showing effect of changing oncoming flow durection to the power output');
    disp('_________________________________________________________________________________________');
    
    fprintf('\nThis part of the assessment will generate a plot of how power output varies with wind direction');
    fprintf('\nUser will be asked to specify a month for this part of the assessment to be carried out in');
    fprintf('\nAfter the plot is generated, they will then be asked if they would like to generate a new plot for a different month');
    % Resetting Hellman Exponent and Height above sea level of location
    if GEOLOC == "deltling";
        HEIGHT = 200;
        ALPHA = 0.11;
    elseif GEOLOC == "collinfirth"
        HEIGHT = 100;
        ALPHA = 0.27;
    elseif GEOLOC == "nestling"
        HEIGHT = 450;
        ALPHA = 0.06;
    elseif GEOLOC == "kerrygold"
        HEIGHT = 100;
        ALPHA = 0.4;
    end

    NewMonth = 1;
    while NewMonth == 1; % Loop to ask user if they would like to repeat PART 4 for a different month

        % Calculating wind speed at hub height of turbine at each site
        WIND_SPEED = WIND_SPEED10.*((HEIGHT/10)^ALPHA);

        MONTH = input('\n\n     Please input the number of the month in which the assessment will be carried out(i.e, January = 1): ','s');
        MONTH = str2double(MONTH);
        while isempty(MONTH) | MONTH < 1 | MONTH > 12 | MONTH ~= floor(MONTH);
            MONTH = input('     Please enter an integer between 1 and 12: ','s');
            MONTH = str2double(MONTH);
        end

        % Choosing air density and Wind speed for different months
        for i = MONTH
            DENSITY = Density(i);
            WIND_SPEED = WIND_SPEED(i,1);
        end

        % Defining the constant variables
        TURBINE_CENTRES= [xCENTRE, yCENTRE, zCENTRE];
        YAW_ANGLES = zeros(length(xCENTRE),1);
        DIAMETERS = 110*ones(length(xCENTRE),1);
        POWER_CURVE = readmatrix("farmpowercurve.xlsx"); % Reads power curve data
        LOCATION = [0,0,0];

        % Wind Direction
        WIND_DIRECTION = [0:360]; % Creating array for 360 degree wind direction
        FINALPOWER = zeros(1,length(WIND_DIRECTION)); % Creating array for power for each value of wind direction
        
        % If statement for Cp limit warning
        if WIND_SPEED < CUT_IN | WIND_SPEED > CUT_OFF;
            disp('     WARNING!');
            disp('     Wind Speed is either below cut in or above cut off speed. ');
        else % Looping FLORIS if no error
            for i = 1:length(WIND_DIRECTION);
                [POWER,SPEED]=floris(WIND_SPEED, DENSITY, WIND_DIRECTION(1,i), TURBINE_CENTRES, YAW_ANGLES, DIAMETERS, POWER_CURVE, LOCATION);
                for j = 1:length(POWER)
                    if POWER(j) > 4500000;
                        POWER(j) = 4500000;
                    end
                end
                POWERSUM = sum(POWER); % Summing up all powers from each turbine
                FINALPOWER(1,i) = POWERSUM./1000000; % Converting to MW and storing in FINALPOWER array
            end
            
            % Plotting total power against wind direction
            figure(4)
            plot(WIND_DIRECTION,FINALPOWER);
            title('Power Output for Varying Wind Direction');
            xlabel('Wind Direction In Degrees');
            ylabel('Power Output (MW)');
        end
        % Asking user if they want to perform assessment for different
        % month
        disp('     Would you like to perform a new assessment for a different month?');
        NewMonth = input('     Enter 1 for YES and 0 for NO: ','s');
        NewMonth = str2double(NewMonth);
        NewMonthnan = isnan(NewMonth);
        while NewMonth > 1 | NewMonth < 0 |  NewMonth ~= floor(NewMonth) | NewMonthnan == 1
            disp('     Invalid input');
            NewMonth = input('     Please enter an 1 for YES and 0 for NO: ','s');
            NewMonth = str2double(NewMonth);
            NewMonthnan = isnan(NewMonth);
        end
            
    end
    
    % Asking user if they would like to repeat assessment for different
    % location
    disp('     Assessment complete.');
    disp('     Would you like to perform another assessment? ');
    Repeat = input('     Enter 1 for YES and 0 for NO: ','s');
    Repeat = str2double(Repeat);
    Repeatnan = isnan(Repeat);
    while isempty(Repeat) | Repeat > 1 | Repeat < 0 | Repeatnan == 1
        disp('     Invalid input');
        Repeat = input('     Please enter an 1 for YES and 0 for NO: ','s');
        Repeat = str2double(Repeat);
        Repeatnan = isnan(Repeat);
    end
end