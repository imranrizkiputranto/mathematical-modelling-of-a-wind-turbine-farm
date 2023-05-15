Design and Computing Summative - Orkland Islands Wind Farm Assessment Tool
Author: Imran Rizki Putranto

This code generates a wind farm assessment for the Orkland Islands wind
farm based on the user's input of the specific site and month the user
would like to conduct the assessment.

The assessment is split into 4 parts: 
	Part 1: Generating a layout of the wind farm based on the desired site
	Part 2: Plotting the wind speed throughout the farm based on the desired month

	Part 3: Plotting a graph of power output throughout the farm of the farm 
        	for the desired site from the average, upper and lower quartile 
        	wind speeds
	Part 4: Plotting a graph of varying power outputs for varying wind
         	directions across a specific site during a specified month


Guidelines:
PART 1: Part 1 of the code produces a grid layout of the wind farm for the user-specified location.
	This layout (number of turbines) and location information will be used throughout the rest of the assessment. 
	The user can choose one of three location (Deltling, Collinfirth and Nestling; KerryGold is not used due to important flora and fauna in the area)
	The best possible layout for each location has been input into the code and users only need to specify which location they would like.
	The maximum usable area has already been accounted for. (Deltling: 3200 hectares, Collinfirth: 1100 hectares, Nestling: 2100 hectares)


PART 2: Part 2 of the code produces a visualisation of the wind speeds across the wind farm for one given oncoming wind speed and direction
	based on the previously chosen location.
	The user will be asked to input the month they would like to conduct the assessment in (to get wind speed data from "task3data.xlsx" file
	and the user can also input any wind direction from 0 to 360 degrees (0 being North, 90 being East, etc.)
	A coloured plot will then be generated, with the x axis and y axis showing the distances along and across the farm and the colour representing
	the speeds at that point.
	The user can then choose to repeat the assessment for a different month/direction if they would like to.

PART 3: Part 3 of the code automatically produces a graph of the power outputs throughout the year for each month, and the overall power output for the year,
	for the mean, 25th percentile and 75th percentile wind speeds for the previously specified location. 
	No user input is required for this part.

Part 4: Part 4 of the code shows how power output varies with a varying wind direction from 0 to 360 degrees for a given month.
	The user will be asked to input the month and the code automatically generates the plot.
	The user can then choose to repeat the assessment for a different month if they would like to.
	
Repetition: After PART 4 of the code is executed, the user will then be asked if they would like to repeat the assessment for a new location if they would like,
	    which repeats the whole code from the start, unless the user decides to terminate the script.