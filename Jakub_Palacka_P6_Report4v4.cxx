#include <iostream>
#include <string>
#include <stdio.h>
//#include <unistd.h>
#include <math.h>
using namespace std;
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include <iomanip>
#include <fstream>
#include "new_gnuplot.cxx"

const int x_length = 202;	//length of container in x-direction including 2 extra columns for entropy calculations
const int y_length = 202;	//length of container in y-direction including 2 extra rows for entropy calculations
const int max_x = 201;		//maximum value of x for x-coordinates
const int min_x = 1;		//minimun value of x for x-coordinates
const int max_y = 201;		//maximun value of y for y-coordinates
const int min_y = 1;		//minimun value of y for y-coordinates
const int x_origin = 101;	//origin x-coordinate
const int y_origin = 101;	//origin y-coordinate
const double pi = 3.14159;	//approximation for pi
const int number_particles = 400;	//total number of particles
//const int number_steps = 10000;	//number of steps to be taken

int main()
{
	int x_coord[number_particles];	//array holding x-coordinate for each particle
	int y_coord[number_particles];	//array holding y-coordinate for each particle
	double x[number_particles];		//array holding x-coordinate for each particle as a double
	double y[number_particles];		//array holding y-coordinate for each particle as a double
	double(*site)[x_length] = new double[x_length][y_length];				//pointer declaration for site grid
	double(*probability)[x_length] = new double[x_length][y_length];		//pointer declaration for probability values
	double(*probability_sum)[x_length] = new double[x_length][y_length];	//pointer declaration for Pi*lnPi

	double x_range[2] = { min_x, max_x*1.0};	//setting x-axis range for graphs of particle distribution
	double y_range[2] = { min_y, max_y*1.0 };	//setting y-axis range for graphs of particle distribution

	double* timestep, * meansqr, * rms, * entropy_t, * entropy_prime, * sqrt_t, * drop_size_RMS, *drop_size_disp, *disp;	//setting pointers for large arrays
	double meansqr_displacement[number_particles];	//array holding mean-square displacement values for each particle
	double radius[number_particles];	//array holding values of displacement of each particle, will be transferred to disp later

	double mean_square, root_mean_square, entropy;	//variables used in finding total mean-square displacement, RMS displacement and entropy of the system at each step k
	int i, j, k, l, step_x, x_previous, y_previous, dt, min_steps, number_steps, input;	//variables used for loops, container barriers and calculation of rate of change in entropy

	//initialising variables to 0
	mean_square = 0;
	root_mean_square = 0;
	entropy = 0;

	cout << "Please enter the number of steps you would like the program to carry out (must be greater than 500). \n";
	cout << "The program prints graphs at step numbers 200, 500, 1000, 3000, 7000, 30000, 50000 and 100000.\n";
	cout << "At each of these intervals 7 graphs are printed. \n";

	cout << "\n";

	cout << "Please be aware that if the number of steps to be executed is 100000 or greater then a total of 56 graphs will be printed, this is intentional and not the result of a bug. \n";
	cout << "Number of steps: ";
	cin >> number_steps;

	cout << "\n";

	cout << "The program will now run for " << number_steps << " steps.\n";
	cout << "If you would like the program to continue please enter the number 1, if you would like to exit please enter the number 2.\n";
	cout << "Input: ";
	cin >> input;

	timestep = new double[number_steps];	//array holding timestep values
	entropy_t = new double[number_steps];	//array holding entropy over time values
	sqrt_t = new double[number_steps];		//array holding sqrt(t) values
	meansqr = new double[number_steps];		//array holding mean-square displacement values for each step k
	rms = new double[number_steps];			//array holding rms displacement values
	disp = new double[number_steps];
	drop_size_RMS = new double[number_steps];	//aray holding drop size values for r = RMS
	drop_size_disp = new double[number_steps];	//array holding drop size values for r = max. displacement
	entropy_prime = new double[number_steps];	//array holding values of rate of change of entropy

	if (input == 1)
	{

	
	srand(time(NULL));	//generate rand seed for random number generation
	
	//printing values to .csv file
	ofstream outdata;


	//fill arrays with zeros
	fill(site[0], site[0] + x_length * y_length, 0);
	fill(probability[0], probability[0] + x_length * y_length, 0);
	fill(probability_sum[0], probability_sum[0] + x_length * y_length, 0);

	//initialise particles in 5x5 square about origin (101, 101)
	//first row
	l = 0;
	for (j = 0; j < 5; j++)
	{
		for (i = 16 * j; i < 16 * (j + 1); i++)
		{
			x_coord[i] = (x_origin -2) + l;
			y_coord[i] = y_origin -2;
			x[i] = x_coord[i];
			y[i] = y_coord[i];
		}
		l = l + 1;
	}

	//second row
	l = 0;
	for (j = 5; j < 10; j++)
	{
		for (i = 16 * j; i < 16 * (j + 1); i++)
		{
			x_coord[i] = (x_origin - 2) + l;
			y_coord[i] = y_origin - 1;
			x[i] = x_coord[i];
			y[i] = y_coord[i];
		}
		l = l + 1;
	}

	//third row
	l = 0;
	for (j = 10; j < 15; j++)
	{
		for (i = 16 * j; i < 16 * (j + 1); i++)
		{
			x_coord[i] = (x_origin - 2) + l;
			y_coord[i] = y_origin;
			x[i] = x_coord[i];
			y[i] = y_coord[i];
		}
		l = l + 1;
	}

	//fourth row
	l = 0;
	for (j = 15; j < 20; j++)
	{
		for (i = 16 * j; i < 16 * (j + 1); i++)
		{
			x_coord[i] = (x_origin - 2) + l;
			y_coord[i] = y_origin + 1;
			x[i] = x_coord[i];
			y[i] = y_coord[i];
		}
		l = l + 1;
	}

	//fifth row
	l = 0;
	for (j = 20; j < 25; j++)
	{
		for (i = 16 * j; i < 16 * (j + 1); i++)
		{
			x_coord[i] = (x_origin - 2) + l;
			y_coord[i] = y_origin + 2;
			x[i] = x_coord[i];
			y[i] = y_coord[i];
		}
		l = l + 1;
	}

	gnuplot_one_function_square_jpg("Particle Positions - 0 steps", "points", "X", "Y", x, y, number_particles, "Positions_0.pdf", x_range, y_range);

	//Initialise occupied sites
	for (i = 0; i < number_particles; i++)
	{
		site[x_coord[i]][y_coord[i]] = 1; //if a particle is at coordinates (i, j) make site[i][j] = 1

	}


	//calculation of initial entropy of system
	for (i = min_x; i < max_x; i++)
	{
		for (j = min_y; j < max_y; j++)
		{
			probability[i][j] = 1 / 8.0 * (site[i][j - 1] + site[i + 1][j] + site[i][j + 1] + site[i - 1][j] + site[i + 1][j + 1] + site[i - 1][j + 1] + site[i + 1][j - 1] + site[i - 1][j - 1]); //probability of a site being occupied

			//eliminating undefined values
			if (probability[i][j] != 0)
			{
				probability_sum[i][j] = (-1) * (probability[i][j] * log(probability[i][j]));

			}
			else
			{
				probability_sum[i][j] = 0;
			}
			entropy = entropy + probability_sum[i][j];
		}
	}

	entropy_t[0] = entropy;	//initial entropy of the system
	entropy = 0;	//resetting variable for later use

	timestep[0] = 0;	//initial time value of system


	//calcution of initial mean-sqaure displacement
	for (i = 0; i < number_particles; i++)
	{
		meansqr_displacement[i] = ((1.0 * (x_coord[i]) - (x_origin*1.0)) * (1.0 * (x_coord[i]) - (x_origin*1.0)) + (1.0 * (y_coord[i]) - (y_origin*1.0)) * (1.0 * (y_coord[i]) - (y_origin*1.0))) / (number_particles);		//mean square displacement per particle - [(x2-x1)^2 + (y2-y1)^2]/n																																	//root-mean-square displacement of each particle from ( 100, 100 ) - origin		rms displacement per particle
		mean_square = mean_square + meansqr_displacement[i];	//total mean square displacement

	}

	//calculation of initial rms displacement
	root_mean_square = sqrt(mean_square);

	//calculating displacement of each particle and finding maximum
	for (i = 0; i < number_particles; i++)
	{
		radius[i] = sqrt((1.0 * (x_coord[i]) - (x_origin * 1.0)) * (1.0 * (x_coord[i]) - (x_origin * 1.0)) + (1.0 * (y_coord[i]) - (y_origin * 1.0)) * (1.0 * (y_coord[i]) - (y_origin * 1.0)));	//calculation of displacement

		//finding maximum initial displacement
		if (radius[0] < radius[i])
		{
			radius[0] = radius[i];
		}

		else
		{
			radius[0] = radius[0];
		}
	}

	disp[0] = radius[0];			//maximum displacement
	meansqr[0] = mean_square;		//mean-square displacement at t = 0
	rms[0] = root_mean_square;		//RMS displacement at t = 0
	timestep[0] = 0;				//initial value of t
	root_mean_square = 0;			//reseting root_mean_square variable to 0 for next use
	mean_square = 0;				//reseting mean_sqaure variable to 0 for next use

	//calculation of initial drop size (r = max. displacment)
	drop_size_disp[0] = pi * disp[0] * disp[0];

	//calculation of initial drop size  (r = RMS)
	drop_size_RMS[0] = pi * rms[0] * rms[0];



	//main position change loop
	for (k = 1; k < number_steps; k++)
	{
		cout << "Step: " << k << "\n";

		//reset all occupied arrays for new inputs
		fill(site[0], site[0] + x_length * y_length, 0);
		fill(probability[0], probability[0] + x_length * y_length, 0);
		fill(probability_sum[0], probability_sum[0] + x_length * y_length, 0);



		//changing position of 400 particles
		for (i = 0; i < number_particles; i++)
		{
			x_previous = x_coord[i];
			y_previous = y_coord[i];

			//Generate random number between 1 and 8 for x-coordinate and y-coordinate
			step_x = rand() % 8 + 1;

			//If random number = 1 then x + 1
			if (step_x == 1)
			{
				x_coord[i] = x_coord[i] + 1;

			}

			//If random number = 2 then x - 1)
			else if (step_x == 2)
			{
				x_coord[i] = x_coord[i] - 1;

			}


			//If random number = 3 then y + 1
			else if (step_x == 3)
			{
				y_coord[i] = y_coord[i] + 1;

			}

			//If random number = 4 then y - 1
			else if (step_x == 4)
			{
				y_coord[i] = y_coord[i] - 1;

			}

			//If random number = 5 then x + 1 y + 1
			else if (step_x == 5)
			{
				x_coord[i] = x_coord[i] + 1;
				y_coord[i] = y_coord[i] + 1;
			}

			//If random number = 6 then x + 1 y -1
			else if (step_x == 6)
			{
				x_coord[i] = x_coord[i] + 1;
				y_coord[i] = y_coord[i] - 1;
			}

			//If random number = 7 then x - 1 y + 1
			else if (step_x == 7)
			{
				x_coord[i] = x_coord[i] - 1;
				y_coord[i] = y_coord[i] + 1;
			}

			//If random number = 8 then x - 1 y - 1
			else
			{
				x_coord[i] = x_coord[i] - 1;
				y_coord[i] = y_coord[i] - 1;
			}

			//if change in position results in particle outside the container return particle to previous position
			if (x_coord[i] > max_x || x_coord[i] < min_y)
			{
				x_coord[i] = x_previous;
				y_coord[i] = y_previous;
			}

			//if change in position results in particle outside the container return particle to previous position
			else if (y_coord[i] > max_y || y_coord[i] < min_y)
			{
				x_coord[i] = x_previous;
				y_coord[i] = y_previous;
			}

			y[i] = y_coord[i];		//converting integer x-coordinates to double to use in gnuplot
			x[i] = x_coord[i];		//converting integer y-coordinates to double to use in gnuplot

			//set occupied lattice locations to 1
			site[x_coord[i]][y_coord[i]] = 1;

		}

		//calculation of entropy
		for (i = min_x; i < max_x; i++)
		{
			for (j = min_y; j < max_y; j++)
			{
				probability[i][j] = 1 / 8.0 * (site[i][j - 1] + site[i + 1][j] + site[i][j + 1] + site[i - 1][j] + site[i + 1][j + 1] + site[i - 1][j + 1] + site[i + 1][j - 1] + site[i - 1][j - 1]);	//probability of a site being occupied

				//eliminating undefined values
				if (probability[i][j] != 0)
				{
					probability_sum[i][j] = (-1) * (probability[i][j] * log(probability[i][j]));

				}
				else
				{
					probability_sum[i][j] = 0;
				}
				entropy = entropy + probability_sum[i][j];
			}
		}

		entropy_t[k] = entropy;	//saving entropy of system at step k
		entropy = 0;	//resetting variable to 0 for later use



		//calcution of mean-sqaure displacement
		for (i = 0; i < number_particles; i++)
		{
			meansqr_displacement[i] = ((1.0 * (x_coord[i]) - (x_origin * 1.0)) * (1.0 * (x_coord[i]) - (x_origin * 1.0)) + (1.0 * (y_coord[i]) - (y_origin * 1.0)) * (1.0 * (y_coord[i]) - (y_origin * 1.0))) / (number_particles);		//mean square displacement per particle - [(x2-x1)^2 + (y2-y1)^2]/n																																	//root-mean-square displacement of each particle from ( 100, 100 ) - origin		rms displacement per particle
			mean_square = mean_square + meansqr_displacement[i];	//total mean square displacement

		}


		//calculation of rms displacement
		root_mean_square = sqrt(mean_square);

		meansqr[k] = mean_square; //mean-square displacement at step k
		rms[k] = root_mean_square;	//RMS displacement at step k
		timestep[k] = 1;
		root_mean_square = 0;	//resetting variable to 0 for later use
		mean_square = 0;		//resetting variable to 0 for later use



		timestep[k] = k;	//saving current step number
		sqrt_t[k] = sqrt(k);	//calculating sqrt(time)

		//calculation of maximum displacement from origin at step k
		for (i = 0; i < number_particles; i++)
		{
			radius[i] = sqrt((1.0 * (x_coord[i]) - (x_origin * 1.0)) * (1.0 * (x_coord[i]) - (x_origin * 1.0)) + (1.0 * (y_coord[i]) - (y_origin * 1.0)) * (1.0 * (y_coord[i]) - (y_origin * 1.0)));	//calculating displacement

			//finding max. displacement
			if (radius[0] < radius[i])
			{
				radius[0] = radius[i];
			}

			else
			{
				radius[0] = radius[0];
			}
		}

		disp[k] = radius[0];	//saving max. displacment at step k

		//calculation of drop size at step k (r = displacement)
		drop_size_disp[k] = pi * disp[k] * disp[k];

		//calculation of drop size at step k (r = RMS)
		drop_size_RMS[k] = pi * rms[k] * rms[k];

		
				//printing graphs at 200 steps
				if (k == 199)
				{
					gnuplot_one_function_jpg("Change in Entropy over Time - 200 steps", "linespoints", "Steps", "Entropy", timestep, entropy_t, number_steps, "Entropy_Versus_Time_200.jpg");
					gnuplot_one_function_jpg("Mean Square Displacement - 200 steps", "linespoints", "Steps", "Mean Square Displacement", timestep, meansqr, number_steps, "MSQ_Versus_Time_200.jpg");
					gnuplot_one_function_jpg("RMS Dsiplacement - 200 steps", "linespoints", "Steps", "RMS Dsiplacement", timestep, rms, number_steps, "RMS_Versus_Time_200.jpg");
					gnuplot_one_function_jpg("Dsiplacement - 200 steps", "linespoints", "Steps", "Dsiplacement", timestep, disp, number_steps, "Displacement_Versus_Time_200.jpg");
					gnuplot_one_function_square_jpg("Particle Positions - 200 steps", "points", "X", "Y", x, y, number_particles, "Positions_200.pdf", x_range, y_range);
					gnuplot_one_function_jpg("Drop Size (r = RMS) - 200 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_RMS, number_steps, "Drop_Size(RMS)_Versus_sqrt(Time)_200.jpg");
					gnuplot_one_function_jpg("Drop Size (r = Displacement)- 200 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_disp, number_steps, "Drop_Size(disp)_Versus_sqrt(Time)_200.jpg");
				}

				//printing graphs at 500 steps
				if (k == 499)
				{
					gnuplot_one_function_jpg("Change in Entropy over Time - 500 steps", "linespoints", "Steps", "Entropy", timestep, entropy_t, number_steps, "Entropy_Versus_Time_500.jpg");
					gnuplot_one_function_jpg("Mean Square Displacement - 500 steps", "linespoints", "Steps", "Mean Square Displacement", timestep, meansqr, number_steps, "MSQ_Versus_Time_500.jpg");
					gnuplot_one_function_jpg("RMS Dsiplacement - 500 steps", "linespoints", "Steps", "RMS Dsiplacement", timestep, rms, number_steps, "RMS_Versus_Time_500.jpg");
					gnuplot_one_function_jpg("Dsiplacement - 500 steps", "linespoints", "Steps", "Dsiplacement", timestep, disp, number_steps, "Displacement_Versus_Time_500.jpg");
					gnuplot_one_function_square_jpg("Particle Positions - 500 steps", "points", "X", "Y", x, y, number_particles, "Positions_500.pdf", x_range, y_range);
					gnuplot_one_function_jpg("Drop Size (r = RMS) - 500 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_RMS, number_steps, "Drop_Size(RMS)_Versus_sqrt(Time)_500.jpg");
					gnuplot_one_function_jpg("Drop Size (r = Displacement) - 500 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_disp, number_steps, "Drop_Size(disp)_Versus_sqrt(Time)_500.jpg");
				}

				//printing graphs at 1000 steps
				if (k == 999)
				{
					gnuplot_one_function_jpg("Change in Entropy over Time - 1000 steps", "linespoints", "Steps", "Entropy", timestep, entropy_t, number_steps, "Entropy_Versus_Time_1000.jpg");
					gnuplot_one_function_jpg("Mean Square Displacement - 1000 steps", "linespoints", "Steps", "Mean Square Displacement", timestep, meansqr, number_steps, "MSQ_Versus_Time_1000.jpg");
					gnuplot_one_function_jpg("RMS Dsiplacement - 1000 steps", "linespoints", "Steps", "RMS Dsiplacement", timestep, rms, number_steps, "RMS_Versus_Time_1000.jpg");
					gnuplot_one_function_jpg("Dsiplacement - 1000 steps", "linespoints", "Steps", "Dsiplacement", timestep, disp, number_steps, "Displacement_Versus_Time_1000.jpg");
					gnuplot_one_function_square_jpg("Particle Positions - 1000 steps", "points", "X", "Y", x, y, number_particles, "Positions_1000.pdf", x_range, y_range);
					gnuplot_one_function_jpg("Drop Size (r = RMS) - 1000 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_RMS, number_steps, "Drop_Size(RMS)_Versus_sqrt(Time)_1000.jpg");
					gnuplot_one_function_jpg("Drop Size (r = Displacement) - 1000 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_disp, number_steps, "Drop_Size(disp)_Versus_sqrt(Time)_1000.jpg");
				}

				//printing graphs at 3000 steps
				if (k == 2999)
				{
					gnuplot_one_function_jpg("Change in Entropy over Time - 3000 steps", "linespoints", "Steps", "Entropy", timestep, entropy_t, number_steps, "Entropy_Versus_Time_3000.jpg");
					gnuplot_one_function_jpg("Mean Square Displacement - 3000 steps", "linespoints", "Steps", "Mean Square Displacement", timestep, meansqr, number_steps, "MSQ_Versus_Time_3000.jpg");
					gnuplot_one_function_jpg("RMS Dsiplacement - 3000 steps", "linespoints", "Steps", "RMS Dsiplacement", timestep, rms, number_steps, "RMS_Versus_Time_3000.jpg");
					gnuplot_one_function_jpg("Dsiplacement - 3000 steps", "linespoints", "Steps", "Dsiplacement", timestep, disp, number_steps, "Displacement_Versus_Time_3000.jpg");
					gnuplot_one_function_square_jpg("Particle Positions - 3000 steps", "points", "X", "Y", x, y, number_particles, "Positions_3000.pdf", x_range, y_range);
					gnuplot_one_function_jpg("Drop Size (r = RMS) - 3000 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_RMS, number_steps, "Drop_Size(RMS)_Versus_sqrt(Time)_3000.jpg");
					gnuplot_one_function_jpg("Drop Size (r = Displacement) - 3000 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_disp, number_steps, "Drop_Size(disp)_Versus_sqrt(Time)_3000.jpg");
				}

				
				//printing graphs at 7000 steps
				if (k == 6999)
				{
					gnuplot_one_function_jpg("Change in Entropy over Time - 7000 steps", "linespoints", "Steps", "Entropy", timestep, entropy_t, number_steps, "Entropy_Versus_Time_7000.jpg");
					gnuplot_one_function_jpg("Mean Square Displacement - 7000 steps", "linespoints", "Steps", "Mean Square Displacement", timestep, meansqr, number_steps, "MSQ_Versus_Time_7000.jpg");
					gnuplot_one_function_jpg("RMS Dsiplacement - 7000 steps", "linespoints", "Steps", "RMS Dsiplacement", timestep, rms, number_steps, "RMS_Versus_Time_7000.jpg");
					gnuplot_one_function_jpg("Dsiplacement - 7000 steps", "linespoints", "Steps", "Dsiplacement", timestep, disp, number_steps, "Displacement_Versus_Time_7000.jpg");
					gnuplot_one_function_square_jpg("Particle Positions - 7000 steps", "points", "X", "Y", x, y, number_particles, "Positions_7000.pdf", x_range, y_range);
					gnuplot_one_function_jpg("Drop Size (r = RMS) - 7000 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_RMS, number_steps, "Drop_Size(RMS)_Versus_sqrt(Time)_7000.jpg");
					gnuplot_one_function_jpg("Drop Size (r = Displacement) - 7000 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_disp, number_steps, "Drop_Size(disp)_Versus_sqrt(Time)_7000.jpg");
				}
		
				//printing graphs at 30000 steps
				if (k == 29999)
				{
					gnuplot_one_function_jpg("Change in Entropy over Time - 30000 steps", "linespoints", "Steps", "Entropy", timestep, entropy_t, number_steps, "Entropy_Versus_Time_30000.jpg");
					gnuplot_one_function_jpg("Mean Square Displacement - 30000 steps", "linespoints", "Steps", "Mean Square Displacement", timestep, meansqr, number_steps, "MSQ_Versus_Time_30000.jpg");
					gnuplot_one_function_jpg("RMS Dsiplacement - 30000 steps", "linespoints", "Steps", "RMS Dsiplacement", timestep, rms, number_steps, "RMS_Versus_Time_30000.jpg");
					gnuplot_one_function_jpg("Dsiplacement - 30000 steps", "linespoints", "Steps", "Dsiplacement", timestep, disp, number_steps, "Displacement_Versus_Time_30000.jpg");
					gnuplot_one_function_square_jpg("Particle Positions - 30000 steps", "points", "X", "Y", x, y, number_particles, "Positions_30000.pdf", x_range, y_range);
					gnuplot_one_function_jpg("Drop Size (r = RMS) - 30000 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_RMS, number_steps, "Drop_Size(RMS)_Versus_sqrt(Time)_30000.jpg");
					gnuplot_one_function_jpg("Drop Size (r = Displacement) - 30000 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_disp, number_steps, "Drop_Size(disp)_Versus_sqrt(Time)_30000.jpg");
				}


				//printing graphs at 50000 steps
				if (k == 49999)
				{
					gnuplot_one_function_jpg("Change in Entropy over Time - 50000 steps", "linespoints", "Steps", "Entropy", timestep, entropy_t, number_steps, "Entropy_Versus_Time_50000.jpg");
					gnuplot_one_function_jpg("Mean Square Displacement - 50000 steps", "linespoints", "Steps", "Mean Square Displacement", timestep, meansqr, number_steps, "MSQ_Versus_Time_50000.jpg");
					gnuplot_one_function_jpg("RMS Dsiplacement - 50000 steps", "linespoints", "Steps", "RMS Dsiplacement", timestep, rms, number_steps, "RMS_Versus_Time_50000.jpg");
					gnuplot_one_function_jpg("Dsiplacement - 50000 steps", "linespoints", "Steps", "Dsiplacement", timestep, disp, number_steps, "Displacement_Versus_Time_50000.jpg");
					gnuplot_one_function_square_jpg("Particle Positions - 50000 steps", "points", "X", "Y", x, y, number_particles, "Positions_50000.pdf", x_range, y_range);
					gnuplot_one_function_jpg("Drop Size (r = RMS) - 50000 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_RMS, number_steps, "Drop_Size(RMS)_Versus_sqrt(Time)_50000.jpg");
					gnuplot_one_function_jpg("Drop Size (r = Displacement) - 50000 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_disp, number_steps, "Drop_Size(disp)_Versus_sqrt(Time)_50000.jpg");
				}

				if (k == 99999)
				{
					gnuplot_one_function_jpg("Change in Entropy over Time - 100000 steps", "linespoints", "Steps", "Entropy", timestep, entropy_t, number_steps, "Entropy_Versus_Time_100000.jpg");
					gnuplot_one_function_jpg("Mean Square Displacement - 100000 steps", "linespoints", "Steps", "Mean Square Displacement", timestep, meansqr, number_steps, "MSQ_Versus_Time_100000.jpg");
					gnuplot_one_function_jpg("RMS Dsiplacement - 100000 steps", "linespoints", "Steps", "RMS Dsiplacement", timestep, rms, number_steps, "RMS_Versus_Time_100000.jpg");
					gnuplot_one_function_jpg("Dsiplacement - 100000 steps", "linespoints", "Steps", "Dsiplacement", timestep, disp, number_steps, "Displacement_Versus_Time_100000.jpg");
					gnuplot_one_function_square_jpg("Particle Positions - 100000 steps", "points", "X", "Y", x, y, number_particles, "Positions_100000.pdf", x_range, y_range);
					gnuplot_one_function_jpg("Drop Size (r = RMS) - 100000 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_RMS, number_steps, "Drop_Size(RMS)_Versus_sqrt(Time)_100000.jpg");
					gnuplot_one_function_jpg("Drop Size (r = Displacement) - 100000 steps", "linespoints", "Steps", "Drop Size", timestep, drop_size_disp, number_steps, "Drop_Size(disp)_Versus_sqrt(Time)_100000.jpg");
				}

	}



	
/**/
	dt = 22;
	for ( i = dt; i < number_steps - dt; i++)
	{
		entropy_prime[i] = (entropy_t[i - dt] - entropy_t[i + dt]) / (2.0 * dt);
	}
	
	min_steps = 0;
	
	for ( i = dt; i < number_steps - dt; i++)
	{
		if (entropy_prime[dt] < entropy_prime[i])
		{
			entropy_prime[dt] = entropy_prime[i];
			min_steps = i;
		}

		else
		{
			entropy_prime[dt] = entropy_prime[dt];
			min_steps = min_steps;
		}
	}

	outdata.open("output_entropy.csv", ios::app);
	outdata << "Rate of Change of Entropy" << "," << "Side Length" << "," << "No. of Steps" << endl;
	outdata << entropy_prime[20] << "," << x_length << "," << min_steps << endl;

	//gnuplot_one_function_jpg("Square-root t versus t", "linespoints", "Steps", "Sqrt(t)", timestep, sqrt_t, number_steps, "Sqrt_t_Versus_Time.jpg");
	
	cout << "Graphs have been printed in the folder containing the .exe of this program. \n";
	
	system("pause");
	
	return 0;
}

else if (input == 2 )
	{
		return 0;

	}
	
	return 0;

}
