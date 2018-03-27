/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 *  Updated on: Mar 26, 2018
 *   Co-Author: Aldo Arriaga
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	//define the number of particles
	num_particles = 100;
	//particles = std::vector<Particle>(num_particles);
	default_random_engine gen;

	//create Gaussian distributions for x,y and theta based on input parameters for sampling with noise
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_q(theta, std[2]);

	//initialize each particle by randomly sampling the corresponding Gaussian distribution
	for (int i = 0; i < num_particles; ++i) {
		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_q(gen);
		particle.weight = 1.0f/num_particles;
		particles.push_back(particle);
		weights.push_back(1/num_particles);
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine

	default_random_engine gen;

	double v = velocity;
	double yawd = yaw_rate;

	for (int i = 0; i < num_particles; i++) {
		//extract i_th particle's state variables
		double x = particles[i].x;
		double y = particles[i].y;
		double yaw = particles[i].theta;
		double x_p, y_p, yaw_p;

		//intermediate variables to avoid calling trig functions multiple times
		double cos_yaw = cos(yaw);
		double sin_yaw = sin(yaw);

		//CTRV model
		yaw_p = yaw + yawd * delta_t;

		//check division by zero
		if (fabs(yawd)>0.0001) {
			x_p = x + (v / yawd) * (sin(yaw + yawd * delta_t) - sin_yaw);
			y_p = y + (v / yawd) * (cos_yaw - cos(yaw + yawd * delta_t));
		}
		else {
			x_p = x + v * cos_yaw * delta_t;
			y_p = y + v * sin_yaw * delta_t;
		}

		normal_distribution<double> noise_x(x_p, std_pos[0]);
		normal_distribution<double> noise_y(y_p, std_pos[1]);
		normal_distribution<double> noise_q(yaw_p, std_pos[2]);

		particles[i].x = noise_x(gen);
		particles[i].y = noise_y(gen);
		particles[i].theta = noise_q(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (int i = 0; i < observations.size(); i++) {

		double nearest_d = 10000; // initializing with a large value

		for (int j = 0; j < predicted.size(); j++) {

			double x_diff = observations[i].x - predicted[j].x;
			double y_diff = observations[i].y - predicted[j].y;
			double magnitude = sqrt((x_diff*x_diff) + (y_diff*y_diff));	

			if (magnitude < nearest_d) {
				nearest_d = magnitude;
				observations[i].id = predicted[j].id;	
			}		
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	double sigma_x = std_landmark[0];
	double sigma_y = std_landmark[1];

	double norm_constant = (1 / (2 * M_PI * sigma_x * sigma_y));
	double sum_weights = 0;

	// for each particle in the filter
	for (int i = 0; i < num_particles; i++) {
		// weight
		double weight_p=1.0;

		// step 1. Check for landmarks within sensor range and put them in a LandmarkObs vector "predicted_lm"
		std::vector<LandmarkObs> predicted_lm;

		for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
			
			double x_diff = particles[i].x - map_landmarks.landmark_list[j].x_f;
			double y_diff = particles[i].y - map_landmarks.landmark_list[j].y_f;
			double magnitude = sqrt((x_diff*x_diff) + (y_diff*y_diff));

			if (magnitude < sensor_range) {
				LandmarkObs predicted;
				predicted.id = map_landmarks.landmark_list[j].id_i;
				predicted.x = map_landmarks.landmark_list[j].x_f;
				predicted.y = map_landmarks.landmark_list[j].y_f;
				predicted_lm.push_back(predicted);
				
			}
		}

		// step 2. transform observation to map refefence frame and put them in a new vector "transformed_meas"

		std::vector < LandmarkObs> transformed_meas;
		for (int j = 0; j < observations.size(); j++) {
			LandmarkObs transformed;
			double sin_q = sin(particles[i].theta);
			double cos_q = cos(particles[i].theta);

			transformed.x = particles[i].x + cos_q * observations[j].x - sin_q * observations[j].y;
			transformed.y = particles[i].y + sin_q * observations[j].x + cos_q * observations[j].y;
			transformed_meas.push_back(transformed);
		}

		// step 3. Associate observations (measurements) to landmarks within sensor range
		dataAssociation(predicted_lm, transformed_meas);

		//step 4. Compute probabilities for each transformed measurement
		for (int j = 0; j<transformed_meas.size();j++){

			//find landmark associated to measurement
			LandmarkObs nearest_lm;
			for (int k = 0; k < predicted_lm.size(); k++) {
				if (predicted_lm[k].id == transformed_meas[j].id) {
					nearest_lm = predicted_lm[k];
					
					break;
				}
			}
			//compute Gaussian
			double x_diff = transformed_meas[j].x - nearest_lm.x;
			double y_diff = transformed_meas[j].y - nearest_lm.y;
			double weight = exp(-((x_diff*x_diff / (2 * sigma_x*sigma_x)) + (y_diff*y_diff / (2 * sigma_y*sigma_y))));
			
			weight_p *= weight;
		}

		//multiply weight by Gaussian normalization constant
		weight_p *= norm_constant;
		sum_weights += weight_p;
		particles[i].weight = weight_p;
		weights[i] = weight_p;
	}

	//Overall weights normalization "alphas"
	for (int i = 0; i < num_particles; i++) {
		particles[i].weight = particles[i].weight/sum_weights;
		weights[i] = weights[i]/sum_weights;
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	//Distribution for index
	uniform_int_distribution<unsigned int> ind_dist(0, num_particles - 1);
	unsigned int ind = ind_dist(gen);

	// Declare and initialize max weight (pointer to the location of the max element)
	double max_weight = *max_element(weights.begin(), weights.end());
	double beta = 0.0;

	uniform_real_distribution<double> weight_dist(0, max_weight*2);
	vector<Particle> resampled_particles;

	//Resampling wheel from the lesson
	for (int i = 0; i < num_particles; i++) {
		beta += weight_dist(gen);

		while (weights[ind] < beta) {

			beta -= weights[ind];
			ind = (ind + 1) % num_particles;
		}

		// Add selected sample to new particles
		resampled_particles.push_back(particles[ind]);
	}

	// new set of particles to particles vector
	particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
