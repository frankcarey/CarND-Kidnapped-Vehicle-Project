/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
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
#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// Setup a normal distribution to use for particles.
	normal_distribution <double> dist_x (x, std[0]);
	normal_distribution <double> dist_y (y, std[1]);
	normal_distribution <double> dist_theta (theta, std[2]);

	// Set the number of particles to use.
	num_particles = 100;

	// Set length of weights and particles to num_particles.
	weights.resize(num_particles);
	particles.resize(num_particles);

	// Generate num_particles particles.
	default_random_engine gen;
	for (int i=0; i<num_particles; ++i) {
		particles.at(i).x = i;
		particles.at(i).x = dist_x(gen);
		particles.at(i).y = dist_y(gen);
		particles.at(i).theta = dist_theta(gen);
		particles.at(i).weight = 1.0;
	}

	for (int i=0; i<particles.size(); i++) {
		Particle p = particles[i];
	}

	cout << "Initialized.." << endl;
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// This predicts the new locations of the particles given the model of the car, previous measurements, and change in time
	// since the last prediction.

	// Use a psedo random number generator
	default_random_engine gen;

	// Depending if yaw rate is near zero, we need to use different set of equations to predict the x,y positions
	// and heading angle which is yaw angle theta to avoid a division by zero error.
	for (int i = 0; i < num_particles; ++i) {
		if (fabs(yaw_rate) < 0.0000001)
		{
			// Add random gaussian noise for x position
			double x_pred = particles.at(i).x + (velocity * delta_t * cos(particles.at(i).theta));
			normal_distribution <double> dist_x_pred(0, std_pos[0]);
			particles.at(i).x = x_pred + dist_x_pred(gen);

			// Add random gaussian noise for y position
			double y_pred = particles.at(i).y + (velocity * delta_t * sin(particles.at(i).theta));
			normal_distribution <double> dist_y_pred(0, std_pos[1]);
			particles.at(i).y = y_pred +  dist_y_pred(gen);

			// No need to update theta since yaw rate is zero.
		} else {
			// Add random gaussian noise for x position.
			double x_pred = particles.at(i).x + (velocity/yaw_rate)*(sin((particles.at(i).theta) + (yaw_rate * delta_t)) - (sin(particles.at(i).theta)));
			normal_distribution <double> dist_x_pred(0, std_pos[0]);
			particles.at(i).x = x_pred + dist_x_pred(gen);

			// Add random gaussian noise for y position.
			double y_pred = particles.at(i).y + (velocity / yaw_rate)*((cos(particles.at(i).theta)) - (cos(particles.at(i).theta + (yaw_rate * delta_t))));
			normal_distribution <double> dist_y_pred(0, std_pos[1]);
			particles.at(i).y = y_pred +  dist_y_pred(gen);
			double theta_pred = particles.at(i).theta + ((yaw_rate * delta_t));

			// Add random gaussian noise for theta
			normal_distribution <double> dist_theta_pred(0, std_pos[2]);
			particles.at(i).theta = theta_pred + dist_theta_pred(gen);

		}
	}


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	// NOTE: Observations MUST be previously converted to mapCordinates.

	for (int sensorCount = 0; sensorCount < observations.size(); ++sensorCount)
	{
		double min_dist = INFINITY;
		for (int predictedCount = 0; predictedCount < predicted.size(); ++predictedCount)
		{
			double dist_euclid = dist(predicted[predictedCount].x, predicted[predictedCount].y,
																observations[sensorCount].x, observations[sensorCount].y);

			// if this distance is less than min distance,associate current particle with that landmark index
			if (dist_euclid < min_dist)
			{
				// change the min distance to calculated euclidean distance
				min_dist = dist_euclid;

				// add the landmark index to predicted vector id
				observations[sensorCount].id = predicted[predictedCount].id;
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


	// This updates the weights of each particle by generating what the measurement would be to each particles nearby
	// landmarks, comparing that to the actual measurements. All weights should add up one one, with the closest matches
	// having higher weights and poorly fitting particles getting lower weights.

	std::vector<LandmarkObs> map_obs;
	std::vector<LandmarkObs> predicted_obs;

	// loop through each particle to calculate its weight
	for (int i = 0; i < num_particles; ++i)
	{
		double weight_for_each_particle = 1.0;
		std::vector<int> associations;
		std::vector<double> sense_x;
		std::vector<double> sense_y;

		// Set size of transformed observations vector
		map_obs.resize(observations.size());

		// Clear any previously predicted observations.
		predicted_obs.clear();

		// For each observation from sensor,convert it to map coordinates
		for(int sensorCount = 0; sensorCount < observations.size(); ++sensorCount) {
			// Transform each sensor measurements from local(car) coordinates to global(map) coordinates.
			map_obs[sensorCount].x = particles.at(i).x + (cos(particles.at(i).theta) * (observations[sensorCount].x)) - (sin(particles.at(i).theta) * (observations[sensorCount].y));
			map_obs[sensorCount].y = particles.at(i).y + (sin(particles.at(i).theta) * (observations[sensorCount].x)) + (cos(particles.at(i).theta) * (observations[sensorCount].y));
		}

		// Iterate through the landmarks and create a list of predicted landmark observations from this particle to that landmark
		// for use in associating measurements to landmarks.
		for(int landmarkCount = 0; landmarkCount < map_landmarks.landmark_list.size(); ++landmarkCount) {
			// Clear existing associations.

			double dist_euclid = dist(map_landmarks.landmark_list[landmarkCount].x_f, map_landmarks.landmark_list[landmarkCount].y_f,
																particles.at(i).x, particles.at(i).y);
			// If this landmark wouldn't be outside the range of the sensor, then add it as a
			// predicted observation.
			if (dist_euclid <= sensor_range) {
				LandmarkObs lm = LandmarkObs();
				lm.id = map_landmarks.landmark_list[landmarkCount].id_i;
				lm.x = double(map_landmarks.landmark_list[landmarkCount].x_f);
				lm.y = double(map_landmarks.landmark_list[landmarkCount].y_f);
				predicted_obs.push_back(lm);

				associations.push_back(lm.id);
				sense_x.push_back(lm.x);
				sense_y.push_back(lm.y);

				//cout << "distance_to_landmark " << map_landmarks.landmark_list[landmarkCount].id_i <<  ": " << dist_euclid << endl;
			}
		}
		SetAssociations(particles.at(i), associations, sense_x, sense_y);

		// apply data association for each sensor measurement and create a predicted vector for each sensor measurement
		dataAssociation(predicted_obs, map_obs);
		//predicted.push_back(observations[sensorCount]);

		// calculate weight using multivariate gaussian probability distribution
		for (int obs = 0; obs < map_obs.size(); ++obs)
		{
			/* start process of assigning weights to each particle
			using multivariate gaussian probability distribution*/

			// retrieve landmark x and y position associated with ith sensor
			// measurement based on the index that is stored during dataAssociation().
			// We need to subtract one here because the index of the landmark is one less than the id of the landmark :/
			int landmark_index = map_obs[obs].id - 1;
			//cout << "observation: " << obs << " likely landmark (index not id):" << landmark_index << endl;
			// retrieve the x and y positions of the landmark
			double x_landmark_diff =  map_obs.at(obs).x - map_landmarks.landmark_list.at(landmark_index).x_f;
			double y_landmark_diff =  map_obs.at(obs).y - map_landmarks.landmark_list.at(landmark_index).y_f;
			//cout << "x_landmark_diff: " << x_landmark_diff << " y_landmark_diff: " << y_landmark_diff << endl;


			// calculate normalizer
			double gauss_norm = (1 / (2 * M_PI * std_landmark[0] * std_landmark[1]));

			//calculate the exponent
			double exponent = ((pow(x_landmark_diff,2))/(2 * pow(std_landmark[0],2))) +
												((pow(y_landmark_diff, 2)) / (2 * pow(std_landmark[1], 2)));

			weight_for_each_particle = weight_for_each_particle * gauss_norm * exp(-exponent);
			//cout << "weight_for_each_particle: " << weight_for_each_particle << endl;

		}

		particles.at(i).weight = weight_for_each_particle;
		weights.at(i) = weight_for_each_particle;
		//	cout << "Weights:" << endl;
		//cout << "Weight of particle :" << weights.at(ii) << endl;
	}




}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::vector<Particle> resampled_particles;
	resampled_particles.clear();
	std::random_device rd;
	std::mt19937 gen(rd());
	std::discrete_distribution<> d(weights.begin(), weights.end());
	std::map<int, int> m;
	for (int n = 0; n<num_particles; ++n) {
		Particle p = particles[d(gen)];
		resampled_particles.push_back(p);
		//printf("\nPOST Init particle: id: %d | x: %f | y: %f | theta: %f |\n\n", p.id, p.x, p.y, p.theta);
	}
	particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
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
