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
#include <vector>
#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	printf("\n Initialize.. \n");
	num_particles = 100;
	weights.resize(100, 1.0);

	default_random_engine gen;

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for (int i=0; i<num_particles; i++) {
		Particle p;
		p.id=i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;
		particles.push_back(p);
	}

	for (int i=0; i<particles.size(); i++) {
		Particle p = particles[i];
		printf("\nPOST Init particle: id: %d | x: %f | y: %f | theta: %f |\n\n", p.id, p.x, p.y, p.theta);
	}

	cout << "Initialized";
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	//printf("\n Prediction.. \n");

	// Setup generators to add gaussian noise.
	default_random_engine gen;
	normal_distribution<double> dist_noise_x(0, std_pos[0]);
	normal_distribution<double> dist_noise_y(0, std_pos[1]);
	normal_distribution<double> dist_noise_theta(0, std_pos[2]);

	for (int i=0; i<num_particles; i++) {
		if (fabs(yaw_rate) < 0.00001){ //Yawrate zero case Handle
			particles[i].x += velocity * delta_t*cos(particles[i].theta);
			particles[i].y += velocity * delta_t*sin(particles[i].theta);
		}
		particles[i].x += ((velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta))) + dist_noise_x(gen);
		particles[i].y += ((velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta) + yaw_rate * delta_t)) + dist_noise_y(gen);
		particles[i].theta += ((yaw_rate * delta_t) + dist_noise_theta(gen));
		//printf("\nprediciton particle: id: %d | x: %f | y: %f | theta: %f |\n\n", particles[i].id, particles[i].x, particles[i].y, particles[i].theta);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (int obs_idx=0; obs_idx < observations.size(); obs_idx++) {
		double min_dist = -1;
		for(int pred_idx=0; pred_idx < predicted.size(); pred_idx++) {
			double curr_dist = dist(predicted[pred_idx].x, predicted[pred_idx].y, observations[obs_idx].x, observations[obs_idx].y);
			if (min_dist == -1 || min_dist > curr_dist) {
				min_dist = curr_dist;
				observations[obs_idx].id = predicted[pred_idx].id;
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		 std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
	//printf("\n updateWeights.. \n");
	std::vector<LandmarkObs> landmarks;
	std::vector<LandmarkObs> txObservations;

	// Convert observations to the car's coordinate system.
	for(int p_idx=0; p_idx< particles.size(); p_idx++) {

		int    p_id = particles[p_idx].id;
		double p_x = particles[p_idx].x;
		double p_y = particles[p_idx].y;
		double p_theta = particles[p_idx].theta;


		//printf("\nParticle: id: %d | x: %f | y: %f | theta: %f | \n\n", p_id, p_x, p_y, p_theta);


		txObservations.clear();
		for(int obs_idx=0; obs_idx < observations.size(); obs_idx++) {
			//printf("\n\tObservations (car space): x: %f | y: %f | \n", observations[obs_idx].x, observations[obs_idx].y);

			LandmarkObs tx_obs {};
			// Convert observations back from Particle space to Global Space.
			tx_obs.x = cos(p_theta) * observations[obs_idx].x - sin(p_theta) * observations[obs_idx].y + p_x;
			tx_obs.y = sin(p_theta) * observations[obs_idx].x + cos(p_theta) * observations[obs_idx].y + p_y;
			tx_obs.id = observations[obs_idx].id;
			//printf("\t             (world space): x: %f | y: %f | \n", tx_obs.x, tx_obs.y);
			txObservations.push_back(tx_obs);
		}

		// Get all landmarks that should be within range of the sensor.
		landmarks.clear();
		for ( int lm_idx=0; lm_idx< map_landmarks.landmark_list.size();lm_idx++) {
			Map::single_landmark_s lm = map_landmarks.landmark_list[lm_idx];
			if (dist(lm.x_f, lm.y_f, p_x, p_y) <= sensor_range) {
				landmarks.push_back(LandmarkObs{lm.id_i, lm.x_f, lm.y_f});
			}
		}

		dataAssociation(landmarks, txObservations);
		vector<double> sense_x;
    vector<double> sense_y;
		vector<int> associations;

		double weight = 1.;
		for(int obs_idx=0; obs_idx < txObservations.size(); obs_idx++) {
			int lm_id {};
			double lm_x {};
			double lm_y {};

			// Find the landmark associated with the Observation.
			for (int lm_idx=0; lm_idx < landmarks.size(); lm_idx++) {
				if (landmarks[lm_idx].id == txObservations[obs_idx].id ) {
					lm_id = landmarks[lm_idx].id;
					lm_x = landmarks[lm_idx].x;
					lm_y = landmarks[lm_idx].y;
					sense_x.push_back(lm_x);
					sense_y.push_back(lm_y);
					associations.push_back(lm_id);
					//printf("\tClosest LM (world space): closest_lm_id: %d | x: %f | y: %f | \n", lm_id, lm_x, lm_y);
					break;
				}
			}

			double d_x = lm_x - txObservations[obs_idx].x;
			double d_y = lm_y - txObservations[obs_idx].y;

			double gauss_norm = (1 / (2 * M_PI * std_landmark[0] * std_landmark[1]));
			double exponent= (d_x*d_x)/(2 * std_landmark[0]*std_landmark[0]) + (d_y*d_y)/(2 * std_landmark[1] * std_landmark[1]);
			weight *= gauss_norm * exp(-exponent);
			//cout << "\t       weight: " << weight << "\n\n";

		}

		//particles[p_idx] = SetAssociations(particles[p_idx], associations, sense_x, sense_y);
		// Update this particles weight.
		particles[p_idx].weight = weight;
		// Add to the list of all of the particle weights for use later.
		weights[p_idx] = weight;

	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::default_random_engine gen;
	std::discrete_distribution<> resample_dist(weights.begin(), weights.end());
	std::vector<Particle> new_particles;
	new_particles.clear();
	for(int n=0; n<num_particles; ++n) {
		Particle p = particles[resample_dist(gen)];
		p.id=n;
		new_particles.push_back(p);
	}
	// Override the old particles with the new.
	particles = new_particles;
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

	particle.associations= std::move(associations);
 	particle.sense_x = std::move(sense_x);
 	particle.sense_y = std::move(sense_y);

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
