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
		particles.at(i).x = dist_x(gen);
		particles.at(i).y = dist_y(gen);
		particles.at(i).theta = dist_theta(gen);
		particles.at(i).weight = 1.0;
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
	printf("\n Prediction.. \n");

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
		printf("\nprediciton particle: id: %d | x: %f | y: %f | theta: %f |\n\n", particles[i].id, particles[i].x, particles[i].y, particles[i].theta);
	}
}

//void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
//	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
//	//   observed measurement to this particular landmark.
//	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
//	//   implement this method and use it as a helper during the updateWeights phase.
//
//}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs>& predicted, std::vector<LandmarkObs>& observations, Map map_landmarks)
{

	for (int sensorCount = 0; sensorCount < observations.size(); ++sensorCount)
	{
		double min_dist = INFINITY;
		// transform each sensor measurement ,which is in local(car) coordinates to global(map) coordinates
		/*
		Equation to transform car coordinates to map coordinates
		|Xmap| = |cos(theta)  -sin(theta)   Xparticle |      | Xcar |
		|Ymap| = |sin(theta)   cos(theta)   Yparticle |   X  | Ycar |
		|  1 | = |    0             0            1    |      |   1  |

		*/
		//// use above equation for each sensor
		//observations[sensorCount].x = particles.at(ii).x + (cos(particles.at(ii).theta) * (observations[sensorCount].x)) - (sin(particles.at(ii).theta) * (observations[sensorCount].y));
		//observations[sensorCount].y = particles.at(ii).y + (sin(particles.at(ii).theta) * (observations[sensorCount].x)) + (cos(particles.at(ii).theta) * (observations[sensorCount].y));

		for (int landmarkCount = 0; landmarkCount < map_landmarks.landmark_list.size(); ++landmarkCount)
		{
			// TODO## restrict far away landmarks suing sensor range

			// Calculate euclidean distance between each landmark pos and sensor pos
			//double dist_euclid = sqrt(pow(map_landmarks.landmark_list[landmarkCount].x_f - observations[sensorCount].x, 2) +
			//					 pow(map_landmarks.landmark_list[landmarkCount].y_f - observations[sensorCount].y, 2));
			double dist_euclid = dist(map_landmarks.landmark_list[landmarkCount].x_f, map_landmarks.landmark_list[landmarkCount].y_f,
																observations[sensorCount].x, observations[sensorCount].y);

			// if this distance is less than min distance,associate current particle with that landmark index
			if (dist_euclid < min_dist)
			{
				// change the min distance to calculated euclidean distance
				min_dist = dist_euclid;

				// add the landmark index to predicted vector id
				observations[sensorCount].id = map_landmarks.landmark_list[landmarkCount].id_i;
			}

		}
		//cout << "transf_obs ( " << observations[sensorCount].x << " " << observations[sensorCount].y << "--->"
		//	<< "landmark position : ( " << map_landmarks.landmark_list[observations[sensorCount].id-1].x_f << " " <<
		//	map_landmarks.landmark_list[observations[sensorCount].id-1].y_f << " )" << " " << "ID: " <<
		//observations[sensorCount].id << endl;
		// now add the sensor measurement with its associated id to predicted vector
		predicted.push_back(observations[sensorCount]);
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

	std::vector<LandmarkObs> transformed_obs;

	// loop through each particle to calculate its weight
	for (int ii = 0; ii < num_particles; ++ii)
	{
		double weight_for_each_particle = 1.0;

		// set size of transformed observations vector
		transformed_obs.resize(observations.size());
		//cout << "Transformations :" << endl;
		// for each observation from sensor,convert it to map coordinates
		for (int sensorCount = 0; sensorCount < observations.size(); ++sensorCount)
		{
			// transform each sensor measurement ,which is in local(car) coordinates to global(map) coordinates
			/*
			//	Equation to transform car coordinates to map coordinates
			//	|Xmap| = |cos(theta)  -sin(theta)   Xparticle |      | Xcar |
			//	|Ymap| = |sin(theta)   cos(theta)   Yparticle |   X  | Ycar |
			//	|  1 | = |    0             0            1    |      |   1  |
			//
			//	*/
			// use above equation for each sensor to transform observations
			transformed_obs[sensorCount].x = particles.at(ii).x + (cos(particles.at(ii).theta) * (observations[sensorCount].x)) - (sin(particles.at(ii).theta) * (observations[sensorCount].y));
			transformed_obs[sensorCount].y = particles.at(ii).y + (sin(particles.at(ii).theta) * (observations[sensorCount].x)) + (cos(particles.at(ii).theta) * (observations[sensorCount].y));

			//cout << "observation: (" << observations[sensorCount].x <<" " << observations[sensorCount].y << ") --->" <<
			//	"transfomed_obs: (" << transformed_obs[sensorCount].x << " " << transformed_obs[sensorCount].y << ")" << endl;
		}

		//cout << "Data associations :" << endl;
		// apply data association for each sensor measurement and create a predicted vector for each sensor measurement
		dataAssociation(predicted, transformed_obs, map_landmarks);




		// calculate weight using multivariate gaussian probability distribution
		for (int ii = 0; ii < observations.size(); ++ii)
		{
			/* start process of assigning weights to each particle
			using multivariate gaussian probability distribution*/

			// retrieve landmark x and y position associated with ith sensor
			// measurement based on the index that is stored in predicted vector
			int associated_index = predicted[ii].id;

			// retrieve the x and y positions of the landmark
			double x_landmark = map_landmarks.landmark_list[associated_index-1].x_f;
			double y_landmark = map_landmarks.landmark_list[associated_index-1].y_f;


			// calculate normalizer
			double gauss_norm = (1 / (2 * M_PI * std_landmark[0] * std_landmark[1]));

			//calculate the exponent
			double exponent = ((pow(predicted[ii].x - x_landmark,2))/(2 * pow(std_landmark[0],2))) +
												((pow(predicted[ii].y - y_landmark, 2)) / (2 * pow(std_landmark[1], 2)));

			weight_for_each_particle = weight_for_each_particle * gauss_norm * exp(-exponent);

		}

		particles.at(ii).weight = weight_for_each_particle;
		weights.at(ii) = weight_for_each_particle;
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
		resampled_particles.push_back(particles[d(gen)]);
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
