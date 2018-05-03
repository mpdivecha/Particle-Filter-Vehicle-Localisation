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

#include "particle_filter.h"

using namespace std;

/**
 * Initializes particle filter by initializing particles to Gaussian
 *   distribution around first position and all the weights to 1.
 * @param {double} x  Initial x position [m] (simulated estimate from GPS)
 * @param {double} y  Initial y position [m]
 * @param {double} theta   Initial orientation [rad]
 * @param {double[]} std   Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
 *   standard deviation of yaw [rad]]
 */ 
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
    /* We create 1000 particles to begin with. We initialize three distributions
     *  for x, y and theta respectively. The mean and std devs are as supplied.
     *  We sample the initial configurlations of each particle from these distributions
     */
    default_random_engine gen;

    num_particles = 500;

    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    for (int i = 0; i < num_particles; i++) {
        
        Particle p;
        p.id = i;
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_theta(gen);
        particles.push_back(p);
    }
    weights.resize(num_particles);

    is_initialized = true;
}

/**
 * prediction Predicts the state for the next time step
 *   using the process model.
 * @param delta_t Time between time step t and t+1 in measurements [s]
 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
 *   standard deviation of yaw [rad]]
 * @param velocity Velocity of car from t to t+1 [m/s]
 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
 */
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    
    // Some short names
    double dt = delta_t;
    double v = velocity;
    double yr = yaw_rate;

    // lightweight random generator
	default_random_engine gen;

    for (int i = 0; i < num_particles; i++) {
        double x_0 = particles[i].x;
        double y_0 = particles[i].y;
        double theta_0 = particles[i].theta;

        // The CTRV model
        double x_1, y_1, theta_1;
        if (yr != 0)
        {
            x_1 = x_0 + v/yr*(sin(theta_0 + yr*dt) - sin(theta_0));
            y_1 = y_0 + v/yr*(cos(theta_0) - cos(theta_0 + yr*dt));
            theta_1 = theta_0 + yr*dt;
        }
        else
        {
            x_1 = x_0 + v*dt*cos(theta_0);
            y_1 = y_0 + v*dt*sin(theta_0);
            theta_1 = theta_0;
        }

        // Gaussian distributions for x, y and theta
        normal_distribution<double> dist_x(x_1, std_pos[0]);
        normal_distribution<double> dist_y(y_1, std_pos[1]);
        normal_distribution<double> dist_theta(theta_1, std_pos[2]);

        // Update the particle state with predictions
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
    }

}

/**
 * dataAssociation Finds which observations correspond to which landmarks (likely by using
 *   a nearest-neighbors data association).
 * @param predicted Vector of predicted landmark observations
 * @param observations Vector of landmark observations
 */
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

/**
 * updateWeights Updates the weights for each particle based on the likelihood of the 
 *   observed measurements. 
 * @param sensor_range Range [m] of sensor
 * @param std_landmark[] Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]]
 * @param observations Vector of landmark observations
 * @param map Map class containing map landmarks
 */
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

    /**
     * First transform the observations to map coordinates.
     * Associate each observation with the nearest map landmark using the dist() function
     *   (this can be done in the above function)
     * Assign weights to each of the particles
     */
    
    double& sig_x = std_landmark[0];    // reference instead of copy for speed
    double& sig_y = std_landmark[1];    // reference instead of copy for speed
    double gauss_norm= (1/(2 * M_PI * sig_x * sig_y));
    //weights.clear();
    for(int i = 0; i < num_particles; i++)
    {
        Particle& p = particles[i];
        double P_weight = 1.0f;
        // Clear the associations
        p.associations.clear();
        p.sense_x.clear();
        p.sense_y.clear();
        // Transform from vehicle to map coordinates
        for (LandmarkObs o : observations)
        {
            double x_m = p.x + cos(p.theta)*o.x - sin(p.theta)*o.y;
            double y_m = p.y + sin(p.theta)*o.x + cos(p.theta)*o.y;

            // Find the nearest landmark to be associated with this observation
            double best_dist = 9999999;
            Map::single_landmark_s best_landmark;
            for (Map::single_landmark_s l : map_landmarks.landmark_list)
            {
                double range = dist(x_m, y_m, l.x_f, l.y_f);
                if (range <= sensor_range && range < best_dist)
                {
                    best_dist = range;
                    best_landmark = l;
                }
            }
            o.id = best_landmark.id_i;

            // Set-up associations and coords
            p.associations.push_back(o.id);
            p.sense_x.push_back(x_m);
            p.sense_y.push_back(y_m);

            // The exponent
            double exponent = (pow((x_m - best_landmark.x_f),2))/(2*pow(sig_x,2)) + (pow((y_m - best_landmark.y_f),2))/(2*pow(sig_y,2));

            // The weight
            double weight = gauss_norm * exp(-exponent);

            P_weight *= weight;
        }
        // Update the particle weight
        p.weight = P_weight;
        weights[i] = P_weight;
    }
}

/**
 * resample Resamples from the updated set of particles to form
 *   the new set of particles.
 */
void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    default_random_engine gen;
    discrete_distribution<int> sampler(weights.begin(), weights.end());

    // Create a new vector. The old vector should not be overwritten until after
    //  the sampling is done.
    vector<Particle> newParticles;
    for (int i = 0; i < num_particles; i++)
    {
        newParticles.push_back(particles[sampler(gen)]);
    }

    particles.clear();
    particles = newParticles;
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
