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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    is_initialized = true ;
    num_particles = 100 ;

	default_random_engine gen;

	// TODO: Set standard deviations for x, y, and theta
	 double std_x = std[0];
	 double std_y = std[1];
	 double std_theta = std[2];


	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);


	
	for (int i = 0; i < num_particles; ++i) {
		 Particle pt;
		 pt.id = i;
		 pt.x = dist_x(gen);
		 pt.y = dist_y(gen);
		 pt.theta = dist_theta(gen);
		 pt.weight = 1;
		 particles.push_back(pt);
		
	//	 cout << "Sample " << i + 1 << " " << particles[i].x << " " << particles[i].y << " " << particles[i].theta << endl;
	}
	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    default_random_engine gen;

	for (int i = 0; i < num_particles; ++i) {
		
		double xpost, ypost,thetapost ;
		
		if (abs(yaw_rate)>1e-8){
			xpost = particles[i].x + velocity/yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			ypost = particles[i].y + velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			thetapost = particles[i].theta + yaw_rate*delta_t;
		}else{
			thetapost = particles[i].theta;
			xpost = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			ypost = particles[i].y + velocity*delta_t*sin(particles[i].theta);
		}
		 

         normal_distribution<double> dist_x2(xpost, std_pos[0]);
         normal_distribution<double> dist_y2(ypost, std_pos[1]);
         normal_distribution<double> dist_theta2(thetapost, std_pos[2]);
         particles[i].x = dist_x2(gen);
		 particles[i].y = dist_y2(gen);
		 particles[i].theta = dist_theta2(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {

    // Gather std values for readability
        // Iterate over all particles
    for (int i = 0; i < num_particles; ++i) {   
        // List all landmarks within sensor range
        vector<LandmarkObs> predicted;

		for (const auto& map_landmark: map_landmarks.landmark_list){
			int land_id = map_landmark.id_i;
			double land_x =(double) map_landmark.x_f;
			double land_y =(double) map_landmark.y_f;
			//double distance = sqrt((particles[i].x-land_x)*(particles[i].x-land_x) +(particles[i].y-land_y)*(particles[i].y-land_y) );
            double distance = dist(particles[i].x,particles[i].y,land_x,land_y);													
			if (distance<sensor_range){
				LandmarkObs MarksInrange;
				MarksInrange.id = land_id;
				MarksInrange.x = land_x;
				MarksInrange.y = land_y;
				predicted.push_back(MarksInrange);
				//cout<< "123456789";
			}
		}
		
        // List all observations in map coordinates
        vector<LandmarkObs> Observed;
        for (size_t j = 0; j < observations.size(); ++j) {

            // Convert observation from particle(vehicle) to map coordinate system
            LandmarkObs rotate;
            rotate.x = cos(particles[i].theta) * observations[j].x - sin(particles[i].theta) * observations[j].y + particles[i].x;
            rotate.y = sin(particles[i].theta) * observations[j].x + cos(particles[i].theta) * observations[j].y + particles[i].y;

            Observed.push_back(rotate);
        }

        // Find which observations correspond to which landmarks (associate ids)
        //dataAssociation(predicted, Observed);

		
		for (auto& obs: Observed) {
			double min_distance = 1e16;
			//cout << min_distance;
			//cout << "\n";
			for (auto& pred: predicted){
				double distance = sqrt((obs.x-pred.x)*(obs.x-pred.x) +(obs.y-pred.y)*(obs.y-pred.y) );
				if (distance < min_distance){
					obs.id = pred.id ;
					min_distance = distance;
				}
			}
		}
        // Compute the likelihood for each particle, that is the probablity of obtaining
        // current observations being in state (particle_x, particle_y, particle_theta)
			double likelyhood = 1;
							  
			double mu_x, mu_y;
			for (const auto& obs:Observed){
																						 
														   
				for (const auto& land:predicted){
					if (obs.id == land.id){
						mu_x = land.x;
						mu_y = land.y;
						double norm_factor = 2*M_PI*std_landmark[0]*std_landmark[1];
						double prob = exp( - (pow(obs.x-mu_x,2)/(2*std_landmark[0]*std_landmark[0]) + pow(obs.y-mu_y,2)/(2*std_landmark[1]*std_landmark[1])) ) ;								  
						likelyhood *= prob/norm_factor;	
					}
															  
					  		
				}
			}
			particles[i].weight = likelyhood;
							   
		}
			
    double norm_fact = 0;
    for (const auto& part : particles){
        norm_fact += part.weight;
	}
    // Normalize weights s.t. they sum to one
    for (auto& part : particles){
        //part.weight /= (norm_fact + numeric_limits<double>::epsilon());
		part.weight /= norm_fact;
	}
   
  
  
 
 
}

																					 
								 
void ParticleFilter::resample() {
	default_random_engine gen;
	vector<double> pw;
	for (const auto& pt: particles){
		pw.push_back(pt.weight);
	}
								  
	discrete_distribution<int> weighted_distribution(pw.begin(),pw.end()) ;													

																									   

										 
	vector<Particle> resample;
	for (int i = 0; i < num_particles; ++i) {
		int k = weighted_distribution(gen);
		resample.push_back(particles[k]);
	}
									
	particles = resample;
	
									
	for(auto& part:particles){
		part.weight=1;
	}
	
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
