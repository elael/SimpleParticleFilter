/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

static std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, std::array<double,3> std) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 20;  // TODO: Set the number of particles

  // Set of current particles
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  particles.reserve(num_particles);
  std::generate_n(std::back_inserter(particles), num_particles,
    [&,n=-1]() mutable {return Particle{++n, dist_x(gen), dist_y(gen), dist_theta(gen), 0, {}, {}, {}};}
  );
  
  // Vector of weights of all particles
  weights.reserve(num_particles);
  std::fill_n(std::back_inserter(weights), num_particles, 1);
  
  // Flag, if filter is initialized
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, std::array<double,3> std_pos, 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  if (fabs(yaw_rate) < 1e-5){
    auto& distance = velocity;
    distance *= delta_t;
    for (auto & particle: particles)
    {    
      // Straight Line when yaw_rate is null
      particle.x += distance * cos(particle.theta);
      particle.y += distance * sin(particle.theta);

      // Add GPS-like noises, as process noise
      particle.x = normal_distribution(particle.x, std_pos[0])(gen);
      particle.y = normal_distribution(particle.y, std_pos[1])(gen);
      particle.theta = normal_distribution(particle.theta, std_pos[2])(gen);    
    }
  }
  else{
    auto& distance = velocity;
    distance /= yaw_rate;

    auto& theta_displacement = yaw_rate;
    theta_displacement *= delta_t;

    for (auto & particle: particles)
    {    
      // Coordinated Turn Model
      particle.x += distance * (sin(particle.theta + theta_displacement) - sin(particle.theta));
      particle.y += distance * (cos(particle.theta) - cos(particle.theta + theta_displacement));
      particle.theta += theta_displacement;

      // Add GPS-like noises, as process noise
      particle.x = normal_distribution(particle.x, std_pos[0])(gen);
      particle.y = normal_distribution(particle.y, std_pos[1])(gen);
      particle.theta = normal_distribution(particle.theta, std_pos[2])(gen);    
    }
  }
  


}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  for (auto& obs: observations)
  {
    double nn_dist = std::numeric_limits<double>::infinity(); 

    for (const auto& mark: predicted)
    if(double mark_dist = dist(obs.x, obs.y, mark.x, mark.y); mark_dist < nn_dist){
      obs.id = mark.id;
      nn_dist = mark_dist;
    }

  }
  

}

void ParticleFilter::updateWeights(double sensor_range, std::array<double,2> std_landmark, 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  double max_weight = -std::numeric_limits<double>::infinity();

  for (auto & particle: particles){

    // Filter inrange Map landmarks
    vector<LandmarkObs> landmarks_inrange;

    for(const auto& land_mark: map_landmarks.landmark_list)
    if (fabs(land_mark.x_f - particle.x) < sensor_range && fabs(land_mark.y_f - particle.y) < sensor_range)
    landmarks_inrange.push_back({land_mark.id_i, land_mark.x_f, land_mark.y_f});

    // Transform observations to Map coordinates
    vector<LandmarkObs> observations_inmap;
    observations_inmap.reserve(observations.size());

    for(const auto& obs: observations)
    observations_inmap.push_back({
      obs.id, 
      particle.x + cos(particle.theta)*obs.x - sin(particle.theta)*obs.y,
      particle.y + sin(particle.theta)*obs.x + cos(particle.theta)*obs.y});

    // Associate correct id to observations
    // dataAssociation(landmarks_inrange, observations_inmap);
    
    // clear stored associations
    particle.associations.clear();
    particle.associations.reserve(observations_inmap.size());
    particle.sense_x.clear();
    particle.sense_x.reserve(observations_inmap.size());
    particle.sense_y.clear();
    particle.sense_y.reserve(observations_inmap.size());

    for (auto& obs: observations_inmap)
    {
      double nn_dist = std::numeric_limits<double>::infinity(); 

      const LandmarkObs* best_mark;
      for (const auto& mark: landmarks_inrange)
      if(double mark_dist = dist(obs.x, obs.y, mark.x, mark.y); mark_dist < nn_dist){
        obs.id = mark.id;
        nn_dist = mark_dist;
        best_mark = &mark;
      }
      
      particle.associations.emplace_back(obs.id);
      particle.sense_x.emplace_back(obs.x);
      particle.sense_y.emplace_back(obs.y);

      const auto x_error = (obs.x - best_mark->x)/std_landmark[0];
      const auto y_error = (obs.y - best_mark->y)/std_landmark[1];
      particle.weight -= x_error*x_error + y_error*y_error;
    }

    if (particle.weight > max_weight)
      max_weight = particle.weight; 
  }

  // normalize max weight to zero
  for (auto & particle: particles)
    particle.weight -= max_weight;

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::vector<double> weights;
  weights.reserve(particles.size());
  for (auto & particle: particles){
    weights.emplace_back(exp(particle.weight/2.0));
    particle.weight = 0;
  }

  std::discrete_distribution particle_index(weights.begin(), weights.end());
  
  std::vector<Particle> new_particles;
  new_particles.reserve(particles.size());
  std::generate_n(std::back_inserter(new_particles), particles.size(), [&](){return particles[particle_index(gen)];});

  particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
