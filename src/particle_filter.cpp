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
   * Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * Add random Gaussian noise to each particle.
   */
  num_particles = 20;  // Set the number of particles

  // Set of current particles
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  particles.reserve(num_particles);
  std::generate_n(std::back_inserter(particles), num_particles,
    [&,n=-1]() mutable {return Particle{++n, dist_x(gen), dist_y(gen), dist_theta(gen), 0, {}, {}, {}};}
  );
  
  // Flag, if filter is initialized
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, std::array<double,3> std_pos, 
                                double velocity, double yaw_rate) {
  /**
   * Add measurements to each particle and add random Gaussian noise.
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

void ParticleFilter::updateWeights(double sensor_range, std::array<double,2> std_landmark, 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
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
    
    // clear stored associations
    particle.associations.clear();
    particle.associations.reserve(observations_inmap.size());
    particle.sense_x.clear();
    particle.sense_x.reserve(observations_inmap.size());
    particle.sense_y.clear();
    particle.sense_y.reserve(observations_inmap.size());

    // add observations information
    for (auto& obs: observations_inmap)
    {
      
      // select best match from landmarks for the given observation
      double nn_dist = std::numeric_limits<double>::infinity(); 
      const LandmarkObs* best_mark;
      for (const auto& mark: landmarks_inrange)
      if(double mark_dist = dist(obs.x, obs.y, mark.x, mark.y); mark_dist < nn_dist){
        obs.id = mark.id;
        nn_dist = mark_dist;
        best_mark = &mark;
      }
      
      // add associations for the observation
      particle.associations.emplace_back(obs.id);
      particle.sense_x.emplace_back(obs.x);
      particle.sense_y.emplace_back(obs.y);

      // add error from observation
      const auto x_error = (obs.x - best_mark->x)/std_landmark[0];
      const auto y_error = (obs.y - best_mark->y)/std_landmark[1];
      particle.weight -= x_error*x_error + y_error*y_error;
    }

    // keep track of best particle
    if (particle.weight > max_weight){
      best_particle_pt = &particle;
      max_weight = particle.weight;
    }
  }

  // normalize max weight to zero
  for (auto & particle: particles)
    particle.weight -= max_weight;

}

const Particle& ParticleFilter::best_particle() const{
  // Just return the first particle if undecided 
  if(best_particle_pt == nullptr)
    return particles[0];

  return *best_particle_pt;
}

void ParticleFilter::resample() {
  /**
   * Resample particles with replacement with probability proportional 
   *   to their weight.
   */
  // create weight vector and clear original weight
  std::vector<double> weights;
  weights.reserve(particles.size());
  for (auto & particle: particles){
    weights.emplace_back(exp(particle.weight/2.0));
    particle.weight = 0;
  }

  // create new set of particles based on particle weights
  std::vector<Particle> new_particles;
  new_particles.reserve(particles.size());
  std::discrete_distribution particle_index(weights.begin(), weights.end());
  std::generate_n(std::back_inserter(new_particles), particles.size(), [&](){return particles[particle_index(gen)];});

  // set new particles and nullify best_particle_pt
  best_particle_pt = nullptr;
  particles = std::move(new_particles);
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
