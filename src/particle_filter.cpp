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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  std::default_random_engine gen;
  num_particles = 100;  // TODO: Set the number of particles
  std::normal_distribution<double> gx(x, std[0]);
  std::normal_distribution<double> gy(y, std[1]);
  std::normal_distribution<double> gtheta(theta, std[2]);
  Particle particle;
  for(int i = 0; i < num_particles; i++){
    particle.x = gx(gen);
    particle.y = gy(gen);
    particle.theta = gtheta(gen);
    particle.id = i;
    particle.weight = 1.0;
    particles.push_back(particle);
    weights.push_back(1.0);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;
  double xf, yf, thetaf;
  for(int i = 0; i < particles.size(); ++i){
    if(yaw_rate == 0){
      thetaf = particles[i].theta;
      xf = particles[i].x + velocity * delta_t * cos(particles[i].theta);
      yf = particles[i].y + velocity * delta_t * sin(particles[i].theta);
    }
    else{
      thetaf = particles[i].theta + yaw_rate*delta_t;
      xf = particles[i].x + velocity/yaw_rate*(sin(thetaf) - sin(particles[i].theta));
      yf = particles[i].y + velocity/yaw_rate*(cos(particles[i].theta) - cos(thetaf));
    }
    std::normal_distribution<double> dist_x(xf, std_pos[0]);
    std::normal_distribution<double> dist_y(yf, std_pos[1]);
    std::normal_distribution<double> dist_theta(thetaf, std_pos[2]);
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);

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
  
  double distant;
  //std::cout << "predicted:" << predicted;
  //std::cout << "observations:" << observations;
  for(int i = 0; i < observations.size(); ++i){
    double min = std::numeric_limits<double>::max();
    int min_index = -1;
    for(int j = 0; j < predicted.size(); ++j){
      distant = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
      if(distant < min){
        min = distant;
        min_index = predicted[j].id ;
      }
    }
    observations[i].id = min_index;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
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
  double x_map, y_map;
  vector<LandmarkObs> predicted, landmarks;
  LandmarkObs single_predicted, single_landmark;
  double gauss_norm, exponent, d, lx, ly;
  //dataAssociation(predicted, observations);
  for(int i = 0; i < particles.size(); ++i){
    predicted.clear();
    for(int j = 0; j < observations.size(); ++j){
      x_map = particles[i].x + (cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y);
      y_map = particles[i].y + (sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y);
      single_predicted.x = x_map;
      single_predicted.y = y_map;
      single_predicted.id = observations[j].id;
      predicted.push_back(single_predicted);
    }
    landmarks.clear();
    for(int k = 0; k < map_landmarks.landmark_list.size(); ++k){
      single_landmark.x = map_landmarks.landmark_list[k].x_f;
      single_landmark.y = map_landmarks.landmark_list[k].y_f;
      single_landmark.id = map_landmarks.landmark_list[k].id_i;
      d = dist(single_landmark.x, single_landmark.y, particles[i].x, particles[i].y);
      if(d <= sensor_range){
        landmarks.push_back(single_landmark);
      }
    }
    dataAssociation(landmarks, predicted);
    particles[i].weight = 1.0;
    for(int k = 0; k < predicted.size(); ++k){
      for(int m = 0; m < landmarks.size(); ++m){
        if(landmarks[m].id == predicted[k].id){
          lx = landmarks[m].x;
          ly = landmarks[m].y;
          break;
        }
      }
      //gauss_norm = ;
      //std::cout<<"gauss_norm = "<<gauss_norm;
      //exponent = (pow(predicted[k].x - lx, 2) / (2 * pow(std_landmark[0], 2))) + (pow(predicted[k].y - ly, 2) / (2 * pow(std_landmark[1], 2)));
      //std::cout<<"exponent = "<<exponent;
      particles[i].weight *= 1.0 / (2 * M_PI * std_landmark[0] * std_landmark[1]) * exp(-(pow(predicted[k].x - lx, 2) / (2 * pow(std_landmark[0], 2)) + (pow(predicted[k].y - ly, 2) / (2 * pow(std_landmark[1], 2)))));
      //std::cout<<particles[i].weight<<std::endl; 
    }
    weights[i] = particles[i].weight;
    //std::cout<<"particles "<<i<<" weights: "<<weights[i]<<std::endl;
  }
  //std::cout<<"weights:";
  //for(int n = 0; n < weights.size(); ++n){
    //std::cout<<" "<<weights[n];
  //}
}

void ParticleFilter::resample(){
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;
  vector<Particle> particles_temp;
  std::discrete_distribution<> d(weights.begin(), weights.end());
  int index;
  particles_temp.clear();
  for(int n = 0; n < num_particles; ++n){
    index = d(gen);
    particles_temp.push_back(particles[index]);
  }
  particles = particles_temp;
  
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