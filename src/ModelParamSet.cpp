/**
 * @file ModelParamSet.cpp
 * @author Robin Dowell 
 * @brief Routines necessary for loading data and segments.
 * @version 0.1
 * @date 2022-02-08
 * 
 */
#include "ModelParamSet.h"

#include <math.h>   
#include <stdio.h>
#include <unistd.h>

#include <string>
#include <vector>
#include <iostream>  // for cout?

#include "../thirdparty/nlohmann/json.hpp"
#include "split.h"
#include "helper.h"
#include "Models.h"

ModelParams::ModelParams() {
  // Empty constructor.
  mu = sigma = lambda = pi = footprint = omega[0] = 0.0;
  omega[1] = -1; omega[2] = -1;
}

/**
 * @brief Constructor
 * 
 * @param mu        The center of the bidirectional, position of loading
 * @param sigma     The variance on the loading position (mu)
 * @param lambda    The length of the exponential.  EMG = N(mu,sigma) + Exp(lambda)
 * @param pi        The strand bias
 * @param footprint     The offset parameter (aka Beta)
 * @param omega     The weights (array as individual params)
 * @param v_b     negL
 * @param v_a     posL
 */
ModelParams::ModelParams(double v_mu, double v_sigma, double v_lambda, 
            double v_pi, double v_footprint, double v_omega0, double v_omega1, 
            double v_omega2, double v_b, double v_a) {
    mu = v_mu;
    sigma = v_sigma;
    lambda = v_lambda;
    pi = v_pi;
    footprint = v_footprint;
    omega[0] = v_omega0;
    omega[1] = v_omega1;
    omega[2] = v_omega2;
    negL = v_b;
    posL = v_a;
}

/**
 * @brief Given a bidirectional model, transfer parameters
 * 
 * @param model 
 */
void ModelParams::SetfromBidirectional(Bidirectional model) {
  mu = model.getMu();
  sigma = model.getSigma();
  lambda = model.getLambda();
  footprint = model.getFootprint();
  pi = model.getPi();
  omega[0] = model.getWeight();
}

void ModelParams::SetfromFullModel(FullModel model) {
   SetfromBidirectional(model.bidir);
   omega[1] = model.forwardElongation.getWeight();
   omega[2] = model.reverseElongation.getWeight();
   posL = model.forwardElongation.uni.upper;
   negL = model.reverseElongation.uni.lower;
}

/**
 * @brief Convert the parameters into an ouput string.
 * 
 * @return std::string 
 */
std::string ModelParams::write() {
  std::string output =    std::to_string(mu) +"\t" + std::to_string(sigma) +"\t" + 
    std::to_string(lambda) + "\t" + std::to_string(pi) +"\t" + std::to_string(footprint) 
    +"\t" + std::to_string(omega[0]) + "\t" + std::to_string(omega[1])
    + "\t" + std::to_string(omega[2]) + "\t" + std::to_string(negL)
    + "\t" + std::to_string(posL) + "\n";
   return output;
}

/**
 * @brief This is the opposite of the write() function.  Expects
 * the write() function's output as its input.
 * 
 * @param single  Tab delimited list of parameters, 
 * in order: mu, sigma, lambda, pi, footprint, omegas (3), negL, posL.
 */
void ModelParams::read(std::string single) {
   vector<std::string> split_tab = split_by_tab(single);

   mu = std::stod(split_tab[0]);
   sigma = std::stod(split_tab[1]);
   lambda = std::stod(split_tab[2]);
   pi = std::stod(split_tab[3]);
   footprint = std::stod(split_tab[4]);
   omega[0] = std::stod(split_tab[5]);
   omega[1] = std::stod(split_tab[6]);
   omega[2] = std::stod(split_tab[7]);
   negL = std::stod(split_tab[8]);
   posL = std::stod(split_tab[9]);

   // Should we check that we don't have more parameters?
}

/**
 * @brief The start coordinate for a bed file with this call.
 * 
 * @return std::string 
 */
double ModelParams::getBedStart() {
    return max(mu-(sigma+lambda), 0.0);
}

/**
 * @brief The end coordinate for a bed file with this call.
 * 
 * @return std::string 
 */
double ModelParams::getBedEnd() {
   return (mu + (sigma+lambda));
}

/**
 * @brief Return the parameters as vector of strings.
 * 
 * @remark Unclear if this should be concerned with precision in conversion to string
 * 
 * @return std::vector<std::string> 
 */
std::vector<std::string> ModelParams::fetch_as_strings() {
   std::vector<std::string> asStrings(10);
   // cout << std::to_string(mu) << std::endl;

   asStrings[0] = std::to_string(mu);
   asStrings[1] = std::to_string(sigma);
   asStrings[2] = std::to_string(lambda);
   asStrings[3] = std::to_string(pi);
   asStrings[4] = std::to_string(footprint);
   asStrings[5] = std::to_string(omega[0]);
   asStrings[6] = std::to_string(omega[1]);
   asStrings[7] = std::to_string(omega[2]);
   asStrings[8] = std::to_string(negL);
   asStrings[9] = std::to_string(posL);

   return asStrings;
}

std::string ModelParams::writeAsJSON() {
  std::string output;
  output += "{\"mu\":" + tfit::prettyDecimal(mu,-1) + ",";
  output += "\"sigma\":" + tfit::prettyDecimal(sigma,2) + ",";
  output += "\"lambda\":" + tfit::prettyDecimal(lambda,2) + ",";
  output += "\"pi\":" + tfit::prettyDecimal(pi,2) + ",";
  output += "\"fp\":" + tfit::prettyDecimal(footprint,2) + ",";
  output += "\"weights\":{";
  output +=  "\"bidir\":" + tfit::prettyDecimal(omega[0], 2) + 
         ", \"forward\":" + tfit::prettyDecimal(omega[1], 2) +
         ", \"reverse\":" + tfit::prettyDecimal(omega[2], 2) + "},";
  output += "\"length_forward\":" + tfit::prettyDecimal(posL,-1) + ",";
  output += "\"length_reverse\":" + tfit::prettyDecimal(negL,-1) + "}";
  return output;
}

void ModelParams::readFromJSON(std::string entry) {
   nlohmann::json input;
   input = nlohmann::json::parse(entry);
   mu = input.value("mu", 0.);
   sigma = input.value("sigma", 0.);
   lambda = input.value("lambda", 0.);
   pi = input.value("pi", 0.);
   footprint = input.value("fp", 0.);
   omega[0] = input["weights"].value("bidir", 0.);
   omega[1] = input["weights"].value("forward", 0.);
   omega[2] = input["weights"].value("reverse", 0.);
   posL = input.value("length_forward", 0.);
   negL = input.value("length_reverse", 0.);
}


/******************************************************/
/* Class: ModelParamsSet */
/*************************/

ModelParamSet::ModelParamSet() {
   K = 0;
   log_likelihood = 0;
}

/**
 * @brief Allocate a set of ModelParams
 * 
 * @param K Number of items in the collection.
 */
ModelParamSet::ModelParamSet(int k) {
    K = k;
    ModelParams *nodata = NULL;
    for (int i = 0; i < K; i++) {
       collection.push_back(nodata) ;
    }
    log_likelihood = 0;
  // cout << K << "," << log_likelihood << std::endl;
}

/**
 * @brief Destroy the ModelParamSet
 * 
 */
ModelParamSet::~ModelParamSet() {
   for (int i = 0; i < K; i++) {
      delete(collection[i]);
   }
}

/**
 * @brief Convert a set of distinct fits into something like 
 * the K_models output string.  NOTE it's not exactly K_models 
 * as we don't keep the wf and wp parameters here.
 * 
 * UNTESTED
 * 
 * @return std::string 
 */
std::string ModelParamSet::write() {
   std::string asString;
   asString = "toWrite!";
  return asString;
}

/**
 * @brief Convert a string of K model output into a Set_EMGparameters 
 * 
 */
void ModelParamSet::readFromKmodels(std::string line) {
   // if K > 0 then there is already data here!

   std::string data = line.substr(1, line.size() - 1); // removes the ~ first character
   std::vector<std::string> tab_split = split_by_tab(data);

   // These are values for the entire set (size and loglikelihood).
   std::vector<std::string> comma_split = split_by_comma(tab_split[0], "");
   K = stoi(comma_split[0]);
   log_likelihood = stod(comma_split[1]);

   // Now we need to build all the EMGparameter sets, of which there are K.
   // Recall the string is formatted:
   // mus+"\t"+sigmas+"\t"+lambdas+"\t" + pis+"\t" + fps + "\t" + ws 
   //    + "\t" + bs + "\t" + as
   std::vector<std::vector<double>> par_as_double(11, std::vector<double>(K)); // [params][model]
   for (int i = 1; i < 6; i++) { // mu to footprint
      std::vector<std::string> temp = split_by_comma(tab_split[i], "");
      for (int j = 0; j < K; j++) { // each of the K models
         par_as_double[i][j] = stod(temp[j]);
      }
   }
   // omega is tab_split[6] but it's got a bar in it so has to have
   // separate parsing.  But should go into par_as_double[6-8].
   // w, fw, rw|w, fw, rw| ...
   std::vector<std::string> temp = split_by_bar(tab_split[6], "");
   // Note that temp.size should eq K.
   for (int i = 0; i < K; i++) {
      std::vector<std::string> wset = split_by_comma(temp[i], "");
      par_as_double[6][i] = stod(wset[0]);   // omega[0]
      par_as_double[7][i] = stod(wset[1]);   // omega[1]
      par_as_double[8][i] = stod(wset[2]);   // omega[2]
   }

   for (int i = 7; i < 9; i++) { // the last two fields 
     std::vector<std::string> temp = split_by_comma(tab_split[i], "");
     for (int j = 0; j < K; j++) { // each of the K models
       par_as_double[i+2][j] = stod(temp[j]);
     }
   }

   for (int i = 0; i < K; i++) {
      collection.push_back(new ModelParams(par_as_double[1][i], par_as_double[2][i],
               par_as_double[3][i], par_as_double[4][i], par_as_double[5][i],
               par_as_double[6][i], par_as_double[7][i], par_as_double[8][i],
               par_as_double[9][i], par_as_double[10][i]));
   }
}

std::string ModelParamSet::writeAsKmodels() {
  std::string asString = "~" + to_string(K) + "," + to_string(log_likelihood);

  // Note that output is currently all K values of mu (et. al.) 
  // But storage is in K EMGparameters -- this leads to a indexing
  // issues in the read/write functions.
  std::vector<std::vector<std::string>> params(K);
  for (int i = 0; i < K; i++) { // each of the K models
    if (collection[i]) {
      params[i] = collection[i]->fetch_as_strings();
      // cout << params[i][0] << std::endl;
    } else {   // This model is null
      params[i] = std::vector<std::string>(11,""); // Set to empty strings
    }
  }

  std::string temp;
  for (int j = 0; j < 5; j++) { // each of the parameters
     temp = params[0][j];
     // cout << temp << std::endl;
     for (int i = 1; i < K; i++) {
       temp = temp + "," + params[i][j];
     }
     asString = asString + "\t" + temp;
  }
  // Now the omegas
  temp = params[0][5] + "," + params[0][6] + "," + params[0][7];
  for (int i = 1; i < K; i++) {
     temp = temp + "|" + params[i][5] + "," + params[i][6] + "," + params[i][7];
  }
  asString = asString + "\t" + temp;

  for (int j = 8; j < 10; j++) {
     temp = params[0][j];
     // cout << temp << std::endl;
     for (int i = 1; i < K; i++) {
       temp = temp + "," + params[i][j];
     }
     asString = asString + "\t" + temp;
  }
  return asString;
}

std::string ModelParamSet::writeJSON() {
   std::string output;
   output = "{\"num_models\":" + tfit::prettyDecimal(K,-1);
   output += ", \"log-likelihood\":" + tfit::prettyDecimal(log_likelihood,2) + ", ";
   output += "\"model_params\": [";
   if (K > 1) {
     output += collection[0]->writeAsJSON();
     for (int i = 1; i < K; i++) {
       output += "," + collection[i]->writeAsJSON();
     }
   }
   output += "]}";
   return output;
}

void ModelParamSet::readJSON(std::string entry) {
   // if K > 0 then this already has data!
   nlohmann::json input;
   input = nlohmann::json::parse(entry);
   K = input.value("num_models", 0.);
   log_likelihood = input.value("log-likelihood", 0.);
   for (int i = 0; i < K; i++) {
      collection.push_back(new ModelParams());
      // std::cout << input["model_params"][i].dump() << std::endl;
      collection[i]->readFromJSON(input["model_params"][i].dump());
   }
}
