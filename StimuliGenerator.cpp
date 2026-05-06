/*
 * Norepinephrine–dopamine sequential integration in astroglia drives behavioral state 
 * transitions 
 *
 * Sijia Wang1,2†, Shilin Fang1†, Thiago O. Bezerra3†, Zhiyuan Wang1,2, Qimeng Zhao1,2, 
 * Xiaoou Wang1,2, Liang Wang4, Kui Wang1, Lijun Chen1,2, Funing Li1, Yuchen Gong1,2, 
 * Yongxin Yang1,2, Jiesi Feng5, Yu Zheng5, Yulong Li5, Jun Chu4, Gerald M. Pao6, Xufei 
 * Du1,2, Jiu-lin Du1,2, Antonio C. Roque3, Yu Mu1,2*
 * 
 * Affiliations:
 * 
 * 1 Institute of Neuroscience, State Key Laboratory of Brain Cognition and Brain-inspired 
 * Intelligence Technology, Center for Excellence in Brain Science and Intelligence Technology, 
 * Chinese Academy of Sciences, 320 Yue-Yang Road, Shanghai 200031, China.
 * 
 * 2 University of Chinese Academy of Sciences,19A Yu-Quan Road, Beijing 100049, China.
 *
 * 3 Department of Physics, School of Philosophy, Sciences and Letters of Ribeirão Preto, 
 * University of São Paulo, Ribeirão Preto, São Paulo, Brazil.
 *
 * 4 Research Center for Primate Neuromodulation and Neuroimaging, Shenzhen Institutes of 
 * Advanced Technology, Chinese Academy of Sciences, Shenzhen 518055, China.
 *
 * 5 Academy for Advanced Interdisciplinary Studies, School of Psychological and Cognitive 
 * Sciences, Peking University, Beijing 100871, China.
 *
 * 6 Okinawa Institute of Science and Technology Graduate University, Biological Nonlinear 
 * Dynamics Data Science Unit, Okinawa 904-0495, Japan.
 *
 * † These authors contributed equally to this work.
 * * Corresponding author. Email: my@ion.ac.cn
 *
 * This file provides the methods used to generate inputs.
 */
#pragma once

#include "StimuliGenerator.h"
#include <string>
#include <fstream>
#include <istream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <ios>
#include <math.h>
#include <cmath>

StimuliGenerator::StimuliGenerator(std::string parameter_file, int N_stimuli, int N_total_stimulated, double dt): parameter_file(parameter_file),
N_stimuli(N_stimuli), N_total_stimulated(N_total_stimulated), dt(dt){

    ImportParameters(parameter_file);
    ConfigureStimulator(dt);

};// end of StimuliGenerator constructor

void StimuliGenerator::ImportParameters(std::string parameter_file){

    std::string infile_stimulator_file = "Stimuli/" + parameter_file;
    std::ifstream in_file(infile_stimulator_file);

    if(!in_file){

        std::cout << "Can't open file in ImportParameters from PoissonStimulator" << std::endl;
        exit(EXIT_FAILURE);

    }else{

        std::string line;
        //getline(in_file, line);                           // Skip first line, header
        
        while(std::getline(in_file, line)){
    
            // ensure correct format for lines
            line.erase(remove(line.begin(), line.end(), ' '), line.end());        // remove ' ' 
            std::replace(line.begin(), line.end(), ';', ' ');                       // replace ';' by ' '

            std::stringstream ss(line);
            std::string auxs;

            ss >> auxs;
            g_max.push_back(std::stod(auxs));

            ss >> auxs;
            t0.push_back(std::stod(auxs));

            ss >> auxs;
            t_rise.push_back(std::stod(auxs));

            ss >> auxs;
            t_decay.push_back(std::stod(auxs));

            ss >> auxs;                                             // number of comparts stimulated
            ss >> auxs;
            ss >> auxs;
            GetStimulatedUnits(auxs);

        }//looping through rows

    }// end of the condition ensuring no problem in opening the file 

}// end of ImportParameters


void StimuliGenerator::GetStimulatedUnits(std::string stimulated_units_string){

    // Stimulated units separated by '&'
    std::replace(stimulated_units_string.begin(), stimulated_units_string.end(), '&', ' ');   // replace '&' by ' ', separate by stimulated units

    std::stringstream ss(stimulated_units_string);
    std::string auxs;

    int aux_N_stimulated = 0;
    bool first = true;
    
    while(ss >> auxs){

        aux_N_stimulated++;
        _stimulated_units.push_back(std::stoi(auxs) - 1);                   // save the astrocyte compartment index
        
    }// end while looping through AN connections

    N_stimulated.push_back(aux_N_stimulated);

}; //end of GetStimulatedUnits

void StimuliGenerator::ConfigureStimulator(double dt){

    // Configure stimulus concentration and time vector
    s.resize(N_total_stimulated);
    std::fill(s.begin(), s.end(), 0);

    // Configure pointer vector
    int aux_index = 0;
    for (int i_s = 0; i_s < N_stimuli; i_s++) {
        stimulated_units.push_back(&_stimulated_units[aux_index]);
        aux_index += N_stimulated[i_s];
    }
    
};// end of ConfigureStimulator


void StimuliGenerator::DoubleExp(double t, double g_min) {
    // double g_max, double t0, double t_rise, double t_decay
    double t_peak, K;

    for (int i_s = 0; i_s < N_stimuli; i_s++) {

        if(t >= t0[i_s]){

            for (int i_c = 0; i_c < N_stimulated[i_s]; i_c++) {

                t_peak = t0[i_s] + t_rise[i_s] * t_decay[i_s] / (t_decay[i_s] - t_rise[i_s]) * log(t_decay[i_s] / t_rise[i_s]);
                K = 1 / (exp(-(t_peak - t0[i_s]) / t_decay[i_s]) - exp(-(t_peak - t0[i_s]) / t_rise[i_s]));
                s[stimulated_units[i_s][i_c]] = g_max[i_s] * K * (exp(-(t - t0[i_s]) / t_decay[i_s]) - exp(-(t - t0[i_s]) / t_rise[i_s])) + g_min;

            }
        }
        else {

            for (int i_c = 0; i_c < N_stimulated[i_s]; i_c++) {

                s[stimulated_units[i_s][i_c]] = g_min;

            }

        }
    }
}// end of DoubleExp


void StimuliGenerator::AlphaSynapse(double t) {
    // double g_max, double t0, double tau
    for (int i_s = 0; i_s < N_stimuli; i_s++) {

        if (t >= t0[i_s]) {

            for (int i_c = 0; i_c < N_stimulated[i_s]; i_c++) {

                s[stimulated_units[i_s][i_c]] = 2 * g_max[i_s] * (t - t0[i_s]) / t_rise[i_s] * exp(-(t - t0[i_s] - t_rise[i_s]) / t_rise[i_s]);

            }
        }
    }
}// end of AlphaSynapse
