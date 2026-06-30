/*
 * Astroglia integrate sequential norepinephrine-dopamine signals to drive behavioral
 * state transitions
 *
 * Sijia Wang1,2†, Shilin Fang1†, Thiago O. Bezerra3†, Zhiyuan Wang1,2, Xiaoou Wang1,2,
 * Liang Wang4, Kui Wang1, Lijun Chen1,2, Funing Li1, Yuchen Gong1,2, Yongxin Yang1,2,
 * Jiesi Feng5, Yu Zheng5, Yulong Li5, Jun Chu6, Gerald M. Pao7, Xu-fei Du1,2,
 * Antonio C. Roque3* , Yu Mu1,2*
 *
 *
 *
 * Affiliations:
 *
 * 1 Institute of Neuroscience, State Key Laboratory of Brain Cognition and Brain-inspired Intelligence Technology,
 * Center for Excellence in Brain Science and Intelligence Technology, Chinese Academy of Sciences, 320
 * Yue-Yang Road, Shanghai 200031, China.
 *
 * 2 University of Chinese Academy of Sciences, 19A Yu-Quan Road, Beijing 100049, China.
 *
 * 3 Department of Physics, School of Philosophy, Sciences and Letters of Ribeirão Preto, University of São Paulo,
 * Ribeirão Preto, São Paulo, Brazil.
 *
 * 4 Faculty of Biomedical Engineering, Shenzhen University of Advanced Technology, Shenzhen, 518055, China.
 *
 * 5 State Key Laboratory of Membrane Biology, Peking University School of Life Sciences, Beijing 100871, China.
 *
 * 6 Research Center for Primate Neuromodulation and Neuroimaging, Shenzhen Institutes of Advanced Technology,
 * Chinese Academy of Sciences, Shenzhen 518055, China.
 *
 * 7 Okinawa Institute of Science and Technology Graduate University, Biological Nonlinear Dynamics Data Science
 * Unit, Okinawa 904-0495, Japan.
 *
 *
 * †These authors contributed equally to this work.
 * *Corresponding author. Email: my@ion.ac.cn / antonior@usp.br
 * 
 * This header defines the Simulator class.
 */
#pragma once

#include <string>
#include <fstream>
#include <istream>
#include <iostream>
#include "Astrocyte.h"
#include "StimuliGenerator.h"
#include <vector>

class Simulator{

    double t_total = 100;       // in ms
    double dt = 0.01e-3;       // in ms
    int sample_rate = 1000;  // output sample rate
    int n_points;
    int n_points_out;

    std::string simulation_parameter_file;
    std::string output_file;
    std::string outfile_astrocytes;

    std::string ne_ext_parameter_file;
    int N_ne_ext_stim = 0;
    int N_ne_ext_total = 0;

    std::string da_ext_parameter_file;
    int N_da_ext_stim = 0;
    int N_da_ext_total = 0;
    

    public:

        Simulator(std::string simulation_parameter_file);
        void SimulateAstrocyticNetwork();
        void ChangeParameters(std::string parameter, std::string value);

};// End of Neurons:Cell class
