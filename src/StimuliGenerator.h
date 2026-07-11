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
 * This header defines the StimuliGenerator class.
 */
#pragma once

#include <vector>
#include <random>

class StimuliGenerator{

    public:

        std::string parameter_file;
        int N_stimuli = 1;
        int N_total_stimulated = 0;

        std::vector <int> N_stimulated;
        std::vector <int *> stimulated_units;
        std::vector <int> _stimulated_units;  

        std::vector<double> g_max;
        std::vector<double> t0;
        std::vector<double> t_rise;
        std::vector<double> t_decay;

        std::vector<double> s;

        double dt;

        StimuliGenerator(std::string parameter_file, int N_stimuli, int N_total_stimulated, double dt);
        void ImportParameters(std::string parameter_file);
        void GetStimulatedUnits(std::string stimulated_units_string);
        void ConfigureStimulator(double dt);
        void DoubleExp(double t, double g_min);
        void AlphaSynapse(double t);



}; // end of StimuliGenerator class
