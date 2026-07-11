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
 * This file provides the methods for running the simulation, interface for classes 
 * integration, and saving results.
 */
#pragma once

#include <string>
#include <fstream>
#include <istream>
#include <iostream>
#include "Simulator.h"
#include <vector>

Simulator::Simulator(std::string simulation_parameter_file): simulation_parameter_file(simulation_parameter_file){

    //srand(time(NULL));

    // Configure simulation
    std::ifstream spf(simulation_parameter_file);
    std::string mod, param, val;
    while(spf >> mod >> param >> val) {
                
        if(mod == "Simulator") ChangeParameters(param, val);

    }// End while to process the simulation parameter file
    spf.close();

    outfile_astrocytes = "AR_" + output_file;
        
};// End of Simulator constructor

void Simulator::SimulateAstrocyticNetwork(){

    n_points = (int) (t_total/dt);
    n_points_out = (int) (n_points/sample_rate);
    double t = 0.0;

    
    Astrocyte astrocytes = Astrocyte(simulation_parameter_file);
    // --- Stimuli --- //

    // NE
    StimuliGenerator NE(ne_ext_parameter_file, N_ne_ext_stim, N_ne_ext_total, dt);

    // DA
    StimuliGenerator DA(da_ext_parameter_file, N_da_ext_stim, N_da_ext_total, dt);


    for(int it = 0; it < n_points; it++){
        
        t = it*dt;


        // --------------- Update Neurotransmitter Variables --------------- //

        DA.AlphaSynapse(t);
        NE.DoubleExp(t, astrocytes.NErest);

        astrocytes.IncrementNeurotransmitter(&(astrocytes.NE), &(astrocytes.ne_stim_comparts), NE.N_stimuli, &(NE.N_stimulated), &(NE.s));
        astrocytes.IncrementNeurotransmitter(&(astrocytes.DA), &(astrocytes.da_stim_comparts), DA.N_stimuli, &(DA.N_stimulated), &(DA.s));


        // ------------------ Update Astrocyte Variables ------------------ //        
        astrocytes.UpdateVariables(dt);
            
        // ------------------------- Save variables ------------------------ //

        if(it % sample_rate == 0){

            astrocytes.SaveData(t, &astrocytes.Ca_i, "Cai", output_file);
            astrocytes.SaveData(t, &astrocytes.IP3, "IP3", output_file);
            astrocytes.SaveData(t, &astrocytes.cAMP, "cAMP", output_file);            
            astrocytes.SaveData(t, &astrocytes.NE, "NE", output_file);
            astrocytes.SaveData(t, &astrocytes.DA, "DA", output_file);

        }

    }// end for time loop

};// end of SimulateAstrocyticNetwork

void Simulator::ChangeParameters(std::string parameter, std::string value){
    
    if(parameter == "t_total") t_total = std::stof(value);
    else if(parameter == "dt") dt = std::stof(value);
    else if(parameter == "sample_rate") sample_rate = std::stoi(value);
    else if(parameter == "output_file") output_file = value;
    else if(parameter == "ne_ext_parameter_file") ne_ext_parameter_file = value;
    else if(parameter == "N_ne_ext_stim") N_ne_ext_stim = std::stoi(value);
    else if(parameter == "N_ne_ext_total") N_ne_ext_total = std::stoi(value);
    else if(parameter == "da_ext_parameter_file") da_ext_parameter_file = value;
    else if(parameter == "N_da_ext_stim") N_da_ext_stim = std::stoi(value);
    else if(parameter == "N_da_ext_total") N_da_ext_total = std::stoi(value);
    else{
            std::cout << "Parameter " << parameter << " not recognized in the simulation parameters file." << std::endl;
            exit(EXIT_FAILURE);

    }// End of conditional for parameter selection

}//end of ChangeParameters
