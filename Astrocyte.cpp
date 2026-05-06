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
 *
 * This file provides the  model equations, auxiliary methods for loading parameters, 
 * and methods for solving the system of ODEs.
 */
#pragma once

#include <string>
#include <fstream>
#include <istream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <random>
#include <ios>
#include <memory> 
#include <math.h>
#include "Astrocyte.h"
#include <iomanip>

Astrocyte::Astrocyte(std::string parameter_file):parameter_file(parameter_file){

    // Configure simulation
    std::ifstream pf(parameter_file);
    std::string mod, param, val;
    while(pf >> mod >> param >> val) { 
                
        if(mod == "Astrocyte") ChangeParameters(param, val);

    }// End while to process the simulation parameter file
    pf.close();
        
    // Configure neuron and astrocytic networks and set the output files names

    ConfigureAstrocyticVectors();

    LoadAstrocyteModel(astrocytic_network_file);

    ConfigureAstrocyticStimulators(ne_parameter_file, &ne_stim_comparts, &_ne_stim_comparts, &N_stimulated_ne);
    ConfigureAstrocyticStimulators(da_parameter_file, &da_stim_comparts, &_da_stim_comparts, &N_stimulated_da);

};// End of Astrocyte constructor;

void Astrocyte::ConfigureAstrocyticVectors(){

    // Resize vectors containing morphological parameters
    n_c_connections.resize(N_comparts);
    c_connections.resize(N_comparts);
    _c_connections.resize(N_comparts_connections);

    A.resize(N_comparts);
    Vol.resize(N_comparts);
    r_ER.resize(N_comparts);

    // --- Resize and initialize vectors --- //
    Ca_i.resize(N_comparts);
    Ca_ER.resize(N_comparts);
    IP3.resize(N_comparts);
    h.resize(N_comparts);
    Ca_o.resize(N_comparts);
    Na_i.resize(N_comparts);
    Na_o.resize(N_comparts);
    K_i.resize(N_comparts);
    K_o.resize(N_comparts);
    V.resize(N_comparts);
    NE.resize(N_comparts);
    DA.resize(N_comparts);
    Ca_i_prev.resize(N_comparts);
    c_signal.resize(N_comparts);

    AC.resize(N_comparts);
    AC_Gas.resize(N_comparts);
    D1R.resize(N_comparts);
    D1R_DA.resize(N_comparts);
    D1R_DA_Gas.resize(N_comparts);
    Gas.resize(N_comparts);
    GasA.resize(N_comparts);
    AC_Gas_ATP.resize(N_comparts);
    cAMP.resize(N_comparts);
    R2C2.resize(N_comparts);
    R2C2cAMP2.resize(N_comparts);
    R2C2cAMP4.resize(N_comparts);
    PKA.resize(N_comparts);

    std::fill(Ca_i.begin(), Ca_i.end(), Ca_irest);
    std::fill(Ca_ER.begin(), Ca_ER.end(), Ca_ERrest);
    std::fill(IP3.begin(), IP3.end(), IP3_rest);
    std::fill(h.begin(), h.end(), h_rest);
    std::fill(Ca_o.begin(), Ca_o.end(), Ca_orest);
    std::fill(Na_i.begin(), Na_i.end(), Na_irest);
    std::fill(Na_o.begin(), Na_o.end(), Na_orest);
    std::fill(K_i.begin(), K_i.end(), K_irest);
    std::fill(K_o.begin(), K_o.end(), K_orest);
    std::fill(V.begin(), V.end(), v_rest);
    std::fill(NE.begin(), NE.end(), NErest);
    std::fill(DA.begin(), DA.end(), DArest);
    std::fill(Ca_i_prev.begin(), Ca_i_prev.end(), Ca_irest);
    std::fill(c_signal.begin(), c_signal.end(), 0);

    std::fill(AC.begin(), AC.end(), AC_rest);
    std::fill(AC_Gas.begin(), AC_Gas.end(), AC_Gas_rest);
    std::fill(D1R.begin(), D1R.end(), D1R_rest);
    std::fill(D1R_DA.begin(), D1R_DA.end(), D1R_DA_rest);
    std::fill(D1R_DA_Gas.begin(), D1R_DA_Gas.end(), D1R_DA_Gas_rest);
    std::fill(Gas.begin(), Gas.end(), Gas_rest);
    std::fill(GasA.begin(), GasA.end(), GasA_rest);
    std::fill(AC_Gas_ATP.begin(), AC_Gas_ATP.end(), AC_Gas_ATP_rest);
    std::fill(cAMP.begin(), cAMP.end(), cAMP_rest);
    std::fill(R2C2.begin(), R2C2.end(), R2C2_rest);
    std::fill(R2C2cAMP2.begin(), R2C2cAMP2.end(), R2C2cAMP2_rest);
    std::fill(R2C2cAMP4.begin(), R2C2cAMP4.end(), R2C2cAMP4_rest);
    std::fill(PKA.begin(), PKA.end(), PKA_rest);


    // Resize stimuli vectors
    ne_stim_comparts.resize(N_ne_total);
    _ne_stim_comparts.resize(2*N_ne_total);
    N_stimulated_ne.resize(N_ne_stim);

    da_stim_comparts.resize(N_da_total);
    _da_stim_comparts.resize(2*N_da_total);
    N_stimulated_da.resize(N_da_stim);

}// end of ConfigureAstrocyticVectors

// ------------ Support Functions for Morphology ---------- //

void Astrocyte::LoadAstrocyteModel(std::string astrocyte_file){
   
    std::string infile_topology_file = "Networks/" + astrocyte_file;
    std::ifstream in_file(infile_topology_file);

    if(!in_file){

        std::cout << "Can't open file in LoadAstrocyteModel" << std::endl;
        in_file.close();
        exit(EXIT_FAILURE);

    }else{

        int N = 0;

        std::string line;

        while(std::getline(in_file, line)){

            // Get the connection list from file
            line.erase(remove(line.begin(), line.end(), ' '), line.end());        // remove ' ' 
            std::replace(line.begin(), line.end(), ';', ' ');                     // replace ';' by ' '

            std::stringstream ss(line);
            std::string auxs;

            ss >> auxs;                                             // skip fourth col cell ID
            ss >> auxs;                                             // number of compartments
            ss >> auxs;                                             // morphological file
            LoadMorphologicalParams(auxs);  // Load morphological parameters of the astrocyte model

        }//looping through rows

        in_file.close();

    }//end of the condition ensuring no problem in opening the file

}// end of Load Morphological Model;

void Astrocyte::LoadMorphologicalParams(std::string morphology_file){

    std::string infile_morphology_file = "Morphology/" + morphology_file;
    std::ifstream in_file(infile_morphology_file);

    if(!in_file){

        std::cout << "Can't open file in LoadMorphologicalParam in astrocyte module." << std::endl;
        in_file.close();
        exit(EXIT_FAILURE);

    }else{

        // Auxiliary number of connections and compartment connections
        int i_c = 0;
        int aux_compart_index = 0;
            
        // Get file lines
        std::string line;
        
        while(std::getline(in_file, line)){

            // Get the connection list from file
            line.erase(remove(line.begin(), line.end(), ' '), line.end());        // remove ' ' 
            std::replace(line.begin(), line.end(), ',', ' ');                       // replace ',' by ' '

            std::stringstream ss(line);
            std::string auxs;

            ss >> auxs;                                             // first col ID
            ss >> auxs;                                             // second col compartment type
            ss >> auxs;                                             // third col number of connections

            n_c_connections[i_c] = std::stoi(auxs);

            ss >> auxs;                                             // fourth col compartment area           

            A[i_c] = std::stod(auxs);

            ss >> auxs;                                              // fifth col compartment volume
            Vol[i_c] = std::stod(auxs);

            ss >> auxs;                                              // sixth col compartment r_ER
            r_ER[i_c] = std::stod(auxs);

            ss >> auxs;                                              // skip compartment r_m
			ss >> auxs;                                              // skip compartment n_c
            ss >> auxs;                                             // connections
            std::replace(auxs.begin(), auxs.end(), '&', ' ');       // separate compartments

            std::stringstream ss_cc(auxs);
            std::string auxs_cc;

            if (n_c_connections[i_c] != 0) {

                c_connections[i_c] = &_c_connections[aux_compart_index];
                while (ss_cc >> auxs_cc) {
                    _c_connections[aux_compart_index] = std::stoi(auxs_cc);
                    aux_compart_index++;
                }// end of while looping through connections
            }
                        
            i_c++;

        }// end looping through compartments 

        in_file.close();
                
    }//end of the condition ensuring no problem in opening the file

}; // end of LoadMorphologicalParams

void Astrocyte::ConfigureAstrocyticStimulators(std::string stimulator_parameter_file, std::vector <int *> *stim_comparts, 
std::vector <int> *_stim_comparts, std::vector <int> *N_stimulated){

    std::string infile_stimulator_file = "Stimuli/" + stimulator_parameter_file;
    std::ifstream in_file(infile_stimulator_file);

    if(!in_file){

        std::cout << "Can't open " << infile_stimulator_file << " in ConfigureAstrocyticStimulators" << std::endl;
        in_file.close();
        exit(EXIT_FAILURE);

    }else{

        std::string line;
        int i_stim = 0;
        int aux_N_stimulated = 0;

        while(std::getline(in_file, line)){

            // Get the connection list from file
            line.erase(remove(line.begin(), line.end(), ' '), line.end());        // remove ' ' 
            std::replace(line.begin(), line.end(), ';', ' ');                     // replace ';' by ' '

            std::stringstream ss(line);
            std::string auxs;

            ss >> auxs;                                             // skip first col, g_max			
            ss >> auxs;                                             // skip second col, t0
            ss >> auxs;                                             // skip third col, t_rise			
			ss >> auxs;                                             // skip third col, t_decay
            ss >> auxs;                                             // number of comparts stimulated			
            (*N_stimulated)[i_stim] = std::stoi(auxs);

            ss >> auxs;                                             // astrocyte stim. parameters
            // Each astrocyte stimulus is given by a string as 
            // "(astrocyte number, compartment number, position in s vector)&()..."
            
            auxs.erase(remove(auxs.begin(), auxs.end(), '('), auxs.end());        // remove '(' 
            auxs.erase(remove(auxs.begin(), auxs.end(), ')'), auxs.end());        // remove ')' 
            std::replace(auxs.begin(), auxs.end(), '&', ' ');                     // replace '&' by ' '

            // get each astrocytic stim. triple. Now separated as (..., ..., ...) (..., ..., ...)

            std::stringstream sss(auxs);
            std::string auxss;
            
            while(sss >> auxss){
    
                // get each triple components. Now as (... ... ...)
                std::replace(auxss.begin(), auxss.end(), ',', ' ');                     // replace ',' by ' '

                std::stringstream ssss(auxss);
                std::string auxsss;

                ssss >> auxsss;                                                       // get compartment ID
                (*_stim_comparts)[2*aux_N_stimulated+0] = std::stoi(auxsss);

                ssss >> auxsss;                                                       // get position in s vector
                (*_stim_comparts)[2*aux_N_stimulated+1] = std::stoi(auxsss);

                aux_N_stimulated++;
            }// end of while looping through triple elements

            //aux_N_stimulated += (*N_stimulated)[i_stim];
            i_stim++;

        }//looping through triples

        int aux_index = 0;
        for(int i_s = 0; i_s < (*N_stimulated).size(); i_s++){
            (*stim_comparts)[i_s] = &((*_stim_comparts)[2*aux_index]);
            aux_index += (*N_stimulated)[i_s];
        }
        
        //for(int i_s = 0; i_s < 12; i_s++){
        //    std::cout << i_s << " " << (*_stim_comparts)[i_s] << ",";
        //}
        
        in_file.close();

    }//end of the condition ensuring no problem in opening the file

}// end of ConfigureAstrocyticStimulators

// ----------------------- Model Equations ----------------- //

double Astrocyte::f_prod_PLCb_NE(double NE, double Ca_i){

    double NE_alpha = pow(NE, alpha);
    double Ca_i_alpha = pow(K_R + K_p * Ca_i / (Ca_i + K_pi), alpha);

    return v_beta * NE_alpha / (NE_alpha + Ca_i_alpha);
}

double Astrocyte::f_prod_PLCb_DA(double DA, double Ca_i) {
      
    //double K_pDA = 10e-3 * 40.53419055237087;

    double DA_beta = pow(DA, beta); 
    double Ca_i_beta = pow(K_DA + K_pDA * Ca_i / (Ca_i + K_pi), beta);

    return v_DA * DA_beta / (DA_beta + Ca_i_beta);
}

double Astrocyte::f_prod_PLCd(double IP3, double Ca_i){

    double Ca_i_2 = pow(Ca_i, 2);

    return v_delta / (1 + IP3/kappa_delta) * Ca_i_2 / (Ca_i_2 + pow(K_PLCdelta, 2));
}

double Astrocyte::f_degr_IP3_3K(double Ca_i, double IP3){

    double Ca_i_4 = pow(Ca_i, 4);

    return v_3K * Ca_i_4 / (Ca_i_4 + pow(K_D, 4)) * IP3 / (IP3 + K_3);
}

double Astrocyte::f_degr_IP_5P(double IP3){return r_5P * IP3;}

double Astrocyte::f_J_IP3R(double A, double Vol, double IP3, double Ca_i, double h, double Ca_ER, double PKA){

    double r_C = r_C_max * PKA / (PKA + k_PKA) + r_C_min;
    return F * Vol / A * r_C * pow(IP3 / (IP3 + d_1), 3) * pow(Ca_i / (Ca_i + d_5), 3) * pow(h, 3) * (Ca_ER - Ca_i);
}

double Astrocyte::f_J_CERleak(double A, double Vol, double Ca_ER, double Ca_i){return F * Vol / A * r_L * (Ca_ER - Ca_i);}

double Astrocyte::f_J_SERCA(double A, double Vol, double Ca_i){

    double Ca_i_2 = pow(Ca_i, 2);

    return F * Vol / A * v_ER * Ca_i_2 / (Ca_i_2 + pow(K_ER, 2));
}

double Astrocyte::f_J_NKA(double Na_i, double K_o){

    double Na_i_15 = pow(Na_i, 1.5);

    return I_NKAmax * (Na_i_15 / (Na_i_15 + pow(K_NKAmN, 1.5))) * (K_o / (K_o + K_NKAmK));
}

double Astrocyte::f_J_NCX(double Na_o, double Ca_o, double Na_i, double V, double Ca_i){

    double Na_o_3 = pow(Na_o, 3);
    double Na_i_3 = pow(Na_i, 3);

    double alpha = I_NCXmax * (Na_o_3 / (Na_o_3 + pow(K_NCXmN, 3))) * (Ca_o / (Ca_o + K_NCXmC));
    double beta = Na_i_3 / Na_o_3 * exp(eta * V * F / (R * T)) - Ca_i / Ca_o * exp((eta - 1)* V * F / (R * T));
    double gama = 1 + k_sat * exp((eta - 1)* V * F / (R * T));

    return alpha * beta / gama;
}

double Astrocyte::f_J_Naleak(double V){return g_Naleak * (V - E_Na);}

double Astrocyte::f_J_Kleak(double V){return g_Kleak * (V - E_K);}

double Astrocyte::Diffusion(int i_c, double D, std::vector<double> *var){

    double J_diff = 0.0;
    double var_pre = (*var)[i_c]; 
    double var_pos = 0;

    double corr = 0.000375;
    int i_pos = 0;

    if (n_c_connections[i_c] != 0) {
        for (int i_cc = 0; i_cc < n_c_connections[i_c]; i_cc++) {

            i_pos = c_connections[i_c][i_cc];

            if (i_c == 0) { corr = 0.000375; }
            else { corr = 1; }

            var_pos = (*var)[i_pos];
            J_diff += corr * D * (var_pos - var_pre);

        }// end for looping through each compartment's connections
    }

    return J_diff;

};// end CalculateCompartmentCaDiffusion

// ------------------------- ODE ------------------------- //

double Astrocyte::dCa_odt(int i_c){

    double J_NCX = f_J_NCX(Na_o[i_c], Ca_o[i_c], Na_i[i_c], V[i_c], Ca_i[i_c]);
    double J_diff_Cao  = Diffusion(i_c, D_Cao, &Ca_o);
        
    return -A[i_c]/(Vol[i_c]*F) * J_NCX + J_diff_Cao;
}

double Astrocyte::dNa_odt(int i_c){
        
    double J_diff_Nao  = Diffusion(i_c, D_Nao, &Na_o);
    
    double J_NKA    = f_J_NKA(Na_i[i_c], K_o[i_c]);
    double J_NCX    = f_J_NCX(Na_o[i_c], Ca_o[i_c], Na_i[i_c], V[i_c], Ca_i[i_c]);
    double J_Naleak = f_J_Naleak(V[i_c]);
    
    return -A[i_c] / (Vol[i_c] * F) * (- 3 * J_NKA - 3 * J_NCX - J_Naleak) + J_diff_Nao;
}

double Astrocyte::dK_odt(int i_c){ 
    
    double J_diff_Ko = Diffusion(i_c, D_Ko, &K_o);
    
    double J_NKA   = f_J_NKA(Na_i[i_c], K_o[i_c]);
    double J_Kleak = f_J_Kleak(V[i_c]);
        
    return -A[i_c] / (Vol[i_c] * F) * (2 * J_NKA - J_Kleak) + J_diff_Ko;
}

double Astrocyte::dCa_idt(int i_c){
    
    // Diffusio
    double J_diff_Cai  = Diffusion(i_c, D_Ca, &Ca_i);
    
    // Transmembrane Currents
    double J_NCX     = f_J_NCX(Na_o[i_c], Ca_o[i_c], Na_i[i_c], V[i_c], Ca_i[i_c]);

    double J_memb    = A[i_c]/(Vol[i_c] * F) * J_NCX;

    // TransER Currents
    double J_IP3R    = f_J_IP3R (A[i_c], Vol[i_c], IP3[i_c], Ca_i[i_c], h[i_c], Ca_ER[i_c], PKA[i_c]);
    double J_SERCA   = f_J_SERCA(A[i_c], Vol[i_c], Ca_i[i_c]);
    double J_CERleak = f_J_CERleak(A[i_c], Vol[i_c], Ca_ER[i_c], Ca_i[i_c]);

    double J_ER      = A[i_c]*sqrt(r_ER[i_c])/(Vol[i_c]*F) * (J_IP3R - J_SERCA + J_CERleak);

    return J_memb + J_ER + J_diff_Cai;
}

double Astrocyte::dhdt(int i_c){

    return a_2*(d_2 * (IP3[i_c] + d_1) / (IP3[i_c] + d_3) * (1 - h[i_c]) - Ca_i[i_c] * h[i_c]);
}

double Astrocyte::dIP3dt(int i_c){
    
    double J_diff_IP3  = Diffusion(i_c, D_IP3, &IP3);
    
    double prod_PLCb_NE  = f_prod_PLCb_NE(NE[i_c], Ca_i[i_c]);
    double prod_PLCb_DA = f_prod_PLCb_DA(DA[i_c], Ca_i[i_c]);
    double prod_PLCd     = f_prod_PLCd(IP3[i_c], Ca_i[i_c]); 
    double degr_IP3_3K   = f_degr_IP3_3K(Ca_i[i_c], IP3[i_c]);
    double degr_IP_5P    = f_degr_IP_5P(IP3[i_c]);
    
    return prod_PLCb_NE + prod_PLCb_DA + prod_PLCd - degr_IP3_3K - degr_IP_5P + J_diff_IP3;
}

double Astrocyte::dCa_ERdt(int i_c){
    
    double J_diff_CaER = Diffusion(i_c, D_CaER, &Ca_ER);
    
    double J_IP3R    = f_J_IP3R (A[i_c], Vol[i_c], IP3[i_c], Ca_i[i_c], h[i_c], Ca_ER[i_c], PKA[i_c]);
    double J_SERCA   = f_J_SERCA(A[i_c], Vol[i_c], Ca_i[i_c]);
    double J_CERleak = f_J_CERleak(A[i_c], Vol[i_c], Ca_ER[i_c], Ca_i[i_c]);
    
    return A[i_c] * sqrt(r_ER[i_c]) / (Vol[i_c] * F * r_ER[i_c]) * (-J_IP3R + J_SERCA - J_CERleak) + J_diff_CaER;
}

double Astrocyte::dNa_idt(int i_c){

    double J_diff_Nai = Diffusion(i_c, D_Na, &Na_i);
    
    double J_NKA    = f_J_NKA(Na_i[i_c], K_o[i_c]);
    double J_NCX    = f_J_NCX(Na_o[i_c], Ca_o[i_c], Na_i[i_c], V[i_c], Ca_i[i_c]);
    double J_Naleak = f_J_Naleak(V[i_c]);
    
    return A[i_c] / (Vol[i_c] * F) * (- 3 * J_NKA - 3 * J_NCX - J_Naleak) + J_diff_Nai;
}

double Astrocyte::dK_idt(int i_c){

    double J_diff_Ki = Diffusion(i_c, D_K, &K_i);
    
    double J_NKA   = f_J_NKA(Na_i[i_c], K_o[i_c]);
    double J_Kleak = f_J_Kleak(V[i_c]);
        
    return A[i_c] / (Vol[i_c] * F) * (2 * J_NKA - J_Kleak) + J_diff_Ki;
}

double Astrocyte::dVdt(int i_c){
    
    // Transmembrane currents
    double J_NCX    = f_J_NCX(Na_o[i_c], Ca_o[i_c], Na_i[i_c], V[i_c], Ca_i[i_c]);
    double J_NKA    = f_J_NKA(Na_i[i_c], K_o[i_c]);
    double J_Naleak = f_J_Naleak(V[i_c]);
    double J_Kleak  = f_J_Kleak(V[i_c]);

    // TransER currents
    double J_IP3R    = f_J_IP3R (A[i_c], Vol[i_c], IP3[i_c], Ca_i[i_c], h[i_c], Ca_ER[i_c], PKA[i_c]);
    double J_SERCA   = f_J_SERCA(A[i_c], Vol[i_c], Ca_i[i_c]);
    double J_CERleak = f_J_CERleak(A[i_c], Vol[i_c], Ca_ER[i_c], Ca_i[i_c]);
    
    return -1 / C_m * (J_NCX + J_NKA + J_Naleak + J_Kleak -2*J_IP3R + 2*J_SERCA - 2*J_CERleak);
}

double Astrocyte::dACGasdt(int i_c) {return k1_ACGas * AC[i_c] * GasA[i_c] - (k2_ACGas + k3_ACGas * ATP) * AC_Gas[i_c];}
double Astrocyte::dGasAdt(int i_c) {return k1_GasA * D1R_DA_Gas[i_c] - k2_GasA * GasA[i_c];}
double Astrocyte::dD1R_DA_Gasdt(int i_c) {return k1_D1R_DA_Gas * D1R_DA[i_c] * Gas[i_c] - (k2_D1R_DA_Gas + k3_D1R_DA_Gas) * D1R_DA_Gas[i_c];}
double Astrocyte::dGasdt(int i_c) {return k1_Gas * GasA[i_c] - k2_Gas * D1R_DA[i_c] * Gas[i_c];}
double Astrocyte::dD1R_DAdt(int i_c) {return k1_D1R_DA * D1R[i_c] * DA[i_c] - k2_D1R_DA * D1R_DA[i_c];}
double Astrocyte::dAC_Gas_ATPdt(int i_c) {return k1_AC_Gas_ATP * AC_Gas[i_c] * ATP - k2_AC_Gas_ATP * AC_Gas_ATP[i_c];}
double Astrocyte::dcAMPdt(int i_c) {return k1_cAMP * AC_Gas_ATP[i_c] - k2_cAMP * AC_Gas[i_c] * cAMP[i_c] - V_PDE * cAMP[i_c] / (cAMP[i_c] + kcat);}
double Astrocyte::dR2C2dt(int i_c) {return k1_R2C2 * R2C2[i_c] * pow(cAMP[i_c], 2) - k2_R2C2 * R2C2cAMP2[i_c];}
double Astrocyte::dR2C2cAMP2dt(int i_c) {return k1_R2C2cAMP2 * R2C2[i_c] * pow(cAMP[i_c], 2) + k2_R2C2cAMP2 * R2C2cAMP4[i_c] - k3_R2C2cAMP2 * R2C2cAMP2[i_c] * pow(cAMP[i_c], 2);}
double Astrocyte::dR2C2cAMP4dt(int i_c) {return k1_R2C2cAMP4 * R2C2cAMP2[i_c] * pow(cAMP[i_c], 2) - (k2_R2C2cAMP4 + k3_R2C2cAMP4 * PKA[i_c]) * R2C2cAMP4[i_c];}
double Astrocyte::dPKAdt(int i_c) { return k1_PKA * R2C2cAMP4[i_c] - k2_PKA * PKA[i_c] * R2cAMP4; }


void Astrocyte::UpdateVariables(double dt){

    trigger_c_signal = false;

    // Calculate k4 and update variables (k0)
    for(int i_c = 0; i_c < N_comparts; i_c++){

        c_signal[i_c] = 0;
        Ca_i_prev[i_c] = Ca_i[i_c];

        Ca_i[i_c]  += dt * dCa_idt(i_c);
        Ca_ER[i_c] += dt * dCa_ERdt(i_c);
        IP3[i_c]   += dt * dIP3dt(i_c);
        h[i_c]     += dt * dhdt(i_c);
        Ca_o[i_c]  += dt * dCa_odt(i_c);
        Na_i[i_c]  += dt * dNa_idt(i_c);
        Na_o[i_c]  += dt * dNa_odt(i_c);
        K_i[i_c]   += dt * dK_idt(i_c);
        K_o[i_c]   += dt * dK_odt(i_c);
        V[i_c]     += dt * dVdt(i_c);

        AC[i_c]         -= dt * dACGasdt(i_c);
        AC_Gas[i_c]     += dt * dACGasdt(i_c);
        Gas[i_c]        += dt * dGasdt(i_c);
        D1R[i_c]        -= dt * dD1R_DAdt(i_c);
        D1R_DA[i_c]     += dt * dD1R_DAdt(i_c);
        D1R_DA_Gas[i_c] += dt * dD1R_DA_Gasdt(i_c);
        GasA[i_c]       += dt * dGasAdt(i_c);
        AC_Gas_ATP[i_c] += dt * dAC_Gas_ATPdt(i_c);
        cAMP[i_c]       += dt * dcAMPdt(i_c);
        R2C2[i_c]       -= dt * dR2C2dt(i_c);
        R2C2cAMP2[i_c]  += dt * dR2C2cAMP2dt(i_c);
        R2C2cAMP4[i_c]  += dt * dR2C2cAMP4dt(i_c);
        PKA[i_c]        += dt * dPKAdt(i_c);

        CheckCaSignalInAstrocyte(i_c);

    }// end for looping through compartments for k4 and updating variables

}// end of UpdateVariables

// ------------------------- Support Functions ------------------------- //

void Astrocyte::CheckCaSignalInAstrocyte(int i_c){

    if((Ca_i[i_c] > ca_thresh) * (Ca_i_prev[i_c] < ca_thresh)){
        c_signal[i_c] = 1;
        trigger_c_signal = true;
    };

};// end of CheckCaSignal

void Astrocyte::SaveData(float time, std::vector <double> *data, std::string data_name, std::string output_file){

    //std::string output_string = "";
    //output_string = output_string + std::to_string(time);

    std::ostringstream oss;

    oss << std::fixed << std::setprecision(6) << time;  // 15 digits after decimal

    for(int i_c = 0; i_c < N_comparts; i_c++){
            
        //output_string = output_string + ";";
        //output_string = output_string + std::to_string((*data)[i_c]);

        oss << ";";
        oss << std::fixed << std::setprecision(12) << (*data)[i_c];  // 15 digits after decimal
            
    };// End looping through compartments 

    //output_string = output_string + "\n";

    oss << "\n";
    std::string s = oss.str();

    output_file = "Output/" + data_name + output_file;

    std::ofstream FILE(output_file, std::ios::app);
    FILE << s;
    //FILE << output_string;
    FILE.close();

}; // end of CheckCaSignal

void Astrocyte::SaveCaSignalsTime(float time, std::string output_file){

    if(trigger_c_signal){

        std::string output_string = "";
        output_string = output_string + std::to_string(time);

        output_file = "Output/" + output_file;
        std::ofstream FILE(output_file, std::ios::app);

        for(int i_c = 0; i_c < N_comparts; i_c++){

            if(c_signal[i_c] == 1){
                
                output_string = output_string + ";";
                output_string = output_string + std::to_string(i_c);
                
            }// end conditional for spiking

        };// End looping through compartments

        output_string = output_string + "\n";
        
        FILE << output_string;
        FILE.close();

    }

};// end of SaveCaSignalsTime

void Astrocyte::IncrementNeurotransmitter(std::vector <double> *Neut, std::vector <int *> *stimulated_units, int N_stimuli,
std::vector <int> *N_stimulated, std::vector <double> *s){

    for(int i_s = 0; i_s < N_stimuli; i_s++){
        for(int i_c = 0; i_c < (*N_stimulated)[i_s]; i_c++){

            (*Neut)[(*stimulated_units)[i_s][2*i_c + 0] - 1] = (*s)[(*stimulated_units)[i_s][2*i_c + 1] - 1];
                
        }
    }
}// end of IncrementNeurotransmitter

void Astrocyte::ChangeParameters(std::string parameter, std::string value){

    if (parameter == "N_comparts"){N_comparts = std::stoi(value);}
    else if(parameter == "N_comparts_connections"){N_comparts_connections = std::stoi(value);}
    else if(parameter == "astrocytic_network_file"){astrocytic_network_file = value;}
    else if(parameter == "ne_parameter_file"){ne_parameter_file = value;}
    else if(parameter == "N_ne_stim"){N_ne_stim = std::stoi(value);}
    else if(parameter == "N_ne_total"){N_ne_total = std::stoi(value);}
    else if(parameter == "da_parameter_file"){da_parameter_file = value;}
    else if(parameter == "N_da_stim"){N_da_stim = std::stoi(value);}
    else if(parameter == "N_da_total"){N_da_total = std::stoi(value);}
    else if(parameter == "F"){F = std::stod(value);}
    else if(parameter == "Ca_irest"){Ca_irest = std::stod(value);}
    else if(parameter == "Ca_ERrest"){Ca_ERrest = std::stod(value);}
    else if(parameter == "Ca_orest"){Ca_orest = std::stod(value);}
    else if(parameter == "Na_irest"){Na_irest = std::stod(value);}
    else if(parameter == "Na_orest"){Na_orest = std::stod(value);}
    else if(parameter == "K_irest"){K_irest = std::stod(value);}
    else if(parameter == "K_orest"){K_orest = std::stod(value);}
    else if(parameter == "IP3_rest"){IP3_rest = std::stod(value);}
    else if(parameter == "h_rest"){h_rest = std::stod(value);}
    else if(parameter == "v_rest"){v_rest = std::stod(value);}
    else if(parameter == "NErest"){NErest = std::stod(value);}
    else if(parameter == "DArest"){DArest = std::stod(value);}
    else if(parameter == "K_p"){K_p = std::stod(value);}
    else if(parameter == "K_pi"){K_pi = std::stod(value);}
    else if(parameter == "v_delta"){v_delta = std::stod(value);}
    else if(parameter == "kappa_delta"){kappa_delta = std::stod(value);}
    else if(parameter == "K_PLCdelta"){K_PLCdelta = std::stod(value);}
    else if(parameter == "v_3K"){v_3K = std::stod(value);}
    else if(parameter == "K_D"){K_D = std::stod(value);}
    else if(parameter == "K_3"){K_3 = std::stod(value);}
    else if(parameter == "r_5P"){r_5P = std::stod(value);}
    else if(parameter == "K_R"){K_R = std::stod(value);}
    else if(parameter == "v_beta"){v_beta = std::stod(value);}
    else if(parameter == "alpha"){alpha = std::stod(value);}
    else if (parameter == "v_DA") { v_DA = std::stod(value); }
    else if (parameter == "K_DA") { K_DA = std::stod(value); }
    else if(parameter == "r_L"){r_L = std::stod(value);}
    else if(parameter == "v_ER"){v_ER = std::stod(value);}
    else if(parameter == "K_ER"){K_ER = std::stod(value);}
    else if(parameter == "d_1"){d_1 = std::stod(value);}
    else if(parameter == "d_5"){d_5 = std::stod(value);}
    else if(parameter == "d_2"){d_2 = std::stod(value);}
    else if(parameter == "d_3"){d_3 = std::stod(value);}
    else if(parameter == "a_2"){a_2 = std::stod(value);}
    else if (parameter == "r_C_max") { r_C_max = std::stod(value); }
    else if (parameter == "r_C_min") { r_C_min = std::stod(value); }
    else if (parameter == "k1_R2C2") { k1_R2C2 = std::stod(value); }
    else if (parameter == "k2_R2C2") { k2_R2C2 = std::stod(value); }
    else if (parameter == "k1_R2C2cAMP2") { k1_R2C2cAMP2 = std::stod(value); }
    else if (parameter == "k2_R2C2cAMP2") { k2_R2C2cAMP2 = std::stod(value); }
    else if (parameter == "k3_R2C2cAMP2") { k3_R2C2cAMP2 = std::stod(value); }
    else if (parameter == "k1_R2C2cAMP4") { k1_R2C2cAMP4 = std::stod(value); }
    else if (parameter == "k2_R2C2cAMP4") { k2_R2C2cAMP4 = std::stod(value); }
    else if (parameter == "k3_R2C2cAMP4") { k3_R2C2cAMP4 = std::stod(value); }
    else if (parameter == "R2cAMP4") { R2cAMP4 = std::stod(value); }
    else if (parameter == "k1_PKA") { k1_PKA = std::stod(value); }
    else if (parameter == "k2_PKA") { k2_PKA = std::stod(value); }
    else if (parameter == "k_PKA") { k_PKA = std::stod(value); }
    else if (parameter == "k1_D1R_DA") { k1_D1R_DA = std::stod(value); }
    else if (parameter == "k3_R2C2cAMP4") { k3_R2C2cAMP4 = std::stod(value); }
    else if(parameter == "I_NKAmax"){I_NKAmax = std::stod(value);}
    else if(parameter == "K_NKAmN"){K_NKAmN = std::stod(value);}
    else if(parameter == "K_NKAmK"){K_NKAmK = std::stod(value);}
    else if(parameter == "I_NCXmax"){I_NCXmax = std::stod(value);}
    else if(parameter == "K_NCXmN"){K_NCXmN = std::stod(value);}
    else if(parameter == "K_NCXmC"){K_NCXmC = std::stod(value);}
    else if(parameter == "k_sat"){k_sat = std::stod(value);}
    else if(parameter == "eta"){eta = std::stod(value);}
    else if(parameter == "C_m"){C_m = std::stod(value);}
    else if(parameter == "g_Naleak"){g_Naleak = std::stod(value);}
    else if(parameter == "E_Na"){E_Na = std::stod(value);}
    else if(parameter == "g_Kleak"){g_Kleak = std::stod(value);}
    else if(parameter == "E_K"){E_K = std::stod(value);}
    else if(parameter == "D_Ca"){D_Ca = std::stod(value);}
    else if(parameter == "D_CaER"){D_CaER = std::stod(value);}
    else if(parameter == "D_IP3"){D_IP3 = std::stod(value);}
    else if(parameter == "D_Na"){D_Na = std::stod(value);}
    else if(parameter == "D_K"){D_K = std::stod(value);}
    else if(parameter == "D_Cao"){D_Cao = std::stod(value);}
    else if(parameter == "D_Nao"){D_Nao = std::stod(value);}
    else if(parameter == "D_Ko"){D_Ko = std::stod(value);}
    else if(parameter == "ca_thresh"){ca_thresh = std::stod(value);}
    else{
        std::cout << "Parameter " << parameter << " not recognized in the simulation parameters file for" 
        "astrocyte configuration." << std::endl;
        exit(EXIT_FAILURE);}

}; //end of ChangeParameters

void Astrocyte::TriggerInconsistencySizeWarning(std::string error, int expected, int configured){

    std::cout << "The size of the vector " << error << " does not match the configuration." << std::endl;
    std::cout << "It was expected " << expected << ", but it has size " << configured << "." << std::endl;
    exit(EXIT_FAILURE);

};// end of TriggerInconsistencySizeWarning 
