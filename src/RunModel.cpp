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
 * This file provides the code to execute the simulations configured with the text 
 * files in the Input folder.
 */

#include <chrono>
#include <string>
#include <iostream>
#include <filesystem>

#include "Simulator.h"

int main(){

    auto tic = std::chrono::high_resolution_clock::now();

    std::string file_name = "";    

    for (const auto & entry : std::filesystem::directory_iterator("Inputs/")) {

        auto tic_p = std::chrono::high_resolution_clock::now();

        file_name = entry.path().string();
            
        std::cout << file_name << std::endl;
        Simulator simulator(file_name);
        simulator.SimulateAstrocyticNetwork();

        auto toc_p = std::chrono::high_resolution_clock::now();
        auto duration_p = std::chrono::duration_cast<std::chrono::milliseconds>(toc_p - tic_p);

        std::ofstream FILE("Simulated", std::ios::app);
        FILE << file_name << ";" << float (duration_p.count())/1000 << "\n";
        FILE.close();

        
    }//end for looping through Inputs directory

    auto toc = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic);

    std::cout << "Execution time: " << float (duration.count())/1000 << " s" << std::endl;
    return 0;
}
