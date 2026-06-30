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
 * This header defines the Astrocyte class used to build the astrocyte model.
 */
#pragma once

#include <memory>
#include <vector>

class Astrocyte{

    public:

        // Astrocyte Morphology
        int N_comparts;
        int N_comparts_connections;

        std::vector<int> n_c_connections;     // n = compartments
        std::vector<int *> c_connections;      // n = compartments
        std::vector<int> _c_connections;      // n = total compartments connections;

        // Morphological parameters
        std::vector<double> A;
        std::vector<double> Vol;
        std::vector<double> r_ER;
        
        // State variables
        std::vector<double> Ca_i;
        std::vector<double> Ca_ER;

        std::vector<double> Ca_i_prev;
        std::vector<double> IP3;
        std::vector<double> h;
        std::vector<double> Ca_o;
        std::vector<double> Na_i;
        std::vector<double> Na_o;
        std::vector<double> K_i;
        std::vector<double> K_o;
        std::vector<double> V;


        std::vector<double> AC;
        std::vector<double> AC_Gas;
        std::vector<double> D1R;
        std::vector<double> D1R_DA;
        std::vector<double> D1R_DA_Gas;
        std::vector<double> Gas;
        std::vector<double> GasA;
        std::vector<double> AC_Gas_ATP;
        std::vector<double> cAMP;
        std::vector<double> R2C2;
        std::vector<double> R2C2cAMP2;
        std::vector<double> R2C2cAMP4;
        std::vector<double> PKA;
        std::vector<double> PDE4;
        std::vector<double> PDE4A;

        std::vector<int> c_signal;
        bool trigger_c_signal = false;

        // Neurotransmitter variables
        std::vector<double> NE;
        std::vector<double> DA;

        // Simulation configuration attributes
        std::string parameter_file;
        std::string astrocytic_network_file;
        std::string ne_parameter_file;
        std::string da_parameter_file;

        // --- Astrocyte model parameters --- //

        // Physical Constants
        double T = 303.16;              // K
        double F = 96500;               // C/mole
        double R = 8.314;               // J/(mole * K)


        // Resting State
        double Ca_irest  = 0.000069;       // mmole/m**3
        double Ca_ERrest = 2.335e-2;       // mmole/m**3
        double Ca_orest  = 1.8;            // mmole/m**3
        double Na_irest  = 15;             // mmole/m**3
        double Na_orest  = 145;            // mmole/m**3
        double K_irest   = 100;            // mmole/m**3
        double K_orest   = 3;              // mmole/m**3
        double IP3_rest  = 0;              // mmole/m**3
        double h_rest    = 0.8029;         // 1
        double v_rest    = -85e-3;         // volt
        double NErest    = 0;              // mmole/m**3
        double DArest    = 10e-3;          // mmole/m**3

        double AC_rest         = 2.342;     // mmole/m**3
        double AC_Gas_rest     = 2.206e-5;  // mmole/m**3
        double D1R_rest        = 0.618e-1;  // mmole/m**3
        double D1R_DA_rest     = 0.206e-3;  // mmole/m**3
        double D1R_DA_Gas_rest = 3.743e-05; // mmole/m**3
        double Gas_rest        = 0.728e-2;  // mmole/m**3
        double GasA_rest       = 7.487e-05; // mmole/m**3
        double AC_Gas_ATP_rest = 1.867e-14; // mmole/m**3
        double cAMP_rest       = 0;         // mmole/m**3
        double R2C2_rest       = 3.590;     // mmole/m**3
        double R2C2cAMP2_rest  = 0;         // mmole/m**3
        double R2C2cAMP4_rest  = 0;         // mmole/m**3
        double PKA_rest        = 0;         // mmole/m**3
        double PDE4_rest       = 1; 
        double PDE4A_rest      = 0; 


        // IP3 Dynamics
        // PLC B Inhibition
        double K_p  = 0.461e-1; // mmole/m**3
        double K_pi = 0.6e-3;   // mmole/m**3
        double K_R = 0.599e-2;  // mM
        double v_beta = 8e-05;  // mM * s
        double alpha = 0.7;

        double v_DA = 0;      // mM * s
        double K_DA = 0.203;  // mM * s
        double K_pDA = 0.405; // mM * s

        double beta = 0.5;

        // IP3 PLC Delta Synthesis
        double v_delta = 0.025e-3;  // mM * s
        double kappa_delta = 1.5e-3;    // mM
        double K_PLCdelta = 0.1e-3;    // mM

        // IP3-3K Degradation
        double v_3K = 2e-3;      // mM / s
        double K_D = 0.7e-3;    // mM
        double K_3 = 1e-3;      // mM

        // IP-5P Degradation
        double r_5P = 0.45;            // 1/s

        // Ca ER Leak Current
        double r_L = 0.11;             // 1/s


        // SERCA Current
        double v_ER = 0.008;            // mM * s
        double K_ER = 0.1e-3;           // mM


        // h Dynamics
        double d_1 = 0.13e-3;           // mM
        double d_5 = 0.10351e-3;        // mM
        double d_2 = 1.049e-3;          // mM
        double d_3 = 0.9434e-3;         // mM
        double a_2 = 0.2;               // m**3/(mmole * s)
        double r_C_max = 5;                 // 1/s
        double r_C_min = 5;                 // 1/s

        double k_PKA = 5e-7;

        // NKA
        double I_NKAmax = 1.52;         // amp / metre ** 2
        double K_NKAmN = 10;            // mM
        double K_NKAmK = 1.5;           // mM

        // NCX
        double I_NCXmax = 0.0001;       // amp/metre**2
        double K_NCXmN = 87.5;          // mM
        double K_NCXmC = 1.380;         // mM
        double k_sat = 0.1;             // 1
        double eta = 0.35;              // 1


        // Voltage Parameter
        double C_m = 1.0e-2;            //F/m**2 


        // Leakage Currents
        double g_Naleak = 13.482807180178595;     // S/m**2
        double E_Na = 61e-3;         // V
        double g_Kleak = 145.81406273541853;    // S/m**2
        double E_K = -94e-3;        // V


        // cAMP-PKA Pathway
        double k1_ACGas = 38.5;  // 1/uMs
        double k2_ACGas = 50;    // 1/s
        double k3_ACGas = 0.128; // 1/uMs
        double k4_ACGas = 56.92; // 1/s
        double ATP = 2000;       // uM

        double k1_GasA = 20;  // 1/s
        double k2_GasA = 10;  // 1/s

        double k1_D1R_DA_Gas = 10;     // 1/uMs
        double k2_D1R_DA_Gas = 0.001;  // 1/s
        double k3_D1R_DA_Gas = 20;     // 1/s

        double k1_Gas = 10;            // 1/s
        double k2_Gas = 10;            // 1/uMs

        double k1_D1R_DA = 3.33;       // 1/uMs
        double k2_D1R_DA = 10;         // 1/s
        double k3_D1R_DA = 10;         // 1/s

        double k1_AC_Gas_ATP = 0.128;  // 1/uMs
        double k2_AC_Gas_ATP = 56.92;  // 1/s
        double k3_AC_Gas_ATP = 0.259;  // 1/uMs

        double k1_cAMP = 56.92;      // 1/s
        double k2_cAMP = 0.259;      // 1/uMs
        double kf_R2C2 = 0.035;      // 1/uM2 s
        double kf_R2C2cAMP2 = 0.27;  // 1/uM2 s
        double kb_R2C2 = 0.06;       // 1/s
        double kb_R2C2cAMP4 = 0.28;  // 1/s

        double V_PDE = 8.25;         // uM/s
        double kcat = 10;            // uM

        double V_PDE4 = 5;           // 1/s
        double kcat4 = 5e-5;         // 1
        double kcat5 = 0.3;          // uM

        double k1_R2C2 = 0.035;      // 1/uM2 s
        double k2_R2C2 = 0.06;       // 1/s

        double k1_R2C2cAMP2 = 0.035; // 1/uM2 s
        double k2_R2C2cAMP2 = 0.28;  // 1/s 
        double k3_R2C2cAMP2 = 0.27;  // 1/uM2 s

        double k1_R2C2cAMP4 = 0.27; // 1/uM2 s
        double k2_R2C2cAMP4 = 0.28; // 1/s
        double R2cAMP4 = 0.69;      // uM

        double k1_PKA = 1;        // 1/s
        double k2_PKA = 1;        // 1/uM2 s
        double k_pep = 0;         // 1/s
        double k_dgr = 0.25;      // 1/s
        double k_degr = 0.1;      // 1/uM

        double k1_PDE4A = 11;     // 1/uM s
        double k2_PDE4A = 0.01;   // 1/s
        
        
        // Diffusion Constants
        double D_Ca   = 0.2;          // 1/s
        double D_CaER = 0.001;        // 1/s
        double D_IP3  = 0.2;          // 1/s
        double D_Na   = 0.316;        // 1/s
        double D_K    = 0.938;        // 1/s
        double D_Cao  = 4.52;         // 1/s
        double D_Nao  = 26.6;         // 1/s
        double D_Ko   = 1.732;        // 1/s

        // Ca Signal Threshold
        double ca_thresh = 0.15e-3;
 
        int N_ne_stim;
        int N_ne_total;
        std::vector<int *> ne_stim_comparts;
        std::vector<int> _ne_stim_comparts;
        std::vector<int> N_stimulated_ne;

        int N_da_stim;
        int N_da_total; 
        std::vector<int *> da_stim_comparts;
        std::vector<int> _da_stim_comparts;
        std::vector<int> N_stimulated_da;
        
        Astrocyte(std::string parameter_file);

        // Support Functions for Morphology

        void ConfigureAstrocyticVectors();
        void LoadAstrocyteModel(std::string astrocyte_file);
        void LoadMorphologicalParams(std::string morphology_file);

        // Model equations
        double f_prod_PLCb_NE(double NE, double Ca_i);
        double f_prod_PLCb_DA(double DA_o, double Ca_i);
        double f_prod_PLCd(double IP3, double Ca_i);
        double f_degr_IP3_3K(double Ca_i, double IP3);
        double f_degr_IP_5P(double IP3);
        double f_J_IP3R (double A, double Vol, double IP3, double Ca_i, double h, double Ca_ER, double PKA);
        double f_J_CERleak(double A, double Vol, double Ca_ER, double Ca_i);
        double f_J_SERCA(double A, double Vol, double Ca_i);
        double f_J_NKA(double Na_i, double K_o);
        double f_J_NCX(double Na_o, double Ca_o, double Na_i, double V, double Ca_i);
        double f_J_Naleak(double V);
        double f_J_Kleak(double V);

        double Diffusion(int i_c, double D, std::vector<double> *var);

        // -------- ODEs -------- //

        double dCa_odt(int i_c);
        double dNa_odt(int i_c);
        double dK_odt(int i_c);
        double dCa_idt(int i_c);
        double dCa_ERdt(int i_c);
        double dhdt(int i_c);
        double dIP3dt(int i_c);
        double dNa_idt(int i_c);
        double dK_idt(int i_c);
        double dVdt(int i_c);

        double dACGasdt(int i_c);
        double dGasAdt(int i_c);
        double dD1R_DA_Gasdt(int i_c);
        double dGasdt(int i_c);
        double dD1R_DAdt(int i_c);
        double dAC_Gas_ATPdt(int i_c);
        double dcAMPdt(int i_c);
        double dR2C2dt(int i_c);
        double dR2C2cAMP2dt(int i_c);
        double dR2C2cAMP4dt(int i_c);
        double dPKAdt(int i_c);
        double dPDE4Adt(int i_c);

        void UpdateVariables(double dt);

        // Support Functions

        void ConfigureAstrocyticStimulators(std::string stimulator_parameter_file, std::vector <int *> *stim_comparts, 
        std::vector <int> *_stim_comparts, std::vector <int> *N_stimulated);
        void CheckCaSignalInAstrocyte(int i_c);
        void SaveData(float time, std::vector <double> *data, std::string data_name, std::string output_file);
        void SaveCaSignalsTime(float time, std::string output_file);
        void ChangeParameters(std::string parameter, std::string value);
        void TriggerInconsistencySizeWarning(std::string error, int expected, int configured);

        void IncrementNeurotransmitter(std::vector <double> *Neut, std::vector <int *> *stimulated_units, int N_stimuli,
        std::vector <int> *N_stimulated, std::vector <double> *s);


};// end of class Astrocyte