# Norepinephrine–dopamine sequential integration in astroglia drives behavioral state transitions 

Sijia Wang<sup>1,2†</sup>, Shilin Fang<sup>1†</sup>, Thiago O. Bezerra<sup>3†</sup>, Zhiyuan Wang<sup>1,2</sup>, Qimeng Zhao<sup>1,2</sup>, Xiaoou Wang<sup>1,2</sup>, Liang Wang<sup>4</sup>, Kui Wang<sup>1</sup>, Lijun Chen<sup>1,2</sup>, Funing Li<sup>1</sup>, Yuchen Gong<sup>1,2</sup>, Yongxin Yang<sup>1,2</sup>, Jiesi Feng<sup>5</sup>, Yu Zheng<sup>5</sup>, Yulong Li<sup>5</sup>, Jun Chu<sup>4</sup>, Gerald M. Pao<sup>6</sup>, Xufei Du<sup>1,2</sup>, Jiu-lin Du<sup>1,2</sup>, Antonio C. Roque<sup>3</sup>, Yu Mu<sup>1,2*</sup>
Affiliations:

<sup>1</sup> Institute of Neuroscience, State Key Laboratory of Brain Cognition and Brain-inspired Intelligence Technology, Center for Excellence in Brain Science and Intelligence Technology, Chinese Academy of Sciences, 320 Yue-Yang Road, Shanghai 200031, China.

<sup>2</sup> University of Chinese Academy of Sciences,19A Yu-Quan Road, Beijing 100049, China.

<sup>3</sup> Department of Physics, School of Philosophy, Sciences and Letters of Ribeirão Preto, University of São Paulo, Ribeirão Preto, São Paulo, Brazil.

<sup>4</sup> Research Center for Primate Neuromodulation and Neuroimaging, Shenzhen Institutes of Advanced Technology, Chinese Academy of Sciences, Shenzhen 518055, China.

<sup>5</sup> Academy for Advanced Interdisciplinary Studies, School of Psychological and Cognitive Sciences, Peking University, Beijing 100871, China.

<sup>6</sup> Okinawa Institute of Science and Technology Graduate University, Biological Nonlinear Dynamics Data Science Unit, Okinawa 904-0495, Japan.

<sup>†</sup> These authors contributed equally to this work.

<sup>*</sup> Corresponding author. Email: my@ion.ac.cn

The astrocyte cAMP model was implemented in C++. The codes are organized in four source files — `Astrocyte.cpp`, `Simulator.cpp`, `StimuliGenerator.cpp`, and `RunModel.cpp` and the corresponding header files (`Astrocyte.h`, `Simulator.h`, and `StimuliGenerator.h`).

- **Astrocyte.cpp:** Model equations, auxiliary methods for loading parameters, and methods for solving the system of ODEs;
- **StimuliGenerator.cpp:** Defines the class used to generate inputs;
- **Simulator.cpp:** Contains methods for running the simulation, interface for classes integration, and saving results;
- **RunModel.cpp:** Executes the simulations configured with the text files in the Input folder;

**Input** folder contains text files that configure the simulations. **Morphology** and **Networks** folders contain the astrocyte morphological model used in the simulations.
After compiling the code using the provided **Makefile** (for Windows), run `RunModel` to start the simulations. Simulation results are saved in the **Output** folder (create this folder if it does not exist in the working directory).

## Input File

Each input file defines the configuration for one simulation. Each line in these files consists of three space-separated columns:

| Column | Description                                         |
| :----- | :-------------------------------------------------- |
| 1      | target object, either Simulator (for simulation parameters) or Astrocyte (for model parameters) |
| 2      | parameter name |
| 3      | parameter value |

All available parameters are defined in the ChangeParameters methods in `Astrocyte.cpp` and `Simulator.cpp`.

The following entries are mandatory in every input file:

- **t_total**: Total simulation time (s);

- **dt**: Integration time step (s);

- **sample_rate**: Sampling frequency for saving data points;

- **ne_ext_parameter_file, ne_parameter_file, da_ext_parameter_file, da_parameter_file**: Text files defining noradrenergic and dopaminergic stimuli (see below);

- **N_ne_ext_stim, N_ne_stim, N_da_ext_stim, N_da_stim**: Number of noradrenergic and dopaminergic stimuli;

- **N_ne_ext_total, N_ne_total, N_da_ext_total, N_da_total**: Number of stimulated compartments;

- **astrocytic_network_file**: File describing the astrocyte morphological model;

- **N_comparts**: Number of compartments;

- **N_comparts_connections**: Number of connections between compartments;
- 
- **output_file**: Output file name; 

In this implementation:
N_ne_ext_stim = N_ne_stim, N_ne_ext_total = N_ne_total,
N_da_ext_stim = N_da_stim, and N_da_ext_total = N_da_total.

**Example:**
```
Simulator t_total 100
Simulator dt 0.01e-3
Simulator sample_rate 1000
Simulator output_file NADA_Example.txt
Simulator ne_ext_parameter_file NEStim20s
Simulator N_ne_ext_stim 1
Simulator N_ne_ext_total 1
Simulator da_ext_parameter_file DAStim0s
Simulator N_da_ext_stim 1
Simulator N_da_ext_total 1
Astrocyte astrocytic_network_file Lin1C_MorphFile.txt
Astrocyte N_comparts 1
Astrocyte ne_parameter_file NEStim20s
Astrocyte N_ne_stim 1
Astrocyte N_ne_total 1
Astrocyte da_parameter_file DAStim0s
Astrocyte N_da_stim 1
Astrocyte N_da_total 1
Astrocyte N_comparts_connections 0
```

## Astrocyte Network Files

Files in the Network folder specify the morphological model of each astrocyte. In this work, only a single astrocyte is simulated. Each file therefore contains a single line with seven semicolon-separated columns:

| Column | Description                                         |
| :----- | :-------------------------------------------------- |
| 1      | Astrocyte index                                     |
| 2      | Number of compartments                              |
| 3      | Filename of the morphology model                    |
| 4 - 7  | Unused in this version of the model                 |

**Example:**
```
1;1;Lin1C.txt;0;none;0;none
```

## Astrocyte Morphology Files

Astrocyte morphology files define the morphological parameters of each compartment. The somatic compartment is modeled as a sphere, while all other compartments are cylindrical. Each line contains nine comma-separated entries:

| Column | Description                                         |
| :----- | :-------------------------------------------------- |
| 1      | Compartment index                                   |
| 2      | Compartment type (1 = soma, 7 = astrocytic process) |
| 3      | Number of connections                               |
| 4      | Surface area                                        |
| 5      | Volume                                              |
| 6      | $r_{ER}$                                          |
| 7      | $r_m$                                             |
| 8      | Presence of mitochondria (1 = yes, 0 = no)          |
| 9      | Connected compartments (separated by `&`)           |


**Example:**
```
1,7,0,1.1890055529211375e-11,1.1250139346215622e-17,0.14962679563205702,0.0016080207548116077,0.0,-1
```

Since in this work we are simulating a single compartment, morphological files contain only one entry.

## Stimuli Files

Stimuli files specify the parameters of the noradrenergic and dopaminergic inputs. Each line configures one stimulus and contains seven semicolon-separated columns:

| Column | Description                                                                                                |
| :----- | :--------------------------------------------------------------------------------------------------------- |
| 1      | Input maxi amplitude (g_max)                                                                        |
| 2      | Input start time (t0)                                                                                     |
| 3      | Input rising time constant (t_rise)                                                                                |
| 4      | Input decay time constant (t_decay)                                                                       |
| 5      | Number of stimulated compartments |
| 6    | Ordered pairs `(a, b)` where `a` is the compartment index and `b` is the stimulus index (separated by `&`) |
| 7      | Stimulus indexes (separated by `&`)                                                                        |

These parameters configure the exponential function describing noradrenergic input and the alpha synapse simumlating dopaminergic input.

**Example (NE):**
```
0.079;20;1;7;1;(1,1);1
```

**Example (DA):**
```
0.77329358;20;7.36388329;1;1;(1,1);1
```
