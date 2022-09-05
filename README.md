# DeepCNN_AAPredict

*Work in progress*

A deep CNN model to predict amino acid identity based on the chemical environment.

Chemical environment is represented as a box with densities of atom types surrounding the C-alpha
protein backbone. The box dimensions are 20 x 20 x 20. The atom types are 'N', 'C', 'O', and 'S'.
Use gen_res_centric_channels_v1.py to generate the chemical environment from pdb files.

Use res_centric_deeplearning_v1.py to train a NN model that uses the chemical environment as input
data.

The model is not final. The parameters need tuning and model architecture may need changes.

