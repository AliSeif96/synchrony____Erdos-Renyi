# Measuring quantify the level of synchrony in a large population of neurons within an Erdos-Renyi network.

In these programs, we created an Erdos-Renyi network for 100 neurons with a probability of connecting 0.3. then We performed these neurons with a Huxkin Huxley model with an external current of 2 µA/cm^2 for 2,500 milliseconds and stored the voltage for the last 500 milliseconds. Then, to better understand the production Erdos-Renyi network, we drew   how neurons connect to each other in the network in three dimensions and binary. And after calculating the spike train, we calculated the synchronization rate.

### For measure the level of synchrony we must pass 5 steps.


Step number  | work done | language used | outputs
:-------------: | :-------------: | :-------------: | :-------------:
first step  | Making an Erdos-Renyi network with Hodgkin-Huxley model | C++| temp.txt - datas.txt - bineryprint.txt
second step  |3D visualization of Erdos-Renyi Network | python| 3D Network.html
Third step  | plotting how neurons connect to each other in the network in a binary way | python| Binary.html
fourth step  | Calculating the spike train and draw the raster plot | MATLAB| train.png - spike_Phi_I.txt
fifth step  | Calculating the overall synchronization of the network | C++| level of synchrony
