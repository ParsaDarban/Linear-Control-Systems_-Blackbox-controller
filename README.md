# Linear-Control-Systems_-Blackbox-controller

In the analysis and design of controllers for systems, the first step is having an appropriate mathematical model of the system. This allows us to thoroughly examine the key characteristics of the real system through the model. For different systems, a specific model does not necessarily exist, and various models can be defined based on the understanding of the system and its application. Moreover, in some cases, due to insufficient knowledge in that domain or lack of access to system details, it is not possible to define a model based on the physics of the system. In such situations, the system must be identified using its input and output data.

In industry, a model for the system is not always available. Therefore, it is necessary to first identify the system using various methods based on the conditions and estimate a model, which can then be used for required analyses and designs. For a system, its model can be estimated by examining its impulse response, step response, response to sinusoidal inputs at different frequencies, and so on. Depending on the system, only some of these methods may be applicable.

In cases where a model is estimated solely based on input and output data, the resulting model is referred to as a "black-box" model of the system. This method has diverse applications, including in machine learning techniques that operate solely based on data. It is also used in control systems for system identification.

In this project, we are dealing with a system for which internal details are not accessible. Initially, we aim to derive a model for this system using our knowledge of linear systems in the time and frequency domains. Then, we will analyze this model, design an appropriate controller based on it, and finally, connect the controller to the original system and evaluate its performance.

You can find the complete information on data preprocessing, models, and their evaluation in the file LCS_CA.pdf, which is written in English.
