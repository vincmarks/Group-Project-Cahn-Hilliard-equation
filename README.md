# Group-Project-Swift-Hohenberg model

This is a repo to work on the group project. 

 
For the Swift-Hohenberg model:
$\partial_t \phi=M \vec{\nabla}^2 \frac{\delta F}{\delta \phi}$
, $F=\int \mathrm{d}^d x \frac{1}{2} \phi\left(-\epsilon+\left(q_0^2+\vec{\nabla}^2\right)\right)^2 \phi+\frac{\phi^4}{4}$

We have found the following update formula:
$\hat{\phi}_{i+1,m} = \frac{\hat{\phi}_{i,m} - \Delta t  M  k_m^2  \sum_{j} (\phi_{i,j})^3 e^{-ik_m x_j} }{1 + \Delta t M  k_m^2 ((q_0^2 - k_m^2)^2 - \epsilon)}$
