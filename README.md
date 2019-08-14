# NSS

Neutron Star solver.

This code is a solver for the Tolman–Oppenheimer–Volkoff (TOV) equation and is designed to reproduce the density profile of neutron Stars. To do so a polytropic equation-of-state is assumed.

Visit [this website](http://caad-group.ddnss.de/neutronstar_new/) for an interactive version with preview of your data.

To use the code download the files "solver.cpp" as well as its config file "cf".
Then, compile "solver.cpp" and edit the config file:
```shell
`g++ -std=c++11 -o solver solver.cpp`
```
Afterwards, simply run `./solver`.

The contents of the config file are as follows:

| Parameter       | Value        | Option
| --------------- | ------------ | --------
|mode             |Multiple      | Single,Multiple
|method           |Rk4           | Euler,Rk2,Rk3,Rk4
|to_vary          |M             | M,R
|kappa            |1.98183e-6    | 
|gamma            |2.75          | 
|maximal_radius   |20.0          | 
|total_mass       |1.6           | 
|mass_resolution  |1e-3          | 
|radius_min       |16.0          | 
|radius_max       |20.0          | 
|mass_min         |1.1           | 
|mass_max         |2.2           | 
|nns              |100           | 
|                 |              |
|file_nameS       |Single.csv    |
|file_nameM       |Multiple.csv  |
|frac             |1e-4          |
|dr               |400.0         |
|stop             |1.0           |
|dr_div           |350.0         |
|max_steps        |1e5           |
|lower_boundary   |5.0e13        |
|upper_boundary   |5.0e17        |
|initial_guess    |3.0e14        |

It is recommended that the user only changes the parameters up to `frac`.
What follows is a short description of the parameters:

#### mode
If `Single` is selected than the code will simulate a single Neutron Star of given mass and radius and will save key parameters in `file_nameS`.
In this case all iteration steps are saved.
If `Multiple` is selected the code will simulate `nns` Neutron Stars. Depending on the option `to_vary` either the mass of the Neutron Star or the radius of the Neutron Star will vary.
Here only the central density on global parameters of the Neutron Stars are saved.

#### method
Here the method of simulation is selected. Currently the Euler-method as well as Runge-Kutta of second, third and fourth order are supported.

#### to_vary
If `M` is selected, the code will vary the mass. Likewise, if `R` is selected, the code will vary the radius of the Neutron Star.
The individual masses or radii of each simulated Neutron Star are linearly distributed over the minimal and maximal values given (i.e. `radius_min`, `radius_max` or `mass_min`, `mass_max`).

#### kappa
`kappa` is the proportionality constant between the density and the pressure in the equation-of-state.

#### gamma
`gamma` is the exponent of the density in the equation-of-state.

#### maximal_radius
This is the maximal radius in km assigned to the Neutron Star. It is best to leave this at a high value, otherwise th code will give unphysical results, if any.

#### total_mass
`total_mass` gives the mass of the Neutron Star in units of solar masses.

#### mass_resolution
This value determines when the simulation is stopped.
It compares the accumulated mass over the simulation to the total mass (`total_mass`) with 1.0.

#### nns
`nns` given the total number of Neutron Star to be simulated when option `Multiple` is selected.
Keep in mind that the computation time increases linearly with increasing values of `nns`!

#### file_nameS
This value gives the file name of the output file of simulation of a single Neutron Star.

#### file_nameM
This value gives the file name of the output file of simulation of multiple Neutron Stars.

