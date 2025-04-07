[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/x_pJNjA_)
# ECE231_Assignment8
# Lab Assignment 8

Create a class called "Field" that contains a `double *value`.  Modify Lab Assignment 7 such that the Electric and Magnetic Field classes created now inherit from the base class Field.  Create a default constructor in the base class Field that sets the components to 0.  Create another constructor that can assign the values when the variable is created, e.g.

```
Electric Field(0,1e5,1e3);
```
This constructor should call the base class constructor to assign the component values. 

In the base class Field, make a function called printMagnitude that prints the components of the electric and magnetic fields.

# Class specific functions
In electric field class create a private variable that stores a value for a calculated electric field. Write a member method/function that calculates the electric field.  The electric field can be computed via Gauss' Law given by

```math
\int \mathbf{E} \cdot d\mathbf{S} = \frac{Q}{\epsilon_0}
```

where $\mathbf{E}$ is the electric field, $Q$ is the total charge, and $\epsilon_0$ is the permittivity of free space.  For a spherical charge, this equation is solved and found to be 

```math
\text{E} = \frac{Q}{4 \pi r^2 \epsilon_0}
```

Write a function that calculates the electric field at a given distance and stores it in the private variable, $r$, away from the spherical charge and total spherical charge $Q$.  Charge is in units of Coulomb's.  

For the magnetic field create a private variable that stores the value of a calculated magnetic field. 

From Ampere's law, the magnetic field can be found via

```math
\int \mathbf{B} \cdot dl = \mu_0 I
```

where $\mathbf{B}$ is the magnetic field, $I$ is the total current flowing through a surface, and $\mu_0$ is the permeability of free space.  For a current-carrying wire, the magnetic field is found to be

```math
\text{B} = \frac{I}{2 \pi r \mu_0}
```

Write a function that solves for the magnetic field at a given distance, $r$, away from the wire and for a given current, $I$ and stores it in the private variable you created.  The current is in units of Amps.

# Operator overload
Create two overloaded operators, one that overloads the '+' operator for two electric or magnetic fields (you will need to implement the overloaded operator in each derived class).  The overloaded operator will enable the functionality of

```
Electric_Field e1, e2, e3;
e1 = e2 + e3;
```
You will need to implement a copy constructor as well for this to work properly.  Finally, overload the '<<' operator such that the functionality of
```
Electric_Field e1;
std::cout << e1 << std::endl;
```
The '<<' operator should output the components of e1 in the above example.  You will need to implement the '<<' overloaded operator in each derived class.

# Code Usage
Create an electric and magnetic field object and use the constructor to set values for each three components.  Use the printMagnitude function to print the three components of each.  Next, use made-up values (you can come up with them) to calculate the electric and magnetic fields based upon Gauss' law and Ampere's law.  Print the resulting values to the screen.  Demonstrate the use of the overloaded operator functions in both derived classes.

# Grading
- 20% base class implemented and inheritance
- 20% class specific functions implemented
- 20% + operator overloaded
- 20% << operator overloaded
- 20% expected code output
