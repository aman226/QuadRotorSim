# QuadRotorSim

## Table of Contents

- [About](#about)
- [Getting Started](#getting_started)
- [Usage](#usage)
- [Contributing](../CONTRIBUTING.md)

## About <a name = "about"></a>
This is a quadrotor simulator that I am working on. It is a work in progress. I am using this project to learn about quadrotor dynamics and control. Right now, LQR (Linear Quadratic Regulator) control is implemented and is working. I am working on implementing MPC (Model Predictive Control) and other control methods. The current model is a simple point mass model.
## Getting Started <a name = "getting_started"></a>

These instructions will get you a copy of the project up and running on your local machine for development and testing. See [deployment](#deployment) for notes on how to deploy the project on a live system.

### Prerequisites

For this project, you need to have Armadillo Library installed. It is a state-of-the-art Linear Algebra that supports vectorization and parallelization.

### Installing

In order to run this, follow these steps.

```
git clone --recursive https://github.com/aman226/QuadRotorSim
cd QuadRotorSim
mkdir build
cd build
cmake ..
make
cd ../bin
/QuadRotorSim
```

## Usage <a name = "usage"></a>

Add notes about how to use the system.
