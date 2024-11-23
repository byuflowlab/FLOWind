# FLOWFarm.jl

**Summary:** Wind farm simulation tool for gradient-based optimization.

**Authors:** Jared J. Thomas, Andrew P.J. Stanley

## Features
- Swap out models without changing anything else in the simulation setup
- Smooth/continous model implementations
- Runs on a single core, across multiple cores (threaded), or on multiple machines (distributed).
- Designed so that new model implementations can be included by adding a single method
- Allows for Wake Expansion Continuation (WEC) as described here

## Installation

### Install FLOWFarm
```
julia
(v1.x) pkg> dev https://github.com/byuflowlab/FLOWFarm.jl.git
```

## Testing

To test FLOWFarm, run the following from the top directory:

```
julia
julia
]
activate .
test
```

## Documentation

* Begin with the [quick start tutorial](Tutorial.md).
* More advanced topics are covered in the [how-to guide](How_to.md).
* Theory details, and links, can be found in the [theory](Explanation.md) page.
* Doc strings can be found in the [references](Reference.md) page.

## How FLOWFarm is structured
FLOWFarm was designed to be highly modular. In a wind farm simulation we are trying to predict the flow properties (wind speed, turbulence intensity, etc) at points in our space. The predictions at individual points can then be combined to predict the average inflow conditions at a wind turbine's rotor swept area, or to create a flow field across the whole wind farm. These predictions come from a combination of models: wake deficit models, wake deflection models, turbulence intensity models, wake combination models, and wind sheer models.

FLOWFarm uses these simple engineering-level models to predict the flow properties at each desired point. However, the published models for flows in wind farms typically combine various aspects of the modeling into a single model predicting wake deficit, wake deflection, turbulence intensity, and wind sheer in what seems to be a single cohesive model. The monolithic approach makes it difficult to research improvements to inidividual aspects of the modeling or to combine portions of different models. To solve the monolithic model problem we have broken up the existing models and implemented them in their individual components of deficit, deflection, sheer, turbulence, and wake combination. If you are interested in adding models to this framework you will need to identify each of the model components and implement them separately in the appropriate places within flowfarm (e.g. `wake_deficit_models.jl` and `wake_deflection_models.jl`). While it may take some effort to understand the model(s) of interest sufficiently to implement, this framework allows a wide range of studies once a set of models is included properly.

**Citing:**
Thomas, McOmber, and Ning "Wake Expansion Continuation: Multi-Modality Reduction in the Wind Farm Layout Optimization Problem" *Wind Energy* (in review), -->

