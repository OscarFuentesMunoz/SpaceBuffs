# Installation
Run the following:  
`conda create -n theorigin -c conda-forge python=3.8 heyoka.py numpy jupyter matplotlib pykep`

The package also requires kepler equation solver by https://github.com/dfm/kepler.py  
Run the following:  
`python -m pip install kepler.py`

# About the package

Basic functions that will be helpful for the project are put together in the "origin" package. The package contains the following building blocks.

### Ephemeris
An `Ephemeris` object is defined to easily handle time history of a satellite's/debris' states. The time history of `Ephemeis` can be accessed by
* time: `Ephemeris.time`
* Cartesian state: `Ephemeris.get_cartesian()`
* Keplerian state: `Ephemeris.get_keplerian()`  

To build an ephemeris object, you can call
`Ephemeris(time:np.ndarray, state:np.ndarray, state_type:str)`

### Propagator
The propagator object can be constructed by `Propagator`. There are three methods defined for a `Propagator` object.
* `Propagator.time`: Set the initial time of the propagation job.
* `Propagator.state`: Set the initial state of the propagation job.
* `Propagator.cram`: Update the Cr(A/m) of the dynamics.

The propagation object can be called as `Propagator(time_array)`, which returns the `Ephemeris` object. See `validation.py` for example propagation.


### Data parser
Data parsers are defined for both debris and satellites. 
* Debris: `debris_parser(deb_id:int, is_training:bool)`  
Returns a dict with keys: `"deb_id", "path", "time", "ephemeris"`.  
If `is_training=True`, then it also has the following keys: `"sat_id", "cram", "event_ephemeris"`.
* Satellite: `satellite_parse(sat_id:int)`  
Returns a dict with keys: `"sat_id", "ephemeris"`
