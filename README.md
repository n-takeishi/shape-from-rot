# Shape From Rotation

(c) Naoya Takeishi 2018

## Prerequisite

- lapack
- blas
- tbb
- boost 1.43 or later
- OpenCV 3.0.0 or later
- gtsam 4.0.0-alpha1 or later

## How to build

Use cmake with the CMakeLists.txt.

## How to execute

Let `ROOTDIR` be the relative path to the root directory, in which experiment's input (config/data) and output files are located. Basically, the main program can be run by

```
./main -d ROOTDIR
```

## What should be in the input directory

In `ROOTDIR`, input files must be prepared. At least, `config.txt` and `camera.txt` must be provided. Optionally, observations of spacecraft's movement can be also specified.

### config.txt

Configuration listed below must be specified.

- USE_SCMOVER         for using observations of spacecraft's rotation movement, set this to true (default: false)
- USE_SCMOVET         for using observations of spacecraft's translation movement, set this to true (default: false)
- CAMERA_SKEW         camera intrinsic parameter (pixel skewness, default: 0)
- CAMERA_FOCAL_X      camera intrinsic parameter (focal length in [px])
- CAMERA_FOCAL_Y      camera intrinsic parameter (focal length in [px])
- CAMERA_CENTER_X     camera intrinsic parameter (focal center in [px])
- CAMERA_CENTER_Y     camera intrinsic parameter (focal center in [px])
- NOISESTD_CAMERA     noise standard deviation of camera observations, i.e., landmarks, in [px]
- NOISESTD_SCMOVER_X  noise standard deviation of spacecraft's rotation movement in [rad]
- NOISESTD_SCMOVER_Y  noise standard deviation of spacecraft's rotation movement in [rad]
- NOISESTD_SCMOVER_Z  noise standard deviation of spacecraft's rotation movement in [rad]
- NOISESTD_SCMOVET_X  noise standard deviation of spacecraft's translation movement in [km]
- NOISESTD_SCMOVET_Y  noise standard deviation of spacecraft's translation movement in [km]
- NOISESTD_SCMOVET_Z  noise standard deviation of spacecraft's translation movement in [km]
- PRIORMEAN_XT_X      prior for initial position of spacecraft; if asteroid is located in the center of the initial image, set this to 0 [km]
- PRIORMEAN_XT_Y      prior for initial position of spacecraft; if asteroid is located in the center of the initial image, set this to 0 [km]
- PRIORMEAN_XT_Z      prior for initial position of spacecraft; set this according to the initial altitude of spacecraft in [km]

### camera.txt

This file should contain observations of camera, i.e., positions of landmarks in images. See `expt/toy/camera.txt` for example of the format.