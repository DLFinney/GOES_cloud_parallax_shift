# GOES Cloud Parallax Shift

Code which adjusts GOES cloud fields for parallax shift, and regrids them.

Developed as part of the DCMEX project: https://cloudsense.ac.uk/dcmex/

## Purpose

These scripts are shared as inspiration for others coding similar problems, and for others to check what has been used for research work. The scripts weren't developed specifically for sharing so you may find they are not simply plug and play, but the steps have been made as clear as possible.

## Installation

### Requirements

Create a virtual environment and install the required dependencies:

```bash
pip install -r requirements.txt
```

### Dependencies

The main dependencies are:
- `numpy` - numerical computing
- `xarray` - multi-dimensional arrays
- `matplotlib` - plotting
- `goes2go` - GOES satellite data access
- `xesmf` - regridding
- `pandas` - data manipulation

## Usage

### Basic Workflow

1. **Download GOES data**: Use `download_goes_subregion_regrid.py` to download satellite data
2. **Apply parallax correction**: Use `parallax_latlons_goes_cloud.py` to correct for parallax shift

### Example

See the `example_input_files/` and `example_output_files/` directories for sample data and expected outputs.

## About Parallax Shift

In a satellite image, the apparent position of a cloud can have an error of several kilometres if the parallax effect is not accounted for. This effect arises because, when calculating the latitude and longitude coordinates of pixels in a satellite image, an assumption must be made for the height of that point above the centre of the earth. A simple assumption of mean sea level can introduce errors on the order of kilometres. 

Errors are largest for higher clouds towards the edge of the satellite image. Since DCMEX looks at deep convective clouds, and because New Mexico is towards the edge of the GOES 16 satellite image, the parallax error could not be missed when trying to relate aircraft tracks to cloud images from the satellite. Since DCMEX studies growing congestus with highly variable cloud heights, an adaptive approach for parallax correction was developed. This method individually corrected each pixel in each timestep by utilizing the GOES cloud-top height product.

## Alternative Solutions

**SatPy**: The [SatPy library](https://satpy.readthedocs.io/en/stable/api/satpy.modifiers.parallax.html) also provides parallax correction capabilities. See [Issue #1](https://github.com/DLFinney/GOES_cloud_parallax_shift/issues/1) for discussion on comparing approaches.

## Scripts

### `download_goes_subregion_regrid.py`

Downloads GOES satellite data for a specified region and time period. 

**Important**: The lat/lons calculated in this script assume the image is of mean sea level height. There is NO parallax adjustment in this script.

### `parallax_latlons_goes_cloud.py`

Applies parallax correction using GOES cloud-top height data and outputs both corrected coordinate files and regridded data.

## Contributing

This code was developed for research purposes. If you find errors or have suggestions for improvements, please open an issue or submit a pull request.

## Citation

If you use this code in your research, please cite the DCMEX project and acknowledge the original author.

