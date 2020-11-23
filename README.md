# Overlap Save CUDA


#### Technologies
Project was created and tested with:
* Windows 10
* Visual Studio 2017
* CUDA 10.1.243


#### Description
Project used to test capabilities of CUDA. This program compares speed of convolution two signals via overlap-save algorithm between CPU and GPU.


#### Setup
- Download and install CUDA from following link:
https://developer.nvidia.com/cuda-10.1-download-archive-update2?target_os=Windows&target_arch=x86_64&target_version=10&target_type=exelocal 
- Prepare Samples.csv file that contain samples of (audio) signal (put each channel in separate column, and set delimiter as comma (,)) 
- Place Samples.csv file in overlap_save_cuda\overlap_save_cuda\ catalogue
- Set parameters in "DANE DO EDYCJI" section


#### Run
- Run overlap_save_cuda\overlap_save_cuda.sln