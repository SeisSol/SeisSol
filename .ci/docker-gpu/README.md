## SeisSol test image
The final image contains a standard SeisSol software stack. Additionall, it
contains a GPU programming stack, namely: CUDA 10, ROCM 4.2, opensycl 0.9.5

### Building
```
docker-compose build
```

Ideally, CI pipeline should automatically rebuild the image if
the following directories are changed:

- `.ci/docker`
- `preprocessing/meshing/cube_c`
- `preprocessing/meshing/gmsh2gambit`
- `preprocessing/science/rconv`



### Running locally
Make sure that you have an Nvidia driver installed that supports CUDA 10, and [nvidia-docker2](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html).
```
docker run --rm -it --gpus all ravilmobile/seissol-base:gcc-8_cuda-10
```
