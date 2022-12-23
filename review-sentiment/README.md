# Customer Review Sentiment Analysis

| Build    | ![Review Sentiment](https://github.com/XAI-Demonstrator/xai-demonstrator/workflows/Review%20Sentiment/badge.svg?branch=master)                                                                                                     |
|----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Backend  | [![Coverage Status](https://coveralls.io/repos/github/XAI-Demonstrator/xai-demonstrator/badge.svg?branch=x-cov-sentiment-backend)](https://coveralls.io/github/XAI-Demonstrator/xai-demonstrator?branch=x-cov-sentiment-backend)   |
| Frontend | [![Coverage Status](https://coveralls.io/repos/github/XAI-Demonstrator/xai-demonstrator/badge.svg?branch=x-cov-sentiment-frontend)](https://coveralls.io/github/XAI-Demonstrator/xai-demonstrator?branch=x-cov-sentiment-frontend) |

## Components

- Multilingual sentiment analysis provided by [NLPTown's bert-base-multilingual-uncased-sentiment](https://huggingface.co/nlptown/bert-base-multilingual-uncased-sentiment) model
- Explanations generated via Facebook's [Captum](https://captum.ai/) library
- Backend microservice built with [FastAPI](https://fastapi.tiangolo.com/)
- Frontend created with [VueJS](https://vuejs.org/)

## Running the service locally (with Docker)

To launch a standalone instance of the service:

```shell
cd review-sentiment
docker-compose up
```

## Running the service locally (without Docker)

During development, it can be desirable to launch the service directly and not as a Docker container.
For instance, this allows to get rapid feedback on changes to the backend code.

```shell
cd review-sentiment
./build_frontend.sh
cd sentiment-backend
uvicorn sentiment.main:app
```

## GPU Acceleration with CUDA
- The sentiment service can be run with GPU acceleration using Pytorch's [CUDA package](https://pytorch.org/docs/stable/cuda.html).
- In order to enable GPU support use `docker build -f Dockerfile_CUDA .` to build the Docker image from CUDA_Dockerfile. Launch a container from this image with `docker run -p 80:8000 --gpus all *image ID*`. The flags are required to expose the container port 8000 (where the service is listening) to machine port 80 (where HTTP requests are coming in) and to grant the container access to the GPU.
- There are compatibility conditions which need to be satisfied for CUDA to work.
  1. CUDA requires the installation of NVIDIA drivers. The drivers have to be supported by the utilized CUDA version. An overview of the minimum driver required for a specific version of CUDA can be found [here](https://docs.nvidia.com/deploy/cuda-compatibility/index.html#binary-compatibility__table-toolkit-driver).
  2. The PyTorch version has to be compatible with the CUDA version. A list of previous PyTorch versions, their supported CUDA versions and download instructions can be found [here](https://pytorch.org/get-started/previous-versions/).
- The configuration has been tested using GCP's "GPU Optimized Debian m32 (with CUDA 10.0)" image (specifications: NVIDIA driver version 410.104, CUDA 10.0 and PyTorch 1.4.0).
