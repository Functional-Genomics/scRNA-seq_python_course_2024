# Course materials for Single Cell Analysis in Python 2024, EBI

Trainers: Iris Yu, Jiawei Wang, Anna Vathrakokili, Yuyao Song, Andrian Yang, Nadav Yayon

Contributor: Hugo Tavares


## Setting up the environment

To run the scripts for generating counts from raw fastq files `Demonstrations/01_*.sh`, please install CellRanger using [these instructions](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in).

To run all the other demo notebooks, follow the instructions below. It walks you through creating and using a container that has all the packages needed to run the demo notebooks.


1. Install [Docker](https://docs.docker.com/engine/install/) or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)  
	- Singularity is best for an HPC environment as a regular HPC user.  Check if your institute's HPC already has it.  
	- Docker might be easier to install in your local machine, depending on your machine.  
2. Pull the prebuilt training image that Andrian Yang built, [here](https://github.com/andr-kun/scRNAseq2024-container/pkgs/container/scrnaseq2024-container).  This image is based on the official scanpy image [gcfntnu:/scanpy:latest](https://hub.docker.com/layers/gcfntnu/scanpy/latest/images/sha256-368ed6c468c13d8f205a2831c7815777f2f51179f5ea1c1c78800f6b3e04c475?context=explore).
	- If using Signularity, do: 
	```
	singularity pull scrnaseq2024.sif docker://ghcr.io/andr-kun/scrnaseq2024-container:latest
	```
	where `scrnaseq2024.sif` is the name you would like to save the image file as.
	- If using Docker, do:
	```
	docker pull ghcr.io/andr-kun/scrnaseq2024-container:latest
	```  
	- Alternatively, build the image using [Andrian’s Dockerfile](https://github.com/andr-kun/scRNAseq2024-container/blob/main/Dockerfile).
3. Run Jupyter lab using the image you built, either using Docker or Singularity. If Singularity, you run it the same way we call it during the course sessions.
	- Singularity
	```
	singularity exec scrnaseq2024.sif jupyter lab
	```
	- Docker
	```
	docker run <image> jupyter lab
	```
	Note: Run the command while standing on this repo's root directory, or whichever directory you would like to be the root of Jupyter when it initialises.