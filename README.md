# ProteinViz

ProteinViz is a repository developed using Python and Flask for the 3D representation of proteins, facilitating advanced internal studies. This tool enables researchers and scientists to visualize protein structures in a three-dimensional space, aiding in the analysis and understanding of their properties and interactions.

## Features

- **3D Representation:** Utilizes advanced rendering techniques to generate three-dimensional representations of protein structures.
- **Interactive Interface:** Provides an interactive interface for users to explore and manipulate protein structures.
- **Customization:** Allows users to customize the visualization settings to suit their specific research needs.
- **Integration:** Seamlessly integrates with Python and Flask, providing a user-friendly web application for protein visualization.

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/Narla-Venkata-Anand-Sai-Kumar/ProteinViz.git
   ```

2. Navigate to the project directory:
   ```bash
   cd ProteinViz
   ```

3. Install the necessary dependencies:
   ```bash
   pip install -r requirements.txt
   ```

4. Run the Flask application:
   ```bash
   python app.py
   ```

5. Access the application through a web browser at `http://localhost:5000`.

## Alternative with docker 
# Dockerized Python Application with Gunicorn

This guide explains how to build and run the Python application inside a Docker container using **Gunicorn** as the application server.

## Prerequisites
- Docker installed on your system.
- A `Dockerfile` in the project directory.
- A `requirements.txt` file specifying the dependencies.

## Building the Docker Image
Run the following command in the project directory (where the `Dockerfile` is located):
```sh
docker build -t proteinviz .
```
This command:
- Uses the `Dockerfile` to create an image.
- Tags the image as `proteinviz`.

## Running the Docker Container
After building the image, you can run the container using:
```sh
docker run -p 80:80 proteinviz
```
This command:
- Maps port `80` of the container to port `80` of the host.
- Runs the container based on `proteinviz` image.

## Running in Detached Mode
To run the container in the background (detached mode):
```sh
docker run -d -p 80:80 proteinviz
```

open [http://localhost:80](http://localhost:80) in your browser.

## Contributions

Contributions to ProteinViz are welcome! Whether you want to enhance the visualization capabilities, improve the user interface, or add new features, feel free to fork the repository and submit a pull request.

## License

This project is licensed under the [specify your license type here, e.g., MIT License].

## Contact

For any questions or feedback regarding ProteinViz, please reach out to venkatnarla0@gmail.com.
