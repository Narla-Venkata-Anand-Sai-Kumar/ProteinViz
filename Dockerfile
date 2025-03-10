FROM python:3.8-slim

# Set the working directory
WORKDIR /app

# Install system dependencies for OpenCV and TensorFlow
RUN apt-get update && apt-get install -y \
    libgl1-mesa-glx \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*  

# Copy application files
COPY . /app

# Install Python dependencies
RUN pip3 install --no-cache-dir --default-timeout=1000 -r requirements.txt

#Epose the port ( user based )
EXPOSE 80

# Run the application
ENTRYPOINT ["gunicorn", "-w", "4", "-b", "0.0.0.0:80", "app:app"]
