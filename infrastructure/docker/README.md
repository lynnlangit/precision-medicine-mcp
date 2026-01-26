# Docker Infrastructure for Precision Medicine Platform

This directory contains Docker-based infrastructure for the Precision Medicine platform, including containerized bioinformatics tools and analysis environments.

## Contents

### [base-images/](base-images/)
Docker base images for bioinformatics and machine learning workflows.

**Images Available:**
- **python-base** (~500MB) - Python 3.12 with bioinformatics libraries
- **r-base** (~2-3GB) - R with Bioconductor packages for genomics analysis
- **tensorflow-base** (~8-10GB GPU / ~4-5GB CPU) - TensorFlow with ML/DL tools

**Features:**
- Size-optimized images
- Non-root user execution
- Local-only (not pushed to registries)

**Quick Start:**
```bash
cd base-images
docker build -t precision-medicine/python-base:latest ./python-base
docker build -t precision-medicine/r-base:latest ./r-base
docker build -t precision-medicine/tensorflow-base:latest ./tensorflow-base
```

## Architecture

All Docker images follow these principles:

### Security
- Run as non-root user (`biouser`)
- Minimal attack surface (slim base images)
- No unnecessary tools or dependencies

### Size Optimization
- Multi-stage builds where applicable
- Aggressive layer cleanup
- Combined RUN statements
- Removal of build dependencies post-installation
- No cache directories (`--no-cache-dir`)

### Local Development
- **IMPORTANT:** All images are for LOCAL USE ONLY
- Do NOT push to Docker Hub or any remote registry
- Tagged with `precision-medicine/` prefix for organization
- Designed for development and testing environments

## Requirements

- **Docker Desktop** or **Docker Engine**
- **Disk Space:**
  - Minimum: 15GB free
  - Recommended: 30GB+ (for all images)
- **Memory:** 8GB+ RAM (16GB+ recommended for TensorFlow)

## Installation & Setup

### 1. Install Docker

**Windows:**
- Download from: https://www.docker.com/products/docker-desktop
- Install and restart
- Ensure WSL 2 backend is enabled (Settings → General)

**Linux:**
```bash
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
sudo usermod -aG docker $USER
```

**macOS:**
```bash
brew install --cask docker
```

### 2. Build Images

```bash
cd infrastructure/docker/base-images
docker build -t precision-medicine/python-base:latest ./python-base
docker build -t precision-medicine/r-base:latest ./r-base
docker build -t precision-medicine/tensorflow-base:latest ./tensorflow-base
```

## Usage Examples

### Running Containers

**Python Environment:**
```bash
docker run -it --rm -v "$(pwd):/data" precision-medicine/python-base:latest
```

**R Environment:**
```bash
docker run -it --rm -v "$(pwd):/data" precision-medicine/r-base:latest
```

**TensorFlow (GPU):**
```bash
docker run -it --rm --gpus all -v "$(pwd):/data" precision-medicine/tensorflow-base:latest
```

**TensorFlow (CPU):**
```bash
docker run -it --rm -v "$(pwd):/data" precision-medicine/tensorflow-base:latest
```

### Running Scripts

**Python Script:**
```bash
docker run --rm -v "$(pwd):/data" precision-medicine/python-base:latest python /data/script.py
```

**R Script:**
```bash
docker run --rm -v "$(pwd):/data" precision-medicine/r-base:latest Rscript /data/analysis.R
```

### Interactive Sessions

**Jupyter Notebook:**
```bash
docker run -it --rm -p 8888:8888 -v "$(pwd):/data" precision-medicine/python-base:latest \
  bash -c "pip install jupyter && jupyter notebook --ip=0.0.0.0 --allow-root --no-browser"
```

## Security Best Practices

### Regular Updates

1. **Rebuild periodically** when base images update
2. **Use vulnerability scanners** of your choice (Trivy, Grype, Snyk, etc.)
3. **Review scan reports** after each build

### Container Runtime Security

- Run as non-root user (built into images)
- Use read-only mounts when possible: `-v "$(pwd):/data:ro"`
- Limit resources: `--memory=4g --cpus=2`
- Drop capabilities: `--cap-drop=ALL`

## Troubleshooting

### Docker Issues

**Docker daemon not running:**
- Windows: Check system tray for Docker icon, restart Docker Desktop
- Linux: `sudo systemctl start docker`

**Permission denied:**
- Linux: Add user to docker group: `sudo usermod -aG docker $USER`
- Windows: Run Docker Desktop as Administrator

**Out of disk space:**
```bash
# Clean up unused containers, images, and volumes
docker system prune -a --volumes

# Check disk usage
docker system df
```

### Build Issues

**Build fails with network errors:**
- Check internet connection
- Configure Docker proxy if behind corporate firewall
- Try again (sometimes mirrors are temporarily down)

**Out of memory during build:**
- Increase Docker Desktop memory limit (Settings → Resources)
- Build images one at a time

**Slow builds:**
- First build downloads base images (one-time)
- Subsequent builds use cache (much faster)
- R packages compilation is slow (use pre-built binaries when possible)

## Contributing

When adding new Docker images or tools:

1. **Follow naming conventions** - Use descriptive, lowercase names with hyphens
2. **Document packages** - List all installed packages in README
3. **Optimize size** - Use slim bases, multi-stage builds, cleanup
4. **Security first** - Non-root users, minimal dependencies
5. **Test thoroughly** - Build and test before committing
6. **Update docs** - Add to this README and relevant documentation

## Additional Resources

- **[Base Images README](base-images/README.md)** - Detailed documentation
- **Docker Documentation:** https://docs.docker.com/

## License

These Docker configurations are part of the Precision Medicine Platform.
