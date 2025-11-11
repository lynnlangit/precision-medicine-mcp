"""Configuration settings for mcp-multiomics server."""

from pathlib import Path
from typing import Optional

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict


class MultiOmicsConfig(BaseSettings):
    """Configuration for multi-omics analysis server.

    Environment Variables:
        MULTIOMICS_DATA_DIR: Directory for multi-omics data files
        MULTIOMICS_CACHE_DIR: Directory for cached results
        MULTIOMICS_DRY_RUN: Enable mock mode without R dependencies
        MULTIOMICS_R_HOME: Path to R installation (optional)
        MULTIOMICS_LOG_LEVEL: Logging level (DEBUG, INFO, WARNING, ERROR)
        MULTIOMICS_MAX_FEATURES: Maximum features per modality
        MULTIOMICS_MIN_SAMPLES: Minimum samples required
    """

    model_config = SettingsConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        case_sensitive=False,
    )

    # Data directories
    data_dir: Path = Field(
        default=Path("/workspace/data/multiomics"),
        description="Directory for multi-omics data files",
    )

    cache_dir: Path = Field(
        default=Path("/workspace/cache/multiomics"),
        description="Directory for cached results and intermediate files",
    )

    # R configuration
    r_home: Optional[Path] = Field(
        default=None,
        description="Path to R installation (auto-detected if not set)",
    )

    dry_run: bool = Field(
        default=True,
        description="Enable mock mode without R dependencies",
    )

    # Analysis parameters
    max_features: int = Field(
        default=5000,
        description="Maximum number of features per modality",
        ge=10,
        le=100000,
    )

    min_samples: int = Field(
        default=3,
        description="Minimum number of samples required for analysis",
        ge=3,
    )

    # Statistical thresholds
    fdr_threshold: float = Field(
        default=0.05,
        description="FDR threshold for significance",
        ge=0.0,
        le=1.0,
    )

    stouffer_weights: bool = Field(
        default=True,
        description="Use weighted Stouffer's method (by sample size)",
    )

    # Logging
    log_level: str = Field(
        default="INFO",
        description="Logging level",
    )

    # Performance
    n_jobs: int = Field(
        default=-1,
        description="Number of parallel jobs (-1 for all CPUs)",
    )

    timeout_seconds: int = Field(
        default=600,
        description="Timeout for long-running operations",
        ge=60,
    )

    def __init__(self, **kwargs):
        """Initialize config and create directories if needed."""
        super().__init__(**kwargs)

        # Create directories if in real mode
        if not self.dry_run:
            self.data_dir.mkdir(parents=True, exist_ok=True)
            self.cache_dir.mkdir(parents=True, exist_ok=True)


# Global config instance
config = MultiOmicsConfig()
