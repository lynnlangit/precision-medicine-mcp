"""Cost tracking and monitoring utilities for precision medicine workflows.

This module provides utilities to track, estimate, and report costs for:
- External API calls (TCGA, HuggingFace, Seqera)
- Cloud compute resources (AWS Batch, Azure, GCP)
- Storage costs (data caching, intermediate files)
- AI model inference costs

Use these tools to:
1. Estimate costs before running analysis
2. Track actual costs during execution
3. Generate cost reports per patient/workflow
4. Set budget alerts and limits
"""

import json
import logging
import os
import time
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Callable
from functools import wraps
from contextlib import contextmanager

logger = logging.getLogger(__name__)


# ============================================================================
# Cost Pricing Data (Updated 2025)
# ============================================================================

PRICING = {
    # Cloud Compute (per hour)
    "compute": {
        "aws_batch_small": 0.096,      # t3.large (2 vCPU, 8GB RAM)
        "aws_batch_medium": 0.192,     # t3.xlarge (4 vCPU, 16GB RAM)
        "aws_batch_large": 0.384,      # t3.2xlarge (8 vCPU, 32GB RAM)
        "azure_standard_d4": 0.192,    # 4 vCPU, 16GB RAM
        "gcp_n1_standard_4": 0.190,    # 4 vCPU, 15GB RAM
    },

    # Storage (per GB-month)
    "storage": {
        "aws_s3_standard": 0.023,
        "aws_s3_ia": 0.0125,           # Infrequent Access
        "azure_blob_hot": 0.018,
        "gcp_standard": 0.020,
    },

    # Data Transfer (per GB)
    "data_transfer": {
        "aws_internet_out": 0.09,
        "aws_internet_in": 0.00,
        "azure_internet_out": 0.087,
        "gcp_internet_out": 0.12,
    },

    # External APIs (estimated)
    "api": {
        "tcga_query": 0.00,            # Free (NCI funded)
        "tcga_download_gb": 0.00,      # Free
        "hf_inference_small": 0.06,    # Per 1M tokens (~100M param model)
        "hf_inference_large": 0.60,    # Per 1M tokens (~1B param model)
        "seqera_platform_hour": 0.02,  # Platform fee per compute hour
    },

    # Common Analysis Costs (estimated)
    "analysis": {
        "fastq_qc_per_gb": 0.01,       # FastQC, quality trimming
        "rna_seq_per_sample": 0.50,    # Alignment + quantification
        "variant_calling_per_sample": 1.20,  # BWA + GATK
        "multiomics_integration": 0.25,  # Per sample across modalities
        "spatial_analysis_per_slide": 2.00,  # 10X Visium analysis
    }
}


@dataclass
class CostItem:
    """Individual cost item record."""
    timestamp: str
    category: str  # "compute", "storage", "api", "analysis"
    operation: str
    amount_usd: float
    quantity: float = 0.0
    unit: str = ""
    metadata: Dict[str, Any] = None

    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}


@dataclass
class CostSummary:
    """Cost summary for a workflow or patient."""
    total_usd: float
    by_category: Dict[str, float]
    items: List[CostItem]
    start_time: str
    end_time: str
    duration_seconds: float
    patient_id: Optional[str] = None
    workflow_id: Optional[str] = None


class CostTracker:
    """Track costs for precision medicine workflows.

    Usage:
        tracker = CostTracker(patient_id="patient_001")

        # Track individual costs
        tracker.add_cost("compute", "rna_seq_alignment", 1.20)
        tracker.add_cost("storage", "cache_results", 0.05)

        # Get summary
        summary = tracker.get_summary()
        print(f"Total cost: ${summary.total_usd:.2f}")
    """

    def __init__(
        self,
        patient_id: Optional[str] = None,
        workflow_id: Optional[str] = None,
        budget_limit_usd: Optional[float] = None,
        cost_log_path: Optional[Path] = None
    ):
        """Initialize cost tracker.

        Args:
            patient_id: Patient identifier for cost attribution
            workflow_id: Workflow identifier
            budget_limit_usd: Optional budget limit (triggers warnings)
            cost_log_path: Path to save cost logs (default: ./cost_logs/)
        """
        self.patient_id = patient_id
        self.workflow_id = workflow_id
        self.budget_limit_usd = budget_limit_usd

        self.items: List[CostItem] = []
        self.start_time = datetime.now().isoformat()

        # Setup cost logging
        if cost_log_path is None:
            cost_log_path = Path.cwd() / "cost_logs"
        self.cost_log_path = Path(cost_log_path)
        self.cost_log_path.mkdir(parents=True, exist_ok=True)

        logger.info(f"CostTracker initialized for patient={patient_id}, workflow={workflow_id}")
        if budget_limit_usd:
            logger.info(f"Budget limit: ${budget_limit_usd:.2f}")

    def add_cost(
        self,
        category: str,
        operation: str,
        amount_usd: float,
        quantity: float = 0.0,
        unit: str = "",
        metadata: Optional[Dict[str, Any]] = None
    ):
        """Add a cost item.

        Args:
            category: Cost category ("compute", "storage", "api", "analysis")
            operation: Operation name (e.g., "rna_seq_alignment")
            amount_usd: Cost in USD
            quantity: Optional quantity (e.g., 10.5 for 10.5 GB)
            unit: Optional unit (e.g., "GB", "hours", "samples")
            metadata: Optional additional metadata
        """
        item = CostItem(
            timestamp=datetime.now().isoformat(),
            category=category,
            operation=operation,
            amount_usd=amount_usd,
            quantity=quantity,
            unit=unit,
            metadata=metadata or {}
        )

        self.items.append(item)

        # Log cost
        logger.info(
            f"💰 Cost added: {operation} = ${amount_usd:.4f} "
            f"({quantity} {unit})" if unit else f"💰 Cost added: {operation} = ${amount_usd:.4f}"
        )

        # Check budget
        current_total = self.get_total()
        if self.budget_limit_usd and current_total > self.budget_limit_usd:
            logger.warning(
                f"⚠️  BUDGET EXCEEDED: ${current_total:.2f} > ${self.budget_limit_usd:.2f}"
            )

    def get_total(self) -> float:
        """Get total cost so far."""
        return sum(item.amount_usd for item in self.items)

    def get_by_category(self) -> Dict[str, float]:
        """Get costs grouped by category."""
        by_category = {}
        for item in self.items:
            by_category[item.category] = by_category.get(item.category, 0.0) + item.amount_usd
        return by_category

    def get_summary(self) -> CostSummary:
        """Get cost summary."""
        end_time = datetime.now().isoformat()
        start_dt = datetime.fromisoformat(self.start_time)
        end_dt = datetime.fromisoformat(end_time)
        duration = (end_dt - start_dt).total_seconds()

        return CostSummary(
            total_usd=self.get_total(),
            by_category=self.get_by_category(),
            items=self.items,
            start_time=self.start_time,
            end_time=end_time,
            duration_seconds=duration,
            patient_id=self.patient_id,
            workflow_id=self.workflow_id
        )

    def save_log(self, filename: Optional[str] = None):
        """Save cost log to file.

        Args:
            filename: Optional filename (default: auto-generated)
        """
        if filename is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            patient_prefix = f"{self.patient_id}_" if self.patient_id else ""
            filename = f"{patient_prefix}cost_log_{timestamp}.json"

        log_path = self.cost_log_path / filename

        summary = self.get_summary()

        log_data = {
            "summary": {
                "total_usd": summary.total_usd,
                "by_category": summary.by_category,
                "start_time": summary.start_time,
                "end_time": summary.end_time,
                "duration_seconds": summary.duration_seconds,
                "patient_id": summary.patient_id,
                "workflow_id": summary.workflow_id,
            },
            "items": [asdict(item) for item in self.items]
        }

        with open(log_path, 'w') as f:
            json.dump(log_data, f, indent=2)

        logger.info(f"Cost log saved to: {log_path}")
        return log_path

    def print_summary(self):
        """Print formatted cost summary."""
        summary = self.get_summary()

        print("\n" + "=" * 80)
        print("COST SUMMARY")
        print("=" * 80)

        if summary.patient_id:
            print(f"Patient: {summary.patient_id}")
        if summary.workflow_id:
            print(f"Workflow: {summary.workflow_id}")

        print(f"Duration: {summary.duration_seconds:.1f} seconds")
        print(f"\nTotal Cost: ${summary.total_usd:.2f}")

        if self.budget_limit_usd:
            remaining = self.budget_limit_usd - summary.total_usd
            pct_used = (summary.total_usd / self.budget_limit_usd) * 100
            print(f"Budget: ${self.budget_limit_usd:.2f} ({pct_used:.1f}% used)")
            print(f"Remaining: ${remaining:.2f}")

        print("\nCosts by Category:")
        for category, amount in sorted(summary.by_category.items()):
            pct = (amount / summary.total_usd) * 100 if summary.total_usd > 0 else 0
            print(f"  {category:20s} ${amount:8.2f} ({pct:5.1f}%)")

        print("\nTop 5 Operations by Cost:")
        sorted_items = sorted(self.items, key=lambda x: x.amount_usd, reverse=True)
        for item in sorted_items[:5]:
            print(f"  {item.operation:30s} ${item.amount_usd:8.4f}")

        print("=" * 80 + "\n")


@contextmanager
def track_cost(
    tracker: CostTracker,
    category: str,
    operation: str,
    cost_estimator: Optional[Callable[[], float]] = None
):
    """Context manager for tracking operation costs.

    Usage:
        tracker = CostTracker(patient_id="patient_001")

        with track_cost(tracker, "compute", "rna_seq_alignment"):
            # Run expensive operation
            run_rna_seq_pipeline()

        # Cost automatically tracked
    """
    start_time = time.time()

    logger.info(f"Starting operation: {operation}")

    try:
        yield
    finally:
        duration = time.time() - start_time

        # Estimate cost if estimator provided
        if cost_estimator:
            cost = cost_estimator()
        else:
            # Default: assume some cost based on duration
            cost = duration * 0.01  # $0.01 per second (placeholder)

        tracker.add_cost(
            category=category,
            operation=operation,
            amount_usd=cost,
            quantity=duration,
            unit="seconds"
        )

        logger.info(f"Completed operation: {operation} (${cost:.4f}, {duration:.1f}s)")


def track_operation_cost(
    category: str,
    operation: str = None,
    cost_estimator: Optional[Callable[..., float]] = None
):
    """Decorator to track operation costs.

    Usage:
        @track_operation_cost("compute", cost_estimator=lambda: 1.20)
        def run_variant_calling(sample_id: str):
            # Expensive operation
            pass

    Note: Requires a CostTracker instance to be passed as 'cost_tracker' kwarg
    or set as global tracker.
    """
    def decorator(func: Callable) -> Callable:
        op_name = operation or func.__name__

        @wraps(func)
        def sync_wrapper(*args, **kwargs):
            # Get tracker from kwargs or global
            tracker = kwargs.pop('cost_tracker', None)
            if tracker is None:
                tracker = globals().get('_global_cost_tracker')

            if tracker is None:
                # No tracker available, just run function
                logger.warning(
                    f"No cost tracker available for {op_name}, costs not tracked"
                )
                return func(*args, **kwargs)

            start_time = time.time()

            try:
                result = func(*args, **kwargs)
                return result
            finally:
                duration = time.time() - start_time

                # Estimate cost
                if cost_estimator:
                    cost = cost_estimator(*args, **kwargs)
                else:
                    cost = duration * 0.01

                tracker.add_cost(
                    category=category,
                    operation=op_name,
                    amount_usd=cost,
                    quantity=duration,
                    unit="seconds"
                )

        @wraps(func)
        async def async_wrapper(*args, **kwargs):
            tracker = kwargs.pop('cost_tracker', None)
            if tracker is None:
                tracker = globals().get('_global_cost_tracker')

            if tracker is None:
                logger.warning(
                    f"No cost tracker available for {op_name}, costs not tracked"
                )
                return await func(*args, **kwargs)

            start_time = time.time()

            try:
                result = await func(*args, **kwargs)
                return result
            finally:
                duration = time.time() - start_time

                if cost_estimator:
                    cost = cost_estimator(*args, **kwargs)
                else:
                    cost = duration * 0.01

                tracker.add_cost(
                    category=category,
                    operation=op_name,
                    amount_usd=cost,
                    quantity=duration,
                    unit="seconds"
                )

        # Detect if function is async
        import asyncio
        if asyncio.iscoroutinefunction(func):
            return async_wrapper
        else:
            return sync_wrapper

    return decorator


# ============================================================================
# Cost Estimators
# ============================================================================

class CostEstimator:
    """Estimate costs for common operations before execution."""

    @staticmethod
    def compute_cost(
        instance_type: str,
        duration_hours: float
    ) -> float:
        """Estimate cloud compute cost.

        Args:
            instance_type: Instance type (e.g., "aws_batch_medium")
            duration_hours: Estimated duration in hours

        Returns:
            Estimated cost in USD
        """
        hourly_rate = PRICING["compute"].get(instance_type, 0.20)
        return hourly_rate * duration_hours

    @staticmethod
    def storage_cost(
        size_gb: float,
        duration_days: float = 30,
        storage_class: str = "aws_s3_standard"
    ) -> float:
        """Estimate storage cost.

        Args:
            size_gb: Data size in GB
            duration_days: Storage duration in days (default: 30)
            storage_class: Storage class (e.g., "aws_s3_standard")

        Returns:
            Estimated cost in USD
        """
        monthly_rate = PRICING["storage"].get(storage_class, 0.023)
        months = duration_days / 30.0
        return size_gb * monthly_rate * months

    @staticmethod
    def rna_seq_cost(num_samples: int) -> float:
        """Estimate RNA-seq analysis cost.

        Args:
            num_samples: Number of samples

        Returns:
            Estimated cost in USD
        """
        per_sample = PRICING["analysis"]["rna_seq_per_sample"]
        return num_samples * per_sample

    @staticmethod
    def variant_calling_cost(num_samples: int) -> float:
        """Estimate variant calling cost.

        Args:
            num_samples: Number of samples

        Returns:
            Estimated cost in USD
        """
        per_sample = PRICING["analysis"]["variant_calling_per_sample"]
        return num_samples * per_sample

    @staticmethod
    def multiomics_integration_cost(
        num_samples: int,
        num_modalities: int = 3
    ) -> float:
        """Estimate multi-omics integration cost.

        Args:
            num_samples: Number of samples
            num_modalities: Number of omics modalities (RNA, protein, phospho, etc.)

        Returns:
            Estimated cost in USD
        """
        base_cost = PRICING["analysis"]["multiomics_integration"]
        # Cost scales with modalities
        return num_samples * base_cost * (num_modalities / 3.0)

    @staticmethod
    def spatial_analysis_cost(num_slides: int) -> float:
        """Estimate spatial transcriptomics analysis cost.

        Args:
            num_slides: Number of tissue slides

        Returns:
            Estimated cost in USD
        """
        per_slide = PRICING["analysis"]["spatial_analysis_per_slide"]
        return num_slides * per_slide

    @staticmethod
    def patient_one_workflow_cost() -> Dict[str, float]:
        """Estimate full PatientOne workflow cost.

        Returns:
            Dictionary with cost breakdown
        """
        costs = {
            "genomic_qc": 0.15,  # FASTQ QC + validation
            "variant_calling": 1.20,  # BWA + GATK
            "rna_seq": 0.50,  # Alignment + quantification
            "multiomics_integration": 0.75,  # 3 modalities
            "spatial_analysis": 2.00,  # 1 slide
            "imaging_analysis": 1.50,  # DeepCell + processing
            "tcga_comparison": 0.10,  # TCGA data queries
            "compute_overhead": 0.80,  # Misc compute
            "storage": 0.20,  # 30-day caching
        }

        costs["total"] = sum(v for k, v in costs.items() if k != "total")

        return costs


# ============================================================================
# Budget Alerts
# ============================================================================

class BudgetAlert:
    """Monitor costs and trigger alerts when thresholds are exceeded."""

    def __init__(
        self,
        tracker: CostTracker,
        thresholds: List[float],
        alert_callback: Optional[Callable[[float, float], None]] = None
    ):
        """Initialize budget alert monitor.

        Args:
            tracker: CostTracker instance to monitor
            thresholds: List of threshold amounts (USD) that trigger alerts
            alert_callback: Optional callback function(current_cost, threshold)
        """
        self.tracker = tracker
        self.thresholds = sorted(thresholds)
        self.alert_callback = alert_callback or self._default_alert
        self.triggered_thresholds = set()

    def check(self):
        """Check current cost against thresholds."""
        current_cost = self.tracker.get_total()

        for threshold in self.thresholds:
            if (current_cost >= threshold and
                threshold not in self.triggered_thresholds):

                self.triggered_thresholds.add(threshold)
                self.alert_callback(current_cost, threshold)

    def _default_alert(self, current_cost: float, threshold: float):
        """Default alert handler."""
        logger.warning(
            f"🚨 BUDGET ALERT: Cost ${current_cost:.2f} exceeded threshold ${threshold:.2f}"
        )
        print(
            f"\n🚨 BUDGET ALERT 🚨\n"
            f"Current cost: ${current_cost:.2f}\n"
            f"Threshold: ${threshold:.2f}\n"
        )


# ============================================================================
# Example Usage
# ============================================================================

if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    print("=" * 80)
    print("COST TRACKING EXAMPLES")
    print("=" * 80 + "\n")

    # Example 1: Basic cost tracking
    print("Example 1: Basic Cost Tracking\n")

    tracker = CostTracker(
        patient_id="patient_001",
        workflow_id="ovarian_cancer_analysis",
        budget_limit_usd=50.0
    )

    # Track various operations
    tracker.add_cost("analysis", "fastq_qc", 0.15, quantity=10, unit="GB")
    tracker.add_cost("analysis", "variant_calling", 1.20, quantity=1, unit="sample")
    tracker.add_cost("analysis", "rna_seq", 0.50, quantity=1, unit="sample")
    tracker.add_cost("analysis", "multiomics_integration", 0.75, quantity=3, unit="modalities")
    tracker.add_cost("analysis", "spatial_analysis", 2.00, quantity=1, unit="slide")
    tracker.add_cost("compute", "cloud_compute", 0.80, quantity=4.2, unit="hours")
    tracker.add_cost("storage", "data_caching", 0.20, quantity=8.7, unit="GB-month")

    tracker.print_summary()

    # Save log
    log_path = tracker.save_log()
    print(f"Cost log saved to: {log_path}\n")

    # Example 2: Cost estimation
    print("\nExample 2: Cost Estimation\n")

    estimator = CostEstimator()

    print("PatientOne Workflow Estimated Costs:")
    patient_one_costs = estimator.patient_one_workflow_cost()
    for component, cost in patient_one_costs.items():
        print(f"  {component:30s} ${cost:6.2f}")

    print(f"\n  {'TOTAL':30s} ${patient_one_costs['total']:6.2f}")

    # Example 3: Budget alerts
    print("\n\nExample 3: Budget Alerts\n")

    alert_tracker = CostTracker(patient_id="patient_002")

    alerts = BudgetAlert(
        tracker=alert_tracker,
        thresholds=[5.0, 10.0, 15.0]
    )

    # Simulate operations
    alert_tracker.add_cost("analysis", "operation_1", 3.0)
    alerts.check()

    alert_tracker.add_cost("analysis", "operation_2", 4.0)  # Triggers $5 threshold
    alerts.check()

    alert_tracker.add_cost("analysis", "operation_3", 5.0)  # Triggers $10 threshold
    alerts.check()

    print("\n" + "=" * 80)
