"""
Cost Calculator for MCP Server Token Usage

Converts token usage to USD costs based on Anthropic API pricing.
Supports Claude Sonnet 4.5, Opus 4.5, and Haiku 4 models.
"""

from typing import Dict, Tuple
from dataclasses import dataclass


@dataclass
class ModelPricing:
    """Pricing per million tokens for a Claude model."""
    input_price: float  # USD per million input tokens
    output_price: float  # USD per million output tokens
    name: str


# Anthropic API Pricing (as of January 2025)
# Source: https://www.anthropic.com/pricing
PRICING = {
    "claude-sonnet-4-5-20250929": ModelPricing(
        input_price=3.00,
        output_price=15.00,
        name="Claude Sonnet 4.5"
    ),
    "claude-opus-4-5-20251101": ModelPricing(
        input_price=15.00,
        output_price=75.00,
        name="Claude Opus 4.5"
    ),
    "claude-haiku-4-20250924": ModelPricing(
        input_price=0.25,
        output_price=1.25,
        name="Claude Haiku 4"
    ),
    # Aliases for convenience
    "sonnet": ModelPricing(
        input_price=3.00,
        output_price=15.00,
        name="Claude Sonnet 4.5"
    ),
    "opus": ModelPricing(
        input_price=15.00,
        output_price=75.00,
        name="Claude Opus 4.5"
    ),
    "haiku": ModelPricing(
        input_price=0.25,
        output_price=1.25,
        name="Claude Haiku 4"
    ),
}

# Default model used by MCP servers
DEFAULT_MODEL = "claude-sonnet-4-5-20250929"


def calculate_cost(
    input_tokens: int,
    output_tokens: int,
    model: str = DEFAULT_MODEL
) -> Dict[str, float]:
    """
    Calculate cost in USD for given token usage.

    Args:
        input_tokens: Number of input tokens
        output_tokens: Number of output tokens
        model: Model name (default: claude-sonnet-4-5-20250929)

    Returns:
        Dictionary with:
        - input_cost: Cost of input tokens (USD)
        - output_cost: Cost of output tokens (USD)
        - total_cost: Total cost (USD)
        - model_name: Human-readable model name

    Example:
        >>> calculate_cost(10000, 2000)
        {
            'input_cost': 0.03,
            'output_cost': 0.03,
            'total_cost': 0.06,
            'model_name': 'Claude Sonnet 4.5'
        }
    """
    # Get pricing for model (default to Sonnet if not found)
    pricing = PRICING.get(model, PRICING[DEFAULT_MODEL])

    # Calculate costs (price is per million tokens)
    input_cost = (input_tokens / 1_000_000) * pricing.input_price
    output_cost = (output_tokens / 1_000_000) * pricing.output_price
    total_cost = input_cost + output_cost

    return {
        "input_cost": round(input_cost, 6),
        "output_cost": round(output_cost, 6),
        "total_cost": round(total_cost, 6),
        "model_name": pricing.name,
    }


def calculate_cost_per_insight(
    total_cost: float,
    insights_generated: int
) -> float:
    """
    Calculate cost per analytical insight.

    Args:
        total_cost: Total USD cost
        insights_generated: Number of insights generated (e.g., DEG analyses, pathway reports)

    Returns:
        Cost per insight in USD

    Example:
        >>> calculate_cost_per_insight(0.50, 5)
        0.1  # $0.10 per insight
    """
    if insights_generated == 0:
        return 0.0
    return round(total_cost / insights_generated, 4)


def estimate_monthly_cost(
    daily_tool_calls: int,
    avg_input_tokens: int,
    avg_output_tokens: int,
    model: str = DEFAULT_MODEL,
    days_per_month: int = 30
) -> Dict[str, float]:
    """
    Estimate monthly costs based on usage patterns.

    Args:
        daily_tool_calls: Average number of tool calls per day
        avg_input_tokens: Average input tokens per call
        avg_output_tokens: Average output tokens per call
        model: Model name
        days_per_month: Days per month (default: 30)

    Returns:
        Dictionary with monthly projections

    Example:
        >>> estimate_monthly_cost(100, 5000, 1000)
        {
            'daily_cost': 0.21,
            'monthly_cost': 6.30,
            'annual_cost': 75.60,
            'total_calls_per_month': 3000
        }
    """
    # Calculate single call cost
    single_call = calculate_cost(avg_input_tokens, avg_output_tokens, model)

    # Project costs
    daily_cost = single_call["total_cost"] * daily_tool_calls
    monthly_cost = daily_cost * days_per_month
    annual_cost = monthly_cost * 12

    return {
        "daily_cost": round(daily_cost, 2),
        "monthly_cost": round(monthly_cost, 2),
        "annual_cost": round(annual_cost, 2),
        "total_calls_per_month": daily_tool_calls * days_per_month,
        "model_name": single_call["model_name"],
    }


def compare_model_costs(
    input_tokens: int,
    output_tokens: int
) -> Dict[str, Dict[str, float]]:
    """
    Compare costs across all available Claude models.

    Args:
        input_tokens: Number of input tokens
        output_tokens: Number of output tokens

    Returns:
        Dictionary mapping model names to cost breakdowns

    Example:
        >>> compare_model_costs(10000, 2000)
        {
            'Claude Sonnet 4.5': {'total_cost': 0.06, ...},
            'Claude Opus 4.5': {'total_cost': 0.30, ...},
            'Claude Haiku 4': {'total_cost': 0.0050, ...}
        }
    """
    results = {}

    # Compare main models only (not aliases)
    for model_id in ["claude-sonnet-4-5-20250929", "claude-opus-4-5-20251101", "claude-haiku-4-20250924"]:
        costs = calculate_cost(input_tokens, output_tokens, model_id)
        results[costs["model_name"]] = costs

    return results


def get_optimization_recommendations(
    server_metrics: Dict[str, Dict]
) -> list[str]:
    """
    Generate cost optimization recommendations based on usage patterns.

    Args:
        server_metrics: Dictionary of server metrics with token usage

    Returns:
        List of optimization recommendations

    Example:
        >>> get_optimization_recommendations({
        ...     'mcp-multiomics': {'total_cost': 5.00, 'avg_output_tokens': 10000},
        ...     'mcp-spatialtools': {'total_cost': 0.50, 'avg_output_tokens': 1000}
        ... })
        [
            'Consider using Haiku for mcp-spatialtools (low token usage)',
            'Review mcp-multiomics output verbosity (high output tokens)'
        ]
    """
    recommendations = []

    for server_name, metrics in server_metrics.items():
        # High cost servers
        if metrics.get("total_cost", 0) > 2.0:
            recommendations.append(
                f"⚠️ {server_name}: High cost (${metrics['total_cost']:.2f}). "
                f"Review tool complexity or consider caching frequent queries."
            )

        # High output token servers (verbose responses)
        if metrics.get("avg_output_tokens", 0) > 5000:
            recommendations.append(
                f"📝 {server_name}: High output tokens ({metrics['avg_output_tokens']:,}). "
                f"Consider more concise tool responses or use Haiku for simple queries."
            )

        # Low usage servers that could use Haiku
        if (metrics.get("total_cost", 0) < 0.10 and
            metrics.get("avg_output_tokens", 0) < 1000):
            recommendations.append(
                f"💡 {server_name}: Low token usage. Consider switching to Haiku "
                f"to reduce costs by ~92% (Haiku: $0.25/$1.25 vs Sonnet: $3/$15)."
            )

    # General recommendations
    if not recommendations:
        recommendations.append(
            "✅ Token usage is optimized. Continue monitoring for changes."
        )

    return recommendations


def format_cost_summary(cost_data: Dict[str, float]) -> str:
    """
    Format cost data as human-readable summary.

    Args:
        cost_data: Cost breakdown from calculate_cost()

    Returns:
        Formatted string summary

    Example:
        >>> format_cost_summary({'input_cost': 0.03, 'output_cost': 0.03, 'total_cost': 0.06})
        'Total: $0.060 (Input: $0.030 | Output: $0.030)'
    """
    return (
        f"Total: ${cost_data['total_cost']:.3f} "
        f"(Input: ${cost_data['input_cost']:.3f} | "
        f"Output: ${cost_data['output_cost']:.3f})"
    )


# Token estimation helpers for common operations
OPERATION_ESTIMATES = {
    "differential_expression": {
        "input_tokens": 8000,  # Gene expression matrix + parameters
        "output_tokens": 2000,  # DEG results + statistics
        "description": "Differential expression analysis"
    },
    "spatial_autocorrelation": {
        "input_tokens": 12000,  # Spatial coordinates + expression data
        "output_tokens": 3000,  # Moran's I results + interpretation
        "description": "Spatial autocorrelation (Moran's I)"
    },
    "pathway_enrichment": {
        "input_tokens": 5000,  # Gene list + pathway database query
        "output_tokens": 4000,  # Enriched pathways + statistics
        "description": "Pathway enrichment analysis"
    },
    "cell_deconvolution": {
        "input_tokens": 10000,  # Expression matrix + signatures
        "output_tokens": 2500,  # Cell type scores by region
        "description": "Cell type deconvolution"
    },
    "fhir_query": {
        "input_tokens": 3000,  # Patient ID + query parameters
        "output_tokens": 1500,  # Clinical data (demographics, meds, labs)
        "description": "FHIR clinical data query"
    },
    "multiomics_integration": {
        "input_tokens": 15000,  # Multiple omics layers
        "output_tokens": 5000,  # Integrated results + correlations
        "description": "Multi-omics data integration"
    },
}


def estimate_operation_cost(operation: str, model: str = DEFAULT_MODEL) -> Dict[str, any]:
    """
    Estimate cost for a common MCP operation.

    Args:
        operation: Operation name (e.g., 'differential_expression')
        model: Model name

    Returns:
        Dictionary with token estimates and costs

    Example:
        >>> estimate_operation_cost('differential_expression')
        {
            'operation': 'differential_expression',
            'description': 'Differential expression analysis',
            'input_tokens': 8000,
            'output_tokens': 2000,
            'total_cost': 0.054,
            ...
        }
    """
    if operation not in OPERATION_ESTIMATES:
        raise ValueError(f"Unknown operation: {operation}. Available: {list(OPERATION_ESTIMATES.keys())}")

    estimates = OPERATION_ESTIMATES[operation]
    costs = calculate_cost(estimates["input_tokens"], estimates["output_tokens"], model)

    return {
        "operation": operation,
        "description": estimates["description"],
        "input_tokens": estimates["input_tokens"],
        "output_tokens": estimates["output_tokens"],
        **costs,
    }


if __name__ == "__main__":
    # Demo usage
    print("=== MCP Cost Calculator Demo ===\n")

    # Example 1: Single tool call
    print("Example 1: Differential Expression Analysis")
    de_cost = estimate_operation_cost("differential_expression")
    print(f"  Tokens: {de_cost['input_tokens']:,} in / {de_cost['output_tokens']:,} out")
    print(f"  Cost: ${de_cost['total_cost']:.4f}\n")

    # Example 2: PatientOne complete workflow (10 tool calls)
    print("Example 2: PatientOne Complete Workflow (10 tools)")
    total_input = 80000
    total_output = 25000
    workflow_cost = calculate_cost(total_input, total_output)
    print(f"  Tokens: {total_input:,} in / {total_output:,} out")
    print(f"  {format_cost_summary(workflow_cost)}\n")

    # Example 3: Model comparison
    print("Example 3: Model Cost Comparison (10K in / 2K out)")
    comparison = compare_model_costs(10000, 2000)
    for model_name, costs in comparison.items():
        print(f"  {model_name}: ${costs['total_cost']:.4f}")
    print()

    # Example 4: Monthly projection
    print("Example 4: Monthly Cost Projection (50 analyses/day)")
    monthly = estimate_monthly_cost(50, 8000, 2000)
    print(f"  Daily: ${monthly['daily_cost']:.2f}")
    print(f"  Monthly: ${monthly['monthly_cost']:.2f}")
    print(f"  Annual: ${monthly['annual_cost']:.2f}")
