"""
Metrics Aggregator for MCP Server Observability

Loads and processes metrics from workflow logs (YAML format).
Calculates costs, aggregates by server/tool, and generates insights.
"""

import yaml
from pathlib import Path
from typing import Dict, List, Any
from datetime import datetime
import pandas as pd

from cost_calculator import calculate_cost, calculate_cost_per_insight, get_optimization_recommendations


class MetricsAggregator:
    """Aggregates and analyzes MCP server metrics."""

    def __init__(self, metrics_file: str):
        """
        Initialize aggregator with metrics file.

        Args:
            metrics_file: Path to YAML metrics file
        """
        self.metrics_file = Path(metrics_file)
        self.raw_data = None
        self.tool_calls_df = None
        self.server_summary_df = None
        self.cost_breakdown = None

        self._load_metrics()
        self._process_metrics()

    def _load_metrics(self):
        """Load raw metrics from YAML file."""
        with open(self.metrics_file, 'r') as f:
            self.raw_data = yaml.safe_load(f)

    def _process_metrics(self):
        """Process raw metrics into DataFrames and cost breakdowns."""
        # Convert tool calls to DataFrame
        tool_calls = self.raw_data.get('tool_calls', [])
        self.tool_calls_df = pd.DataFrame(tool_calls)

        # Add cost calculations to each tool call
        if not self.tool_calls_df.empty:
            self.tool_calls_df['cost'] = self.tool_calls_df.apply(
                lambda row: calculate_cost(
                    row['input_tokens'],
                    row['output_tokens'],
                    self.raw_data['workflow_metadata'].get('model', 'claude-sonnet-4-5-20250929')
                )['total_cost'],
                axis=1
            )

            # Parse timestamp if present
            if 'timestamp' in self.tool_calls_df.columns:
                self.tool_calls_df['timestamp'] = pd.to_datetime(self.tool_calls_df['timestamp'])

        # Create server summary
        self._create_server_summary()

        # Calculate cost breakdown
        self._calculate_cost_breakdown()

    def _create_server_summary(self):
        """Aggregate metrics by server."""
        if self.tool_calls_df.empty:
            self.server_summary_df = pd.DataFrame()
            return

        # Group by server
        server_groups = self.tool_calls_df.groupby('server').agg({
            'tool': 'count',
            'input_tokens': 'sum',
            'output_tokens': 'sum',
            'latency_ms': ['sum', 'mean'],
            'cost': 'sum'
        }).reset_index()

        # Flatten column names
        server_groups.columns = [
            'server', 'tool_calls', 'input_tokens', 'output_tokens',
            'total_latency_ms', 'avg_latency_ms', 'total_cost'
        ]

        # Add status from raw data
        server_metrics = self.raw_data.get('server_metrics', {})
        server_groups['status'] = server_groups['server'].apply(
            lambda s: server_metrics.get(s, {}).get('status', 'unknown')
        )

        # Calculate average tokens per call
        server_groups['avg_input_tokens'] = (
            server_groups['input_tokens'] / server_groups['tool_calls']
        ).round(0).astype(int)

        server_groups['avg_output_tokens'] = (
            server_groups['output_tokens'] / server_groups['tool_calls']
        ).round(0).astype(int)

        # Sort by total cost descending
        server_groups = server_groups.sort_values('total_cost', ascending=False)

        self.server_summary_df = server_groups

    def _calculate_cost_breakdown(self):
        """Calculate detailed cost breakdown."""
        if self.tool_calls_df.empty:
            self.cost_breakdown = {}
            return

        model = self.raw_data['workflow_metadata'].get('model', 'claude-sonnet-4-5-20250929')

        total_input = self.tool_calls_df['input_tokens'].sum()
        total_output = self.tool_calls_df['output_tokens'].sum()
        costs = calculate_cost(total_input, total_output, model)

        # Get insights count
        insights = self.raw_data.get('insights_generated', [])
        cost_per_insight_val = calculate_cost_per_insight(
            costs['total_cost'],
            len(insights)
        ) if insights else 0

        self.cost_breakdown = {
            'total_input_tokens': total_input,
            'total_output_tokens': total_output,
            'total_tokens': total_input + total_output,
            'input_cost': costs['input_cost'],
            'output_cost': costs['output_cost'],
            'total_cost': costs['total_cost'],
            'model': costs['model_name'],
            'cost_per_tool_call': costs['total_cost'] / len(self.tool_calls_df) if len(self.tool_calls_df) > 0 else 0,
            'cost_per_insight': cost_per_insight_val,
            'insights_count': len(insights),
        }

    def get_workflow_metadata(self) -> Dict[str, Any]:
        """Get workflow metadata."""
        return self.raw_data.get('workflow_metadata', {})

    def get_tool_calls(self) -> pd.DataFrame:
        """Get all tool calls as DataFrame."""
        return self.tool_calls_df.copy()

    def get_server_summary(self) -> pd.DataFrame:
        """Get per-server summary as DataFrame."""
        return self.server_summary_df.copy()

    def get_cost_breakdown(self) -> Dict[str, float]:
        """Get cost breakdown."""
        return self.cost_breakdown.copy()

    def get_summary_stats(self) -> Dict[str, Any]:
        """Get high-level summary statistics."""
        summary = self.raw_data.get('summary', {})
        workflow = self.get_workflow_metadata()

        return {
            'workflow_name': workflow.get('name', 'Unknown'),
            'patient_id': workflow.get('patient_id', 'Unknown'),
            'analysis_date': workflow.get('analysis_date', 'Unknown'),
            'total_duration_seconds': workflow.get('total_duration_seconds', 0),
            'total_tool_calls': summary.get('total_tool_calls', 0),
            'successful_calls': summary.get('successful_calls', 0),
            'dry_run_calls': summary.get('dry_run_calls', 0),
            'servers_used': summary.get('servers_used', 0),
            **self.cost_breakdown
        }

    def get_top_cost_tools(self, n: int = 10) -> pd.DataFrame:
        """
        Get top N most expensive tool calls.

        Args:
            n: Number of top tools to return

        Returns:
            DataFrame with top N tools by cost
        """
        if self.tool_calls_df.empty:
            return pd.DataFrame()

        return self.tool_calls_df.nlargest(n, 'cost')[
            ['server', 'tool', 'input_tokens', 'output_tokens', 'cost', 'latency_ms']
        ].copy()

    def get_insights_summary(self) -> pd.DataFrame:
        """Get insights generated with metadata."""
        insights = self.raw_data.get('insights_generated', [])
        if not insights:
            return pd.DataFrame()

        df = pd.DataFrame(insights)

        # Calculate cost per insight
        if self.cost_breakdown:
            cost_per_insight = self.cost_breakdown['cost_per_insight']
            df['estimated_cost'] = cost_per_insight

        return df

    def get_optimization_recommendations(self) -> List[str]:
        """Get cost optimization recommendations."""
        if self.server_summary_df.empty:
            return ["No metrics available for optimization analysis."]

        # Convert server summary to dict format for recommendations
        server_metrics = {}
        for _, row in self.server_summary_df.iterrows():
            server_metrics[row['server']] = {
                'total_cost': row['total_cost'],
                'avg_output_tokens': row['avg_output_tokens'],
                'status': row['status']
            }

        return get_optimization_recommendations(server_metrics)

    def get_timeline_data(self) -> pd.DataFrame:
        """
        Get timeline of tool calls for visualization.

        Returns:
            DataFrame with timestamp, cumulative_cost, cumulative_tokens
        """
        if self.tool_calls_df.empty or 'timestamp' not in self.tool_calls_df.columns:
            return pd.DataFrame()

        timeline = self.tool_calls_df[['timestamp', 'cost', 'input_tokens', 'output_tokens']].copy()
        timeline = timeline.sort_values('timestamp')

        # Calculate cumulative metrics
        timeline['cumulative_cost'] = timeline['cost'].cumsum()
        timeline['cumulative_tokens'] = (
            timeline['input_tokens'] + timeline['output_tokens']
        ).cumsum()

        return timeline

    def get_server_comparison(self) -> pd.DataFrame:
        """
        Get server comparison for visualization.

        Returns:
            DataFrame with server, metric, value for faceted charts
        """
        if self.server_summary_df.empty:
            return pd.DataFrame()

        # Create long-form data for faceted visualization
        data = []
        for _, row in self.server_summary_df.iterrows():
            server = row['server']
            data.extend([
                {'server': server, 'metric': 'Total Cost ($)', 'value': row['total_cost']},
                {'server': server, 'metric': 'Tool Calls', 'value': row['tool_calls']},
                {'server': server, 'metric': 'Avg Latency (ms)', 'value': row['avg_latency_ms']},
                {'server': server, 'metric': 'Input Tokens (K)', 'value': row['input_tokens'] / 1000},
                {'server': server, 'metric': 'Output Tokens (K)', 'value': row['output_tokens'] / 1000},
            ])

        return pd.DataFrame(data)

    def export_summary_report(self, output_path: str):
        """
        Export summary report as text file.

        Args:
            output_path: Path to save report
        """
        stats = self.get_summary_stats()
        top_tools = self.get_top_cost_tools(5)
        recommendations = self.get_optimization_recommendations()

        with open(output_path, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("MCP SERVER COST & PERFORMANCE REPORT\n")
            f.write("=" * 80 + "\n\n")

            f.write(f"Workflow: {stats['workflow_name']}\n")
            f.write(f"Patient: {stats['patient_id']}\n")
            f.write(f"Date: {stats['analysis_date']}\n")
            f.write(f"Duration: {stats['total_duration_seconds']}s\n\n")

            f.write("-" * 80 + "\n")
            f.write("SUMMARY STATISTICS\n")
            f.write("-" * 80 + "\n")
            f.write(f"Total Tool Calls: {stats['total_tool_calls']}\n")
            f.write(f"  - Successful: {stats['successful_calls']}\n")
            f.write(f"  - Dry Run: {stats['dry_run_calls']}\n")
            f.write(f"Servers Used: {stats['servers_used']}\n")
            f.write(f"Insights Generated: {stats['insights_count']}\n\n")

            f.write("-" * 80 + "\n")
            f.write("COST BREAKDOWN\n")
            f.write("-" * 80 + "\n")
            f.write(f"Model: {stats['model']}\n")
            f.write(f"Total Tokens: {stats['total_tokens']:,}\n")
            f.write(f"  - Input: {stats['total_input_tokens']:,}\n")
            f.write(f"  - Output: {stats['total_output_tokens']:,}\n")
            f.write(f"Total Cost: ${stats['total_cost']:.4f}\n")
            f.write(f"  - Input Cost: ${stats['input_cost']:.4f}\n")
            f.write(f"  - Output Cost: ${stats['output_cost']:.4f}\n")
            f.write(f"Cost per Tool Call: ${stats['cost_per_tool_call']:.4f}\n")
            f.write(f"Cost per Insight: ${stats['cost_per_insight']:.4f}\n\n")

            f.write("-" * 80 + "\n")
            f.write("TOP 5 MOST EXPENSIVE TOOL CALLS\n")
            f.write("-" * 80 + "\n")
            for idx, row in top_tools.iterrows():
                f.write(f"{row['server']}.{row['tool']}\n")
                f.write(f"  Tokens: {row['input_tokens']:,} in / {row['output_tokens']:,} out\n")
                f.write(f"  Cost: ${row['cost']:.4f}\n")
                f.write(f"  Latency: {row['latency_ms']:,}ms\n\n")

            f.write("-" * 80 + "\n")
            f.write("OPTIMIZATION RECOMMENDATIONS\n")
            f.write("-" * 80 + "\n")
            for rec in recommendations:
                f.write(f"• {rec}\n")

            f.write("\n" + "=" * 80 + "\n")


def load_sample_metrics(sample_name: str = "patientone_workflow") -> MetricsAggregator:
    """
    Load sample metrics from ui/dashboard/sample_data directory.

    Args:
        sample_name: Name of sample file (without .yaml extension)

    Returns:
        MetricsAggregator instance

    Example:
        >>> aggregator = load_sample_metrics("patientone_workflow")
        >>> print(aggregator.get_cost_breakdown())
    """
    dashboard_dir = Path(__file__).parent
    sample_file = dashboard_dir / "sample_data" / f"{sample_name}.yaml"

    if not sample_file.exists():
        raise FileNotFoundError(f"Sample metrics file not found: {sample_file}")

    return MetricsAggregator(str(sample_file))


if __name__ == "__main__":
    # Demo usage
    print("=== Metrics Aggregator Demo ===\n")

    # Load sample metrics
    aggregator = load_sample_metrics("patientone_workflow")

    # Print summary stats
    stats = aggregator.get_summary_stats()
    print(f"Workflow: {stats['workflow_name']}")
    print(f"Total Cost: ${stats['total_cost']:.4f}")
    print(f"Cost per Insight: ${stats['cost_per_insight']:.4f}")
    print(f"Insights Generated: {stats['insights_count']}\n")

    # Show top 3 expensive tools
    print("Top 3 Most Expensive Tool Calls:")
    top_tools = aggregator.get_top_cost_tools(3)
    for idx, row in top_tools.iterrows():
        print(f"  {row['server']}.{row['tool']}: ${row['cost']:.4f}")
    print()

    # Show optimization recommendations
    print("Optimization Recommendations:")
    for rec in aggregator.get_optimization_recommendations():
        print(f"  • {rec}")
